from typing import NamedTuple, Optional

FLAG_REVERSE = 0x10
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800

_COMPLEMENT = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


class SamRecord(NamedTuple):
    read_id: str
    # basecalled sequence in the same orientation as the raw signal
    sequence: str
    moves: Optional[list[int]]
    stride: Optional[int]
    trimmed_samples: int


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def iter_sam_records(sam_filename, read_ids=None, with_moves=True):
    """
    Yield one SamRecord per read from a basecalled SAM file.

    Works with both unaligned and aligned SAM files: secondary and
    supplementary alignments are skipped, and reads aligned to the reverse
    strand are reverse complemented back, so the sequence always matches
    the order of the move table and the raw signal.
    """
    with open(sam_filename, "r") as sam:
        for line in sam:
            if line.startswith("@"):
                continue

            fields = line.rstrip("\n").split("\t")
            read_id = fields[0]

            if read_ids is not None and read_id not in read_ids:
                continue

            flag = int(fields[1])
            if flag & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY):
                continue

            sequence = fields[9]
            if flag & FLAG_REVERSE:
                sequence = reverse_complement(sequence)

            moves = None
            stride = None
            trimmed_samples = 0

            if with_moves:
                for field in fields[11:]:
                    if field.startswith("mv:B:c,"):
                        values = field[7:].split(",")
                        stride = int(values[0])
                        moves = [int(value) for value in values[1:]]
                    elif field.startswith("ts:i:"):
                        trimmed_samples = int(field[5:])

            yield SamRecord(read_id, sequence, moves, stride, trimmed_samples)


def find_sam_record(sam_filename, read_id: str) -> Optional[SamRecord]:
    for record in iter_sam_records(sam_filename, {read_id}):
        return record
    return None


def base_boundaries(record: SamRecord) -> list[int]:
    """
    Signal sample index at which each base starts, followed by the index
    where the last base ends (len = number of bases + 1).

    Each move table element covers `stride` samples of the signal after the
    first `trimmed_samples` samples; a 1 marks the block where a new base starts.
    """
    starts = [
        record.trimmed_samples + i * record.stride
        for i, move in enumerate(record.moves)
        if move == 1
    ]
    starts.append(record.trimmed_samples + len(record.moves) * record.stride)
    return starts


def signal_range(record: SamRecord, start_base: int, stop_base: int):
    """
    Signal sample range [start, stop) of bases start_base..stop_base
    (1-based and inclusive, as in BLAST output).
    Returns (None, None) if the range is outside the read.
    """
    if not record.moves:
        return None, None

    boundaries = base_boundaries(record)
    if not 1 <= start_base <= stop_base < len(boundaries):
        return None, None

    return boundaries[start_base - 1], boundaries[stop_base]
