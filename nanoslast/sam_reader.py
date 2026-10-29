import shutil
import subprocess
from contextlib import contextmanager
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
    # Dorado splits some reads into several subreads. A subread's signal is the part
    # [signal_start, signal_start + num_samples) of the parent read's signal in POD5/FAST5,
    # and its move table refers to that part. For other reads parent_id == read_id.
    parent_id: str = ""
    signal_start: int = 0
    num_samples: Optional[int] = None
    barcode: str = ""

    @property
    def signal_read_id(self):
        """ID under which the read's signal is stored in POD5/FAST5."""
        return self.parent_id or self.read_id


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def is_bam(filename):
    return str(filename).lower().endswith(".bam")


@contextmanager
def open_alignments(filename, threads=1):
    """Lines of a SAM file, or of a BAM file decoded by samtools."""
    if not is_bam(filename):
        with open(filename, "r") as f:
            yield f
        return

    if shutil.which("samtools") is None:
        raise RuntimeError(
            "Reading BAM files needs samtools on PATH. Install samtools, or convert the BAM "
            "to SAM with 'samtools view -h reads.bam > reads.sam' and use the SAM file."
        )
    process = subprocess.Popen(
        ["samtools", "view", "-h", "-@", str(threads), filename],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    try:
        yield process.stdout
    finally:
        process.stdout.close()
        # samtools is stopped when we leave early (e.g. after the reads we needed)
        if process.poll() is None:
            process.terminate()
        _, stderr = process.communicate()
        if process.returncode not in (0, -13, -15) and stderr.strip():  # 0 or stopped by us
            raise RuntimeError(f"samtools view failed:\n{stderr.strip()}")


def iter_sam_records(sam_filename, read_ids=None, with_moves=True):
    """
    Yield one SamRecord per read from a basecalled SAM or BAM file.

    Works with both unaligned and aligned files: secondary and supplementary
    alignments are skipped, and reads aligned to the reverse strand are reverse
    complemented back, so the sequence always matches the order of the move
    table and the raw signal.
    """
    with open_alignments(sam_filename) as sam:
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
            parent_id = ""
            signal_start = 0
            num_samples = None
            barcode = ""

            for field in fields[11:]:
                tag = field[:5]
                if tag == "mv:B:" and with_moves and field.startswith("mv:B:c,"):
                    values = field[7:].split(",")
                    stride = int(values[0])
                    moves = [int(value) for value in values[1:]]
                elif tag == "ts:i:":
                    trimmed_samples = int(field[5:])
                elif tag == "pi:Z:":
                    parent_id = field[5:]
                elif tag == "sp:i:":
                    signal_start = int(field[5:])
                elif tag == "ns:i:":
                    num_samples = int(field[5:])
                elif tag == "BC:Z:":
                    barcode = field[5:]

            if not parent_id:
                # sp/ns only describe the position of a subread in its parent
                signal_start, num_samples = 0, None

            yield SamRecord(read_id, sequence, moves, stride, trimmed_samples,
                            parent_id, signal_start, num_samples, barcode)


def find_sam_record(sam_filename, read_id: str) -> Optional[SamRecord]:
    for record in iter_sam_records(sam_filename, {read_id}):
        return record
    return None


def read_signal(record: SamRecord, parent_signal):
    """The signal of a read (for a subread, its part of the parent's signal)."""
    if not record.parent_id:
        return parent_signal
    end = None if record.num_samples is None else record.signal_start + record.num_samples
    return parent_signal[record.signal_start:end]


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
    (1-based and inclusive, as in BLAST output), in the read's own signal.
    Returns (None, None) if the range is outside the read.
    """
    if not record.moves:
        return None, None

    boundaries = base_boundaries(record)
    if not 1 <= start_base <= stop_base < len(boundaries):
        return None, None

    return boundaries[start_base - 1], boundaries[stop_base]
