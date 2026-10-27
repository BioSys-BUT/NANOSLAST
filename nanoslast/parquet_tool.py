import csv
import pyarrow as pa
import pyarrow.parquet as pq

from .sam_reader import iter_sam_records, signal_range
from .signal_reader import load_signals


def extract_info_from_tsv(filename):
    """
    Read a BLAST TSV file and keep the best hit (highest bitscore) for each read,
    across all query sequences.
    Returns {read_id: (query_id, sstart, send)}.
    """
    best_hits = {}

    with open(filename, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            read_id = row["sseqid"]
            bitscore = float(row["bitscore"])

            if read_id not in best_hits or bitscore > best_hits[read_id][0]:
                best_hits[read_id] = (
                    bitscore, row["qseqid"], int(row["sstart"]), int(row["send"])
                )

    return {read_id: hit[1:] for read_id, hit in best_hits.items()}


def sam_file_signal_search(sam_filename, hits: dict[str, tuple[str, int, int]]):
    """
    For each read with a BLAST hit, find the part of the raw signal that
    corresponds to the hit, using the move table from the SAM file.
    Reads that are missing from the SAM file or have no move table are skipped.
    """
    result = {}

    for record in iter_sam_records(sam_filename, hits.keys()):
        query_id, start_idx, stop_idx = hits[record.read_id]
        direction = 1

        if start_idx > stop_idx:
            start_idx, stop_idx = stop_idx, start_idx
            direction = -1

        start_signal, stop_signal = signal_range(record, start_idx, stop_idx)
        if start_signal is None:
            continue

        result[record.read_id] = {
            "query_id": query_id,
            "sequence": record.sequence[start_idx - 1 : stop_idx],
            "start_idx": start_idx,
            "stop_idx": stop_idx,
            "start_signal": start_signal,
            "stop_signal": stop_signal,
            "direction": direction,
        }

        if len(result) == len(hits):
            break

    return result


def dump_signals_csv(signals, signal_ranges, filename) -> int:
    """Write the whole signal of each read together with the position of the hit."""
    written = 0

    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "read_id",
                "query_id",
                "sequence",
                "start_idx",
                "stop_idx",
                "signal",
                "start_signal",
                "stop_signal",
                "direction",
            ]
        )

        for read_id, hit in signal_ranges.items():
            if read_id not in signals:
                continue

            signal, _ = signals[read_id]

            writer.writerow(
                [
                    read_id,
                    hit["query_id"],
                    hit["sequence"],
                    hit["start_idx"],
                    hit["stop_idx"],
                    ";".join(map(str, signal)),
                    hit["start_signal"],
                    hit["stop_signal"],
                    hit["direction"],
                ]
            )
            written += 1

    return written


def dump_signals_parquet(signals, signal_ranges, filename) -> int:
    """Write only the part of the signal that corresponds to the hit."""
    read_ids = [read_id for read_id in signal_ranges if read_id in signals]

    table = pa.table(
        {
            "read_id": read_ids,
            "query_id": [signal_ranges[read_id]["query_id"] for read_id in read_ids],
            "signal": [
                signals[read_id][0][
                    signal_ranges[read_id]["start_signal"] : signal_ranges[read_id]["stop_signal"]
                ]
                for read_id in read_ids
            ],
        }
    )
    pq.write_table(table, filename)

    return len(read_ids)


def export_signals(tsv_file, sam_file, signals_dir, output, mode="parquet", progress=print):
    """Whole Export workflow. Returns (number of exported reads, number of reads in TSV)."""
    progress("Loading TSV...")
    hits = extract_info_from_tsv(tsv_file)

    progress("Searching SAM file...")
    signal_ranges = sam_file_signal_search(sam_file, hits)

    progress("Loading signals...")
    signals = load_signals(signals_dir, signal_ranges.keys())

    progress(f"Saving to {mode.upper()}...")
    dump = dump_signals_parquet if mode == "parquet" else dump_signals_csv
    written = dump(signals, signal_ranges, output)

    return written, len(hits)
