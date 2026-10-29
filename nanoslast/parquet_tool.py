import csv
import os
import time
from collections import Counter
from dataclasses import dataclass, field

import pyarrow as pa
import pyarrow.parquet as pq

from .checks import check_sam, check_signals_dir
from .hits import best_hit_per_read, passed_hits, read_hits
from .report import write_report
from .sam_reader import iter_sam_records, read_signal, signal_range
from .signal_reader import SignalIndex

SKIP_REASONS = {
    "read_not_in_sam": "read not found in the SAM/BAM file",
    "no_move_table": "read has no move table (mv tag)",
    "outside_read": "hit positions outside of the basecalled read",
    "signal_not_found": "signal not found in the POD5/FAST5 files",
}


@dataclass
class ExportResult:
    hits_total: int
    reads_total: int
    hits_written: int = 0
    skipped: Counter = field(default_factory=Counter)  # reason -> number of hits
    hits_per_query: Counter = field(default_factory=Counter)  # written hits
    negatives_written: int = 0
    report_file: str = ""

    def summary(self):
        lines = [f"{self.hits_written} of {self.hits_total} hits "
                 f"(in {self.reads_total} reads) exported"]
        if self.negatives_written:
            lines.append(f"{self.negatives_written} reads without hits added as negative samples")
        for reason, n in self.skipped.items():
            lines.append(f"  skipped {n}: {SKIP_REASONS[reason]}")
        return "\n".join(lines)

    def as_report(self):
        return {
            "hits to export": self.hits_total,
            "reads with these hits": self.reads_total,
            "hits exported": self.hits_written,
            "hits exported per query": dict(sorted(self.hits_per_query.items())),
            "hits skipped": {SKIP_REASONS[k]: v for k, v in self.skipped.items()},
        }


def load_export_hits(hits_file, best_only=False):
    hits = passed_hits(read_hits(hits_file))
    return best_hit_per_read(hits) if best_only else hits


def sam_file_signal_search(sam_filename, hits, skipped):
    """
    For each hit, find the part of the read's signal that corresponds to it,
    using the move table from the SAM/BAM file. Hits that cannot be located are
    counted in `skipped` by reason.
    Returns (rows, records): one row per located hit, and the SAM record of each
    of their reads (without sequence and moves) to find the read's signal.
    """
    hits_by_read = {}
    for hit in hits:
        hits_by_read.setdefault(hit.read_id, []).append(hit)

    records = {}
    result = []

    for record in iter_sam_records(sam_filename, hits_by_read.keys()):
        records[record.read_id] = record._replace(sequence="", moves=None)

        for hit in hits_by_read[record.read_id]:
            if not record.moves:
                skipped["no_move_table"] += 1
                continue

            start_signal, stop_signal = signal_range(record, hit.read_start, hit.read_end)
            if start_signal is None:
                skipped["outside_read"] += 1
                continue

            result.append({
                "read_id": hit.read_id,
                "query_id": hit.query_id,
                "sequence": record.sequence[hit.read_start - 1 : hit.read_end],
                "start_idx": hit.read_start,
                "stop_idx": hit.read_end,
                "start_signal": start_signal,
                "stop_signal": stop_signal,
                "direction": 1 if hit.strand == "+" else -1,
                "coverage": hit.coverage,
                "score": hit.score,
                "bitscore": hit.bitscore,
                "at_read_end": hit.at_read_end,
            })

        if len(records) == len(hits_by_read):
            break

    for read_id in hits_by_read.keys() - records.keys():
        skipped["read_not_in_sam"] += len(hits_by_read[read_id])

    return result, records


def load_read_signals(index, records, read_ids):
    """{read_id: signal in pA} of the reads (for subreads, their part of the parent's signal)."""
    parents = index.load({records[r].signal_read_id for r in read_ids})
    return {
        read_id: read_signal(records[read_id], parents[records[read_id].signal_read_id][0])
        for read_id in read_ids
        if records[read_id].signal_read_id in parents
    }


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
                "coverage",
                "score",
                "at_read_end",
            ]
        )

        for hit in signal_ranges:
            if hit["read_id"] not in signals:
                continue

            writer.writerow(
                [
                    hit["read_id"],
                    hit["query_id"],
                    hit["sequence"],
                    hit["start_idx"],
                    hit["stop_idx"],
                    ";".join(map(str, signals[hit["read_id"]])),
                    hit["start_signal"],
                    hit["stop_signal"],
                    hit["direction"],
                    round(hit["coverage"], 4),
                    round(hit["score"], 4),
                    hit["at_read_end"],
                ]
            )
            written += 1

    return written


def dump_signals_parquet(signals, signal_ranges, filename) -> int:
    """Write only the part of the signal that corresponds to the hit."""
    rows = [hit for hit in signal_ranges if hit["read_id"] in signals]

    columns = ["read_id", "query_id", "start_idx", "stop_idx", "start_signal",
               "stop_signal", "direction", "coverage", "score", "at_read_end"]
    table = {column: [hit[column] for hit in rows] for column in columns}
    table["signal"] = [
        signals[hit["read_id"]][hit["start_signal"] : hit["stop_signal"]] for hit in rows
    ]
    pq.write_table(pa.table(table), filename)

    return len(rows)


def export_signals(hits_file, sam_file, signals_dir, output, mode="parquet",
                   best_only=False, progress=print, **dataset_options) -> ExportResult:
    """
    Whole Export workflow: hits → signal positions → signals → CSV/Parquet/HDF5.
    dataset_options are passed to the HDF5 dataset writer (see dataset.write_dataset).
    A report is written next to the output (<output>.report.txt).
    """
    started = time.time()
    progress("Checking inputs...")
    check_sam(sam_file)
    check_signals_dir(signals_dir)

    progress("Loading hits...")
    hits = load_export_hits(hits_file, best_only)
    result = ExportResult(hits_total=len(hits), reads_total=len({hit.read_id for hit in hits}))

    progress("Searching SAM file...")
    signal_ranges, records = sam_file_signal_search(sam_file, hits, result.skipped)

    progress("Indexing POD5/FAST5 files...")
    index = SignalIndex(signals_dir)

    report = {
        "Inputs": {"hits": os.path.abspath(hits_file), "reads": os.path.abspath(sam_file),
                   "signals": os.path.abspath(signals_dir)},
        "Settings": {"format": mode, "hits": "best hit of each read" if best_only else "all hits"},
    }

    if mode == "hdf5":
        from .dataset import write_dataset

        progress("Writing HDF5 dataset...")
        report["HDF5 dataset"] = write_dataset(
            output, hits_file, sam_file, index, signal_ranges, records, result,
            progress=progress, **dataset_options,
        )
    else:
        progress("Loading signals...")
        read_ids = {hit["read_id"] for hit in signal_ranges}
        signals = load_read_signals(index, records, read_ids)

        progress(f"Saving to {mode.upper()}...")
        dump = dump_signals_parquet if mode == "parquet" else dump_signals_csv
        result.hits_written = dump(signals, signal_ranges, output)
        result.hits_per_query.update(h["query_id"] for h in signal_ranges if h["read_id"] in signals)
        result.skipped["signal_not_found"] += sum(1 for h in signal_ranges if h["read_id"] not in signals)

    result.skipped = +result.skipped  # drop zero counts
    report["Results"] = {**result.as_report(), "run time (s)": round(time.time() - started, 1)}
    result.report_file = write_report(output, "export", report)
    return result
