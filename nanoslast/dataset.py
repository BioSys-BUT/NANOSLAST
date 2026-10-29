"""
HDF5 training dataset in the format used for NanoResistNet
(e.g. dataset_19G_<run>.h5): one row per read with its whole signal, all genes
(hits) found in it, and reads without any hit as negative samples.

Root attributes: run_id, n_signals, creation_date, gene_names (JSON list),
plus NANOSLAST metadata (nanoslast_version, normalization, neg_ratio, seed,
align_ratios_meaning, source_*).

Datasets, one element per read; several hits of a read are joined by ";":
    signals        vlen float32   whole signal of the read (z-score normalized by default)
    read_ids       S64
    run_ids        S64
    barcodes       S32            e.g. "barcode13", "unclassified"
    gene_labels    S500           query names, e.g. "blaTEM-1A_1_HM749966;sul2_1_AF542061"; "no_gene"
    start_signals  S200           first sample of each hit in `signals`; "0" for no_gene
    stop_signals   S200           end of each hit (exclusive); len(signal) for no_gene
    bitscores      S200           summed BLAST bitscore of each hit; "-1.0" for no_gene
    align_ratios   S200           fraction of the gene covered by each hit; "-1.0" for no_gene
    strands        S200           "+"/"-" for each hit (not in older datasets); "" for no_gene
    is_negative    int8           1 for no_gene reads
"""
import json
import os
import random
from collections import Counter, defaultdict
from datetime import datetime

import h5py
import numpy as np

from . import __version__
from .hits import read_hits, read_query_names
from .sam_reader import iter_sam_records, read_signal

STRING_WIDTHS = {
    "read_ids": 64, "run_ids": 64, "barcodes": 32, "gene_labels": 500,
    "start_signals": 200, "stop_signals": 200, "bitscores": 200,
    "align_ratios": 200, "strands": 200,
}
BATCH_SIZE = 500


def normalize_barcode(tag):
    """Dorado BC tag 'SQK-RBK114-24_barcode13' -> 'barcode13'; no tag -> 'unclassified'."""
    if not tag:
        return "unclassified"
    return tag.rsplit("_", 1)[-1] if "barcode" in tag else tag


def normalize_signal(signal, normalization):
    signal = np.asarray(signal, dtype=np.float64)
    if normalization == "zscore":
        std = signal.std()
        signal = (signal - signal.mean()) / std if std > 0 else signal - signal.mean()
    elif normalization != "pa":
        raise ValueError(f"Unknown normalization: {normalization}")
    return signal.astype(np.float32)


def default_run_id(sam_file):
    name = os.path.basename(sam_file)
    for extension in (".sam", ".bam"):
        if name.lower().endswith(extension):
            return name[: -len(extension)]
    return name


class DatasetWriter:
    """Appends rows to the HDF5 file in batches, so the whole dataset never has to be in memory."""

    def __init__(self, path, widths):
        self.file = h5py.File(path, "w")
        self.n = 0
        self.widths = widths
        self.file.create_dataset("signals", shape=(0,), maxshape=(None,),
                                 dtype=h5py.vlen_dtype(np.float32), chunks=(256,))
        for name, width in widths.items():
            self.file.create_dataset(name, shape=(0,), maxshape=(None,), dtype=f"S{width}", chunks=(4096,))
        self.file.create_dataset("is_negative", shape=(0,), maxshape=(None,), dtype=np.int8, chunks=(4096,))

    def append(self, rows):
        if not rows:
            return
        start, self.n = self.n, self.n + len(rows)
        for name in ["signals", "is_negative", *self.widths]:
            self.file[name].resize((self.n,))
        signals = self.file["signals"]
        for i, row in enumerate(rows):
            signals[start + i] = row["signals"]
        self.file["is_negative"][start:] = [row["is_negative"] for row in rows]
        for name, width in self.widths.items():
            self.file[name][start:] = np.array([row[name].encode() for row in rows], dtype=f"S{width}")

    def close(self, attrs):
        attrs["n_signals"] = self.n
        for key, value in attrs.items():
            self.file.attrs[key] = value
        self.file.close()


def _positive_meta(read_id, hits, run_id, barcode):
    hits = sorted(hits, key=lambda h: h["start_signal"])
    return {
        "read_ids": read_id,
        "run_ids": run_id,
        "barcodes": barcode,
        "gene_labels": ";".join(h["query_id"] for h in hits),
        "start_signals": ";".join(str(h["start_signal"]) for h in hits),
        "stop_signals": ";".join(str(h["stop_signal"]) for h in hits),
        "bitscores": ";".join(f"{h['bitscore']:.2f}" for h in hits),
        "align_ratios": ";".join(f"{h['coverage']:.4f}" for h in hits),
        "strands": ";".join("+" if h["direction"] == 1 else "-" for h in hits),
        "is_negative": 0,
    }


def _negative_meta(read_id, run_id, barcode):
    return {
        "read_ids": read_id, "run_ids": run_id, "barcodes": barcode,
        "gene_labels": "no_gene", "start_signals": "0", "stop_signals": None,
        "bitscores": "-1.0", "align_ratios": "-1.0", "strands": "", "is_negative": 1,
    }


def _negative_candidates(sam_file, reads_with_hits, index):
    """Reads without any BLAST hit (even below the filters) whose signal is available, by barcode."""
    by_barcode = defaultdict(list)
    for record in iter_sam_records(sam_file, with_moves=False):
        if record.read_id in reads_with_hits or record.signal_read_id not in index:
            continue
        by_barcode[normalize_barcode(record.barcode)].append(record._replace(sequence=""))
    return by_barcode


def write_dataset(output, hits_file, sam_file, index, signal_ranges, records, result,
                  run_id=None, gene_names=None, neg_ratio=0.3, seed=42,
                  normalization="zscore", progress=print):
    """
    Write the HDF5 dataset. signal_ranges/records come from
    parquet_tool.sam_file_signal_search; `result` (ExportResult) is updated.
    Returns information for the run report.
    """
    run_id = run_id or default_run_id(sam_file)
    all_hits = read_hits(hits_file)
    gene_names = gene_names or read_query_names(hits_file) or sorted({hit.query_id for hit in all_hits})

    # positive reads: all located hits of each read, if its signal is available
    hits_by_read = defaultdict(list)
    for hit in signal_ranges:
        hits_by_read[hit["read_id"]].append(hit)
    for read_id in [r for r in hits_by_read if records[r].signal_read_id not in index]:
        result.skipped["signal_not_found"] += len(hits_by_read.pop(read_id))

    positives = {
        read_id: _positive_meta(read_id, hits, run_id, normalize_barcode(records[read_id].barcode))
        for read_id, hits in hits_by_read.items()
    }

    # negative reads: per barcode, neg_ratio x the number of positive reads, sampled reproducibly
    progress("Selecting reads without hits as negative samples...")
    candidates = _negative_candidates(sam_file, {hit.read_id for hit in all_hits}, index)
    positives_per_barcode = Counter(meta["barcodes"] for meta in positives.values())
    rng = random.Random(seed)
    negatives = {}
    for barcode, n_positive in sorted(positives_per_barcode.items()):
        pool = sorted(candidates.get(barcode, []))
        for record in rng.sample(pool, min(len(pool), round(neg_ratio * n_positive))):
            negatives[record.read_id] = record
            records[record.read_id] = record

    metas = list(positives.values()) + [
        _negative_meta(read_id, run_id, normalize_barcode(record.barcode))
        for read_id, record in negatives.items()
    ]

    # fixed-width strings as in the existing datasets; wider only if a value would not fit
    widths = {}
    for name, width in STRING_WIDTHS.items():
        longest = max((len((m[name] or "").encode()) for m in metas), default=0)
        widths[name] = max(width, longest)
        if longest > width:
            progress(f"Warning: {name} values up to {longest} bytes, column widened from S{width}")

    writer = DatasetWriter(output, widths)
    try:
        for start in range(0, len(metas), BATCH_SIZE):
            batch = metas[start: start + BATCH_SIZE]
            parents = index.load({records[m["read_ids"]].signal_read_id for m in batch})
            rows = []
            for meta in batch:
                record = records[meta["read_ids"]]
                signal = read_signal(record, parents[record.signal_read_id][0])
                row = dict(meta, signals=normalize_signal(signal, normalization))
                if row["is_negative"]:
                    row["stop_signals"] = str(len(signal))
                rows.append(row)
            writer.append(rows)
            progress(f"  written {min(start + BATCH_SIZE, len(metas))}/{len(metas)} reads")
    finally:
        writer.close({
            "run_id": run_id,
            "creation_date": datetime.now().isoformat(),
            "gene_names": json.dumps(gene_names),
            "nanoslast_version": __version__,
            "normalization": normalization,
            "neg_ratio": neg_ratio,
            "seed": seed,
            "align_ratios_meaning": "fraction of the gene (query) covered by the hit",
            "source_hits": os.path.abspath(hits_file),
            "source_reads": os.path.abspath(sam_file),
        })

    result.hits_written = sum(len(h) for h in hits_by_read.values())
    result.hits_per_query.update(h["query_id"] for hits in hits_by_read.values() for h in hits)
    result.negatives_written = len(negatives)

    negatives_per_barcode = Counter(m["barcodes"] for m in metas if m["is_negative"])
    barcodes = sorted(set(positives_per_barcode) | set(candidates))
    return {
        "run ID": run_id,
        "rows (reads)": len(metas),
        "reads with hits": len(positives),
        "negative samples (reads without any hit)": len(negatives),
        "signal": "z-score of the whole read" if normalization == "zscore" else "pA",
        "negative samples per read with hits (neg_ratio)": neg_ratio,
        "random seed": seed,
        "per barcode (reads with hits / negative samples / reads without hits available)": {
            barcode: f"{positives_per_barcode[barcode]} / {negatives_per_barcode[barcode]} / "
                     f"{len(candidates.get(barcode, []))}"
            for barcode in barcodes
        },
        "gene names": gene_names,
    }
