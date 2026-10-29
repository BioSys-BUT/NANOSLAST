"""
End-to-end test of the NANOSLAST workflow on synthetic data with known answers.

Each base of each synthetic read has its own signal level (1000 + base index),
so the test can check that every exported signal segment starts exactly at the
first matched base and ends exactly at the last one.

Covers: unaligned and aligned SAM (reverse strand, secondary and supplementary
records, read without a move table), multi-read and single-read FAST5, POD5,
BLAST with several query sequences, joining HSPs into hits, several hits per
read (different queries and two copies of one query), hit filters, partial and
truncated copies, split reads (subreads with an offset in the parent's signal),
barcodes, CSV, Parquet and HDF5 dataset export (negative samples, normalization),
skipped hits and their reasons, run reports, BAM input (when samtools is
available), input checks, hits files of NANOSLAST 1.1.0, plotting, command line
interface, BLAST settings.

Requires BLAST+ on PATH. Run from the repository root (or with nanoslast installed):

    python tests/test_workflow.py            # temporary directory, deleted afterwards
    python tests/test_workflow.py OUT_DIR    # keep data, exports and plots in OUT_DIR
"""
import csv
import json
import os
import random
import re
import subprocess
import shutil
import sys
import tempfile
import uuid
from datetime import datetime, timezone

import h5py
import numpy as np
import pod5 as p5
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from nanoslast.blasting import read_queries, search_reads  # noqa: E402
from nanoslast.checks import (  # noqa: E402
    InputError, check_blast_installed, check_queries, check_sam, check_signals_dir,
)
from nanoslast.hits import HitFilters, hsps_to_hits, passed_hits, read_hits, write_hits  # noqa: E402
from nanoslast.parquet_tool import export_signals  # noqa: E402
from nanoslast.plotting import create_signal_figure, load_read  # noqa: E402
from nanoslast.sam_reader import (  # noqa: E402
    SamRecord, base_boundaries, reverse_complement, signal_range,
)

STRIDE = 6
SAMPLE_RATE = 5000
N_READS = 12  # 7 and 8: subreads of one split read; 9-11: no query
csv.field_size_limit(10**9)


def quiet(message):
    pass


def expect_error(error, function, *args, text=""):
    try:
        function(*args)
    except error as e:
        assert text in str(e), (text, str(e))
        return
    raise AssertionError(f"{function.__name__} did not raise {error.__name__}")


# ---------------------------------------------------------------- unit checks

def check_move_table_conversion():
    record = SamRecord("r", "ACG", [1, 0, 0, 1, 0, 1], 5, 10)
    assert base_boundaries(record) == [10, 25, 35, 40]
    assert signal_range(record, 1, 1) == (10, 25)
    assert signal_range(record, 2, 3) == (25, 40)
    assert signal_range(record, 1, 3) == (10, 40)
    assert signal_range(record, 0, 2) == (None, None)
    assert signal_range(record, 2, 4) == (None, None)
    assert reverse_complement("AACGT") == "ACGTT"


def check_hsp_joining():
    def hsp(qstart, qend, sstart, send, bitscore=100.0):
        return {"qseqid": "g", "sseqid": "r", "pident": "100", "length": str(qend - qstart + 1),
                "qstart": qstart, "qend": qend, "sstart": sstart, "send": send,
                "bitscore": bitscore, "qlen": 100, "slen": 1000}

    # one copy split into two HSPs by an indel -> one hit
    hits = hsps_to_hits([hsp(1, 50, 101, 150), hsp(52, 100, 153, 201)], {"g": 200.0})
    assert len(hits) == 1
    assert (hits[0].read_start, hits[0].read_end, hits[0].n_hsps) == (101, 201, 2)
    assert abs(hits[0].coverage - 0.99) < 1e-9 and hits[0].score == 1.0

    # parts of one copy interrupted by a 300-base insertion: joined only with a larger max_gap
    parts = [hsp(1, 50, 101, 150), hsp(51, 100, 451, 500)]
    assert len(hsps_to_hits(parts, {"g": 200.0})) == 2
    joined = hsps_to_hits(parts, {"g": 200.0}, max_gap=400)
    assert len(joined) == 1 and joined[0].coverage == 1 and (joined[0].read_start, joined[0].read_end) == (101, 500)

    # two complete copies next to each other -> two hits
    hits = hsps_to_hits([hsp(1, 100, 101, 200), hsp(1, 100, 211, 310)], {"g": 100.0})
    assert [(h.read_start, h.read_end) for h in hits] == [(101, 200), (211, 310)]

    # reverse strand, at the read start
    hits = hsps_to_hits([hsp(1, 100, 100, 5)], {"g": 100.0})
    assert (hits[0].strand, hits[0].read_start, hits[0].at_read_end) == ("-", 5, True)


def check_input_checks(out):
    check_queries([("a", "ACGT"), ("b", "acgtn")])
    expect_error(InputError, check_queries, [("a", "ACGT"), ("a", "ACGT")], text="Duplicate")
    expect_error(InputError, check_queries, [("a", "ACGT-X")], text="not nucleotides")
    expect_error(InputError, check_queries, [], text="No query")

    expect_error(InputError, check_sam, f"{out}/no_moves.sam", text="--emit-moves")
    check_sam(f"{out}/no_moves.sam", False)
    expect_error(InputError, check_sam, f"{out}/missing.sam", text="not found")

    os.makedirs(f"{out}/empty_dir")
    expect_error(InputError, check_signals_dir, f"{out}/empty_dir", text="No POD5 or FAST5")

    path = os.environ["PATH"]
    os.environ["PATH"] = ""
    try:
        expect_error(InputError, check_blast_installed, text="BLAST+")
    finally:
        os.environ["PATH"] = path


# ---------------------------------------------------------------- synthetic data

def make_reads(rng):
    reads = []
    for _ in range(N_READS):
        n_bases = rng.randint(450, 600)
        moves = []
        for _ in range(n_bases):
            moves += [1] + [0] * rng.randint(0, 3)
        ts = rng.randint(50, 200)

        # 0 in the trimmed part, base k has level 1000 + k, then a -1 tail
        signal = [0] * ts
        for k, base_moves in enumerate(np.split(np.array(moves), np.flatnonzero(moves)[1:])):
            signal += [1000 + k] * (len(base_moves) * STRIDE)
        signal += [-1] * rng.randint(0, 50)

        reads.append({
            "id": str(uuid.UUID(int=rng.getrandbits(128))),
            "seq": "".join(rng.choice("ACGT") for _ in range(n_bases)),
            "moves": moves,
            "ts": ts,
            "signal": np.array(signal, dtype=np.int16),
            "barcode": f"SQK-RBK114-24_barcode0{1 + len(reads) % 2}",
        })
    reads[11]["barcode"] = ""  # unclassified

    # reads 7 and 8 are subreads of one parent read: parent signal = 7 + 37 junk samples + 8 + tail
    parent_id = str(uuid.UUID(int=rng.getrandbits(128)))
    reads[7].update(parent=parent_id, sp=0)
    reads[8].update(parent=parent_id, sp=len(reads[7]["signal"]) + 37)
    parent_signal = np.concatenate([reads[7]["signal"], np.full(37, -5, np.int16),
                                    reads[8]["signal"], np.full(20, -5, np.int16)])
    parent = {"id": parent_id, "signal": parent_signal}
    return reads, parent


def insert(read, pos, sequence):
    read["seq"] = read["seq"][:pos] + sequence + read["seq"][pos + len(sequence):]
    return (pos + 1, pos + len(sequence))


def insert_motifs(rng, reads):
    """
    Returns queries, hits expected with the default filters and hits expected
    only with lowered filters: {(read_id, first base, last base): (query, strand)}.
    """
    motif1 = "".join(rng.choice("ACGT") for _ in range(120))
    motif2 = "".join(rng.choice("ACGT") for _ in range(90))
    expected, expected_low = {}, {}

    def add(table, i, pos, query, sequence, strand="+"):
        table[(reads[i]["id"], *insert(reads[i], pos, sequence))] = (query, strand)

    add(expected, 0, 60, "motif1", motif1)                       # two different queries
    add(expected, 0, 300, "motif2", motif2)
    add(expected, 1, 100, "motif1", motif1)
    add(expected, 2, 50, "motif1", motif1)                       # two copies of one query
    add(expected, 2, 210, "motif1", motif1)
    add(expected_low, 3, 150, "motif1", motif1[:84])             # 70 % of the query
    add(expected, 4, 200, "motif1", reverse_complement(motif1), "-")
    add(expected, 5, 200, "motif2", motif2)                      # no move table in unaligned SAM
    add(expected_low, 6, 0, "motif1", motif1[30:])               # cut off by the read start
    add(expected, 8, 120, "motif1", motif1)                      # in the second subread of a split read

    return {"motif1": motif1, "motif2": motif2}, expected, expected_low


def sam_line(read, flag=4, seq=None, with_moves=True):
    tags = []
    if with_moves:
        tags = [f"mv:B:c,{STRIDE}," + ",".join(map(str, read["moves"])), f"ts:i:{read['ts']}"]
    if read["barcode"]:
        tags.append(f"BC:Z:{read['barcode']}")
    if "parent" in read:
        tags += [f"pi:Z:{read['parent']}", f"sp:i:{read['sp']}", f"ns:i:{len(read['signal'])}"]
    seq = read["seq"] if seq is None else seq
    return "\t".join([read["id"], str(flag), "*", "0", "0", "*", "*", "0", "0", seq, "*"] + tags) + "\n"


def write_sam_files(out, reads):
    # unaligned; read 5 has no move table
    with open(f"{out}/unaligned.sam", "w") as f:
        f.write("@HD\tVN:1.6\n")
        for i, read in enumerate(reads):
            f.write(sam_line(read, with_moves=i != 5))

    # aligned; reads 1 and 4 on the reverse strand, plus secondary/supplementary records
    with open(f"{out}/aligned.sam", "w") as f:
        f.write("@HD\tVN:1.6\n@SQ\tSN:chr\tLN:10000\n")
        for i, read in enumerate(reads):
            if i in (1, 4):
                f.write(sam_line(read, flag=16, seq=reverse_complement(read["seq"])))
                f.write(sam_line(read, flag=2048 + 16, seq=reverse_complement(read["seq"][:100]), with_moves=False))
            else:
                f.write(sam_line(read, flag=0))
                f.write(sam_line(read, flag=256, seq="*", with_moves=False))

    # basecalled without --emit-moves
    with open(f"{out}/no_moves.sam", "w") as f:
        for read in reads:
            f.write(sam_line(read, with_moves=False))


def write_channel_id(group):
    group.attrs.update(digitisation=8192.0, offset=0.0, range=8192.0, sampling_rate=float(SAMPLE_RATE))


def write_fast5_multi(path, reads):
    with h5py.File(path, "w") as f:
        for read in reads:
            group = f.create_group(f"read_{read['id']}")
            raw = group.create_group("Raw")
            raw.attrs["read_id"] = read["id"].encode()
            raw.create_dataset("Signal", data=read["signal"])
            write_channel_id(group.create_group("channel_id"))


def write_fast5_single(path, read):
    with h5py.File(path, "w") as f:
        group = f.create_group("Raw/Reads/Read_1")
        group.attrs["read_id"] = read["id"].encode()
        group.create_dataset("Signal", data=read["signal"])
        write_channel_id(f.create_group("UniqueGlobalKey/channel_id"))


def write_pod5(path, reads):
    now = datetime.now(timezone.utc)
    run_info = p5.RunInfo(
        acquisition_id="acq", acquisition_start_time=now, adc_max=4095, adc_min=-4096,
        context_tags={}, experiment_name="", flow_cell_id="", flow_cell_product_code="",
        protocol_name="", protocol_run_id="", protocol_start_time=now, sample_id="",
        sample_rate=SAMPLE_RATE, sequencing_kit="", sequencer_position="",
        sequencer_position_type="", software="", system_name="", system_type="", tracking_id={},
    )
    with p5.Writer(path) as writer:
        for n, read in enumerate(reads):
            writer.add_read(p5.Read(
                read_id=uuid.UUID(read["id"]), read_number=n, start_sample=0, median_before=100.0,
                end_reason=p5.EndReason.from_reason_with_default_forced(p5.EndReasonEnum.SIGNAL_POSITIVE),
                calibration=p5.Calibration(offset=0.0, scale=1.0),
                pore=p5.Pore(channel=1, well=1, pore_type="x"), run_info=run_info,
                signal=read["signal"],
            ))


def write_signal_files(signals_dir, reads, parent):
    """Read 6 has no signal file; reads 7 and 8 are stored as their parent read."""
    os.makedirs(f"{signals_dir}/fast5_multi")
    os.makedirs(f"{signals_dir}/fast5_single")
    os.makedirs(f"{signals_dir}/pod5/nested")
    write_fast5_multi(f"{signals_dir}/fast5_multi/batch0.fast5", reads[0:2])
    write_fast5_single(f"{signals_dir}/fast5_single/read2.fast5", reads[2])
    write_fast5_single(f"{signals_dir}/fast5_single/read3.fast5", reads[3])
    write_pod5(f"{signals_dir}/pod5/nested/batch1.pod5", reads[4:6])
    write_pod5(f"{signals_dir}/pod5/batch2.pod5", [parent] + reads[9:12])


# ---------------------------------------------------------------- workflow checks

def check_segment(segment, first_base, last_base):
    levels = np.unique(segment)
    assert levels[0] == 1000 + first_base - 1, (levels[0], first_base)
    assert levels[-1] == 1000 + last_base - 1, (levels[-1], last_base)
    assert len(levels) == last_base - first_base + 1


def found(hits):
    return {(h.read_id, h.read_start, h.read_end): (h.query_id, h.strand) for h in hits}


def check_rows(rows, expected, queries, reads, signal_column):
    for row in rows:
        key = (row["read_id"], int(row["start_idx"]), int(row["stop_idx"]))
        query_id, strand = expected[key]
        assert row["query_id"] == query_id
        assert int(row["direction"]) == (1 if strand == "+" else -1)
        if "sequence" in row:
            read = next(r for r in reads if r["id"] == row["read_id"])
            assert row["sequence"] == read["seq"][key[1] - 1: key[2]]
        segment = signal_column(row)
        check_segment(segment, key[1], key[2])


def check_workflow(out, sam_kind, reads, queries, expected, expected_low, signals_dir):
    sam = f"{out}/{sam_kind}.sam"
    work = f"{out}/{sam_kind}"
    queries = list(queries.items())

    # default filters: complete copies only
    hits_file = f"{work}/hits.tsv"
    n_hits, n_reads = search_reads(sam, queries, hits_file, db_dir=f"{work}/db", progress=quiet)
    hits = passed_hits(read_hits(hits_file))
    assert found(hits) == expected, (found(hits), expected)
    assert open(hits_file).readline() == "# queries: motif1;motif2\n"
    # hits below the filters are kept in the file, marked as not passed
    assert found([h for h in read_hits(hits_file) if not h.passed]) == expected_low
    assert (n_hits, n_reads) == (len(expected), len({k[0] for k in expected}))
    assert all(h.coverage == 1 and h.score == 1 and not h.at_read_end for h in hits)

    # lowered filters: also the partial copy and the copy cut off by the read start
    low_file = f"{work}/hits_low.tsv"
    search_reads(sam, queries, low_file, db_dir=f"{work}/db",
                 filters=HitFilters(min_coverage=0.6, min_score=0.6), progress=quiet)
    low = passed_hits(read_hits(low_file))
    assert found(low) == {**expected, **expected_low}
    assert [h.at_read_end for h in low if h.read_id == reads[6]["id"]] == [True]

    # export: all hits, one row per hit; read 5 has no move table in the unaligned SAM
    no_moves = 1 if sam_kind == "unaligned" else 0
    out_csv = f"{work}/export.csv"
    result = export_signals(hits_file, sam, signals_dir, out_csv, "csv", progress=quiet)
    assert result.hits_written == len(expected) - no_moves
    assert dict(result.skipped) == ({"no_move_table": 1} if no_moves else {})
    rows = list(csv.DictReader(open(out_csv)))
    check_rows(rows, expected, queries, reads,
               lambda r: np.array(r["signal"].split(";"), dtype=float)[int(r["start_signal"]):int(r["stop_signal"])])

    out_parquet = f"{work}/export.parquet"
    result = export_signals(hits_file, sam, signals_dir, out_parquet, "parquet", progress=quiet)
    table = pq.read_table(out_parquet).to_pylist()
    assert len(table) == len(expected) - no_moves
    check_rows(table, expected, queries, reads, lambda r: r["signal"])

    # only the best hit of each read
    result = export_signals(hits_file, sam, signals_dir, f"{work}/best.parquet", "parquet",
                            best_only=True, progress=quiet)
    best = pq.read_table(f"{work}/best.parquet").to_pylist()
    assert len(best) == len({k[0] for k in expected}) - no_moves
    assert len({r["read_id"] for r in best}) == len(best)
    assert next(r for r in best if r["read_id"] == reads[0]["id"])["query_id"] == "motif1"

    # skipped hits: no signal (read 6), not in SAM (made-up read)
    extra = low + [low[0].__class__(**{**low[0].__dict__, "read_id": "not-a-read"})]
    write_hits(extra, f"{work}/hits_extra.tsv")
    result = export_signals(f"{work}/hits_extra.tsv", sam, signals_dir, f"{work}/extra.parquet",
                            "parquet", progress=quiet)
    assert result.skipped["signal_not_found"] == 1
    assert result.skipped["read_not_in_sam"] == 1
    assert result.hits_written == len(low) - 1 - no_moves

    # hits file written by NANOSLAST 1.1.0 (one BLAST HSP per row)
    read_id, first, last = next(iter(expected))
    with open(f"{work}/hits_v1.1.tsv", "w") as f:
        f.write("qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n")
        f.write(f"motif1\t{read_id}\t100.0\t120\t0\t0\t1\t120\t{first}\t{last}\t0.0\t222.0\n")
    result = export_signals(f"{work}/hits_v1.1.tsv", sam, signals_dir, f"{work}/v1.1.parquet",
                            "parquet", progress=quiet)
    assert result.hits_written == 1

    print(f"{sam_kind} SAM: {len(hits)} hits ({len(low)} with lowered filters), "
          f"CSV and Parquet segments exact, skipped hits reported")


def check_dataset(out, reads, expected, signals_dir):
    """HDF5 training dataset from the aligned SAM (all reads have move tables)."""
    sam = f"{out}/aligned.sam"
    hits_file = f"{out}/aligned/hits.tsv"
    read = {r["id"]: r for r in reads}
    positives = {k[0] for k in expected}            # reads 0, 1, 2, 4, 5, 8
    # reads without any hit (3 and 6 have hits below the filters, so they are not negatives):
    # barcode01 (even reads): 10; barcode02 (odd): 7 (a subread), 9; unclassified: 11
    # neg_ratio 1 -> barcode01: min(1, 4 positives), barcode02: min(2, 2), unclassified: 0 positives
    expected_negatives = {reads[i]["id"] for i in (7, 9, 10)}

    for normalization in ("pa", "zscore"):
        output = f"{out}/aligned/dataset_{normalization}.h5"
        # gene names come from the "# queries" line of the hits file
        result = export_signals(hits_file, sam, signals_dir, output, "hdf5", progress=quiet,
                                neg_ratio=1.0, seed=7, normalization=normalization)
        assert result.hits_written == len(expected) and result.negatives_written == 3

        with h5py.File(output) as f:
            assert f.attrs["run_id"] == "aligned" and f.attrs["normalization"] == normalization
            assert json.loads(f.attrs["gene_names"]) == ["motif1", "motif2"]
            assert f.attrs["n_signals"] == len(positives) + 3
            for name, dtype in [("read_ids", "S64"), ("run_ids", "S64"), ("barcodes", "S32"),
                                ("gene_labels", "S500"), ("start_signals", "S200"), ("stop_signals", "S200"),
                                ("bitscores", "S200"), ("align_ratios", "S200"), ("strands", "S200")]:
                assert f[name].dtype == np.dtype(dtype), (name, f[name].dtype)
            assert f["is_negative"].dtype == np.int8 and f["signals"].dtype.metadata["vlen"] == np.float32

            rows = {f["read_ids"][i].decode(): i for i in range(f.attrs["n_signals"])}
            assert set(rows) == positives | expected_negatives

            for read_id, i in rows.items():
                signal = f["signals"][i]
                labels = f["gene_labels"][i].decode().split(";")
                starts = [int(x) for x in f["start_signals"][i].decode().split(";")]
                stops = [int(x) for x in f["stop_signals"][i].decode().split(";")]
                assert len(signal) == len(read[read_id]["signal"])
                assert f["barcodes"][i].decode() == (read[read_id]["barcode"].split("_")[-1] or "unclassified")
                if normalization == "pa":
                    assert np.array_equal(signal, read[read_id]["signal"].astype(np.float32))
                else:
                    assert abs(signal.mean()) < 1e-4 and abs(signal.std() - 1) < 1e-4

                if read_id in expected_negatives:
                    assert f["is_negative"][i] == 1 and labels == ["no_gene"]
                    assert (starts, stops) == ([0], [len(signal)])
                    assert f["bitscores"][i] == b"-1.0" and f["align_ratios"][i] == b"-1.0"
                    continue

                assert f["is_negative"][i] == 0
                hits = sorted((first, last, q, strand) for (r, first, last), (q, strand) in expected.items()
                              if r == read_id)
                assert labels == [h[2] for h in hits]
                assert f["strands"][i].decode().split(";") == [h[3] for h in hits]
                assert [float(x) for x in f["align_ratios"][i].decode().split(";")] == [1.0] * len(hits)
                if normalization == "pa":
                    for (first, last, _, _), start, stop in zip(hits, starts, stops):
                        check_segment(signal[start:stop], first, last)

        report = open(result.report_file).read()
        assert result.report_file == output + ".report.txt"
        # per barcode: reads with hits / negative samples / reads without hits available
        assert re.search(r"barcode01:\s+4 / 1 / 1", report) and re.search(r"barcode02:\s+2 / 2 / 2", report)
        assert re.search(r"unclassified:\s+0 / 0 / 1", report)
        assert re.search(rf"hits exported:\s+{len(expected)}\n", report)

    print("HDF5 dataset OK (format, positions, split read, negatives, normalization, report)")


def check_bam(out, expected, signals_dir):
    if shutil.which("samtools") is None:
        print("BAM input: samtools not found, skipped")
        return
    bam = f"{out}/aligned.bam"
    subprocess.run(["samtools", "view", "-b", "-o", bam, f"{out}/aligned.sam"], check=True)
    search_reads(bam, read_queries(f"{out}/queries.fasta"), f"{out}/bam_hits.tsv", progress=quiet)
    assert found(passed_hits(read_hits(f"{out}/bam_hits.tsv"))) == expected
    result = export_signals(f"{out}/bam_hits.tsv", bam, signals_dir, f"{out}/bam.parquet", "parquet", progress=quiet)
    assert result.hits_written == len(expected)
    print("BAM input OK")


def check_plotting(out, reads, signals_dir, save):
    for i in (0, 2, 4):
        record, signal, sampling_rate = load_read(f"{out}/aligned.sam", signals_dir, reads[i]["id"])
        assert sampling_rate == SAMPLE_RATE
        # reverse strand read must be turned back to the order of the move table
        assert record.sequence == reads[i]["seq"]
        fig = create_signal_figure(signal, sampling_rate, record, reads[i]["id"])
        if save:
            fig.savefig(f"{out}/plot_read{i}.png", dpi=80)

    # second subread of a split read: its own part of the parent's signal
    record, signal, _ = load_read(f"{out}/aligned.sam", signals_dir, reads[8]["id"])
    assert record.parent_id and np.array_equal(signal, reads[8]["signal"])

    # first 12 bases of read 0, labelled, to check by eye that colours match the signal steps
    record, signal, sampling_rate = load_read(f"{out}/unaligned.sam", signals_dir, reads[0]["id"])
    fig = create_signal_figure(signal, sampling_rate, record, reads[0]["id"], signal_range(record, 1, 12))
    if save:
        ax = fig.axes[0]
        for start, nucleotide in zip(base_boundaries(record), record.sequence[:12]):
            ax.text(start / sampling_rate * 1000, 1.01, nucleotide, transform=ax.get_xaxis_transform())
        fig.savefig(f"{out}/plot_zoom.png", dpi=80)

    print("plotting OK")


def check_cli(out, reads, expected, signals_dir):
    from nanoslast.cli import main as cli

    work = f"{out}/cli"
    os.makedirs(work)
    sam = f"{out}/aligned.sam"
    queries = f"{out}/queries.fasta"

    # run = search + export; the database is kept and reused by the second search
    assert cli(["run", "--sam", sam, "--signals", signals_dir, "--query", queries,
                "--db-dir", f"{work}/db", "--hits-output", f"{work}/hits.tsv", "-o", f"{work}/run.parquet"]) == 0
    assert len(pq.read_table(f"{work}/run.parquet")) == len(expected)

    assert cli(["search", "--sam", sam, "--query", queries, "--db-dir", f"{work}/db",
                "--min-coverage", "0.6", "--min-score", "0.6", "-o", f"{work}/hits_low.tsv"]) == 0
    assert len(passed_hits(read_hits(f"{work}/hits_low.tsv"))) == len(expected) + 2

    assert cli(["export", "--sam", sam, "--signals", signals_dir, "--hits", f"{work}/hits.tsv",
                "--best-only", "-o", f"{work}/best.csv"]) == 0
    assert len(list(csv.DictReader(open(f"{work}/best.csv")))) == len({k[0] for k in expected})

    assert cli(["plot", "--sam", sam, "--signals", signals_dir, "--read-id", reads[0]["id"],
                "-o", f"{work}/plot.png"]) == 0
    assert os.path.getsize(f"{work}/plot.png") > 0

    assert cli(["run", "--sam", sam, "--signals", signals_dir, "--query", queries, "--db-dir", f"{work}/db",
                "--neg-ratio", "1", "--run-id", "test_run", "-o", f"{work}/dataset.h5"]) == 0
    with h5py.File(f"{work}/dataset.h5") as f:
        assert f.attrs["run_id"] == "test_run" and json.loads(f.attrs["gene_names"]) == ["motif1", "motif2"]
        assert f.attrs["n_signals"] == len({k[0] for k in expected}) + 3

    # BLAST settings: a 20-nt sequence is too short for megablast (word size 28), found with blastn-short
    short = read_queries(queries)[1][1][:20]          # start of motif2, in reads 0 and 5
    assert cli(["search", "--sam", sam, "--sequence", short, "--db-dir", f"{work}/db", "-o", f"{work}/short.tsv"]) == 0
    assert not passed_hits(read_hits(f"{work}/short.tsv"))
    assert cli(["search", "--sam", sam, "--sequence", short, "--db-dir", f"{work}/db", "--task", "blastn-short",
                "--blast-options=-word_size 7 -reward 1 -penalty -3 -max_target_seqs 10",
                "-o", f"{work}/short.tsv"]) == 0
    assert {h.read_id for h in passed_hits(read_hits(f"{work}/short.tsv"))} == {reads[0]["id"], reads[5]["id"]}
    search_report = open(f"{work}/short.tsv.report.txt").read()
    assert re.search(r"task:\s+blastn-short", search_report) and "-word_size 7 -reward 1" in search_report
    assert re.search(r"max. reads per query:\s+set in the other options", search_report)
    assert cli(["search", "--sam", sam, "--sequence", short, "--blast-options=-outfmt 5", "-o", f"{work}/x.tsv"]) == 1

    # errors are reported with a non-zero exit code, the input checks run before BLAST
    assert cli(["run", "--sam", f"{out}/no_moves.sam", "--signals", signals_dir, "--query", queries,
                "-o", f"{work}/x.parquet"]) == 1
    assert cli(["plot", "--sam", sam, "--signals", signals_dir, "--read-id", "missing", "-o", f"{work}/x.png"]) == 1
    assert cli(["export", "--sam", sam, "--signals", signals_dir, "--hits", f"{work}/hits.tsv",
                "-o", f"{work}/export.txt"]) == 1

    print("command line OK")


def main():
    keep = len(sys.argv) > 1
    out = os.path.abspath(sys.argv[1]) if keep else tempfile.mkdtemp(prefix="nanoslast_test_")
    if keep:
        shutil.rmtree(out, ignore_errors=True)
        os.makedirs(out)

    try:
        check_move_table_conversion()
        check_hsp_joining()
        print("move table conversion and HSP joining OK")

        rng = random.Random(1)
        reads, parent = make_reads(rng)
        queries, expected, expected_low = insert_motifs(rng, reads)
        with open(f"{out}/queries.fasta", "w") as f:
            f.write("".join(f">{name}\n{seq}\n" for name, seq in queries.items()))

        write_sam_files(out, reads)
        signals_dir = f"{out}/signals"
        write_signal_files(signals_dir, reads, parent)

        check_input_checks(out)
        print("input checks OK")

        for sam_kind in ("unaligned", "aligned"):
            check_workflow(out, sam_kind, reads, queries, expected, expected_low, signals_dir)

        check_dataset(out, reads, expected, signals_dir)
        check_bam(out, expected, signals_dir)
        check_plotting(out, reads, signals_dir, save=keep)
        check_cli(out, reads, expected, signals_dir)
    finally:
        if not keep:
            shutil.rmtree(out, ignore_errors=True)

    print("ALL TESTS PASSED")


if __name__ == "__main__":
    main()
