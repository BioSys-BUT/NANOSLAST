"""
End-to-end test of the NANOSLAST workflow on synthetic data with known answers.

Each base of each synthetic read has its own signal level (1000 + base index),
so the test can check that every exported signal segment starts exactly at the
first matched base and ends exactly at the last one.

Covers: unaligned and aligned SAM (reverse strand, secondary and supplementary
records, read without a move table), multi-read and single-read FAST5, POD5,
BLAST with several query sequences, CSV and Parquet export, plotting.

Requires BLAST+ on PATH. Run from the repository root:

    python tests/test_workflow.py            # temporary directory, deleted afterwards
    python tests/test_workflow.py OUT_DIR    # keep data, exports and plots in OUT_DIR
"""
import csv
import os
import random
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

from nanoslast.blasting import blast_search  # noqa: E402
from nanoslast.convert_sam_file import convert_sam_to_fasta  # noqa: E402
from nanoslast.create_blastdb import create_blast_db  # noqa: E402
from nanoslast.parquet_tool import export_signals, extract_info_from_tsv  # noqa: E402
from nanoslast.plotting import create_signal_figure, load_read  # noqa: E402
from nanoslast.sam_reader import (  # noqa: E402
    SamRecord, base_boundaries, reverse_complement, signal_range,
)

STRIDE = 6
SAMPLE_RATE = 5000
N_READS = 6


def check_move_table_conversion():
    record = SamRecord("r", "ACG", [1, 0, 0, 1, 0, 1], 5, 10)
    assert base_boundaries(record) == [10, 25, 35, 40]
    assert signal_range(record, 1, 1) == (10, 25)
    assert signal_range(record, 2, 3) == (25, 40)
    assert signal_range(record, 1, 3) == (10, 40)
    assert signal_range(record, 0, 2) == (None, None)
    assert signal_range(record, 2, 4) == (None, None)
    assert reverse_complement("AACGT") == "ACGTT"


def make_reads(rng):
    reads = []
    for _ in range(N_READS):
        n_bases = rng.randint(400, 600)
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
        })
    return reads


def insert_motifs(rng, reads):
    """Returns queries and {read_id: (first base, last base, query name)}."""
    motif1 = "".join(rng.choice("ACGT") for _ in range(40))
    motif2 = "".join(rng.choice("ACGT") for _ in range(30))
    expected = {}

    # motif1 in reads 0-3 forward, in read 4 reverse complemented
    for i, read in enumerate(reads[:5]):
        pos = rng.randint(50, 300)
        motif = motif1 if i < 4 else reverse_complement(motif1)
        read["seq"] = read["seq"][:pos] + motif + read["seq"][pos + len(motif):]
        expected[read["id"]] = (pos + 1, pos + len(motif), "motif1")

    # motif2 in read 5, and in read 0 where motif1 is the better hit
    for i, pos in ((5, 200), (0, 380)):
        reads[i]["seq"] = reads[i]["seq"][:pos] + motif2 + reads[i]["seq"][pos + 30:]
    expected[reads[5]["id"]] = (201, 230, "motif2")

    return {"motif1": motif1, "motif2": motif2}, expected


def sam_line(read, flag=4, seq=None, with_moves=True):
    tags = []
    if with_moves:
        tags = [f"mv:B:c,{STRIDE}," + ",".join(map(str, read["moves"])), f"ts:i:{read['ts']}"]
    seq = read["seq"] if seq is None else seq
    return "\t".join([read["id"], str(flag), "*", "0", "0", "*", "*", "0", "0", seq, "*"] + tags) + "\n"


def write_sam_files(out, reads):
    # unaligned; read 5 has no move table
    with open(f"{out}/unaligned.sam", "w") as f:
        f.write("@HD\tVN:1.6\n")
        for read in reads[:5]:
            f.write(sam_line(read))
        f.write(sam_line(reads[5], with_moves=False))

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


def write_signal_files(signals_dir, reads):
    os.makedirs(f"{signals_dir}/fast5_multi")
    os.makedirs(f"{signals_dir}/fast5_single")
    os.makedirs(f"{signals_dir}/pod5/nested")
    write_fast5_multi(f"{signals_dir}/fast5_multi/batch0.fast5", reads[0:2])
    write_fast5_single(f"{signals_dir}/fast5_single/read2.fast5", reads[2])
    write_fast5_single(f"{signals_dir}/fast5_single/read3.fast5", reads[3])
    write_pod5(f"{signals_dir}/pod5/nested/batch1.pod5", reads[4:6])


def check_segment(segment, first_base, last_base):
    levels = np.unique(segment)
    assert levels[0] == 1000 + first_base - 1, (levels[0], first_base)
    assert levels[-1] == 1000 + last_base - 1, (levels[-1], last_base)
    assert len(levels) == last_base - first_base + 1


def check_workflow(out, sam_kind, reads, queries, expected, signals_dir):
    sam = f"{out}/{sam_kind}.sam"
    work = f"{out}/{sam_kind}"
    os.makedirs(work)

    fasta = f"{work}/reads.fasta"
    convert_sam_to_fasta(sam, fasta)
    fasta_ids = [line[1:].strip() for line in open(fasta) if line.startswith(">")]
    assert len(fasta_ids) == len(set(fasta_ids)) == N_READS, fasta_ids

    create_blast_db(fasta, work)
    n_hits = blast_search(list(queries.items()), f"{work}/blastdb", f"{work}/hits.tsv")
    hits = extract_info_from_tsv(f"{work}/hits.tsv")
    assert {read_id: hit[0] for read_id, hit in hits.items()} == {
        read_id: exp[2] for read_id, exp in expected.items()
    }, hits

    # read 5 has no move table in the unaligned SAM and must be skipped
    n_exported = 5 if sam_kind == "unaligned" else 6
    reverse_read = reads[4]["id"]

    out_csv = f"{work}/export.csv"
    assert export_signals(f"{work}/hits.tsv", sam, signals_dir, out_csv, "csv", lambda m: None) == (n_exported, 6)
    csv.field_size_limit(10**9)
    rows = list(csv.DictReader(open(out_csv)))
    assert len(rows) == n_exported
    for row in rows:
        first_base, last_base, query_id = expected[row["read_id"]]
        assert row["query_id"] == query_id
        assert (int(row["start_idx"]), int(row["stop_idx"])) == (first_base, last_base)
        assert row["direction"] == ("-1" if row["read_id"] == reverse_read else "1")
        motif = queries[query_id]
        assert row["sequence"] == (reverse_complement(motif) if row["read_id"] == reverse_read else motif)
        signal = np.array(row["signal"].split(";"), dtype=float)
        check_segment(signal[int(row["start_signal"]):int(row["stop_signal"])], first_base, last_base)

    out_parquet = f"{work}/export.parquet"
    assert export_signals(f"{work}/hits.tsv", sam, signals_dir, out_parquet, "parquet", lambda m: None) == (n_exported, 6)
    table = pq.read_table(out_parquet).to_pydict()
    assert len(table["read_id"]) == n_exported
    for read_id, query_id, segment in zip(table["read_id"], table["query_id"], table["signal"]):
        first_base, last_base, expected_query = expected[read_id]
        assert query_id == expected_query
        check_segment(segment, first_base, last_base)

    print(f"{sam_kind} SAM: {n_hits} BLAST hits, {n_exported} CSV and Parquet segments exact")


def check_plotting(out, reads, signals_dir, save):
    for i in (0, 2, 4):
        record, signal, sampling_rate = load_read(f"{out}/aligned.sam", signals_dir, reads[i]["id"])
        assert sampling_rate == SAMPLE_RATE
        # reverse strand read must be turned back to the order of the move table
        assert record.sequence == reads[i]["seq"]
        fig = create_signal_figure(signal, sampling_rate, record, reads[i]["id"])
        if save:
            fig.savefig(f"{out}/plot_read{i}.png", dpi=80)

    # first 12 bases of read 0, labelled, to check by eye that colours match the signal steps
    record, signal, sampling_rate = load_read(f"{out}/unaligned.sam", signals_dir, reads[0]["id"])
    fig = create_signal_figure(signal, sampling_rate, record, reads[0]["id"], signal_range(record, 1, 12))
    if save:
        ax = fig.axes[0]
        for start, nucleotide in zip(base_boundaries(record), record.sequence[:12]):
            ax.text(start / sampling_rate * 1000, 1.01, nucleotide, transform=ax.get_xaxis_transform())
        fig.savefig(f"{out}/plot_zoom.png", dpi=80)

    print("plotting OK")


def main():
    keep = len(sys.argv) > 1
    out = os.path.abspath(sys.argv[1]) if keep else tempfile.mkdtemp(prefix="nanoslast_test_")
    if keep:
        shutil.rmtree(out, ignore_errors=True)
        os.makedirs(out)

    try:
        check_move_table_conversion()
        print("move table conversion OK")

        rng = random.Random(1)
        reads = make_reads(rng)
        queries, expected = insert_motifs(rng, reads)
        with open(f"{out}/queries.fasta", "w") as f:
            f.write("".join(f">{name}\n{seq}\n" for name, seq in queries.items()))

        write_sam_files(out, reads)
        signals_dir = f"{out}/signals"
        write_signal_files(signals_dir, reads)

        for sam_kind in ("unaligned", "aligned"):
            check_workflow(out, sam_kind, reads, queries, expected, signals_dir)

        check_plotting(out, reads, signals_dir, save=keep)
    finally:
        if not keep:
            shutil.rmtree(out, ignore_errors=True)

    print("ALL TESTS PASSED")


if __name__ == "__main__":
    main()
