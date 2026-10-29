import glob
import shlex
import tempfile
import time
from collections import Counter
from dataclasses import dataclass

from Bio import SeqIO
import os

from .checks import check_queries
from .convert_sam_file import convert_sam_to_fasta
from .create_blastdb import create_blast_db, run_blast_program
from .hits import HSP_COLUMNS, MAX_HSP_GAP, HitFilters, hsps_to_hits, read_hsps, write_hits
from .report import write_report

BLAST_TASKS = ["megablast", "dc-megablast", "blastn", "blastn-short"]

# blastn options that NANOSLAST sets itself and that must not be passed in extra_args
RESERVED_OPTIONS = {"-query", "-db", "-out", "-outfmt", "-subject", "-num_threads", "-task", "-evalue"}

# Maximum number of reads reported per query (can be changed with "-max_target_seqs N" in
# the extra options). BLAST's default (500) and the 1000 used before NANOSLAST 1.1.0
# silently dropped hits for sequences present in many reads.
MAX_TARGET_SEQS = 1_000_000


@dataclass
class BlastSettings:
    """
    Settings of blastn. The defaults suit genes and other sequences longer than ~50 nt;
    for short sequences (e.g. barcodes, motifs) use task="blastn-short".
    """
    task: str = "megablast"
    evalue: float = 0.001
    extra_args: str = ""    # other blastn options, e.g. "-word_size 7 -reward 1 -penalty -2"

    def __post_init__(self):
        if self.task not in BLAST_TASKS:
            raise ValueError(f"Unknown BLAST task {self.task}, use one of: {', '.join(BLAST_TASKS)}")
        reserved = RESERVED_OPTIONS & set(self.extra_options())
        if reserved:
            raise ValueError(f"These blastn options are set by NANOSLAST and cannot be in the extra "
                             f"BLAST options: {', '.join(sorted(reserved))} (use the NANOSLAST options instead)")

    def extra_options(self):
        return shlex.split(self.extra_args)

    def scoring_args(self):
        """Options that affect alignments and scores (also used for the self-alignment of queries)."""
        return ["-task", self.task] + self._extra_without_max_target_seqs()

    def search_args(self):
        """Options of the search in the reads."""
        args = ["-evalue", str(self.evalue)] + self.scoring_args()
        options = self.extra_options()
        if "-max_target_seqs" in options:
            i = options.index("-max_target_seqs")
            return args + options[i: i + 2]
        return args + ["-max_target_seqs", str(MAX_TARGET_SEQS)]

    def _extra_without_max_target_seqs(self):
        # -max_target_seqs does not belong to the self-alignment of the queries
        options = self.extra_options()
        if "-max_target_seqs" in options:
            i = options.index("-max_target_seqs")
            del options[i: i + 2]
        return options

    def describe(self):
        return {
            "task": self.task,
            "e-value": self.evalue,
            "other blastn options": self.extra_args or "none",
            "max. reads per query": MAX_TARGET_SEQS if "-max_target_seqs" not in self.extra_options()
                                    else "set in the other options",
        }


def write_queries(queries, fasta_file):
    with open(fasta_file, "w") as f:
        for name, sequence in queries:
            f.write(f">{name}\n{sequence}\n")


def blast_hsps(query_fasta, database_path, output_file, settings, threads=1):
    """Run blastn and save all HSPs (columns HSP_COLUMNS, no header)."""

    # BLAST+ cannot use a database path containing spaces, so blastn runs
    # inside the database directory and gets only the database name
    db_dir, db_name = os.path.split(os.path.abspath(database_path))

    run_blast_program([
        "blastn",
        "-query", os.path.abspath(query_fasta),
        "-db", db_name,
        "-out", os.path.abspath(output_file),
        "-outfmt", "6 " + " ".join(HSP_COLUMNS),
        "-num_threads", str(threads),
        *settings.search_args(),
    ], cwd=db_dir)


def self_scores(query_fasta, tmp_dir, settings):
    """Bitscore of each query aligned to itself (with the same scoring), the maximum a hit can reach."""
    output = os.path.join(tmp_dir, "self.tsv")
    run_blast_program([
        "blastn", "-query", query_fasta, "-subject", query_fasta,
        "-outfmt", "6 qseqid sseqid bitscore", "-out", output,
        *settings.scoring_args(),
    ])

    scores = {}
    with open(output) as f:
        for line in f:
            query, subject, bitscore = line.split("\t")
            if query == subject:
                scores[query] = max(scores.get(query, 0.0), float(bitscore))
    return scores


def find_hits(queries, database_path, output_file, filters=None, blast=None, threads=1,
              progress=print, report_inputs=None):
    """
    Search each query sequence separately in the BLAST database, join HSPs into
    hits (copies of the queries in reads) and save them to a TSV file. All hits
    are saved; the column `passed` says whether a hit passed the filters.
    Reads with hits that did not pass are kept so that they are never used as
    reads without the queries (negative samples).
    queries: [(name, sequence), ...]
    A report is written next to the output (<output>.report.txt).
    Returns (number of hits that passed, number of reads with such a hit).
    """
    started = time.time()
    filters = filters or HitFilters()
    blast = blast or BlastSettings()
    check_queries(queries)

    with tempfile.TemporaryDirectory(prefix="nanoslast_") as tmp_dir:
        query_fasta = os.path.join(tmp_dir, "queries.fasta")
        hsp_file = os.path.join(tmp_dir, "hsps.tsv")
        write_queries(queries, query_fasta)

        progress(f"Running BLAST ({blast.task}) for {len(queries)} query sequence(s)...")
        blast_hsps(query_fasta, database_path, hsp_file, blast, threads)

        hits = hsps_to_hits(read_hsps(hsp_file), self_scores(query_fasta, tmp_dir, blast))

    for hit in hits:
        hit.passed = filters.accepts(hit)
    passed = [hit for hit in hits if hit.passed]
    write_hits(hits, output_file, [name for name, _ in queries])

    progress(f"{len(passed)} of {len(hits)} hits passed the filters "
             f"(coverage >= {filters.min_coverage}, score >= {filters.min_score})")

    passed_per_query = Counter(hit.query_id for hit in passed)
    failed_per_query = Counter(hit.query_id for hit in hits if not hit.passed)
    write_report(output_file, "search", {
        "Inputs": {
            **(report_inputs or {}),
            "BLAST database": os.path.abspath(database_path),
            "query sequences": [name for name, _ in queries],
        },
        "BLAST settings": {**blast.describe(), "threads": threads},
        "Hits": {
            "HSPs joined when at most this far apart (bases)": MAX_HSP_GAP,
            "min. coverage of the query": filters.min_coverage,
            "min. score (relative to the self-alignment)": filters.min_score,
        },
        "Results": {
            "hits found": len(hits),
            "hits passed the filters": len(passed),
            "reads with a passed hit": len({hit.read_id for hit in passed}),
            "passed hits per query": {name: passed_per_query[name] for name, _ in queries},
            "hits below the filters per query": {name: failed_per_query[name] for name, _ in queries},
            "run time (s)": round(time.time() - started, 1),
        },
    })
    return len(passed), len({hit.read_id for hit in passed})


def read_queries(fasta_file) -> list[tuple[str, str]]:
    """Read query sequences from a FASTA file: [(name, sequence), ...]."""
    queries = [(record.id, str(record.seq)) for record in SeqIO.parse(fasta_file, "fasta")]
    if not queries:
        raise ValueError(f"No sequences found in {fasta_file}")
    return queries


def search_reads(sam_file, queries, output_file, db_dir=None, filters=None, blast=None,
                 threads=1, progress=print):
    """
    Whole Search workflow: SAM/BAM → FASTA → BLAST database → blastn → hits.
    If db_dir already contains a database, it is reused; without db_dir the
    database is built in a temporary directory and deleted afterwards.
    Returns (number of hits, number of reads with a hit).
    """
    with tempfile.TemporaryDirectory(prefix="nanoslast_") as tmp_dir:
        db_dir = db_dir or tmp_dir
        db_path = os.path.join(db_dir, "blastdb")

        if glob.glob(db_path + ".n*"):
            progress(f"Using existing BLAST database {db_path}")
        else:
            os.makedirs(db_dir, exist_ok=True)
            fasta_file = os.path.join(db_dir, "reads.fasta")

            progress("Converting SAM to FASTA...")
            convert_sam_to_fasta(sam_file, fasta_file)

            progress("Creating BLAST database...")
            create_blast_db(fasta_file, db_dir)

        return find_hits(queries, db_path, output_file, filters, blast, threads, progress,
                         report_inputs={"reads": os.path.abspath(sam_file)})
