"""Checks of the inputs before a run, with errors that say how to fix the problem."""
import os
import shutil

from .sam_reader import iter_sam_records
from .signal_reader import find_signal_files

SAM_RECORDS_TO_CHECK = 1000
NUCLEOTIDES = set("ACGTUNRYSWKMBDHVacgtunryswkmbdhv")


class InputError(ValueError):
    pass


def check_blast_installed():
    missing = [p for p in ("makeblastdb", "blastn") if shutil.which(p) is None]
    if missing:
        raise InputError(
            f"{' and '.join(missing)} not found. Install NCBI BLAST+ and make sure it is on PATH."
        )


def check_queries(queries):
    if not queries:
        raise InputError("No query sequences given")

    names = [name for name, _ in queries]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        raise InputError(f"Duplicate query names: {', '.join(duplicates)}")

    for name, sequence in queries:
        if not sequence:
            raise InputError(f"Query {name} has an empty sequence")
        invalid = set(sequence) - NUCLEOTIDES
        if invalid:
            raise InputError(f"Query {name} contains characters that are not nucleotides: {''.join(sorted(invalid))}")


def check_sam(sam_file, require_moves=True):
    """Check that the SAM file has reads and, if required, move tables (from the first reads)."""
    if not os.path.isfile(sam_file):
        raise InputError(f"SAM file not found: {sam_file}")

    n_records = n_with_moves = 0
    for record in iter_sam_records(sam_file, with_moves=require_moves):
        n_records += 1
        n_with_moves += bool(record.moves)
        if n_records >= SAM_RECORDS_TO_CHECK:
            break

    if n_records == 0:
        raise InputError(f"No reads found in {sam_file}")

    if require_moves and n_with_moves == 0:
        raise InputError(
            f"The reads in {sam_file} have no move tables (mv tag), so their bases cannot be "
            "linked to the raw signal. Basecall with dorado using --emit-moves "
            "(and --emit-sam, or convert the BAM with samtools view -h)."
        )


def check_signals_dir(path):
    if not os.path.exists(path):
        raise InputError(f"Signal directory not found: {path}")
    if not find_signal_files(path):
        raise InputError(f"No POD5 or FAST5 files found in {path}")
