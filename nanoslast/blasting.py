import tempfile
import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Application import ApplicationError
import os

# Maximum number of reads reported per query. BLAST's default (500) and the
# previous 1000 silently dropped hits for sequences present in many reads.
MAX_TARGET_SEQS = 1_000_000

def blast_search(queries: list[tuple[str, str]], database_path, output_file) -> int:
    """
    Search each query sequence separately in the BLAST database and save
    all hits to one TSV file (qseqid = query name).
    queries: [(name, sequence), ...]
    Returns the number of hits.
    """

    # Write the query sequences to a temporary FASTA file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as query_file:
        for name, sequence in queries:
            query_file.write(f">{name}\n{sequence}\n")

    database = os.path.normpath(database_path)

    # Set up the BLAST command
    blastn_cline = NcbiblastnCommandline(
        query=query_file.name,
        db=database,
        out=output_file,
        outfmt=6,
        max_target_seqs=MAX_TARGET_SEQS,
        evalue=0.001,
    )

    # Execute the BLAST search
    try:
        blastn_cline()
    except ApplicationError as e:
        raise RuntimeError(f"blastn failed:\n{e.stderr}") from e
    finally:
        os.remove(query_file.name)

    if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
        return 0

    blast_table = pd.read_csv(output_file, sep="\t", header=None, dtype={0: str, 1: str})
    blast_table.columns = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ]

    query_lengths = blast_table["qseqid"].map({name: len(sequence) for name, sequence in queries})

    # Filter out hits that are shorter than their query sequence
    filtered_blast_table = blast_table[blast_table["length"] >= query_lengths - 10]

    filtered_blast_table.to_csv(output_file, sep="\t", index=False)

    return len(filtered_blast_table)
