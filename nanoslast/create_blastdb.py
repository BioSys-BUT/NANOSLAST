from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Application import ApplicationError
import os

def create_blast_db(input_file, output_dir):
    db_path = os.path.join(output_dir, "blastdb")
    makeblastdb_cline = NcbimakeblastdbCommandline(input_file=input_file, dbtype="nucl", out=db_path)

    try:
        makeblastdb_cline()
    except ApplicationError as e:
        raise RuntimeError(f"makeblastdb failed:\n{e.stderr}") from e
