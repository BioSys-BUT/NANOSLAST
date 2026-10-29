import os
import subprocess


def run_blast_program(command, cwd=None):
    """Run a BLAST+ program, turning its failures into readable errors."""
    try:
        subprocess.run(command, check=True, capture_output=True, text=True, cwd=cwd)
    except FileNotFoundError as e:
        raise RuntimeError(
            f"{command[0]} not found. Install NCBI BLAST+ and make sure it is on PATH."
        ) from e
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"{command[0]} failed:\n{e.stderr.strip()}") from e


def create_blast_db(input_file, output_dir):
    # BLAST+ cannot use a database path containing spaces (common on Windows),
    # so makeblastdb runs inside output_dir and gets only the database name.
    # -in accepts a path with spaces when it is quoted inside the value.
    run_blast_program(
        ["makeblastdb", "-in", f'"{os.path.abspath(input_file)}"', "-dbtype", "nucl", "-out", "blastdb"],
        cwd=output_dir,
    )
