from .sam_reader import iter_sam_records


def convert_sam_to_fasta(input_sam, output_fasta):
    # one primary record per read, in the orientation of the raw signal
    with open(output_fasta, "w") as fasta_file:
        for record in iter_sam_records(input_sam, with_moves=False):
            fasta_file.write(f">{record.read_id}\n{record.sequence}\n")
