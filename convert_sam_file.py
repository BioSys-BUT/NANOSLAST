def convert_sam_to_fasta(input_sam, output_fasta):

    with open(input_sam, "r") as sam_file, open(output_fasta, "w") as fasta_file:
        for line in sam_file:
            if line.startswith("@"):
                continue

            fields = line.strip().split("\t")
            read_id = fields[0]
            sequence = fields[9]

            fasta_file.write(f">{read_id}\n{sequence}\n")
