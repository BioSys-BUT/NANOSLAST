import customtkinter as ctk
from tkinter import filedialog, messagebox
from tkinter import StringVar, BooleanVar
import os
from PIL import Image
from customtkinter import CTkImage
from Bio import SeqIO

from .background import run_in_background
from .convert_sam_file import convert_sam_to_fasta
from .create_blastdb import create_blast_db
from .blasting import blast_search

def create_search_part(parent_frame):
    main_frame = ctk.CTkFrame(parent_frame)
    main_frame.pack(fill="both", expand=True, padx=20, pady=20)

    main_frame.grid_columnconfigure((0, 1, 2), weight=1, uniform="column")

    # ========== COLUMN 1: SAM → FASTA ==========
    col1 = ctk.CTkFrame(main_frame)
    col1.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)

    ctk.CTkLabel(col1, text="SAM TO FASTA").pack(pady=10)
    sam_path = StringVar(value="")
    sam_path_label = ctk.CTkLabel(col1, text="")
    ctk.CTkButton(col1, text="Select SAM file", command=lambda: select_file(sam_path_label, sam_path, "SAM")).pack()
    sam_path_label.pack()
    ctk.CTkButton(col1, text="Convert SAM file", command=lambda: handle_convert_sam_to_fasta(sam_path_label, sam_path)).pack(pady=5)

    # ========== COLUMN 2: BLAST DB ==========
    col2 = ctk.CTkFrame(main_frame)
    col2.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)

    ctk.CTkLabel(col2, text="BLAST DATABASE").pack(pady=10)
    fasta_path = StringVar(value="")
    fasta_path_label = ctk.CTkLabel(col2, text="")
    ctk.CTkButton(col2, text="Select FASTA file", command=lambda: select_file(fasta_path_label, fasta_path, "FASTA")).pack()
    fasta_path_label.pack()
    ctk.CTkButton(col2, text="Create BLAST database", command=lambda: handle_create_blast_db(fasta_path_label, fasta_path)).pack(pady=5)

    # ========== COLUMN 3: RUN BLAST ==========
    col3 = ctk.CTkFrame(main_frame)
    col3.grid(row=0, column=2, sticky="nsew", padx=10, pady=10)

    ctk.CTkLabel(col3, text="RUN BLAST").pack(pady=10)
    db_dir = StringVar(value="")
    db_dir_label = ctk.CTkLabel(col3, text="")
    ctk.CTkButton(col3, text="Select BLAST database directory", command=lambda: select_directory(db_dir_label, db_dir)).pack()
    db_dir_label.pack()

    # Wrapper for the pasted sequence, hidden when a FASTA file is used
    sequence_frame = ctk.CTkFrame(col3)
    sequence_frame.pack()
    ctk.CTkLabel(sequence_frame, text="Paste sequence").pack()
    sequence_entry = ctk.CTkEntry(sequence_frame, width=200)
    sequence_entry.pack(pady=5)

    run_button = ctk.CTkButton(
        col3,
        text="Run BLAST search",
        command=lambda: run_blast(run_button, sequence_entry, db_dir, use_fasta_var, fasta_sequence_path)
    )
    run_button.pack(pady=5)

    use_fasta_var = BooleanVar(value=False)
    fasta_sequence_path = StringVar(value="")

    def toggle_fasta_input():
        if use_fasta_var.get():
            sequence_frame.pack_forget()
            fasta_button.pack(pady=(5, 0))
            fasta_sequence_label.pack()
        else:
            fasta_button.pack_forget()
            fasta_sequence_label.pack_forget()
            sequence_frame.pack()

    ctk.CTkCheckBox(col3, text="Use FASTA file instead", variable=use_fasta_var, command=toggle_fasta_input).pack(
        pady=(10, 5))

    # Query FASTA file selection, hidden by default
    fasta_sequence_label = ctk.CTkLabel(col3, text="")
    fasta_button = ctk.CTkButton(
        col3,
        text="Select query FASTA file",
        command=lambda: select_sequence_file(fasta_sequence_label, fasta_sequence_path)
    )

    # ========== LOGO ==========
    image_path = os.path.join(os.path.dirname(__file__), "assets", "nanoslast_logo.png")
    image = CTkImage(light_image=Image.open(image_path), size=(200, 78))

    image_label = ctk.CTkLabel(parent_frame, image=image, text="")
    image_label.pack(pady=(10, 20))

    return main_frame

# === Helper functions ===

FASTA_FILETYPES = [("FASTA files", "*.fasta *.fa *.fna *.fas"), ("All files", "*.*")]

# Full paths are kept in StringVars, labels show only the file name

def select_file(label, path_var, filetype):
    types = [("SAM files", "*.sam")] if filetype == "SAM" else FASTA_FILETYPES
    file_path = filedialog.askopenfilename(title=f"Select {filetype} file", filetypes=types)
    if file_path:
        path_var.set(file_path)
        label.configure(text=os.path.basename(file_path))

def select_directory(label, path_var):
    path = filedialog.askdirectory(title="Select directory")
    if path:
        path_var.set(path)
        label.configure(text=os.path.basename(path))

def handle_convert_sam_to_fasta(label, sam_path):
    sam_file = sam_path.get()
    if not sam_file:
        messagebox.showerror("Error", "Please select a SAM file")
        return

    output_file = filedialog.asksaveasfilename(title="Save FASTA file", defaultextension=".fasta",
                                               filetypes=[("FASTA files", "*.fasta")])
    if not output_file:
        return

    run_in_background(
        label,
        lambda progress: convert_sam_to_fasta(sam_file, output_file),
        lambda _: messagebox.showinfo("Done", "FASTA file created successfully"),
    )

def handle_create_blast_db(label, fasta_path):
    fasta_file = fasta_path.get()
    if not fasta_file:
        messagebox.showerror("Error", "Please select a FASTA file")
        return

    output_dir = filedialog.askdirectory(title="Select directory to save BLAST DB")
    if not output_dir:
        return

    run_in_background(
        label,
        lambda progress: create_blast_db(fasta_file, output_dir),
        lambda _: messagebox.showinfo("Done", "BLAST database created"),
    )

def select_sequence_file(label, var):
    file_path = filedialog.askopenfilename(
        title="Select query FASTA file",
        filetypes=FASTA_FILETYPES[:1] + [("Text files", "*.txt"), ("All files", "*.*")]
    )
    if file_path:
        var.set(file_path)
        label.configure(text=os.path.basename(file_path))

def run_blast(widget, entry, db_dir_var, use_fasta_var, fasta_sequence_path):
    db_dir = db_dir_var.get()

    if not db_dir:
        messagebox.showerror("Error", "Please select a BLAST database directory")
        return

    if use_fasta_var.get():
        fasta_file = fasta_sequence_path.get()
        if not fasta_file:
            messagebox.showerror("Error", "Please select a FASTA file with a sequence")
            return
        try:
            # each sequence in the file is searched separately
            queries = [(record.id, str(record.seq)) for record in SeqIO.parse(fasta_file, "fasta")]
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read FASTA file:\n{e}")
            return
        if not queries:
            messagebox.showerror("Error", "No sequences found in the FASTA file")
            return
    else:
        sequence = "".join(entry.get().split())
        if not sequence:
            messagebox.showerror("Error", "Please enter a sequence")
            return
        queries = [("Query", sequence)]

    output_file = filedialog.asksaveasfilename(
        title="Save BLAST output",
        defaultextension=".tsv",
        filetypes=[("TSV files", "*.tsv")]
    )
    if not output_file:
        return

    db_path = os.path.normpath(os.path.join(db_dir, "blastdb"))

    def on_success(n_hits):
        if n_hits == 0:
            messagebox.showinfo("No Hits", "No hits were found.")
        else:
            messagebox.showinfo("Done", f"BLAST search completed: {n_hits} hits")

    run_in_background(
        widget,
        lambda progress: blast_search(queries, db_path, output_file),
        on_success,
        on_error=lambda e: messagebox.showerror("Error", f"BLAST search failed:\n{e}"),
    )

