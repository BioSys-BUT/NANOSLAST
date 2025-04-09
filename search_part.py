import customtkinter as ctk
from tkinter import filedialog, messagebox
from tkinter import StringVar, BooleanVar
import threading
import os
from PIL import Image
from customtkinter import CTkImage

from convert_sam_file import convert_sam_to_fasta
from create_blastdb import create_blast_db
from blasting import blast_search

def create_search_part(parent_frame):
    main_frame = ctk.CTkFrame(parent_frame)
    main_frame.pack(fill="both", expand=True, padx=20, pady=20)

    main_frame.grid_columnconfigure((0, 1, 2), weight=1, uniform="column")

    # ========== COLUMN 1: SAM → FASTA ==========
    col1 = ctk.CTkFrame(main_frame)
    col1.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)

    ctk.CTkLabel(col1, text="SAM TO FASTA").pack(pady=10)
    sam_path_label = ctk.CTkLabel(col1, text="")
    ctk.CTkButton(col1, text="Select SAM file", command=lambda: select_file(sam_path_label, "SAM")).pack()
    sam_path_label.pack()
    ctk.CTkButton(col1, text="Convert SAM file", command=lambda: handle_convert_sam_to_fasta(sam_path_label)).pack(pady=5)

    # ========== COLUMN 2: BLAST DB ==========
    col2 = ctk.CTkFrame(main_frame)
    col2.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)

    ctk.CTkLabel(col2, text="BLAST DATABASE").pack(pady=10)
    fasta_path_label = ctk.CTkLabel(col2, text="")
    ctk.CTkButton(col2, text="Select FASTA file", command=lambda: select_file(fasta_path_label, "FASTA")).pack()
    fasta_path_label.pack()
    ctk.CTkButton(col2, text="Create BLAST database", command=lambda: handle_create_blast_db(fasta_path_label)).pack(pady=5)

    # ========== COLUMN 3: RUN BLAST ==========
    col3 = ctk.CTkFrame(main_frame)
    col3.grid(row=0, column=2, sticky="nsew", padx=10, pady=10)

    ctk.CTkLabel(col3, text="RUN BLAST").pack(pady=10)
    db_dir_label = ctk.CTkLabel(col3, text="")
    ctk.CTkButton(col3, text="Select BLAST database directory", command=lambda: select_directory(db_dir_label)).pack()
    # db_dir_label.pack(pady=5)
    # ctk.CTkLabel(col3, text="Paste sequence").pack(pady=(10, 0))
    db_dir_label.pack()  # <–– bez paddingu

    # ctk.CTkLabel(col3, text="Paste sequence").pack()
    # sequence_entry = ctk.CTkEntry(col3, width=200)
    # sequence_entry.pack(pady=5)
    # # ctk.CTkButton(col3, text="Run BLAST search", command=lambda: run_blast(sequence_entry, db_dir_label)).pack(pady=5)
    # ctk.CTkButton(
    #     col3,
    #     text="Run BLAST search",
    #     command=lambda: run_blast(sequence_entry, db_dir_label, use_fasta_var, fasta_sequence_path)
    # ).pack(pady=5)

    # Obal pro paste-sequence
    sequence_frame = ctk.CTkFrame(col3)
    sequence_frame.pack()
    ctk.CTkLabel(sequence_frame, text="Paste sequence").pack()
    sequence_entry = ctk.CTkEntry(sequence_frame, width=200)
    sequence_entry.pack(pady=5)

    # Tlačítko Run BLAST search
    run_button = ctk.CTkButton(
        col3,
        text="Run BLAST search",
        command=lambda: run_blast(sequence_entry, db_dir_label, use_fasta_var, fasta_sequence_path)
    )
    run_button.pack(pady=5)

    use_fasta_var = BooleanVar(value=False)
    fasta_sequence_path = StringVar(value="")  # cesta k souboru

    def toggle_fasta_input():
        if use_fasta_var.get():
            sequence_frame.pack_forget()
            fasta_button.pack(pady=5)
        else:
            fasta_button.pack_forget()
            sequence_frame.pack()

    ctk.CTkCheckBox(col3, text="Use FASTA file instead", variable=use_fasta_var, command=toggle_fasta_input).pack(
        pady=(10, 5))

    # Button pro výběr FASTA souboru – standardně skrytý
    fasta_button = ctk.CTkButton(
        col3,
        text="Select FASTA file",
        command=lambda: select_sequence_file(fasta_sequence_path)
    )

    # ========== LOGO ==========
    image_path = "nanoslast_logo.png"  # nebo tvá reálná cesta k souboru
    image = CTkImage(light_image=Image.open(image_path), size=(200, 70))

    image_label = ctk.CTkLabel(parent_frame, image=image, text="")
    image_label.pack(pady=(10, 20))  # odsazení od layoutu výš

    return main_frame

# === Podpůrné funkce ===

def select_file(label, filetype):
    types = [("SAM files", "*.sam")] if filetype == "SAM" else [("FASTA files", "*.fasta")]
    file_path = filedialog.askopenfilename(title=f"Select {filetype} file", filetypes=types)
    if file_path:
        label.configure(text=file_path)

def select_directory(label):
    path = filedialog.askdirectory(title="Select directory")
    if path:
        label.configure(text=path)

def handle_convert_sam_to_fasta(label):
    sam_file = label.cget("text")
    if not sam_file:
        messagebox.showerror("Error", "Please select a SAM file")
        return

    output_file = filedialog.asksaveasfilename(title="Save FASTA file", defaultextension=".fasta",
                                               filetypes=[("FASTA files", "*.fasta")])
    if not output_file:
        return

    def worker():
        convert_sam_to_fasta(sam_file, output_file)
        messagebox.showinfo("Done", "FASTA file created successfully")

    threading.Thread(target=worker).start()

def handle_create_blast_db(label):
    fasta_file = label.cget("text")
    if not fasta_file:
        messagebox.showerror("Error", "Please select a FASTA file")
        return

    output_dir = filedialog.askdirectory(title="Select directory to save BLAST DB")
    if not output_dir:
        return

    def worker():
        os.makedirs(output_dir, exist_ok=True)
        create_blast_db(fasta_file, output_dir)
        messagebox.showinfo("Done", "BLAST database created")

    threading.Thread(target=worker).start()

def select_sequence_file(var):
    file_path = filedialog.askopenfilename(
        title="Select FASTA file",
        filetypes=[("FASTA files", "*.fasta"), ("Text files", "*.txt"), ("All files", "*.*")]
    )
    if file_path:
        var.set(file_path)

def run_blast(entry, db_label, use_fasta_var, fasta_sequence_path):
    db_dir = db_label.cget("text")

    if not db_dir:
        messagebox.showerror("Error", "Please select a BLAST database directory")
        return

    if use_fasta_var.get():
        fasta_file = fasta_sequence_path.get()
        if not fasta_file:
            messagebox.showerror("Error", "Please select a FASTA file with a sequence")
            return
        try:
            with open(fasta_file, "r") as f:
                sequence = "".join(line.strip() for line in f if not line.startswith(">"))
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read FASTA file:\n{e}")
            return
    else:
        sequence = entry.get()
        if not sequence:
            messagebox.showerror("Error", "Please enter a sequence")
            return

    output_file = filedialog.asksaveasfilename(
        title="Save BLAST output",
        defaultextension=".tsv",
        filetypes=[("TSV files", "*.tsv")]
    )
    if not output_file:
        return

    db_path = os.path.normpath(os.path.join(db_dir, "blastdb"))

    def worker():
        try:
            blast_search(sequence, db_path, output_file)
            messagebox.showinfo("Done", "BLAST search completed")
        except Exception as e:
            messagebox.showerror("Error", f"BLAST search failed:\n{e}")

    threading.Thread(target=worker).start()

