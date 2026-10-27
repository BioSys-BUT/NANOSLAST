import customtkinter as ctk
from tkinter import filedialog, messagebox
import threading
import os

from .parquet_tool import (
    extract_info_from_tsv,
    sam_file_signal_search,
    load_signals,
)
from .plotting import plot_signal


def create_plot_specific_signal_part(parent_frame):
    frame = ctk.CTkFrame(parent_frame)
    frame.pack(fill="both", expand=True, padx=20, pady=20)
    frame.grid_columnconfigure((0, 1, 2), weight=1, uniform="column")

    state = {'sam_file': None, 'tsv_file': None, 'signal_dir': None, 'read_id': None}

    # === COLUMN 1: SELECT SAM FILE ===
    col1 = ctk.CTkFrame(frame)
    col1.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
    ctk.CTkLabel(col1, text="SELECT SAM FILE").pack(pady=10)
    sam_label = ctk.CTkLabel(col1, text="")
    ctk.CTkButton(col1, text="Select SAM File", command=lambda: select_file(sam_label, state, "sam")).pack()
    sam_label.pack()
    state['sam_label'] = sam_label

    # === COLUMN 2: SELECT SIGNAL DIRECTORY ===
    col2 = ctk.CTkFrame(frame)
    col2.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")
    ctk.CTkLabel(col2, text="SELECT SIGNAL DIRECTORY").pack(pady=10)
    signal_label = ctk.CTkLabel(col2, text="")
    ctk.CTkButton(col2, text="Select Directory", command=lambda: select_directory(signal_label, state)).pack()
    signal_label.pack()
    state['signal_label'] = signal_label

    # === COLUMN 3: SELECT TSV FILE ===
    col3 = ctk.CTkFrame(frame)
    col3.grid(row=0, column=2, padx=10, pady=10, sticky="nsew")
    ctk.CTkLabel(col3, text="SELECT TSV FILE").pack(pady=10)
    tsv_label = ctk.CTkLabel(col3, text="")
    ctk.CTkButton(col3, text="Select TSV File", command=lambda: select_file(tsv_label, state, "tsv")).pack()
    tsv_label.pack()
    state['tsv_label'] = tsv_label

    # === COLUMN 4: PROVIDE READ ID ===
    col4 = ctk.CTkFrame(frame)
    col4.grid(row=1, column=0, padx=10, pady=10, sticky="nsew")
    ctk.CTkLabel(col4, text="PROVIDE READ ID").pack(pady=10)
    read_id_var = ctk.StringVar()
    ctk.CTkEntry(col4, textvariable=read_id_var, width=250).pack(pady=5)
    state['read_id_entry'] = read_id_var

    # === BUTTON TO PLOT ===
    plot_button = ctk.CTkButton(
        col4, text="PLOT SIGNAL", command=lambda: plot_specific_signal(state)
    )
    plot_button.pack(pady=20)

    return frame


def select_file(label, state, filetype):
    ext = {"sam": "*.sam", "tsv": "*.tsv"}.get(filetype, "*")
    file = filedialog.askopenfilename(filetypes=[(f"{filetype.upper()} files", ext)],
                                      title=f"Select {filetype.upper()} file")
    if file:
        label.configure(text=file)
        if filetype == "sam":
            state["sam_file"] = file
        elif filetype == "tsv":
            state["tsv_file"] = file


def select_directory(label, state):
    path = filedialog.askdirectory(title="Select directory")
    if path:
        label.configure(text=path)
        state["signal_dir"] = path


def plot_specific_signal(state):
    sam_file = state['sam_file']
    tsv_file = state['tsv_file']
    signal_dir = state['signal_dir']
    read_id = state['read_id_entry'].get()

    if not all([sam_file, tsv_file, signal_dir, read_id]):
        messagebox.showerror("Error", "Missing required input file(s) or Read ID")
        return

    try:
        # Extract the signal information from the TSV file
        info = extract_info_from_tsv(tsv_file)

        # Find the signal ranges for the given read_id using SAM file
        signal_ranges = {read_id: info.get(read_id)}  # Correct signal_ranges format for SAM search

        # Find signal range for the given read_id (it should be a dictionary)
        signal_ranges = sam_file_signal_search(sam_file, read_id, signal_ranges)  # Fix argument passing

        # Load signals from the specified directory (POD5 or FAST5)
        signals = load_signals(signal_dir, signal_ranges)

        # Get the start and end positions of the signal from the TSV file (BLAST results)
        blast_start, blast_end = info.get(read_id, (None, None))

        if blast_start is None or blast_end is None:
            messagebox.showerror("Error", "Read ID not found in TSV file")
            return

        # Extract the relevant signal portion based on the BLAST result
        signal_in_pA, sampling_rate = signals.get(read_id, (None, None))
        if signal_in_pA is None:
            messagebox.showerror("Error", f"No signal found for {read_id}")
            return

        # Plot the signal for the specified section from BLAST
        plot_signal(signal_in_pA[blast_start:blast_end], signal_ranges, read_id)

        messagebox.showinfo("Done", "Signal plotted successfully")
    except Exception as e:
        messagebox.showerror("Error", str(e))


