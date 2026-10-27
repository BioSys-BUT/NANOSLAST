import customtkinter as ctk
from tkinter import filedialog, messagebox
import threading
import os
from PIL import Image
from customtkinter import CTkImage

from .parquet_tool import (
    extract_info_from_tsv,
    sam_file_signal_search,
    load_signals,
    dump_signals_parquet,
    dump_signals_csv,
)

def create_export_part(parent_frame):
    frame = ctk.CTkFrame(parent_frame)
    frame.pack(fill="both", expand=True, padx=20, pady=20)
    frame.grid_columnconfigure((0, 1, 2), weight=1, uniform="column")

    state = {}

    # Top label
    ctk.CTkLabel(frame, text="SIGNAL EXTRACTOR").pack(pady=10)

    # === INPUTS ===
    input_grid = ctk.CTkFrame(frame)
    input_grid.pack(pady=5)

    def create_button(text, command, col):
        btn = ctk.CTkButton(input_grid, text=text, command=command)
        btn.grid(row=0, column=col, padx=5)
        return btn

    state['fast_dir_label'] = ctk.CTkLabel(frame, text="")
    state['sam_label'] = ctk.CTkLabel(frame, text="")
    state['tsv_label'] = ctk.CTkLabel(frame, text="")

    create_button("Select POD5/FAST5 Directory", lambda: select_directory(state['fast_dir_label']), 0)
    create_button("Select SAM File", lambda: select_file(state['sam_label'], "SAM"), 1)
    create_button("Select TSV File", lambda: select_file(state['tsv_label'], "TSV"), 2)

    state['fast_dir_label'].pack()
    state['sam_label'].pack()
    state['tsv_label'].pack()

    # === OUTPUTS ===
    output_grid = ctk.CTkFrame(frame)
    output_grid.pack(pady=10)

    parquet_btn = ctk.CTkButton(output_grid, text="CREATE PARQUET FILE", command=lambda: export_file(state, mode="parquet"))
    parquet_btn.pack(side="left", padx=5)

    csv_btn = ctk.CTkButton(output_grid, text="CREATE CSV FILE", command=lambda: export_file(state, mode="csv"))
    csv_btn.pack(side="left", padx=5)

    state['parquet_btn'] = parquet_btn
    state['csv_btn'] = csv_btn

    image_path = os.path.join(os.path.dirname(__file__), "assets", "nanoslast_logo.png")
    logo_img = CTkImage(light_image=Image.open(image_path), size=(200, 78))
    image_label = ctk.CTkLabel(frame, image=logo_img, text="")
    image_label.pack(pady=(80, 0))

    return frame

# === Helper functions ===

def select_file(label, filetype):
    ext = {"SAM": "*.sam", "TSV": "*.tsv"}.get(filetype, "*")
    file = filedialog.askopenfilename(filetypes=[(f"{filetype} files", ext)], title=f"Select {filetype} file")
    if file:
        label.configure(text=file)

def select_directory(label):
    path = filedialog.askdirectory(title="Select directory")
    if path:
        label.configure(text=path)

def export_file(state, mode="parquet"):
    sam = state['sam_label'].cget("text")
    signals = state['fast_dir_label'].cget("text")
    tsv = state['tsv_label'].cget("text")

    if not sam or not signals or not tsv:
        messagebox.showerror("Error", "Missing required input file(s)")
        return

    ext = ".parquet" if mode == "parquet" else ".csv"
    filetypes = [("Parquet files", "*.parquet")] if mode == "parquet" else [("CSV files", "*.csv")]
    output = filedialog.asksaveasfilename(title="Save output", defaultextension=ext, filetypes=filetypes)

    if not output:
        return

    btn = state['parquet_btn'] if mode == "parquet" else state['csv_btn']
    btn.configure(state="disabled")

    def worker():
        try:
            btn.configure(text="Loading TSV...")
            info = extract_info_from_tsv(tsv)

            btn.configure(text="Searching SAM file...")
            signal_ranges = sam_file_signal_search(sam, info.keys(), info, mode=mode)

            btn.configure(text="Loading signals...")
            signals_data = load_signals(signals, signal_ranges, mode=mode)

            btn.configure(text=f"Saving to {ext.upper()}...")
            if mode == "parquet":
                dump_signals_parquet(signals_data, signal_ranges, output)
            else:
                dump_signals_csv(signals_data, signal_ranges, output)

            messagebox.showinfo("Done", f"{ext.upper()} file created successfully")
        except Exception as e:
            messagebox.showerror("Error", str(e))
        finally:
            btn.configure(text=f"CREATE {ext.upper()} FILE", state="normal")

    threading.Thread(target=worker).start()
