from collections import defaultdict

import numpy as np
import customtkinter as ctk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from .sam_reader import base_boundaries, find_sam_record
from .signal_reader import load_signals

NUCLEOTIDE_TO_COLOR = {"A": "red", "T": "blue", "U": "blue", "C": "yellow", "G": "green"}
DEFAULT_COLOR = "gray"


def load_read(sam_filename, signals_dir, read_id):
    """Return (SamRecord, signal_in_pA, sampling_rate) for one read."""
    record = find_sam_record(sam_filename, read_id)
    if record is None:
        raise ValueError(f"Read {read_id} not found in the SAM file")
    if not record.moves:
        raise ValueError(f"Read {read_id} has no move table (mv tag) in the SAM file")

    signals = load_signals(signals_dir, {read_id})
    if read_id not in signals:
        raise ValueError(f"Signal of read {read_id} not found in the POD5/FAST5 files")

    signal_in_pA, sampling_rate = signals[read_id]
    return record, signal_in_pA, sampling_rate


def create_signal_figure(signal_in_pA, sampling_rate, record, read_id, sample_range=None):
    """
    Plot the raw signal and highlight each basecalled base with its colour.
    `sample_range=(start, stop)` limits the plot to a part of the signal.
    """
    start, stop = sample_range or (0, len(signal_in_pA))
    time_in_ms = np.arange(start, stop) / sampling_rate * 1000

    fig = Figure(figsize=(14, 6))
    ax = fig.add_subplot()
    ax.plot(time_in_ms, signal_in_pA[start:stop], color="black", linewidth=0.8)
    ax.set_title(f"Signal from read {read_id}")
    ax.set_xlabel("Time [ms]")
    ax.set_ylabel("Signal [pA]")

    # one (start, width) span in ms per base, grouped by colour
    spans = defaultdict(list)
    boundaries = base_boundaries(record)

    for nucleotide, base_start, base_stop in zip(record.sequence, boundaries, boundaries[1:]):
        lo = max(base_start, start)
        hi = min(base_stop, stop)
        if lo < hi:
            color = NUCLEOTIDE_TO_COLOR.get(nucleotide, DEFAULT_COLOR)
            spans[color].append((lo / sampling_rate * 1000, (hi - lo) / sampling_rate * 1000))

    for color, xranges in spans.items():
        # y in axes coordinates (0–1), so the spans cover the full height
        ax.broken_barh(
            xranges, (0, 1), facecolors=color, alpha=0.1, linewidth=0,
            transform=ax.get_xaxis_transform(),
        )

    ax.set_xlim(start / sampling_rate * 1000, stop / sampling_rate * 1000)
    ax.autoscale(axis="y")
    return fig


def show_figure(fig, title):
    """Show a figure in a new window. Must be called from the GUI thread."""
    window = ctk.CTkToplevel()
    window.title(title)

    canvas = FigureCanvasTkAgg(fig, master=window)
    NavigationToolbar2Tk(canvas, window).update()
    canvas.get_tk_widget().pack(fill="both", expand=True)
    canvas.draw()

    # CTkToplevel can open behind the main window
    window.after(200, window.lift)
