import customtkinter as ctk
from tkinter import filedialog, messagebox
import os
from .background import run_in_background
from .plotting import create_signal_figure, load_read, show_figure
from PIL import Image
from customtkinter import CTkImage


def create_plot_part(parent_frame):
    frame = ctk.CTkFrame(parent_frame)
    frame.pack(fill="both", expand=True, padx=20, pady=20)
    frame.grid_columnconfigure((0, 1, 2), weight=1, uniform="column")

    state = {"sam": None, "signals_dir": None}

    # === COLUMN 1 ===
    col1 = ctk.CTkFrame(frame)
    col1.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)
    ctk.CTkLabel(col1, text="SELECT SAM FILE").pack(pady=10)

    sam_label = ctk.CTkLabel(col1, text="")
    ctk.CTkButton(col1, text="Select SAM File", command=lambda: select_file(sam_label, state, "sam")).pack()
    sam_label.pack()

    # === COLUMN 2 ===
    col2 = ctk.CTkFrame(frame)
    col2.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)
    ctk.CTkLabel(col2, text="SELECT SIGNAL DIRECTORY").pack(pady=10)

    signal_label = ctk.CTkLabel(col2, text="")
    ctk.CTkButton(col2, text="Select Directory", command=lambda: select_file(signal_label, state, "signals_dir", directory=True)).pack()
    signal_label.pack()

    # === COLUMN 3 ===
    col3 = ctk.CTkFrame(frame)
    col3.grid(row=0, column=2, sticky="nsew", padx=10, pady=10)
    ctk.CTkLabel(col3, text="PROVIDE A READ ID").pack(pady=10)

    read_id_var = ctk.StringVar()
    ctk.CTkEntry(col3, textvariable=read_id_var, width=250).pack(pady=5)

    ctk.CTkButton(col3, text="PLOT SIGNAL", command=lambda: threaded_plot_signal(frame, state, read_id_var.get())).pack(pady=20)


    # ========== LOGO ==========
    image_path = os.path.join(os.path.dirname(__file__), "assets", "nanoslast_logo.png")
    image = CTkImage(light_image=Image.open(image_path), size=(200, 78))

    image_label = ctk.CTkLabel(parent_frame, image=image, text="")
    image_label.pack(pady=(10, 20))

    return frame


def select_file(label, state, key, directory=False):
    try:
        path = filedialog.askdirectory() if directory else filedialog.askopenfilename()
        if path:
            state[key] = path
            label.configure(text=os.path.basename(path))
    except Exception as e:
        messagebox.showerror("Error", str(e))


def threaded_plot_signal(widget, state, read_id):
    read_id = read_id.strip()

    if not read_id:
        messagebox.showerror("Error", "Please provide a read ID")
        return

    if not state["sam"] or not state["signals_dir"]:
        messagebox.showerror("Error", "Missing SAM file or signal directory")
        return

    def work(progress):
        record, signal_in_pA, sampling_rate = load_read(state["sam"], state["signals_dir"], read_id)
        return create_signal_figure(signal_in_pA, sampling_rate, record, read_id)

    run_in_background(widget, work, lambda fig: show_figure(fig, f"Signal from read {read_id}"))