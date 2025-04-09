import numpy as np
import matplotlib.pyplot as plt
import os
import h5py
import pod5 as p5

def sam_file_signal_search(
    sam_filename, read_id: str
):
    def peek_line(f):
        pos = f.tell()
        line = f.readline()
        f.seek(pos)
        return line

    with open(sam_filename, "r") as sam:
        while peek_line(sam).startswith("@"):
            sam.readline()

        for line in sam:
            if not line.startswith(read_id):
                continue

            line = line.split("\t")
            sequence = line[9]

            for field in line[11:]:
                if field.startswith("mv:B:c"):
                    move_table = field[9:]
                    stride = int(field[7])

                elif field.startswith("ts:i:"):
                    trimmed_samples = int(field[5:])

            if not move_table:
                continue

            moves = list(map(int, move_table.split(",")))

            return sequence, moves, trimmed_samples, stride

    return None


def load_signal_fast5(path: str, read_id: str):
    with h5py.File(path, "r") as file:
        for read in file:
            if file[f"{read}/Raw"].attrs["read_id"].decode("utf-8") != read_id:
                continue

            raw_signal = list(file[f"{read}/Raw/Signal"])
            channel_id = file[f"{read}/channel_id"]

            digitisation = channel_id.attrs["digitisation"]
            offset = channel_id.attrs["offset"]
            range = channel_id.attrs["range"]
            sampling_rate = channel_id.attrs["sampling_rate"]

            signal_in_pA = (raw_signal + offset) * range / digitisation

            return signal_in_pA, sampling_rate

    return None


def load_signal_pod5(path: str, read_id: str):
    with p5.Reader(path) as reader:
        for read in reader.reads([read_id], missing_ok=True):
            signal_in_pA = read.signal_pa
            sampling_rate = read.run_info.sample_rate
            return signal_in_pA, sampling_rate

    return None


def load_signal(dir_path, read_id):
    if isinstance(dir_path, str):
        dir_paths = [dir_path]
    if isinstance(dir_path, list):
        dir_paths = dir_path

    queue = dir_paths.copy()
    files: list[str] = []

    while queue:
        path = queue.pop(0)

        if os.path.isdir(path):
            queue.extend([os.path.join(path, f) for f in os.listdir(path)])
        else:
            files.append(path)

    file_count = len(files)
    for i, path in enumerate(files):
        # progress_callback(1 / file_count)
        if path.endswith(".fast5"):
            signal_tuple = load_signal_fast5(path, read_id)
            if signal_tuple:
                return signal_tuple

        if path.endswith(".pod5"):
            signal_tuple = load_signal_pod5(path, read_id)
            if signal_tuple:
                return signal_tuple

    return None


def plot_signal(
    signal_tuple,
    sam_search_result: tuple[str, list[int], int, int],
    read_id: str,
):

    nucleotide_sequence, moves, trimmed_samples, stride = sam_search_result

    signal_in_pA, sampling_rate = signal_tuple

    # Iterating over the move array to calculate signal_points
    zero_count = 0
    signal_points = []

    for move in moves:
        if move == 1:
            signal_point = (zero_count + 1) * stride
            signal_points.append(signal_point)
            zero_count = 0
        else:
            zero_count += 1

    # Calculate cumulative sum of signal_points
    cumulative_points = np.cumsum(signal_points)

    # Calculate segment boundaries based on cumulative_points starting from trimmed_samples
    segment_boundaries = [trimmed_samples]

    for point in cumulative_points:
        segment_boundaries.append(trimmed_samples + point)

    nucleotide_to_color = {"A": "red", "T": "blue", "C": "yellow", "G": "green"}
    default_color = "gray"

    time_in_ms = np.arange(len(signal_in_pA)) / sampling_rate * 1000

    plt.figure(figsize=(14, 6))
    plt.plot(time_in_ms, signal_in_pA, color="black")
    plt.title(f"Signal from read {read_id}")
    plt.xlabel("Time [ms]")
    plt.ylabel("Signal [pA]")

    segment_boundaries_in_ms = np.array(segment_boundaries) / sampling_rate * 1000

    for i, nucleotide in enumerate(nucleotide_sequence):
        start = segment_boundaries_in_ms[i]
        end = segment_boundaries_in_ms[i + 1]
        color = nucleotide_to_color.get(nucleotide, default_color)
        plt.axvspan(start, end, color=color, alpha=0.1)

    plt.show()
