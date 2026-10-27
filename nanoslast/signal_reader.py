import os
import h5py
import numpy as np
import pod5 as p5


def find_signal_files(dir_path) -> list[str]:
    """Recursively collect .pod5 and .fast5 files from one or more directories."""
    queue = [dir_path] if isinstance(dir_path, str) else list(dir_path)
    files = []

    while queue:
        path = queue.pop(0)

        if os.path.isdir(path):
            queue.extend(os.path.join(path, name) for name in sorted(os.listdir(path)))
        elif path.endswith((".pod5", ".fast5")):
            files.append(path)

    return files


def read_pod5(file_path, read_ids):
    """Yield (read_id, signal_in_pA, sampling_rate) for the requested reads."""
    with p5.Reader(file_path) as reader:
        for read in reader.reads(selection=list(read_ids), missing_ok=True):
            yield str(read.read_id), read.signal_pa, read.run_info.sample_rate


def read_fast5(file_path, read_ids):
    """
    Yield (read_id, signal_in_pA, sampling_rate) for the requested reads.
    Supports both multi-read and single-read FAST5 files.
    """
    with h5py.File(file_path, "r") as f:
        if "Raw/Reads" in f:
            # single-read FAST5: /Raw/Reads/Read_N, shared /UniqueGlobalKey/channel_id
            channel_id = f["UniqueGlobalKey/channel_id"]
            for read_group in f["Raw/Reads"].values():
                read_id = _decode(read_group.attrs["read_id"])
                if read_id in read_ids:
                    yield read_id, *_signal_in_pA(read_group["Signal"], channel_id)
            return

        # multi-read FAST5: /read_<id>/Raw, /read_<id>/channel_id
        for read_group in f.values():
            if "Raw" not in read_group:
                continue

            read_id = _decode(read_group["Raw"].attrs["read_id"])
            if read_id in read_ids:
                yield read_id, *_signal_in_pA(
                    read_group["Raw/Signal"], read_group["channel_id"]
                )


def load_signals(dir_path, read_ids) -> dict:
    """
    Load signals of the given reads from all POD5/FAST5 files in a directory.
    Returns {read_id: (signal_in_pA, sampling_rate)} for the reads that were found.
    """
    remaining = set(read_ids)
    result = {}

    for file_path in find_signal_files(dir_path):
        if not remaining:
            break

        reader = read_pod5 if file_path.endswith(".pod5") else read_fast5

        for read_id, signal, sampling_rate in reader(file_path, remaining):
            result[read_id] = (signal, sampling_rate)

        remaining -= result.keys()

    return result


def _decode(value) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def _signal_in_pA(raw_signal, channel_id):
    digitisation = channel_id.attrs["digitisation"]
    offset = channel_id.attrs["offset"]
    signal_range = channel_id.attrs["range"]
    sampling_rate = channel_id.attrs["sampling_rate"]

    raw = np.asarray(raw_signal[()], dtype=np.float32)
    signal_in_pA = (raw + offset) * signal_range / digitisation

    # float32 like POD5's signal_pa
    return signal_in_pA.astype(np.float32), float(sampling_rate)
