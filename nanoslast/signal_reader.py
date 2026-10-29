import os
from collections import defaultdict

import h5py
import numpy as np
import pod5 as p5


def find_signal_files(dir_path) -> list[str]:
    """Recursively collect .pod5 and .fast5 files from one or more directories (or files)."""
    queue = [dir_path] if isinstance(dir_path, str) else list(dir_path)
    files = []

    while queue:
        path = queue.pop(0)

        if os.path.isdir(path):
            queue.extend(os.path.join(path, name) for name in sorted(os.listdir(path)))
        elif path.endswith((".pod5", ".fast5")):
            files.append(path)

    return files


class SignalIndex:
    """
    Which POD5/FAST5 file contains which read. Built once from the read IDs stored
    in the files (signals are not read), so that signals can then be loaded only
    from the files that contain them.
    """

    def __init__(self, dir_path):
        self.file_of = {}
        for file_path in find_signal_files(dir_path):
            for read_id in _read_ids(file_path):
                self.file_of[read_id] = file_path

    def __contains__(self, read_id):
        return read_id in self.file_of

    def __len__(self):
        return len(self.file_of)

    def load(self, read_ids) -> dict:
        """{read_id: (signal_in_pA, sampling_rate)} for the reads that are in the index."""
        by_file = defaultdict(set)
        for read_id in read_ids:
            if read_id in self.file_of:
                by_file[self.file_of[read_id]].add(read_id)

        result = {}
        for file_path, ids in by_file.items():
            reader = read_pod5 if file_path.endswith(".pod5") else read_fast5
            for read_id, signal, sampling_rate in reader(file_path, ids):
                result[read_id] = (signal, sampling_rate)
        return result


def _read_ids(file_path):
    if file_path.endswith(".pod5"):
        with p5.Reader(file_path) as reader:
            return [str(read_id) for read_id in reader.read_ids]

    with h5py.File(file_path, "r") as f:
        if "Raw/Reads" in f:
            return [_decode(group.attrs["read_id"]) for group in f["Raw/Reads"].values()]
        return [_decode(group["Raw"].attrs["read_id"]) for group in f.values() if "Raw" in group]


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
    return SignalIndex(dir_path).load(read_ids)


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
