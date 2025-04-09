import h5py
import os
import pod5 as p5
import pyarrow as pa
import pyarrow.parquet as pq
import csv

from typing import Literal
from tkinter import messagebox



def main():

    signals = get_signals_from_tsv(tsv_file, sam_file_pod5, [path_dir_pod5])

    # dump signals to csv file
    dump_signals_parquet(signals, parquet_result)
    dump_signals_csv(signals, csv_result)

    print(
        "found",
        sum(1 for signal in signals.values() if signal is not None),
        "/",
        len(signals),
    )


def get_signals_from_tsv(tsv_file, sam_file, signals_dir):
    # extract searching data from tsv file
    info = extract_info_from_tsv(tsv_file)

    # search for signal ranges in sam file
    signal_ranges = sam_file_signal_search(sam_file, info.keys(), info)

    # load signals from pod5 and fast5 files for each read ID and signal range
    signals = load_signals(signals_dir, signal_ranges)

    return signals


def dump_signals_csv(signals, signal_ranges, filename):
    with open(filename, "w", newline="") as f:
        try:
            writer = csv.writer(f)
        except Exception as e:
            messagebox.showerror("Error", e)
            return

        writer.writerow(
            [
                "read_id",
                "sequence",
                "start_idx",
                "stop_idx",
                "signal",
                "start_signal",
                "stop_signal",
                "direction",
            ]
        )
        for read_id, result in signals.items():
            if result is None:
                continue

            (
                start_signal,
                stop_signal,
                sequence,
                start_idx,
                stop_idx,
                direction,
                signal,
            ) = result

            writer.writerow(
                [
                    read_id,
                    sequence[start_idx : stop_idx + 1],
                    start_idx,
                    stop_idx,
                    ";".join(map(str, signal)),
                    start_signal,
                    stop_signal,
                    direction,
                ]
            )


def dump_signals_parquet(signals, _, filename):
    signals = dict(filter(lambda x: x[1] is not None, signals.items()))
    table = pa.table(
        {
            "read_id": signals.keys(),
            "signal": signals.values(),
        }
    )
    pq.write_table(table, filename)


def extract_info_from_tsv(filename):
    with open(filename, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")

        return {row["sseqid"]: (int(row["sstart"]), int(row["send"])) for row in reader}


def signal_searcher(list, start_idx, stop_idx, trimmed_samples, stride=5):
    count_all = 0
    count_ones = 0
    start_signal = None
    stop_signal = None

    for x in list:
        count_all += 1

        if x != 1:
            continue

        count_ones += 1

        if count_ones == start_idx:
            start_signal = count_all * stride + trimmed_samples

        if count_ones == stop_idx:
            stop_signal = count_all * stride + trimmed_samples
            break

    return start_signal, stop_signal


def sam_file_signal_search(
    sam_filename,
    read_ids: set[str] | str,
    signal_ranges: dict[str, tuple[int, int]],
    mode: Literal["parquet", "csv"] = "parquet",
):

    def parquet_result():
        return start_signal, stop_signal

    def csv_result():
        return start_signal, stop_signal, sequence, start_idx, stop_idx, direction

    result_func = {"parquet": parquet_result, "csv": csv_result}[mode]

    def peek_line(f):
        pos = f.tell()
        line = f.readline()
        f.seek(pos)
        return line

    # iterator for converting move table string to integers
    def move_table_parser(move_table: str):
        for char in move_table:
            if char.isdigit():
                yield int(char)

    # if read_ids is a string, convert it to a set
    # so that we can use the same logic for both cases (one or more read ids)
    if isinstance(read_ids, str):
        read_ids: set[str] = {read_ids}

    result: dict[str, any] = {id: None for id in read_ids}

    with open(sam_filename, "r") as sam:
        while peek_line(sam).startswith("@"):
            sam.readline()

        found_ids = 0

        import time

        start = time.time()

        for line in sam:
            line = line.split("\t")
            qname = line[0]
            seq_len = int(line[8])
            sequence = line[9]

            if qname not in read_ids:
                continue

            start_idx, stop_idx = signal_ranges[qname]

            direction = 1

            for field in line[11:]:
                if field.startswith("mv:B:c"):
                    move_table = field[9:]
                    stride = int(field[7])

                elif field.startswith("ts:i:"):
                    trimmed_samples = int(field[5:])

            if start_idx > stop_idx:
                start_idx, stop_idx = stop_idx, start_idx
                direction = -1

            if not move_table:
                continue

            moves = move_table_parser(move_table)

            signal_of_interest = signal_searcher(
                moves, start_idx, stop_idx, trimmed_samples, stride
            )

            start_signal, stop_signal = signal_of_interest

            result[qname] = result_func()

            found_ids += 1

            if all(result.values()):
                break

        print("Elapsed time:", time.time() - start)

    return result


def load_signals(
    dir_path: str | list[str],
    signal_ranges: dict[str, tuple[int, int]],
    mode: Literal["parquet", "csv"] = "parquet",
):
    result = {id: None for id in signal_ranges.keys()}

    if isinstance(dir_path, str):
        dir_paths = [dir_path]
    if isinstance(dir_path, list):
        dir_paths = dir_path

    queue = dir_paths.copy()
    while queue:
        current = queue.pop()

        if os.path.isdir(current):
            queue.extend(os.path.join(current, name) for name in os.listdir(current))
            continue

        file_extension = os.path.splitext(current)[-1]

        if file_extension == ".fast5":
            signals = load_signals_fast5(current, signal_ranges, mode=mode)
        elif file_extension == ".pod5":
            signals = load_signals_pod5(current, signal_ranges, mode=mode)
        else:
            continue

        if signals:
            result.update(signals)

    return result


def load_signals_pod5(
    file_path,
    signal_ranges: dict[str, tuple[int, int]],
    mode: Literal["parquet", "csv"] = "parquet",
):
    def parquet_result():
        return signal_in_pA[start:stop]

    def csv_result():
        return start, stop, *rest, signal_in_pA

    result_func = {"parquet": parquet_result, "csv": csv_result}[mode]
    result = {}

    with p5.Reader(file_path) as reader:

        for read in reader.reads(signal_ranges.keys(), missing_ok=True):

            signal_in_pA = read.signal_pa
            read_id = str(read.read_id)

            start, stop, *rest = signal_ranges[read_id]
            result[read_id] = result_func()

    return result


def load_signals_fast5(
    file_path,
    signal_ranges: dict[str, tuple[int, int]],
    mode: Literal["parquet", "csv"] = "parquet",
):
    def parquet_result():
        return signal_in_pA[start:stop]

    def csv_result():
        return start, stop, *rest, signal_in_pA

    result_func = {"parquet": parquet_result, "csv": csv_result}[mode]
    result = {}

    with h5py.File(file_path, "r") as f:

        for read in f:
            read_id = f[f"{read}/Raw"].attrs["read_id"].decode("utf-8")

            if read_id in signal_ranges:

                raw_signal = list(f[f"{read}/Raw/Signal"])
                channel_id = f[f"{read}/channel_id"]

                digitisation = channel_id.attrs["digitisation"]
                signal_offset = channel_id.attrs["offset"]
                signal_range = channel_id.attrs["range"]

                signal_in_pA = (
                    (raw_signal + signal_offset) * signal_range / digitisation
                )

                start, stop, *rest = signal_ranges[read_id]

                result[read_id] = result_func()

    return result


if __name__ == "__main__":
    main()
