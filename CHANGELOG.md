# Changelog

## 2.0.0 (2026-10-29)

NANOSLAST is now a command line tool. **The graphical interface was removed**; version 1.1.0 is the last one with it.

### Added

- **Command line interface** (`nanoslast search`, `export`, `run`, `plot`). NANOSLAST is now a Python package installable with `pip install git+https://github.com/BioSys-BUT/NANOSLAST.git`. A command without options prints its help with all defaults.
- `--db-dir` option to keep the BLAST database of the reads and reuse it in later searches.
- **All hits of a read are exported**, one row per hit, so a read with several genes or several copies of one gene gives several rows (`--best-only` restores one hit per read).
- **Hits instead of raw BLAST HSPs:** HSPs of one query that lie next to each other in a read are joined into one hit, with its coverage of the query and its score relative to the query's self-alignment.
- **Hit filters** `--min-coverage` and `--min-score` (default 0.8 each). They replace the fixed rule that a hit may be at most 10 nt shorter than the query, which dropped copies with a small gap or cut off by the read end.
- `at_read_end` flag for hits that reach the start or end of the read and may be incomplete.
- **Input checks** before a run (BLAST+ installed, valid query FASTA, SAM with reads and move tables, POD5/FAST5 files present) with errors explaining what to fix.
- The export reports how many hits were skipped and why.
- Parquet export contains the hit positions and scores, not only the signal.
- **HDF5 training dataset** (`-o dataset.h5`) in the format of the NanoResistNet training datasets: one row per read with its whole z-score normalized signal and all its hits, plus reads without any hit as negative samples (`--neg-ratio`, per barcode, reproducible with `--seed`), barcodes, run ID and gene names. Written in batches, so large runs fit in memory.
- **BAM input** (`--bam`, read with samtools).
- **Split reads:** subreads made by Dorado (`pi`, `sp`, `ns` tags) are exported with their part of the parent read's signal; before they were skipped.
- **Text reports** `<output>.report.txt` of `search` (BLAST settings, filters, hits found and passed per query) and of every export (inputs, exported and skipped hits per reason and query, negative samples per barcode, run time).
- **BLAST settings:** `--task` (megablast, dc-megablast, blastn, blastn-short), `--evalue`, and `--blast-options` for any other blastn options (e.g. `-word_size`). Scores of hits are relative to the self-alignment of the query with the same settings.
- POD5/FAST5 files are indexed once and signals are read only from the files that contain the needed reads.
- `--threads` for BLAST.
- Full manual in `docs/manual.md`.

### Changed

- BLAST+ is called directly instead of through `Bio.Blast.Applications`, so Biopython no longer has to be older than 1.86.
- Dependencies are defined in `pyproject.toml`; `requirements.txt` was removed.
- The hits file of `search` has new columns (one row per hit, see the manual), including hits that did not pass the filters (`passed=False`); only passed hits are exported. Hits files from 1.1.0 can still be exported.
- The graphical interface (`python main.py`) was removed, together with the dependencies customtkinter and Pillow. All its functions are in the command line: `search` (Search tab), `export` (Export tab) and `plot` (Plot tab).

### Fixed

- Paths with spaces (e.g. a Windows user name with a space in the temporary directory) no longer break makeblastdb and blastn.

## 1.1.0

### Fixed

- **Signal positions of BLAST hits.** The conversion from base positions to signal samples was shifted: exported segments started one stride (5–6 samples) too late and lost the end of the last matched base. Segments now run exactly from the first sample of the first matched base to the last sample of the last matched base.
- **Base colours in the Plot tab.** Each coloured segment was shifted by about one base, so colours did not match the signal. Thymine (T) was not coloured at all.
- **Reads without a move table** no longer take the move table of the previous read in the SAM file.
- **FASTA files with several query sequences** were joined into one sequence. Each sequence is now searched separately.
- **Hits were silently capped at 1000 reads per query.** `max_target_seqs` is now 1,000,000.
- **Several hits in one read:** the best hit (highest bitscore) is now used instead of the last one in the TSV file.
- The matched sequence in the CSV export was shifted by one base.
- The GUI no longer updates widgets and opens plots from background threads, which Tkinter does not support.

### Added

- Support for **aligned SAM files**: secondary and supplementary alignments are skipped and reads on the reverse strand are reverse complemented back to the orientation of the raw signal.
- Support for **single-read FAST5 files** (in addition to multi-read FAST5 and POD5).
- `query_id` column in the CSV and Parquet exports, with the name of the query sequence found in the read.
- The plot opens in its own window with zoom and pan tools, and long reads are drawn much faster.
- The Export tab reports how many reads were exported.
- The Search tab shows file names instead of truncated paths, and shows the selected query FASTA file.
- Automated test on synthetic data (`tests/test_workflow.py`).
- Version number in the window title.

### Changed

- Source files moved to the `nanoslast/` package. The app is still started with `python main.py`, now from any directory.
- The unfinished *Plot Specific Signal* tab was removed.

## 1.0.0

Version described in the paper (EECSS'25, ICBES 151, doi: [10.11159/icbes25.151](https://doi.org/10.11159/icbes25.151)).
