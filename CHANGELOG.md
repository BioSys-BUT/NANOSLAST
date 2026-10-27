# Changelog

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
