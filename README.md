<p align="center">
  <img src="docs/nanoslast_logo.png" alt="NANOSLAST – Nanopore Signal Local Alignment and Segmentation Tool" width="500">
</p>

# NANOSLAST

NANOSLAST is a Python tool that links local sequence similarity search (BLAST) to raw Oxford Nanopore signals (squiggles). It finds sequences of interest (genes, motifs) in basecalled reads, locates the exact part of the raw signal that corresponds to each copy, and exports these signal segments or a complete HDF5 training dataset for signal-based neural networks. It is a command line tool, runs on servers and clusters without a display, and can also be used from Python.

- Raw signals in **POD5** and **FAST5** (single- and multi-read), basecalled reads in **SAM** or **BAM** (Dorado or Guppy, with move tables), aligned or unaligned, including split reads.
- All copies of all query sequences in each read, with filters on coverage and score.
- Exports to **CSV**, **Parquet** and **HDF5** (training dataset with negative samples), with text reports of settings and results.
- Adjustable BLAST settings (task, e-value, any other `blastn` options).
- Signal positions validated by basecalling the exported segments (99.3 % median identity, boundaries within 1–3 bases).

Typical uses: preparing training data for signal-based gene detection, basecalling or methylation models; extracting signals of specific motifs; comparing signal characteristics between flowcells or sequencers.

**[Full manual](docs/manual.md)**: data preparation, all commands and options, output formats, building a training dataset, troubleshooting.

## What's new in 2.0

- NANOSLAST is a command line tool installable with `pip`; the graphical interface was removed (the last version with it is 1.1.0).
- All copies of all queries in a read are found and exported, instead of one BLAST HSP per read; hits are filtered by coverage and score.
- HDF5 training dataset with negative samples in the format of the NanoResistNet datasets.
- BAM input, split reads, input checks and text reports of every run.

Details in [CHANGELOG.md](CHANGELOG.md).

## Installation

Requires Python 3.10+, [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) on `PATH` and [samtools](https://www.htslib.org/) for BAM input.

```bash
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install git+https://github.com/BioSys-BUT/NANOSLAST.git          # a specific version: ...NANOSLAST.git@v2.0.0

conda install -c bioconda blast samtools   # or: sudo apt install ncbi-blast+ samtools
```

On Windows, install BLAST+ with the installer from NCBI, add its `bin` directory to `PATH` and use SAM files.

## Input data

Reads must be basecalled with **move tables**, which link bases to the raw signal:

```bash
dorado basecaller <model> pod5_dir/ --emit-moves [--kit-name ...] > reads.bam
```

Existing BAM files without move tables cannot be used; NANOSLAST checks this and stops with an explanation. See [Preparing the input data](docs/manual.md#2-preparing-the-input-data).

## Usage

```bash
# signal segments of all copies of the genes from genes.fasta
nanoslast run --bam reads.bam --signals pod5_dir/ --query genes.fasta -o segments.parquet

# HDF5 training dataset: whole reads, all genes, reads without genes as negative samples
nanoslast run --bam reads.bam --signals pod5_dir/ --query genes.fasta --run-id my_run -o dataset.h5

# the same in two steps, keeping the hits
nanoslast search --bam reads.bam --query genes.fasta -o hits.tsv
nanoslast export --bam reads.bam --signals pod5_dir/ --hits hits.tsv -o dataset.h5

# plot one read with its bases
nanoslast plot --bam reads.bam --signals pod5_dir/ --read-id <read_id> -o read.png
```

Hits are copies of the query sequences, made by joining neighbouring BLAST HSPs; by default a hit must cover ≥ 80 % of the query and reach ≥ 80 % of the query's self-alignment score (`--min-coverage`, `--min-score`). BLAST runs as megablast by default; `--task`, `--evalue` and `--blast-options="…"` change it (e.g. `--task blastn-short` for short sequences). `search` and each export write `<output>.report.txt` with the settings and the numbers of found, exported and skipped hits.

| Output | Content |
|---|---|
| hits file (`search`) | one row per hit: read, query, position in the read, strand, coverage, score |
| `.csv` | one row per hit, whole read signal (pA) and the position of the hit in it |
| `.parquet` | one row per hit, signal segment of the hit (pA) |
| `.h5` | one row per read: whole signal (z-score), all hits, negative samples; the format of the NanoResistNet training datasets |

Details of all formats: [Output formats](docs/manual.md#7-output-formats). `nanoslast <command> --help` lists all options.

## Testing

`tests/test_workflow.py` runs the whole workflow on synthetic reads whose signal positions are known exactly, and checks every exported segment, the HDF5 dataset, the input checks and the command line. It needs BLAST+ (and samtools for the BAM test):

```bash
python tests/test_workflow.py            # runs in a temporary directory
python tests/test_workflow.py test_out   # keeps the test data, exports and plots in test_out/
```

## Changelog

See [CHANGELOG.md](CHANGELOG.md). Version 2.0.0 removed the graphical interface (the last version with it is 1.1.0). Version 1.1.0 fixed the conversion of BLAST hit positions to signal positions; signal segments exported with version 1.0.0 are shifted by a few samples.

## Project structure

```
NANOSLAST/
├── pyproject.toml             # package definition and dependencies
├── CHANGELOG.md
├── docs/manual.md             # full manual (and the logo)
├── tests/test_workflow.py     # end-to-end test on synthetic data
└── nanoslast/
    ├── cli.py                 # command line interface
    ├── sam_reader.py          # SAM/BAM reading, split reads, move table → signal positions
    ├── signal_reader.py       # POD5/FAST5 index and signal loading
    ├── checks.py              # input checks
    ├── convert_sam_file.py    # reads → FASTA
    ├── create_blastdb.py      # BLAST database (makeblastdb)
    ├── blasting.py            # BLAST search (blastn), whole search workflow
    ├── hits.py                # joining HSPs into hits, filters, hits file
    ├── parquet_tool.py        # export: signal positions of hits, CSV/Parquet
    ├── dataset.py             # export: HDF5 training dataset
    ├── report.py              # run reports
    └── plotting.py            # signal plot with bases
```

## Citation

If you use NANOSLAST in your work, please cite:

> M. Jakubickova, M. Zbudilova, M. Suriak, H. Vitkova, "NANOSLAST: Nanopore Signal Local Alignment and Segmentation Tool," in *Proceedings of the 11th World Congress on Electrical Engineering and Computer Systems and Sciences (EECSS'25)*, Paris, France, Aug. 2025, Paper No. ICBES 151. doi: [10.11159/icbes25.151](https://doi.org/10.11159/icbes25.151)

## Acknowledgements

This work was supported by GA23-05845S and NW24-09-00126.

Department of Biomedical Engineering, Brno University of Technology.
