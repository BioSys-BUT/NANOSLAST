# NANOSLAST

**Nanopore Signal Local Alignment and Segmentation Tool**

NANOSLAST is a lightweight Python GUI application that links local sequence similarity search (BLAST) to raw Oxford Nanopore signals (squiggles). It lets you find reads containing a sequence of interest, locate the exact part of the raw signal that corresponds to it, and visualize or export that signal segment.

It works with both raw nanopore data formats, **FAST5** and **POD5**, and with basecalled reads in **SAM** format (Dorado or Guppy) that include move tables.

Typical uses:

- extracting signals of specific motifs or functional elements,
- comparing raw signal characteristics between flowcells or sequencers,
- building signal–sequence datasets for training machine learning models (basecalling, methylation detection, etc.).

## Requirements

- Python 3.10 or newer
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html): `makeblastdb` and `blastn` must be on your `PATH`
- Tkinter (bundled with Python on Windows and macOS; on Debian/Ubuntu install `python3-tk`)
- Python packages listed in [requirements.txt](requirements.txt)

> **Note:** NANOSLAST calls BLAST+ through `Bio.Blast.Applications`, which was removed in Biopython 1.86. `requirements.txt` therefore pins `biopython<1.86`.

## Installation

```bash
git clone https://github.com/BioSys-BUT/NANOSLAST.git
cd NANOSLAST

python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

Install BLAST+ separately, for example:

```bash
# Debian/Ubuntu
sudo apt install ncbi-blast+
# conda
conda install -c bioconda blast
```

On Windows, use the installer from NCBI and add its `bin` directory to `PATH`.

## Running

Run the application from the repository directory:

```bash
python main.py
```

The GUI loads `nanoslast_logo.png` relative to the current working directory, so start the app from the repository directory.

## Input data

| Data | Format | Notes |
|---|---|---|
| Basecalled reads | SAM | Must contain move tables (`mv:B:c` tag) and trimmed samples (`ts:i` tag). With Dorado, use `dorado basecaller ... --emit-moves --emit-sam`. |
| Raw signals | FAST5 or POD5 | Pass a directory; it is searched recursively, and FAST5 and POD5 files can be mixed. |
| Query sequence | text or FASTA | Paste a sequence or select a FASTA file. |

## Workflow

The application has three tabs, **Search**, **Export** and **Plot**. Each can be used on its own, or they can be used in sequence.

### 1. Search

Find reads that contain a query sequence.

1. **SAM to FASTA:** select the SAM file and convert it to FASTA.
2. **BLAST database:** select the FASTA file and a directory where the nucleotide BLAST database (named `blastdb`) will be created with `makeblastdb`.
3. **Run BLAST:** select the database directory, paste the query sequence (or tick *Use FASTA file instead* and select a FASTA file) and run `blastn`.

The result is saved as a tab-separated TSV file with a header and standard BLAST columns: `qseqid`, `sseqid` (the read ID), `pident`, `length`, `mismatch`, `gapopen`, `qstart`, `qend`, `sstart`, `send`, `evalue`, `bitscore`.

The search uses `evalue = 0.001` and `max_target_seqs = 1000`. Hits with an alignment shorter than the query length minus 10 nt are discarded.

### 2. Export

Extract the raw signal segments that correspond to the BLAST hits.

1. Select the directory with the POD5/FAST5 files, the SAM file and the TSV file from the Search tab.
2. Click **CREATE PARQUET FILE** or **CREATE CSV FILE** and choose the output path.

For each hit, NANOSLAST combines the alignment coordinates (`sstart`, `send`) with the move table, stride and trimmed samples from the SAM file to calculate the start and end of the signal segment. The signal is converted to picoamperes.

Output columns:

- **CSV:** `read_id`, `sequence` (matched subsequence), `start_idx`, `stop_idx` (positions in the read), `signal` (values separated by `;`), `start_signal`, `stop_signal` (sample indices), `direction` (`1` = forward, `-1` = reverse)
- **Parquet:** `read_id`, `signal` (the extracted signal segment)

If a read has more than one hit in the TSV, only the last one is used.

### 3. Plot

Visualize a whole read's signal together with its basecalled sequence.

1. Select the SAM file and the directory with POD5/FAST5 files.
2. Enter a read ID and click **PLOT SIGNAL**.

The signal is plotted in pA against time in ms. Each base is highlighted with a coloured segment: A red, T blue, C yellow, G green.

## Project structure

| File | Purpose |
|---|---|
| `main.py` | Application entry point, sets up the tabs |
| `search_part.py` | Search tab GUI |
| `convert_sam_file.py` | SAM → FASTA conversion |
| `create_blastdb.py` | BLAST database creation (`makeblastdb`) |
| `blasting.py` | BLAST search (`blastn`) and hit filtering |
| `export_part.py` | Export tab GUI |
| `parquet_tool.py` | Signal segment calculation, FAST5/POD5 loading, CSV/Parquet export |
| `plot_part.py` | Plot tab GUI |
| `plotting.py` | Signal loading and plotting for a single read |
| `plot_specific_part.py` | Experimental tab for plotting a specific segment (currently disabled in `main.py`) |

## Citation

If you use NANOSLAST in your work, please cite:

> M. Jakubickova, M. Zbudilova, M. Suriak, H. Vitkova, "NANOSLAST: Nanopore Signal Local Alignment and Segmentation Tool," in *Proceedings of the 11th World Congress on Electrical Engineering and Computer Systems and Sciences (EECSS'25)*, Paris, France, Aug. 2025, Paper No. ICBES 151. doi: [10.11159/icbes25.151](https://doi.org/10.11159/icbes25.151)

## Acknowledgements

This work was supported by GA23-05845S and NW24-09-00126.

Department of Biomedical Engineering, Brno University of Technology.
