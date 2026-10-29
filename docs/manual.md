# NANOSLAST manual

NANOSLAST finds sequences of interest (for example resistance genes) in basecalled nanopore reads with BLAST and locates the exact part of the raw signal that corresponds to each of them. It exports these signal segments, or a complete HDF5 training dataset for signal-based models such as NanoResistNet, and plots signals with their bases.

Contents:

1. [Installation](#1-installation)
2. [Preparing the input data](#2-preparing-the-input-data)
3. [Quick start](#3-quick-start)
4. [Commands](#4-commands)
5. [How hits are found](#5-how-hits-are-found)
6. [How signal positions are calculated](#6-how-signal-positions-are-calculated)
7. [Output formats](#7-output-formats)
8. [Building a training dataset](#8-building-a-training-dataset)
9. [Large runs and performance](#9-large-runs-and-performance)
10. [Validation](#10-validation)
11. [Troubleshooting](#11-troubleshooting)
12. [Using NANOSLAST from Python](#12-using-nanoslast-from-python)

---

## 1. Installation

Requirements:

- Python 3.10 or newer
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) (`makeblastdb`, `blastn`) on `PATH`
- [samtools](https://www.htslib.org/) on `PATH`, only to read BAM files

```bash
python -m venv venv
source venv/bin/activate                 # Windows: venv\Scripts\activate
pip install git+https://github.com/BioSys-BUT/NANOSLAST.git
nanoslast --version
```

A specific version: `pip install git+https://github.com/BioSys-BUT/NANOSLAST.git@v2.0.0`. For development, clone the repository and run `pip install -e .`.

BLAST+ and samtools, for example with conda: `conda install -c bioconda blast samtools`.

## 2. Preparing the input data

NANOSLAST needs three inputs:

| Input | Format | Notes |
|---|---|---|
| Basecalled reads | SAM or BAM | With **move tables** (`mv` tag) and `ts` tag. Aligned or unaligned. |
| Raw signals | POD5 and/or FAST5 | The files the reads were basecalled from. A directory, searched recursively. Single-read and multi-read FAST5 can be mixed with POD5. |
| Query sequences | FASTA | One or more sequences (genes, motifs); each is searched separately. The FASTA IDs become the gene labels. |

### Basecalling with Dorado

The move table links each basecalled base to its part of the raw signal. Dorado writes it only with `--emit-moves`:

```bash
dorado basecaller <model> pod5_dir/ --emit-moves --kit-name SQK-RBK114-24 > reads.bam
```

- **Existing BAM files without move tables cannot be used**; the reads have to be basecalled again with `--emit-moves` (NANOSLAST stops with an explanation if the moves are missing).
- Use the same POD5 files for basecalling and for NANOSLAST.
- `--kit-name` (barcode trimming and classification) is fine: Dorado keeps the move table consistent with the trimmed sequence, and NANOSLAST stores the barcode (`BC` tag) in the HDF5 dataset.
- BAM can be used directly (needs samtools). A SAM file (`--emit-sam`, or `samtools view -h reads.bam > reads.sam`) works without samtools, e.g. on Windows.
- **Split reads:** Dorado splits some reads into subreads (tags `pi` = parent read, `sp` = start in the parent signal, `ns` = number of samples). NANOSLAST handles them: the signal of a subread is its part of the parent's signal.
- **Aligned files:** secondary and supplementary alignments are skipped and reads on the reverse strand are reverse complemented back to the orientation in which they passed the pore.

## 3. Quick start

```bash
# signal segments of all copies of the genes in the reads
nanoslast run --bam reads.bam --signals pod5_dir/ --query genes.fasta -o segments.parquet

# HDF5 training dataset (whole reads, all genes, negative samples)
nanoslast run --bam reads.bam --signals pod5_dir/ --query genes.fasta --run-id 230720-FNB-mk1c -o dataset.h5

# the same in two steps, keeping the hits
nanoslast search --bam reads.bam --query genes.fasta -o hits.tsv
nanoslast export --bam reads.bam --signals pod5_dir/ --hits hits.tsv -o dataset.h5

# plot one read with its bases
nanoslast plot --bam reads.bam --signals pod5_dir/ --read-id <read_id> -o read.png
```

`search` and every export write a text report next to the output (`<output>.report.txt`) with the settings and results.

## 4. Commands

All commands take the reads either as `--sam FILE` or `--bam FILE` (one of them; BAM is read with samtools). `nanoslast <command> --help` lists all options.

### `nanoslast search`

Finds copies of the query sequences in the reads and saves them to a TSV file ([hits file](#hits-file-search)).

| Option | Default | Meaning |
|---|---|---|
| `--sam FILE` or `--bam FILE` | required | basecalled reads |
| `--query FASTA` or `--sequence SEQ` | required | query sequences |
| `-o FILE` | required | output hits file (TSV) |
| `--db-dir DIR` | temporary | keep the BLAST database of the reads in DIR and reuse it next time |
| `--threads N` | 1 | BLAST threads |
| `--task` | `megablast` | BLAST task: `megablast`, `dc-megablast`, `blastn`, `blastn-short` ([BLAST settings](#blast-settings)) |
| `--evalue` | 0.001 | e-value threshold of HSPs |
| `--blast-options="…"` | none | other `blastn` options passed as they are |
| `--min-coverage`, `--min-score` | 0.8, 0.8 | [hit filters](#5-how-hits-are-found) |

Reads do not need move tables for `search`.

### `nanoslast export`

Exports the signal of the hits from a hits file.

| Option | Default | Meaning |
|---|---|---|
| `--sam FILE` or `--bam FILE` | required | basecalled reads with move tables |
| `--signals DIR` | required | POD5/FAST5 files |
| `--hits FILE` | required | hits file from `search` (or a BLAST TSV from NANOSLAST 1.1.0) |
| `-o FILE` | required | `.csv`, `.parquet` or `.h5` ([formats](#7-output-formats)) |
| `--format` | from extension | `csv`, `parquet`, `hdf5` |
| `--best-only` | off | only the best hit (highest bitscore) of each read |
| `--neg-ratio` | 0.3 | HDF5: negative reads per read with a gene, in each barcode |
| `--seed` | 42 | HDF5: random seed for choosing negative reads |
| `--normalization` | `zscore` | HDF5: `zscore` (whole read) or `pa` |
| `--run-id` | name of the SAM/BAM file | HDF5: run ID stored in the dataset |
| `--gene-names FASTA` | queries listed in the hits file | HDF5: list of all genes stored in the `gene_names` attribute |

### `nanoslast run`

`search` and `export` in one step. Takes the options of both (except `--hits`, `--gene-names`; the gene names are taken from `--query`), plus `--hits-output FILE` to keep the hits file. All inputs are checked before BLAST starts.

### `nanoslast plot`

Plots the raw signal of one read (or subread) with each base highlighted (A red, T blue, C yellow, G green), time in ms on the x-axis and current in pA on the y-axis.

| Option | Meaning |
|---|---|
| `--sam FILE` or `--bam FILE`, `--signals DIR`, `--read-id ID` | the read to plot |
| `-o FILE` | save to a file (`.png`, `.pdf`, `.svg`, …) instead of opening a window |
| `--dpi N` | resolution of the saved plot (default 150) |

## 5. How hits are found

1. The reads are converted to FASTA (one record per read, in signal orientation) and a BLAST nucleotide database is built (`makeblastdb`).
2. All query sequences are searched with `blastn` ([BLAST settings](#blast-settings): by default task megablast, e-value 0.001, at most 1,000,000 reads per query).
3. **HSPs are joined into hits.** BLAST often reports one copy of a gene as several HSPs, split by an indel or a low-quality stretch. HSPs of the same query on the same strand of a read are joined into one hit when they are at most 50 bases apart and cover different parts of the query. An HSP that covers a part of the query already covered starts a new hit (a second copy of the gene).
4. For each hit:
   - **coverage** = the part of the query covered by its HSPs (union) / query length,
   - **score** = sum of the HSP bitscores / bitscore of the query aligned to itself (1 = complete error-free copy),
   - **identity** = % identity weighted by alignment length,
   - **at_read_end** = the hit reaches within 10 bases of the read start or end, so the copy may be cut off.
5. **Filters** decide which hits are used:

   | Filter | Default | Option |
   |---|---|---|
   | coverage | ≥ 0.8 | `--min-coverage` |
   | score | ≥ 0.8 | `--min-score` |

   The identity of each hit is reported in the hits file; the score already decreases with mismatches and gaps, so there is no separate identity filter.

   The defaults correspond to the "mild" criterion of a gene being present in a read used in the NanoResistNet ground truth (coverage ≥ 0.8, score ≥ 0.8); a stricter choice is `--min-coverage 0.95 --min-score 0.9`. Hits that do not pass are still written to the hits file with `passed=False`. They are not exported, but their reads are also never used as negative samples.

A read can therefore have several hits: different genes, or several copies of one gene. All of them are exported unless `--best-only` is used.

**Genes interrupted by an insertion.** When a gene in a read is interrupted by an insertion longer than 50 bases (e.g. an insertion sequence), its parts are separate hits; each covers only part of the gene and usually does not pass the filters, so the gene is not labelled (and the read is not used as a negative sample). This is intended: joining the parts would label the inserted sequence as the gene. The NanoResistNet ground truth sums HSPs over the whole read and counts such a gene. In run 230720 this happened in 1 of 204 reads with genes.

### BLAST settings

| Option | Default | Meaning |
|---|---|---|
| `--task` | `megablast` | `megablast` (fast, for genes and long sequences), `dc-megablast` (more sensitive for divergent sequences), `blastn` (sensitive, slower), `blastn-short` (sequences shorter than ~50 nt, e.g. barcodes and motifs) |
| `--evalue` | 0.001 | HSPs with a higher e-value are not reported |
| `--blast-options="…"` | none | any other `blastn` options, passed as they are, e.g. `--blast-options="-word_size 7 -reward 1 -penalty -2"`. Write them with `=` and in quotes. Options that NANOSLAST sets itself (`-query`, `-db`, `-out`, `-outfmt`, `-subject`, `-num_threads`, `-task`, `-evalue`) are refused. |

Useful `blastn` options for `--blast-options`:

| Option | Default | Meaning |
|---|---|---|
| `-word_size N` | of the task: megablast 28, dc-megablast 11, blastn 11, blastn-short 7 | smaller = more sensitive and slower |
| `-max_target_seqs N` | 1,000,000 (set by NANOSLAST) | maximum number of reads reported for each query |
| `-reward`, `-penalty`, `-gapopen`, `-gapextend` | of the task | scoring of matches, mismatches and gaps |
| `-perc_identity P` | 0 | report only HSPs with at least P % identity |

With megablast (word size 28), queries much shorter than ~30 nt are not found. For short sequences use `--task blastn-short` (word size 7; this is what was used for barcodes in run 230720).

The score of a hit is relative to the self-alignment of the query computed **with the same task, word size and options**, so the filters keep their meaning with any settings. The settings are recorded in the search report.

## 6. How signal positions are calculated

The move table (`mv:B:c,stride,m1,m2,…`) has one element for each block of `stride` signal samples after the first `ts` samples (trimmed by the basecaller). A `1` marks the block in which a new base starts. So base *k* (1-based) starts at sample

    ts + i_k × stride

where *i_k* is the 0-based position of the *k*-th `1` in the move table, and ends where base *k + 1* starts (the last base ends at `ts + len(moves) × stride`).

A hit covering bases *first*–*last* (BLAST coordinates, 1-based, inclusive) is the signal range `[start of first, start of last + 1)`, i.e. from the first sample of the first base to the last sample of the last base. Positions refer to the signal of the read; for a subread of a split read, this is its part `[sp, sp + ns)` of the parent's signal.

Basecalling only the exported segment gives back the sequence of the hit (99.3 % median identity, boundaries within 1–3 bases; see [validation](#10-validation)).

## 7. Output formats

### Hits file (`search`)

Tab-separated, one row per hit. The first line lists all query sequences that were searched, also those without any hit (`# queries: name1;name2;…`), so that a dataset exported from the file has the complete gene list.

| Column | Meaning |
|---|---|
| `read_id` | read (or subread) containing the hit |
| `query_id` | query sequence (FASTA ID; `Query` for `--sequence`) |
| `read_start`, `read_end` | first and last base of the hit in the read (1-based, inclusive, signal orientation) |
| `strand` | `+` if the query matches the read as it is, `-` if its reverse complement |
| `coverage`, `score`, `identity` | see [section 5](#5-how-hits-are-found) |
| `bitscore` | summed bitscore of the HSPs |
| `n_hsps` | number of HSPs joined into the hit |
| `at_read_end` | the copy may be cut off by the read start or end |
| `read_length` | length of the read (bases) |
| `passed` | the hit passed the filters |

### CSV (`.csv`)

One row per hit, with the **whole signal of the read**:

`read_id`, `query_id`, `sequence` (matched part of the read), `start_idx`, `stop_idx` (bases, 1-based inclusive), `signal` (whole read signal in pA, values separated by `;`), `start_signal`, `stop_signal` (0-based, `stop_signal` exclusive: the hit is `signal[start_signal:stop_signal]`), `direction` (`1`/`-1`), `coverage`, `score`, `at_read_end`.

CSV stores signals as text and becomes very large; use Parquet or HDF5 for more than a few thousand reads.

### Parquet (`.parquet`)

One row per hit with **only the signal segment of the hit** (`signal`, list of float32 in pA) and the same columns as CSV except `sequence`.

### HDF5 training dataset (`.h5`)

The format of the NanoResistNet training datasets (e.g. `dataset_19G_<run>.h5`): **one row per read** with its whole signal and all its hits, plus reads without any hit as negative samples.

Root attributes:

| Attribute | Value |
|---|---|
| `run_id` | `--run-id` |
| `n_signals` | number of rows |
| `creation_date` | ISO timestamp |
| `gene_names` | JSON list of all query names |
| `normalization` | `zscore` or `pa` |
| `neg_ratio`, `seed` | settings of the negative samples |
| `nanoslast_version`, `source_hits`, `source_reads`, `align_ratios_meaning` | provenance |

Datasets (one element per row; several hits of a read are joined by `;`, ordered by position):

| Dataset | Type | Content | Negative sample |
|---|---|---|---|
| `signals` | vlen float32 | whole signal of the read: z-score `(x - mean) / std` over the whole read (default), or pA | same |
| `read_ids` | S64 | read ID | |
| `run_ids` | S64 | run ID | |
| `barcodes` | S32 | barcode from the `BC` tag, e.g. `barcode13`; `unclassified` without a barcode | |
| `gene_labels` | S500 | query names, e.g. `blaTEM-1A_1_HM749966;sul2_1_AF542061` | `no_gene` |
| `start_signals` | S200 | first sample of each hit in `signals` | `0` |
| `stop_signals` | S200 | end of each hit (exclusive) | length of the signal |
| `bitscores` | S200 | summed bitscore of each hit (`%.2f`) | `-1.0` |
| `align_ratios` | S200 | coverage of the gene by each hit (`%.4f`) | `-1.0` |
| `strands` | S200 | `+`/`-` of each hit | empty |
| `is_negative` | int8 | 1 for negative samples | 1 |

Reading a row:

```python
import h5py

with h5py.File("dataset.h5") as f:
    i = 0
    signal = f["signals"][i]
    labels = f["gene_labels"][i].decode().split(";")
    starts = [int(x) for x in f["start_signals"][i].decode().split(";")]
    stops = [int(x) for x in f["stop_signals"][i].decode().split(";")]
    for label, start, stop in zip(labels, starts, stops):
        gene_signal = signal[start:stop]
```

Differences from the datasets made by the earlier scripts (`new_pipeline_*.py`):

- **Hit positions are exact.** The earlier scripts started each gene one stride (6 samples) too late and ended it before the end of the last base; in run 230720 the starts of 95 % of genes differ by exactly +6 samples.
- `align_ratios` is the fraction of the gene covered by the hit (≤ 1), not alignment length / gene length.
- `barcodes` is always `barcodeNN` (the earlier scripts mixed `barcodeNN` and `SQK-RBK114-24_barcodeNN`).
- `strands` is added; reads containing hits that did not pass the filters are excluded from the negative samples; split reads are included.
- Gene labels are the FASTA IDs of the queries (with `geny.fa`, fosA is `fosA6_1_KU254579`).

If a value would not fit the fixed string width (e.g. a read with very many genes), the column is widened and a warning is printed, instead of silently cutting the value.

### Reports (`<output>.report.txt`)

Plain text files written next to the outputs:

- **search report** (`hits.tsv.report.txt`): NANOSLAST version, command, reads, BLAST database, query sequences, BLAST settings, joining of HSPs and filters, number of hits found and passed, passed hits and hits below the filters per query, run time;
- **export report** (e.g. `dataset.h5.report.txt`): inputs, format, number of hits to export and exported (total and per query), skipped hits with reasons, run time; for HDF5 also the run ID, number of rows, reads with hits and negative samples, normalization, `neg_ratio`, seed, and per barcode the reads with hits / negative samples / reads without hits available.

Skipped hits:

| Reason | Meaning |
|---|---|
| read not found in the SAM/BAM file | the hits file belongs to other reads |
| read has no move table | basecalled without `--emit-moves` |
| hit positions outside of the basecalled read | the hits file belongs to another basecalling of the reads |
| signal not found in the POD5/FAST5 files | the signal directory does not contain the read (or its parent read) |

## 8. Building a training dataset

One dataset per sequencing run. On a PBS cluster, for example:

```bash
#!/usr/bin/env bash
#PBS -N basecall_moves
#PBS -l select=1:ncpus=8:ngpus=1:mem=32gb
set -euo pipefail
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.2.0 /path/to/run/pod5 \
    --device cuda:0 --kit-name SQK-RBK114-24 --emit-moves > reads.bam
```

```bash
#!/usr/bin/env bash
#PBS -N nanoslast_dataset
#PBS -l select=1:ncpus=8:mem=32gb
set -euo pipefail
source venv/bin/activate
nanoslast run --bam reads.bam --signals /path/to/run/pod5 --query geny.fa \
    --run-id 230720-FNB-mk1c --threads 8 --db-dir blastdb \
    --hits-output hits.tsv -o dataset_19G_230720-FNB-mk1c.h5
```

Notes:

- `--db-dir` keeps the BLAST database, so the dataset can be rebuilt with other filters or `--neg-ratio` without rebuilding the database: `nanoslast search … --db-dir blastdb -o hits.tsv` then `nanoslast export … --hits hits.tsv -o dataset.h5`.
- Negative samples are reads without any BLAST HSP of any query (not even one below the filters), chosen per barcode: `round(neg_ratio × reads with a gene in that barcode)`, reproducibly with `--seed`. Barcodes without reads with a gene get no negatives.
- `--normalization pa` stores the signal in pA instead of z-score.
- Check `dataset.h5.report.txt` for skipped hits and the number of negatives per barcode.

## 9. Large runs and performance

- The SAM/BAM file is read as a stream; it is never loaded into memory.
- POD5/FAST5 files are indexed once (read IDs only) and signals are read only from the files that contain the needed reads.
- The HDF5 dataset is written in batches of 500 reads, so memory does not grow with the dataset size. The negative-sample candidates (read IDs without hits) are kept in memory: roughly 100–200 bytes per read.
- `--threads` speeds up BLAST; `--db-dir` avoids rebuilding the BLAST database.
- Measured on 10 POD5 files of run 230720-FNB-mk1c (40,549 reads, 1.9 GB of POD5, 8 BLAST threads): `nanoslast run` to an HDF5 dataset took 18 s and 0.5 GB of memory. A whole run (about 60× more reads) is therefore expected to take tens of minutes; the BLAST database of all reads is the largest part, so keep it with `--db-dir` when building datasets repeatedly.

## 10. Validation

NANOSLAST was validated on run 230720-FNB-mk1c (Dorado 1.3.1, sup@v5.2.0) with the 19 NanoResistNet genes:

- **Signal segments give back the gene:** each exported segment was basecalled on its own. Median identity with the sequence BLAST found in the read: 99.2–99.3 %; the correct gene in 100 % of 577 segments; boundaries within 1–3 bases (the basecaller loses 1–2 bases at the edges of a segment).
- **Control:** segments shifted by ±250 samples match the part of the read predicted from the move table equally well, so the conversion is correct everywhere in the read.
- **Genes in reads match the ground truth:** for 204 reads, the set of genes found by NANOSLAST equals the NanoResistNet BLAST ground truth (mild criterion) in 204/204 reads.
- **HDF5 dataset:** same format as the existing datasets, identical signals (z-score), positions corrected by one stride, negatives without genes; details in `validace/README.md`.
- Automated test on synthetic data with known positions: `python tests/test_workflow.py`.

## 11. Troubleshooting

| Message | Cause and fix |
|---|---|
| `The reads in … have no move tables (mv tag)` | Basecall again with `dorado basecaller … --emit-moves`. |
| `makeblastdb not found` / `blastn not found` | Install BLAST+ and put it on `PATH`. |
| `Reading BAM files needs samtools on PATH` | Install samtools, or convert the BAM to SAM (`samtools view -h reads.bam > reads.sam`). |
| `No POD5 or FAST5 files found in …` | Wrong signal directory. |
| `Duplicate query names` / `…not nucleotides` | Fix the query FASTA (unique IDs, only nucleotide letters). |
| `Cannot tell the output format from the output file name` | Use `.csv`, `.parquet` or `.h5`, or `--format`. |
| Many hits skipped with *signal not found* | The POD5/FAST5 directory is not the one the reads were basecalled from. |
| Many hits skipped with *outside of the basecalled read* | The hits file comes from another basecalling of the reads; run `search` again on the same SAM/BAM. |
| Fewer hits than expected | Lower `--min-coverage`/`--min-score`; check `passed=False` rows in the hits file. Queries shorter than ~30 nt are not found by megablast. |
| `No display to show the plot on` | `nanoslast plot` without `-o` opens a window; on a server save the plot with `-o read.png`. |

## 12. Using NANOSLAST from Python

```python
from nanoslast.blasting import BlastSettings, read_queries, search_reads
from nanoslast.hits import HitFilters
from nanoslast.parquet_tool import export_signals

queries = read_queries("genes.fasta")
search_reads("reads.bam", queries, "hits.tsv", db_dir="blastdb",
             filters=HitFilters(min_coverage=0.95, min_score=0.9),
             blast=BlastSettings(task="megablast", evalue=1e-10), threads=8)
result = export_signals("hits.tsv", "reads.bam", "pod5_dir", "dataset.h5", mode="hdf5",
                        neg_ratio=0.3, run_id="my_run")
print(result.summary())
```

Signal positions of a single read:

```python
from nanoslast.sam_reader import find_sam_record, signal_range

record = find_sam_record("reads.bam", read_id)
start, stop = signal_range(record, first_base, last_base)   # samples in the read's signal
```
