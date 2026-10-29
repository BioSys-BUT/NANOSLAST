"""
Command line interface:

    nanoslast search  --sam reads.sam --query genes.fasta -o hits.tsv
    nanoslast export  --sam reads.sam --signals pod5_dir --hits hits.tsv -o segments.parquet
    nanoslast run     --sam reads.sam --signals pod5_dir --query genes.fasta -o segments.parquet
    nanoslast run     --bam reads.bam --signals pod5_dir --query genes.fasta -o dataset.h5
    nanoslast plot    --sam reads.sam --signals pod5_dir --read-id READ_ID [-o plot.png]
"""
import argparse
import os
import sys
import tempfile

from . import __version__

FORMATS = {".csv": "csv", ".parquet": "parquet", ".h5": "hdf5", ".hdf5": "hdf5"}


def log(message):
    print(message, file=sys.stderr)


def get_queries(args):
    from .blasting import read_queries

    if args.sequence:
        return [("Query", "".join(args.sequence.split()))]
    return read_queries(args.query)


def export_format(args):
    if args.format:
        return args.format
    extension = os.path.splitext(args.output)[1].lower()
    if extension in FORMATS:
        return FORMATS[extension]
    name = os.path.splitext(args.output)[0]
    raise ValueError(f"Cannot tell the output format from the output file name '{args.output}'. "
                     f"Give the file an extension, e.g. -o {name}.h5 (HDF5 dataset), "
                     f"-o {name}.parquet or -o {name}.csv, or use --format csv, parquet or hdf5")


def get_filters(args):
    from .hits import HitFilters

    return HitFilters(min_coverage=args.min_coverage, min_score=args.min_score)


def get_blast_settings(args):
    from .blasting import BlastSettings

    return BlastSettings(task=args.task, evalue=args.evalue, extra_args=args.blast_options or "")


def dataset_options(args, gene_names=None):
    return {"run_id": args.run_id, "gene_names": gene_names, "neg_ratio": args.neg_ratio,
            "seed": args.seed, "normalization": args.normalization}


def search(args, queries, output):
    from .blasting import search_reads

    n_hits, n_reads = search_reads(args.sam, queries, output, db_dir=args.db_dir,
                                   filters=get_filters(args), blast=get_blast_settings(args),
                                   threads=args.threads, progress=log)
    log(f"{n_hits} hits in {n_reads} reads saved to {output}, report in {output}.report.txt")
    return n_hits


def export(args, hits, mode, gene_names=None):
    from .parquet_tool import export_signals

    options = dataset_options(args, gene_names) if mode == "hdf5" else {}
    result = export_signals(hits, args.sam, args.signals, args.output, mode=mode,
                            best_only=args.best_only, progress=log, **options)
    log(result.summary())
    log(f"Saved to {args.output}, report in {result.report_file}")


def cmd_search(args):
    from .checks import check_blast_installed, check_sam

    queries = get_queries(args)
    get_blast_settings(args)  # invalid settings stop here, before the BLAST database is built
    check_blast_installed()
    check_sam(args.sam, require_moves=False)
    search(args, queries, args.output)


def cmd_export(args):
    from .blasting import read_queries

    mode = export_format(args)
    gene_names = [name for name, _ in read_queries(args.gene_names)] if args.gene_names else None
    export(args, args.hits, mode, gene_names)


def cmd_run(args):
    from .checks import check_blast_installed, check_sam, check_signals_dir

    mode = export_format(args)
    queries = get_queries(args)

    # check everything before the (possibly long) BLAST search
    get_blast_settings(args)
    check_blast_installed()
    check_sam(args.sam)
    check_signals_dir(args.signals)

    with tempfile.TemporaryDirectory(prefix="nanoslast_") as tmp_dir:
        hits = args.hits_output or os.path.join(tmp_dir, "hits.tsv")
        if search(args, queries, hits) == 0 and mode != "hdf5":
            log("No hits, nothing to export.")
            return
        export(args, hits, mode, gene_names=[name for name, _ in queries])


def cmd_plot(args):
    from .plotting import create_signal_figure, load_read

    record, signal_in_pA, sampling_rate = load_read(args.sam, args.signals, args.read_id)

    if args.output:
        fig = create_signal_figure(signal_in_pA, sampling_rate, record, args.read_id)
        fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
        log(f"Plot saved to {args.output}")
    else:
        import matplotlib.pyplot as plt

        if plt.get_backend().lower() == "agg":
            raise ValueError("No display to show the plot on, save it to a file with -o plot.png")
        create_signal_figure(signal_in_pA, sampling_rate, record, args.read_id, fig=plt.figure(figsize=(14, 6)))
        plt.show()


def add_reads_argument(parser):
    reads = parser.add_mutually_exclusive_group(required=True)
    reads.add_argument("--sam", dest="sam", metavar="SAM",
                       help="basecalled reads in SAM with mv and ts tags (dorado --emit-moves); "
                            "aligned or unaligned")
    reads.add_argument("--bam", dest="sam", metavar="BAM",
                       help="the same in BAM (needs samtools)")


def add_query_arguments(parser):
    query = parser.add_mutually_exclusive_group(required=True)
    query.add_argument("--query", metavar="FASTA",
                       help="FASTA file with one or more query sequences, each searched separately")
    query.add_argument("--sequence", metavar="SEQ", help="a single query sequence")
    parser.add_argument("--db-dir", metavar="DIR",
                        help="directory for the BLAST database of the reads; an existing database there "
                             "is reused (default: temporary directory, deleted afterwards)")
    parser.add_argument("--threads", type=int, default=1, help="BLAST threads (default: %(default)s)")

    from .blasting import BLAST_TASKS, BlastSettings
    blast_defaults = BlastSettings()
    blast = parser.add_argument_group("BLAST settings")
    blast.add_argument("--task", choices=BLAST_TASKS, default=blast_defaults.task,
                       help="blastn task: megablast for genes and other sequences longer than ~50 nt, "
                            "blastn-short for short sequences such as barcodes (default: %(default)s)")
    blast.add_argument("--evalue", type=float, default=blast_defaults.evalue,
                       help="e-value threshold of HSPs (default: %(default)s)")
    blast.add_argument("--blast-options", metavar="OPTIONS",
                       help='other blastn options passed as they are; write them with "=" and in quotes, '
                            'e.g. --blast-options="-word_size 7 -reward 1 -penalty -2"')

    from .hits import HitFilters
    defaults = HitFilters()
    filters = parser.add_argument_group(
        "hit filters",
        "BLAST HSPs of one query that lie next to each other in a read are joined into one hit "
        "(one copy of the query). A hit is used when it passes both thresholds; hits that "
        "do not pass are kept in the hits file (passed=False) and their reads are never used as "
        "negative samples.")
    filters.add_argument("--min-coverage", type=float, default=defaults.min_coverage, metavar="FRACTION",
                         help="minimum fraction of the query covered by the hit (default: %(default)s)")
    filters.add_argument("--min-score", type=float, default=defaults.min_score, metavar="FRACTION",
                         help="minimum summed bitscore of the hit relative to the bitscore of the "
                              "query aligned to itself (default: %(default)s)")


def add_export_arguments(parser, gene_names_option):
    parser.add_argument("-o", "--output", required=True,
                        help="output file: .csv, .parquet, or .h5 (HDF5 training dataset)")
    parser.add_argument("--format", choices=["csv", "parquet", "hdf5"],
                        help="output format; needed only if the output file has another extension "
                             "(default: .csv -> csv, .parquet -> parquet, .h5/.hdf5 -> hdf5)")
    parser.add_argument("--best-only", action="store_true",
                        help="export only the best hit (highest bitscore) of each read "
                             "(default: all hits)")

    dataset = parser.add_argument_group(
        "HDF5 dataset",
        "One row per read with its whole signal and all its hits, plus reads without any hit "
        "as negative samples (gene label no_gene).")
    dataset.add_argument("--neg-ratio", type=float, default=0.3,
                         help="negative reads per read with a hit, in each barcode (default: %(default)s)")
    dataset.add_argument("--seed", type=int, default=42,
                         help="random seed for choosing negative reads (default: %(default)s)")
    dataset.add_argument("--normalization", choices=["zscore", "pa"], default="zscore",
                         help="signal stored as z-score of the whole read, or in pA (default: %(default)s)")
    dataset.add_argument("--run-id", help="run ID stored in the dataset (default: name of the SAM/BAM file)")
    if gene_names_option:
        dataset.add_argument("--gene-names", metavar="FASTA",
                             help="FASTA with all query sequences, stored as gene_names "
                                  "(default: the queries listed in the hits file)")


def build_parser():
    parser = argparse.ArgumentParser(
        prog="nanoslast",
        description="NANOSLAST: Nanopore Signal Local Alignment and Segmentation Tool. "
                    "Finds sequences in basecalled nanopore reads with BLAST and extracts "
                    "the corresponding raw signal segments from POD5/FAST5 files.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command", required=True, metavar="COMMAND")

    signals_help = "directory with POD5 and/or FAST5 files (searched recursively)"

    p = commands.add_parser("search", help="find copies of query sequences in the reads (BLAST)",
                            description="Find copies of the query sequences in the reads and save them to a TSV file.")
    add_reads_argument(p)
    p.add_argument("-o", "--output", required=True, help="output TSV file with the hits")
    add_query_arguments(p)
    p.set_defaults(func=cmd_search)

    p = commands.add_parser("export", help="export signal segments of the hits",
                            description="Export the raw signal segment of each hit, or an HDF5 training dataset.")
    add_reads_argument(p)
    p.add_argument("--signals", required=True, help=signals_help)
    p.add_argument("--hits", required=True, help="TSV file with hits from 'nanoslast search'")
    add_export_arguments(p, gene_names_option=True)
    p.set_defaults(func=cmd_export)

    p = commands.add_parser("run", help="search and export in one step",
                            description="Find copies of the query sequences in the reads and export "
                                        "their signal segments or an HDF5 training dataset.")
    add_reads_argument(p)
    p.add_argument("--signals", required=True, help=signals_help)
    add_export_arguments(p, gene_names_option=False)
    p.add_argument("--hits-output", metavar="TSV", help="also keep the hits in this TSV file")
    add_query_arguments(p)
    p.set_defaults(func=cmd_run)

    p = commands.add_parser("plot", help="plot the signal of a read with its bases",
                            description="Plot the raw signal of one read with each basecalled base highlighted.")
    add_reads_argument(p)
    p.add_argument("--signals", required=True, help=signals_help)
    p.add_argument("--read-id", required=True, help="ID of the read to plot")
    p.add_argument("-o", "--output", help="save the plot to a file (png, pdf, svg...) instead of showing it")
    p.add_argument("--dpi", type=int, default=150, help="resolution of the saved plot (default: 150)")
    p.set_defaults(func=cmd_plot)


    return parser


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    # a command without any options (e.g. "nanoslast search") shows its full help
    if not argv or (len(argv) == 1 and argv[0] in ("search", "export", "run", "plot")):
        argv.append("--help")
    args = build_parser().parse_args(argv)
    try:
        args.func(args)
    except KeyboardInterrupt:
        log("Interrupted")
        return 130
    except Exception as e:
        log(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
