"""
BLAST HSPs → copies of query sequences in reads ("hits"), and their filtering.

BLAST often reports one copy of a gene as several HSPs (split by indels or
low-quality stretches). HSPs of the same query on the same strand of a read
that lie next to each other and cover different parts of the query are joined
into one hit. For each hit:

    coverage = union of the query parts covered by its HSPs / query length
    score    = sum of HSP bitscores / bitscore of the query aligned to itself

so a complete, error-free copy has coverage 1 and score 1.
"""
import csv
from dataclasses import asdict, dataclass, fields

HSP_COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]

# Neighbouring HSPs (at most this many bases apart in the read) are joined into one hit
# unless the new HSP mostly covers a part of the query that is already covered,
# which means a second copy of the query follows.
MAX_HSP_GAP = 50
MAX_QUERY_OVERLAP = 0.5

# Hits reaching this close to the start or end of the read may be cut off by the read end
READ_END_DISTANCE = 10


@dataclass
class HitFilters:
    min_coverage: float = 0.8   # fraction of the query covered
    min_score: float = 0.8      # summed bitscore relative to the query's self-alignment

    def accepts(self, hit):
        return hit.coverage >= self.min_coverage and hit.score >= self.min_score


@dataclass
class Hit:
    read_id: str
    query_id: str
    read_start: int       # first base of the hit in the read (1-based, read/signal orientation)
    read_end: int         # last base of the hit (inclusive)
    strand: str           # "+" if the query matches the read as it is, "-" if its reverse complement
    coverage: float
    score: float
    identity: float
    bitscore: float       # sum of HSP bitscores
    n_hsps: int
    at_read_end: bool     # hit reaches the start/end of the read, the copy may be incomplete
    read_length: int
    passed: bool = True   # passed the filters of the search


def read_hsps(filename):
    with open(filename) as f:
        for row in csv.DictReader(f, delimiter="\t", fieldnames=HSP_COLUMNS):
            yield row


def _union_length(intervals):
    total, end = 0, 0
    for lo, hi in sorted(intervals):
        if hi > end:
            total += hi - max(lo, end + 1) + 1
            end = hi
    return total


def _overlap(interval, intervals):
    lo, hi = interval
    return _union_length([(max(lo, a), min(hi, b)) for a, b in intervals if max(lo, a) <= min(hi, b)])


def hsps_to_hits(hsps, self_scores, max_gap=MAX_HSP_GAP):
    """
    Join HSPs into hits. self_scores: {query_id: bitscore of the query against itself}.
    max_gap: largest distance (bases) between HSPs of one hit in the read.
    """
    groups = {}
    for h in hsps:
        sstart, send = int(h["sstart"]), int(h["send"])
        strand = "+" if sstart <= send else "-"
        groups.setdefault((h["sseqid"], h["qseqid"], strand), []).append(
            (min(sstart, send), max(sstart, send), int(h["qstart"]), int(h["qend"]),
             float(h["pident"]), int(h["length"]), float(h["bitscore"]),
             int(h["qlen"]), int(h["slen"]))
        )

    hits = []
    for (read_id, query_id, strand), group in groups.items():
        clusters = []
        for hsp in sorted(group):
            read_lo, read_hi, q_lo, q_hi = hsp[:4]
            current = clusters[-1] if clusters else None
            if (current is not None
                    and read_lo <= max(c[1] for c in current) + max_gap
                    and _overlap((q_lo, q_hi), [(c[2], c[3]) for c in current])
                    <= MAX_QUERY_OVERLAP * (q_hi - q_lo + 1)):
                current.append(hsp)
            else:
                clusters.append([hsp])

        for cluster in clusters:
            query_length, read_length = cluster[0][7], cluster[0][8]
            read_start = min(c[0] for c in cluster)
            read_end = max(c[1] for c in cluster)
            bitscore = sum(c[6] for c in cluster)
            aligned = sum(c[5] for c in cluster)
            hits.append(Hit(
                read_id=read_id,
                query_id=query_id,
                read_start=read_start,
                read_end=read_end,
                strand=strand,
                coverage=_union_length([(c[2], c[3]) for c in cluster]) / query_length,
                score=bitscore / self_scores.get(query_id, max(c[6] for c in cluster)),
                identity=sum(c[4] * c[5] for c in cluster) / aligned,
                bitscore=bitscore,
                n_hsps=len(cluster),
                at_read_end=(read_start <= READ_END_DISTANCE
                             or read_end > read_length - READ_END_DISTANCE),
                read_length=read_length,
            ))

    hits.sort(key=lambda h: (h.read_id, h.read_start))
    return hits


def passed_hits(hits):
    return [hit for hit in hits if hit.passed]


def best_hit_per_read(hits):
    best = {}
    for hit in hits:
        if hit.read_id not in best or hit.bitscore > best[hit.read_id].bitscore:
            best[hit.read_id] = hit
    return [hit for hit in hits if best[hit.read_id] is hit]


HIT_COLUMNS = [f.name for f in fields(Hit)]


QUERIES_LINE = "# queries: "


def write_hits(hits, filename, query_names=None):
    """
    Write hits to a TSV file. The first line lists all query sequences that were
    searched ("# queries: a;b;c"), including those without any hit, so that the
    gene list of a dataset made from the file is complete.
    """
    with open(filename, "w", newline="") as f:
        if query_names:
            f.write(QUERIES_LINE + ";".join(query_names) + "\n")
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(HIT_COLUMNS)
        for hit in hits:
            row = asdict(hit)
            for key in ("coverage", "score"):
                row[key] = round(row[key], 4)
            row["identity"] = round(row["identity"], 3)
            writer.writerow(row[c] for c in HIT_COLUMNS)


def read_hits(filename):
    """
    Read hits saved by write_hits (all hits, see Hit.passed). Also accepts the
    BLAST TSV written by NANOSLAST 1.1.0 (one HSP per row, all passed), where
    each HSP becomes one hit.
    """
    with open(filename) as f:
        reader = csv.DictReader((line for line in f if not line.startswith("#")), delimiter="\t")
        rows = list(reader)
        columns = reader.fieldnames or []

    if "sseqid" in columns:
        hits = []
        for row in rows:
            sstart, send = int(row["sstart"]), int(row["send"])
            hits.append(Hit(
                read_id=row["sseqid"], query_id=row["qseqid"],
                read_start=min(sstart, send), read_end=max(sstart, send),
                strand="+" if sstart <= send else "-",
                coverage=float("nan"), score=float("nan"),
                identity=float(row["pident"]), bitscore=float(row["bitscore"]),
                n_hsps=1, at_read_end=False, read_length=0,
            ))
        return hits

    missing = set(HIT_COLUMNS) - set(columns)
    if missing:
        raise ValueError(f"{filename} is not a NANOSLAST hits file (missing columns: {', '.join(sorted(missing))})")

    return [
        Hit(
            read_id=row["read_id"], query_id=row["query_id"],
            read_start=int(row["read_start"]), read_end=int(row["read_end"]),
            strand=row["strand"], coverage=float(row["coverage"]), score=float(row["score"]),
            identity=float(row["identity"]), bitscore=float(row["bitscore"]),
            n_hsps=int(row["n_hsps"]), at_read_end=row["at_read_end"] == "True",
            read_length=int(row["read_length"]), passed=row["passed"] == "True",
        )
        for row in rows
    ]


def read_query_names(filename):
    """Query sequences listed in the first line of a hits file, or None (older files)."""
    with open(filename) as f:
        first = f.readline().rstrip("\n")
    return first[len(QUERIES_LINE):].split(";") if first.startswith(QUERIES_LINE) else None
