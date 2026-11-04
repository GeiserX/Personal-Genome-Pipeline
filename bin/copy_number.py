#!/usr/bin/env python3
"""copy_number.py: copy number at any locus of bin/cn_loci.tsv from short-read depth.

One method, one catalog row per locus: the depth over the locus's target
segments against the depth over single-copy control segments, as
`samtools bedcov` sums it. With T the target depth and C the control depth
(each combined as the row's target_depth and control_depth say),
    estimate = scale x T / C - offset
then guarded by a minimum control depth, flank agreement and the share of the
target depth from reads of mapping quality 0, and called as the row's `call`
column says. docs/cn-loci.md describes the columns and how to add a locus.

  copy_number.py bed [--catalog TSV] [--locus L ...]
      print the segments of the loci (all bedcov loci by default) as a BED
      (chrom, 0-based start, end, name), sorted by position within a locus
  copy_number.py call [--catalog TSV] --locus L --bedcov T [--bedcov-mapq1 T1]
                      [--json OUT.json] [--table OUT.tsv]
      read `samtools bedcov` over those segments (chrom, start, end, name,
      summed depth; the row's first bedcov pass) and, for a row with a MAPQ 0
      guard, the MAPQ >= 1 pass; write the locus's JSON block and, with
      --table, the step's metric<TAB>value table (bin/hba_copy_number.py and
      bin/lpa_kiv2.py lay those out). Exit 0 whatever the status; 2 on
      unreadable input; 3 (BROKEN_CONTROL) when a control segment has no
      read while the other controls have the depth the row needs: the row's
      coordinates are wrong, or the sample lacks that stretch, and no number
      is given.

  copy_number.py bed --kind interval
      the segments of every interval row (kind interval, one bedcov pass: the
      recurrent CNV regions, the carrier deletion alleles and AZF) in one BED
  copy_number.py intervals [--catalog TSV] --bedcov T --sex male|female|unknown
                           --out OUT.tsv
  copy_number.py intervals --failed --out OUT.tsv
      call every interval row from one `samtools bedcov` over that BED; write
      one line per row (locus, name, chrom, start, end, status, copies,
      expected, change, estimate, reason). A row on chrX or chrY is read
      against the sample's own count there: 1 on a male, 2 on a female's X;
      a chrY row on a female is 'not applicable', and without a sex both are
      'not assessed' (sex_unknown). A row whose control has no read is 'not
      assessed' (broken_control) and the others still run. --failed writes
      every row as 'not assessed' (check_failed), for a step whose depth
      pass failed.

  copy_number.py bed --psv --locus L
      the paralogous sequence variant (PSV) sites of a psv row (call psv:D):
      each site's position on the gene and on its pseudogene, one base each
      (name and name_pseudo), for `samtools mpileup -l`
  copy_number.py call --locus L --bedcov T --pileup P ...
      a psv row (CYP21A2) also reads `samtools mpileup` without a reference
      over those sites (P): the module total from depth as above, then at
      each site the share of reads with the gene's base, pooled over the gene
      and pseudogene positions, times that total gives the copies that carry
      the gene's sequence there. The table names a deletion (one module and
      one gene copy fewer), a possible conversion (fewer gene copies at some
      sites than at the others, or fewer than two) or none; it never says
      which chromosome carries the change.

  copy_number.py check-bam --header H
      refuse a BAM aligned to a reference with ALT or HLA contigs. H is
      `samtools view -H` of the BAM, or its contig names one per line ('-'
      reads stdin). On such a BAM a read that fits a primary locus and its
      ALT copy equally gets MAPQ 0 or lands on the ALT copy, so the depth of
      a duplicated gene drops and a depth-based count reads a deletion that
      is not there. Exit 0 and no output without such a contig; exit 4 with
      one REFUSED line on stderr with one. `call --bam-header H` does the
      same check first.

Two empty tables mean the BAM has no reads on that contig at all (bedcov is
skipped there): status 'not assessed'. Standard library only: it runs in the
Python and Parascopy images and under tests/test_copy_number.py.
"""
import argparse
import csv
import importlib
import json
import math
import os
import sys

CATALOG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cn_loci.tsv")
COLUMNS = ["locus", "kind", "chrom", "targets", "controls", "details", "target_depth", "control_depth",
           "scale", "offset", "bedcov", "mapq0_min", "mapq0_max", "min_control_depth", "max_control_ratio",
           "call", "interval", "readings", "input", "owner", "cross_check", "pmids"]
KINDS = ("gene", "interval", "organelle")
# GRCh38 names of the primary chromosomes (the ClinGen BED spells chrX 'chrx':
# scripts/ci/build-cn-intervals.py normalises it, and a row with any other
# spelling is refused here).
CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]}
SEXES = ("male", "female", "unknown")
INTERVAL_COLUMNS = ["locus", "name", "chrom", "start", "end", "status", "copies", "expected", "change",
                    "estimate", "reason"]
INPUTS = ("bedcov", "cnvpytor-level")
Z95 = 1.96
# The step table of a locus, by the module that lays it out (same folder).
TABLES = {"HBA": "hba_copy_number", "LPA_KIV2": "lpa_kiv2"}
# A psv row: the reads a site needs (gene and pseudogene bases pooled) before
# its copies are read, and the findings a reading must exist for.
PSV_MIN_READS = 30
# A site reads as converted when its copies estimate lies this far under the
# others' count (or under the usual count) or more. At 30x, about 70 reads at
# a site, read sampling alone gives an estimate of 2 copies a standard
# deviation of about 0.25; HG002's six sites read 1.62 to 2.82 (35x), so the
# real spread is wider and 0.6 is a screen setting, not a measured error rate.
PSV_MARGIN = 0.6
PSV_FINDINGS = ("none", "deletion", "conversion", "pseudogene_deletion")
BASES = "ACGT"


class CatalogError(ValueError):
    pass


class AltContigs(ValueError):
    """The BAM was aligned to a reference with ALT or HLA contigs."""


def alt_contigs(lines):
    """The ALT or HLA contigs of `samtools view -H` lines or of contig names
    one per line (the rule of scripts/validate-setup.sh: a name ending in
    _alt, or starting with HLA-)."""
    out = []
    for line in lines:
        line = line.rstrip("\n")
        if line.startswith("@"):
            if not line.startswith("@SQ"):
                continue
            name = next((f[3:] for f in line.split("\t") if f.startswith("SN:")), "")
        else:
            name = line.split("\t")[0].strip()
        if name.endswith("_alt") or name.startswith("HLA-"):
            out.append(name)
    return out


def check_bam(path):
    """Raise AltContigs, with the refusal line, when the header at path
    (or '-' for stdin) names an ALT or HLA contig."""
    if path == "-":
        found = alt_contigs(sys.stdin)
    else:
        with open(path, errors="replace") as f:
            found = alt_contigs(f)
    if found:
        raise AltContigs(f"REFUSED: the BAM was aligned to a reference with {len(found)} ALT or HLA contigs "
                         f"(first: {found[0]}). Reads that fit a primary locus and its ALT copy equally get MAPQ 0 "
                         "or land on the ALT copy, so depth-based copy number reads deletions that are not there. "
                         "Realign to the default no-ALT reference (docs/realignment.md).")


class BrokenControl(ValueError):
    """A control segment with no read while the other controls are covered."""


def _num(row, col, default=None):
    v = row[col]
    if v == ".":
        return default
    try:
        x = float(v)
    except ValueError:
        x = math.nan
    if not math.isfinite(x):
        raise CatalogError(f"{row['locus']}: {col} is {v!r}, not a number or '.'")
    return x


def _segments(row, col):
    """[(name, start, end, label)] of a targets, controls or details cell."""
    if row[col] == ".":
        return []
    out = []
    for part in row[col].split(";"):
        bits = part.split(":", 2)
        try:
            name, (s, e) = bits[0], bits[1].split("-")
            start, end = int(s), int(e)
        except (IndexError, ValueError):
            raise CatalogError(f"{row['locus']}: {col} segment {part!r} is not NAME:start-end[:label]")
        if not name or end <= start or start < 0:
            raise CatalogError(f"{row['locus']}: {col} segment {part!r} has no name or an empty span")
        out.append((name, start, end, bits[2] if len(bits) > 2 else ""))
    return out


def _readings(row):
    out = {}
    if row["readings"] == ".":
        return out
    for part in row["readings"].split("|"):
        k, sep, text = part.partition("=")
        if not sep or not k or not text:
            raise CatalogError(f"{row['locus']}: reading {part!r} is not KEY=text")
        out[k] = text
    return out


def _call(row):
    mode, sep, args = row["call"].partition(":")
    try:
        vals = [float(x) for x in args.split(",")]
    except ValueError:
        vals = []
    if not all(math.isfinite(v) for v in vals):
        vals = []
    if not sep or mode not in ("round", "cut", "class", "psv") or not vals \
            or (mode in ("round", "psv") and len(vals) != 1) or (mode == "class" and len(vals) != 2) \
            or vals != sorted(vals):
        raise CatalogError(f"{row['locus']}: call {row['call']!r} is not round:D, cut:A,B..., class:LOW,HIGH or psv:D")
    return mode, vals


def _psv_sites(row, details):
    """[(name, gene 0-based position, pseudogene 0-based position, gene
    bases, pseudogene bases, text)] of a psv row: every details segment is one
    base, labelled 'PSEUDO_POS GENE>PSEUDO text' (bases among ACGT; In2G's
    gene base is A or C)."""
    out = []
    for name, s, e, label in details:
        bits = label.split(" ", 2)
        try:
            pos = int(bits[0])
            gene, _, pseudo = bits[1].partition(">")
        except (IndexError, ValueError):
            raise CatalogError(f"{row['locus']}: PSV {name} label {label!r} is not 'PSEUDO_POS GENE>PSEUDO text'")
        if e - s != 1 or pos < 0 or not gene or not pseudo or set(gene + pseudo) - set(BASES) \
                or set(gene) & set(pseudo) or len(bits) < 3:
            raise CatalogError(f"{row['locus']}: PSV {name} must be one base with distinct gene and pseudogene "
                               f"bases and a text (got {s}-{e}, {label!r})")
        out.append((name, s, pos, set(gene), set(pseudo), bits[2]))
    return out


def parse_row(row):
    """The catalog row checked, with its parsed fields under '_' keys."""
    loc = row["locus"]
    if row["kind"] not in KINDS:
        raise CatalogError(f"{loc}: kind {row['kind']!r} (one of {', '.join(KINDS)})")
    if row["input"] not in INPUTS:
        raise CatalogError(f"{loc}: input {row['input']!r} (one of {', '.join(INPUTS)})")
    if row["chrom"] not in CHROMS:
        raise CatalogError(f"{loc}: chrom {row['chrom']!r} (GRCh38 names, chr1 to chrM)")
    p = dict(row)
    p["_targets"], p["_controls"], p["_details"] = (_segments(row, c) for c in ("targets", "controls", "details"))
    p["_call"] = _call(row)
    p["_readings"] = _readings(row)
    p["_psv"] = _psv_sites(row, p["_details"]) if p["_call"][0] == "psv" else []
    if p["_call"][0] == "psv" and (not p["_psv"] or row["target_depth"] != "sum"
                                   or not set(PSV_FINDINGS) <= set(p["_readings"])):
        raise CatalogError(f"{loc}: call psv needs PSV sites in details, target_depth sum (gene and pseudogene) "
                           f"and a reading for each of {', '.join(PSV_FINDINGS)}")
    if not p["_targets"]:
        raise CatalogError(f"{loc}: no target segment")
    if row["input"] == "bedcov":
        if not p["_controls"]:
            raise CatalogError(f"{loc}: a bedcov row needs control segments")
        if row["target_depth"] not in ("sum", "pooled") or row["control_depth"] not in ("flanks", "pooled"):
            raise CatalogError(f"{loc}: target_depth sum|pooled and control_depth flanks|pooled")
        p["_scale"], p["_offset"] = _num(row, "scale"), _num(row, "offset", 0.0)
        if not p["_scale"] or p["_scale"] <= 0:
            raise CatalogError(f"{loc}: scale must be a positive number")
        p["_min_control"] = _num(row, "min_control_depth", 0.0)
        p["_max_ratio"] = _num(row, "max_control_ratio")
        p["_mapq0_min"], p["_mapq0_max"] = _num(row, "mapq0_min"), _num(row, "mapq0_max")
        p["_passes"] = row["bedcov"].split(",")
        if (p["_mapq0_min"] is not None or p["_mapq0_max"] is not None) and len(p["_passes"]) != 2:
            raise CatalogError(f"{loc}: a MAPQ 0 guard needs two bedcov passes (all reads, MAPQ >= 1)")
        if p["_call"][0] == "class" and row["interval"] != "cv":
            raise CatalogError(f"{loc}: call class needs interval cv")
        if row["interval"] not in ("cv", "."):
            raise CatalogError(f"{loc}: interval {row['interval']!r} (cv or '.')")
    names = [n for n, *_ in p["_targets"] + p["_controls"] + p["_details"]]
    if len(names) != len(set(names)):
        raise CatalogError(f"{loc}: a segment name is used twice")
    return p


def load_catalog(path=None):
    """{locus: parsed row} of the catalog, in file order."""
    path = path or CATALOG
    with open(path) as f:
        lines = [ln for ln in f if ln.strip() and not ln.startswith("#")]
    rows = list(csv.DictReader(lines, delimiter="\t", quoting=csv.QUOTE_NONE))
    if not rows or list(rows[0].keys()) != COLUMNS:
        raise CatalogError(f"{path}: the header must be {' '.join(COLUMNS)}")
    out, seen = {}, set()
    for r in rows:
        if None in r or None in r.values():
            raise CatalogError(f"{path}: row {r.get('locus')!r} has the wrong number of columns")
        if r["locus"] in out:
            raise CatalogError(f"{path}: locus {r['locus']} is listed twice")
        p = parse_row(r)
        for n, *_ in p["_targets"] + p["_controls"] + p["_details"]:
            if n in seen:
                raise CatalogError(f"{path}: segment name {n} is used by two loci")
            seen.add(n)
        out[r["locus"]] = p
    return out


def regions(row):
    """[(chrom, start, end, name)] of a row's segments, sorted by position."""
    segs = row["_controls"] + row["_targets"] + row["_details"]
    return sorted(((row["chrom"], s, e, n) for n, s, e, _ in segs), key=lambda r: (r[1], r[2]))


def read_bedcov(path, row):
    """{name: mean depth} from `samtools bedcov` output over the row's
    segments, which must all be there and nothing else; {} for an empty table."""
    means, seen = {}, set()
    with open(path, errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 5:
                raise ValueError(f"{path}: expected chrom, start, end, name, depth sum; got {line.strip()!r}")
            chrom, start, end, name = c[0], int(c[1]), int(c[2]), c[3]
            if end <= start:
                raise ValueError(f"{path}: region {name} {chrom}:{start}-{end} is empty or reversed")
            if (chrom, start, end, name) in seen:
                raise ValueError(f"{path}: region {name} {chrom}:{start}-{end} appears twice")
            seen.add((chrom, start, end, name))
            means[name] = float(c[4]) / (end - start)
    if seen and seen != set(regions(row)):
        raise ValueError(f"{path}: its regions differ from the ones `copy_number.py bed --locus {row['locus']}` prints")
    return means


def _combine(mode, segs, m):
    if mode == "sum":
        return sum(m[n] for n, *_ in segs)
    if mode == "flanks":
        return sum(m[n] for n, *_ in segs) / len(segs)
    return sum(m[n] * (e - s) for n, s, e, _ in segs) / sum(e - s for _, s, e, _ in segs)


def _cv(values):
    """Coefficient of variation (sample standard deviation over the mean)."""
    mu = sum(values) / len(values)
    if mu <= 0 or len(values) < 2:
        return 0.0
    return math.sqrt(sum((v - mu) ** 2 for v in values) / (len(values) - 1)) / mu


def cut_copies(value, cuts):
    """The number of cut points a value reaches (cut:0.5,1.5 gives 0, 1 or 2)."""
    return sum(value >= c for c in cuts)


def reading(row, key):
    r = row["_readings"]
    if str(key) in r:
        return r[str(key)]
    if isinstance(key, int) and "more" in r and key > max(int(k) for k in r if k.isdigit()):
        return r["more"].format(n=key)
    return ""


def ploidy(row, sex):
    """The copies of the row's chromosome a sample has: 2, or 1 on a male's
    chrX or chrY; None where it cannot be said (a chrY row on a female or on
    a sample of no known sex, a chrX row without a sex)."""
    if row["chrom"] not in ("chrX", "chrY"):
        return 2
    if sex == "male":
        return 1
    if sex == "female" and row["chrom"] == "chrX":
        return 2
    return None


def assess(row, m, q1=None, sex=None):
    """The locus's block for mean depths m (every read) and q1 (reads of
    MAPQ >= 1, or None for a row without a MAPQ 0 guard; mapq0_share is then
    None).

    status is 'called' or 'not assessed'; reason names the guard that
    stopped it: no_data, low_control_depth, controls_disagree,
    low_mapq0_share, high_mapq0_share, below_floor, between_counts,
    sex_unknown; status 'not applicable' (reason no_y) for a chrY row on a
    female. sex matters only for a row on chrX or chrY: there the controls
    sit on the same chromosome, so the scale follows the sample's count
    there (ploidy). expected is the count with no change."""
    guarded = q1 is not None
    if not guarded and (row["_mapq0_min"] is not None or row["_mapq0_max"] is not None):
        raise ValueError(f"{row['locus']} has a MAPQ 0 guard: the MAPQ >= 1 depths are needed")
    q1 = q1 if guarded else {}
    tg, ct = row["_targets"], row["_controls"]
    b = {"locus": row["locus"], "status": "not assessed", "reason": None, "estimate": None, "total": None,
         "copies": None, "ci95": None, "class": None, "target_depth": None, "control_depth": None,
         "depths": dict(m), "mapq0_share": None, "details": {}, "reading": "", "expected": None}
    # expected: the copies with no change, the ploidy times the reference
    # units (each target segment is one with target_depth sum: AZFc's 4 DAZ)
    units = len(row["_targets"]) if row.get("target_depth") == "sum" else 1
    pl = ploidy(row, sex)
    b["expected"] = None if pl is None else pl * units
    if row["chrom"] == "chrY" and sex == "female":
        b.update(status="not applicable", reason="no_y")
        return b
    if pl is None:
        b["reason"] = "sex_unknown"
        return b
    if not m or any(n not in m for n, *_ in tg + ct):
        b["reason"] = "no_data"
        return b
    for n, *_ in ct:
        rest = [c for c in ct if c[0] != n]
        if m[n] == 0 and rest and _combine(row["control_depth"], rest, m) >= row["_min_control"] > 0:
            raise BrokenControl(f"BROKEN_CONTROL: control segment {n} of {row['locus']} has no read while the "
                                f"other controls have {_combine(row['control_depth'], rest, m):.1f}x; check its "
                                "coordinates in bin/cn_loci.tsv (or the sample lacks that stretch)")
    T, C = _combine(row["target_depth"], tg, m), _combine(row["control_depth"], ct, m)
    tq = _combine(row["target_depth"], tg, {n: q1.get(n, 0.0) for n, *_ in tg})
    # The MAPQ 0 share is known only when the MAPQ >= 1 pass was given.
    share = (1 - tq / T if T > 0 else 0.0) if guarded else None
    b.update(target_depth=T, control_depth=C, mapq0_share=share)
    low = min(m[n] for n, *_ in ct) if row["control_depth"] == "flanks" else C
    if low < row["_min_control"] or C <= 0:
        b["reason"] = "low_control_depth"
        return b
    if row["control_depth"] == "flanks" and row["_max_ratio"] is not None and len(ct) > 1:
        vals = [m[n] for n, *_ in ct]
        if min(vals) <= 0 or max(vals) / min(vals) > row["_max_ratio"]:
            b["reason"] = "controls_disagree"
            return b
    total = row["_scale"] * pl / 2 * T / C
    est = total - row["_offset"]
    b.update(total=total, estimate=est,
             details={n: row["_scale"] * pl / 2 * m[n] / C for n, *_ in row["_details"] if n in m})
    # With no copy left there is no read to judge the aligner by (a row with
    # an offset has a floor of its own below, after this guard).
    no_copy = row["_offset"] == 0 and round(est) == 0
    if row["_mapq0_min"] is not None and not no_copy and share < row["_mapq0_min"]:
        b["reason"] = "low_mapq0_share"
        return b
    if row["_mapq0_max"] is not None and share > row["_mapq0_max"]:
        b["reason"] = "high_mapq0_share"
        return b
    if row["_offset"] > 0 and est <= 0:
        b["reason"] = "below_floor"
        return b
    mode, vals = row["_call"]
    if row["interval"] == "cv":
        rel = math.sqrt(_cv([m[n] for n, *_ in tg]) ** 2 / len(tg) + _cv([m[n] for n, *_ in ct]) ** 2 / len(ct))
        b["ci95"] = (max(0.0, total * math.exp(-Z95 * rel) - row["_offset"]),
                     total * math.exp(Z95 * rel) - row["_offset"])
    if mode == "psv":
        # the module total (gene plus pseudogene copies); psv_assess reads the sites
        n = round(est)
        if abs(est - n) > vals[0]:
            b["reason"] = "between_counts"
            return b
        b["module_total"] = n
    elif mode == "round":
        n = round(est)
        if abs(est - n) > vals[0]:
            b["reason"] = "between_counts"
            return b
        b["copies"] = n
        b["reading"] = reading(row, n)
    elif mode == "cut":
        b["copies"] = cut_copies(est, vals)
        b["reading"] = reading(row, b["copies"])
    else:
        lo, hi = b["ci95"]
        b["class"] = "low" if hi < vals[0] else "high" if lo > vals[1] else "typical"
        b["reading"] = reading(row, b["class"])
    b["status"] = "called"
    return b


def count_bases(bases, ref="N"):
    """{A, C, G, T: reads} of one `samtools mpileup` bases column: read
    starts (^ and its MAPQ character), read ends, indels (+N and -N with
    their N bases), deletions and skips are not bases; '.' and ',' are the
    reference base when there is one."""
    n = dict.fromkeys(BASES, 0)
    i = 0
    while i < len(bases):
        ch = bases[i]
        if ch == "^":
            i += 2
            continue
        if ch in "+-":
            j = i + 1
            while j < len(bases) and bases[j].isdigit():
                j += 1
            i = j + int(bases[i + 1:j] or 0)
            continue
        b = ref.upper() if ch in ".," else ch.upper()
        if b in n:
            n[b] += 1
        i += 1
    return n


def read_pileup(path, row):
    """{site: (gene-base reads, pseudogene-base reads, other reads)} of a psv
    row from `samtools mpileup` over `bed --psv`, each site's two positions
    pooled: a read that carries the gene's base counts for the gene wherever
    the aligner put it. A position missing from the pileup has no read."""
    pos = {}
    for name, g, ps, gb, pb, _ in row["_psv"]:
        pos[g + 1] = pos[ps + 1] = name
    counts = {name: dict.fromkeys(BASES, 0) for name, *_ in row["_psv"]}
    with open(path, errors="replace") as f:
        for line in f:
            if not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 5 or not c[1].isdigit():
                raise ValueError(f"{path}: expected samtools mpileup lines (chrom, pos, ref, depth, bases, ...); "
                                 f"got {line.strip()[:80]!r}")
            if c[0] != row["chrom"] or int(c[1]) not in pos:
                continue
            for k, v in count_bases(c[4], c[2]).items():
                counts[pos[int(c[1])]][k] += v
    out = {}
    for name, g, ps, gb, pb, _ in row["_psv"]:
        n = counts[name]
        gene, pseudo = sum(n[x] for x in gb), sum(n[x] for x in pb)
        out[name] = (gene, pseudo, sum(n.values()) - gene - pseudo)
    return out


def psv_assess(row, b, counts):
    """A psv row's block (from assess) with its sites read. With the module
    total called, each site of PSV_MIN_READS reads or more gets the copies
    carrying the gene's base there (total x gene reads / gene and pseudogene
    reads) and its ratio to the median of the sites. cyp21a2 copies is the
    least of the sites. finding: deletion when the module total and the gene
    copies are both under the usual (the 30 kb chimeric deletion leaves one
    module and one gene copy fewer); conversion at the sites whose estimate
    lies PSV_MARGIN or more under the others' count or under the usual; pseudogene_deletion when the module
    total is under the usual with the gene kept everywhere; else none. Which
    chromosome carries a change is never said."""
    sites = {}
    for name, *_ in row["_psv"]:
        g, p, o = counts.get(name, (0, 0, 0))
        sites[name] = {"gene_reads": g, "pseudo_reads": p, "other_reads": o, "copies": None,
                       "estimate": None, "ratio": None}
    b.update(sites=sites, finding=None, conversion_sites=[])
    if b["status"] != "called" or b.get("module_total") is None:
        b["status"] = "not assessed"
        return b
    usual = b["expected"] // len(row["_targets"])
    read = []
    for name, s in sites.items():
        n = s["gene_reads"] + s["pseudo_reads"]
        if n >= PSV_MIN_READS:
            s["estimate"] = b["total"] * s["gene_reads"] / n
            s["copies"] = round(s["estimate"])
            read.append(s)
    if not read:
        b.update(status="not assessed", reason="psv_low_depth", copies=None)
        return b
    est = sorted(s["estimate"] for s in read)
    mid = len(est) // 2
    median = est[mid] if len(est) % 2 else (est[mid - 1] + est[mid]) / 2
    for s in read:
        s["ratio"] = s["estimate"] / median if median > 0 else None
    b["copies"] = least = min(s["copies"] for s in read)
    floor = max(round(median), usual)
    conv = [name for name, s in sites.items() if s["estimate"] is not None and s["estimate"] < floor - PSV_MARGIN]
    if least < usual and b["module_total"] < b["expected"]:
        b["finding"] = "deletion"
    elif conv:
        b.update(finding="conversion", conversion_sites=conv)
    elif b["module_total"] < b["expected"]:
        b["finding"] = "pseudogene_deletion"
    else:
        b["finding"] = "none"
    text = {name: t for name, *_, t in row["_psv"]}
    b["reading"] = reading(row, b["finding"]).replace(
        "{sites}", ", ".join(f"{n} ({text[n]})" for n in conv))
    return b


def psv_table(row, b):
    """([(metric, value)], one line) of a psv row's step table."""
    f = (lambda x: "." if x is None else f"{x:.2f}")
    d = (lambda x: "." if x is None else str(x))
    rows = [("status", b["status"]), ("cyp21a2_copies", d(b["copies"])), ("module_total", d(b.get("module_total"))),
            ("expected_total", d(b["expected"])), ("estimate", f(b["estimate"])), ("finding", b.get("finding") or "."),
            ("conversion_sites", ",".join(b.get("conversion_sites") or []) or ".")]
    for name, s in (b.get("sites") or {}).items():
        rows.append((f"psv_{name}", f"{s['gene_reads']} gene, {s['pseudo_reads']} pseudogene, {s['other_reads']} other"
                                    f" reads; copies {d(s['copies'])} ({f(s['estimate'])}); ratio {f(s['ratio'])}"))
    rows += [("target_depth", f(b["target_depth"])), ("control_depth", f(b["control_depth"])),
             ("reason", b["reason"] or "."), ("reading", b["reading"] or ".")]
    line = (f"{row['locus']}: {b['status']} ({b['reason'] or b.get('finding')}; {d(b['copies'])} gene copies, "
            f"module total {d(b.get('module_total'))} of {d(b['expected'])})")
    return rows, line


def interval_rows(cat):
    """The rows `intervals` calls: kind interval with one bedcov pass (LCE3BC,
    an interval row with a MAPQ >= 1 pass, stays a gene-dosage row)."""
    return [r for r in cat.values() if r["kind"] == "interval" and r["input"] == "bedcov" and r["_passes"] == ["."]]


def _bed(args):
    cat = load_catalog(args.catalog)
    if args.psv:
        if not args.locus or args.kind:
            raise CatalogError("bed --psv needs --locus (and no --kind)")
        for loc in args.locus:
            if loc not in cat or not cat[loc]["_psv"]:
                raise CatalogError(f"no psv locus {loc!r} in the catalog")
            r = cat[loc]
            pos = [(s, f"{n}") for n, s, *_ in r["_psv"]] + [(p, f"{n}_pseudo") for n, _, p, *_ in r["_psv"]]
            for s, n in sorted(pos):
                print(f"{r['chrom']}\t{s}\t{s + 1}\t{n}")
        return 0
    if args.kind == "interval":
        if args.locus:
            raise CatalogError("bed: --kind interval or --locus, not both")
        loci = [r["locus"] for r in interval_rows(cat)]
    else:
        loci = args.locus or [k for k, r in cat.items() if r["input"] == "bedcov"]
    for loc in loci:
        if loc not in cat or cat[loc]["input"] != "bedcov":
            raise CatalogError(f"no bedcov locus {loc!r} in the catalog")
        for r in regions(cat[loc]):
            print("\t".join(map(str, r)))
    return 0


def _call_cmd(args):
    try:
        if args.bam_header:
            check_bam(args.bam_header)
    except AltContigs as e:
        print(e, file=sys.stderr)
        return 4
    except OSError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    try:
        cat = load_catalog(args.catalog)
        if args.locus not in cat or cat[args.locus]["input"] != "bedcov":
            raise CatalogError(f"no bedcov locus {args.locus!r} in the catalog")
        row = cat[args.locus]
        m = read_bedcov(args.bedcov, row)
        guarded = row["_mapq0_min"] is not None or row["_mapq0_max"] is not None
        if guarded and not args.bedcov_mapq1:
            raise ValueError(f"{args.locus} has a MAPQ 0 guard: --bedcov-mapq1 is needed")
        q1 = read_bedcov(args.bedcov_mapq1, row) if args.bedcov_mapq1 else None
        if q1 is not None and bool(m) != bool(q1):
            raise ValueError(f"{args.bedcov_mapq1 if m else args.bedcov} is empty but the other table is not: "
                             "both must cover the regions, or both be empty (no reads on that contig)")
        block = assess(row, m, q1)
        if row["_psv"]:
            if not args.pileup:
                raise ValueError(f"{args.locus} is a psv row: --pileup (samtools mpileup over `bed --psv`) is needed")
            block = psv_assess(row, block, read_pileup(args.pileup, row) if m else {})
    except BrokenControl as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 3
    except (OSError, ValueError, KeyError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    rows, line = table(row, block)
    try:
        if args.json:
            with open(args.json, "w") as f:
                json.dump(block, f, indent=1, sort_keys=True)
                f.write("\n")
        if args.table:
            with open(args.table, "w") as f:
                f.write("metric\tvalue\n")
                for k, v in rows:
                    f.write(f"{k}\t{v}\n")
    except OSError as e:
        print(f"ERROR: cannot write: {e}", file=sys.stderr)
        return 2
    print(line)
    return 0


def read_bedcov_all(path):
    """[(chrom, start, end, name, depth sum)] of a `samtools bedcov` table."""
    out = []
    with open(path, errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 5:
                raise ValueError(f"{path}: expected chrom, start, end, name, depth sum; got {line.strip()!r}")
            out.append((c[0], int(c[1]), int(c[2]), c[3], float(c[4])))
    return out


def interval_line(row, block):
    """One line of the interval table for a row and its block."""
    tg = row["_targets"]
    f = (lambda x: "." if x is None else f"{x:.2f}")
    copies, exp = block["copies"], block["expected"]
    change = "." if copies is None or exp is None else \
        "deletion" if copies < exp else "duplication" if copies > exp else "no change"
    return [row["locus"], tg[0][3] or row["locus"], row["chrom"], str(min(s for _, s, _, _ in tg)),
            str(max(e for _, _, e, _ in tg)), block["status"], "." if copies is None else str(copies),
            "." if exp is None else str(exp), change, f(block["estimate"]), block["reason"] or "."]


def call_intervals(cat, rows_in, sex):
    """[line] of every interval row from one bedcov table (rows_in, as
    read_bedcov_all gives it); a broken control makes that row 'not
    assessed' (broken_control), never the whole table."""
    by_name = {}
    for c, s, e, n, d in rows_in:
        if n in by_name:
            raise ValueError(f"segment {n} appears twice in the bedcov table")
        by_name[n] = (c, s, e, d)
    out = []
    for row in interval_rows(cat):
        want = regions(row)
        have = [w for w in want if w[3] in by_name]
        if len(have) != len(want):
            raise ValueError(f"the bedcov table lacks some of {row['locus']}'s segments: "
                             f"its regions differ from `copy_number.py bed --kind interval`")
        m = {}
        for c, s, e, n in have:
            if by_name[n][:3] != (c, s, e):
                raise ValueError(f"segment {n} of {row['locus']} has other coordinates in the bedcov table")
            m[n] = by_name[n][3] / (e - s)
        try:
            block = assess(row, m, None, sex)
        except BrokenControl as e:
            print(f"WARNING: {e}", file=sys.stderr)
            block = assess(row, {}, None, sex)
            if block["status"] == "not assessed":
                block["reason"] = "broken_control"
        out.append(interval_line(row, block))
    return out


def _intervals(args):
    try:
        cat = load_catalog(args.catalog)
        if args.failed:
            lines = [interval_line(r, {"status": "not assessed", "copies": None, "estimate": None,
                                       "reason": "check_failed", "expected": None}) for r in interval_rows(cat)]
        else:
            if not args.bedcov or not args.sex:
                raise ValueError("intervals needs --bedcov and --sex (or --failed)")
            lines = call_intervals(cat, read_bedcov_all(args.bedcov), args.sex)
    except (OSError, ValueError, KeyError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    try:
        with open(args.out + ".tmp", "w") as f:
            f.write("\t".join(INTERVAL_COLUMNS) + "\n")
            for ln in lines:
                f.write("\t".join(ln) + "\n")
        os.replace(args.out + ".tmp", args.out)
    except OSError as e:
        print(f"ERROR: cannot write: {e}", file=sys.stderr)
        return 2
    n = {}
    for ln in lines:
        n[ln[5]] = n.get(ln[5], 0) + 1
    print(f"interval rows: {len(lines)} ({', '.join(f'{v} {k}' for k, v in sorted(n.items()))})")
    return 0


def table(row, block):
    """([(metric, value)], one summary line) of a locus's step table."""
    if row["_psv"]:
        return psv_table(row, block)
    mod = TABLES.get(row["locus"])
    if mod:
        return importlib.import_module(mod).table(block)
    f = (lambda x: "." if x is None else f"{x:.2f}")
    rows = [("status", block["status"]), ("copies", "." if block["copies"] is None else block["copies"]),
            ("estimate", f(block["estimate"])), ("target_depth", f(block["target_depth"])),
            ("control_depth", f(block["control_depth"])), ("reason", block["reason"] or "."),
            ("reading", block["reading"] or ".")]
    return rows, f"{row['locus']}: {block['status']} ({block['reason'] or block['reading']})"


def main(argv=None):
    ap = argparse.ArgumentParser(prog="copy_number.py", description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("bed", help="print the segments of the loci")
    b.add_argument("--catalog")
    b.add_argument("--locus", action="append")
    b.add_argument("--kind", choices=["interval"], help="every interval row (as `intervals` calls them)")
    b.add_argument("--psv", action="store_true", help="the PSV sites of a psv row (gene and pseudogene positions)")
    c = sub.add_parser("call", help="the copy number of one locus from samtools bedcov over its segments")
    c.add_argument("--catalog")
    c.add_argument("--locus", required=True)
    c.add_argument("--bedcov", required=True)
    c.add_argument("--bedcov-mapq1", help="samtools bedcov with the row's MAPQ >= 1 pass")
    c.add_argument("--json")
    c.add_argument("--table")
    c.add_argument("--bam-header", help="refuse a BAM with ALT or HLA contigs first (as check-bam)")
    c.add_argument("--pileup", help="samtools mpileup over `bed --psv` (a psv row)")
    iv = sub.add_parser("intervals", help="call every interval row from one samtools bedcov table")
    iv.add_argument("--catalog")
    iv.add_argument("--bedcov", help="samtools bedcov over `copy_number.py bed --kind interval`")
    iv.add_argument("--sex", choices=SEXES)
    iv.add_argument("--failed", action="store_true", help="every row not assessed (check_failed)")
    iv.add_argument("--out", required=True)
    k = sub.add_parser("check-bam", help="refuse a BAM aligned to a reference with ALT or HLA contigs")
    k.add_argument("--header", required=True, help="samtools view -H of the BAM, or contig names; - for stdin")
    args = ap.parse_args(argv)
    if args.cmd == "check-bam":
        try:
            check_bam(args.header)
        except AltContigs as e:
            print(e, file=sys.stderr)
            return 4
        except OSError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 2
        return 0
    if args.cmd == "intervals":
        return _intervals(args)
    if args.cmd == "bed":
        try:
            return _bed(args)
        except (OSError, CatalogError) as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 2
    return _call_cmd(args)


if __name__ == "__main__":
    sys.exit(main())
