#!/usr/bin/env python3
"""build-cn-intervals.py: write the interval rows of bin/cn_loci.tsv.

The interval rows are the copy-number checks step 35 runs with one
`samtools bedcov` pass (bin/copy_number.py intervals):

  - every region row of ClinGen Dosage Sensitivity's recurrent CNV BED,
    version 2.1, GRCh38 (CLINGEN_RECURRENT_CNV_URL and _SHA256 in
    versions.env, CC0): 22q11.2, 1q21.1, 15q11.2, 16p11.2, 17p12, 17p11.2,
    7q11.23, 3q29, Xp22.31 and the rest. The file's breakpoint and LCR rows
    (itemRgb 0,0,0) are the segmental duplications that bound the regions:
    they are cut out of the targets and kept out of the controls. Its chrX
    rows are spelled 'chrx' and are normalised to chrX; any other name that is
    not a GRCh38 chromosome stops the build;
  - two ClinGen curations the BED lacks (HAND_ROWS): the 1p36 terminal
    region (ISCA-37434, haploinsufficiency 3) and the Xq28 region that holds
    MECP2 (ISCA-46304, triplosensitivity 3), from
    ClinGen_region_curation_list_GRCh38.tsv as read on 2026-11-04;
  - three carrier deletion alleles at their ClinVar GRCh38 spans (HAND_ROWS):
    del(GJB6-D13S1830) (VCV 5546), CFTRdele2,3 (VCV 66105) and a 621 bp HBB
    deletion that removes exon 3 (VCV 3769040); and two beta-globin
    deletions at their HbVar GRCh38 spans, checked on 2026-11-04 by placing
    50 bp on each side of the NG_000007.3 breakpoints on NC_000011.10: the
    619 bp HBB deletion common in India and Pakistan (HbVar 979,
    NG_000007.3:g.71609_72227del619, Orkin 1979) and the 13.4 kb Sicilian
    (delta-beta)0 deletion from HBD to 3' of HBB (HbVar 1035,
    NG_000007.3:g.64336_77738del13403, Henthorn 1990). The alpha-globin
    --SEA deletion is read by the HBA row (step 35's HBA table), not here;
  - the Y chromosome AZF regions (HAND_ROWS), read only on a sample with a Y:
    AZFa as USP9Y and DDX3Y, AZFb as KDM5D and EIF1AY (single-copy genes each
    deletion removes), and AZFc as the four DAZ copies added up (a complete
    AZFc deletion leaves none, the partial gr/gr and b2/b3 deletions two;
    which copies are lost, and the breakpoints, are not claimed). Gene spans
    from Ensembl's REST lookup, checked on 2026-11-04. Controls: ZFY and
    NLGN4Y, single-copy genes outside every AZF region.

Targets are cut to pieces of at most 500 kb and lose every stretch inside a
breakpoint or LCR row of the BED, a GIAB v3.6 segmental duplication, an
assembly gap (GIAB's gaps with 15 kb on each side), a satellite or a region of
a model the pinned Parascopy data holds (assets/parascopy_models_v1.7.txt:
Parascopy owns those, CHRNA7 in 15q13.3 among them); the AZF targets keep
their segmental duplications, since the DAZ copies are one. Controls are two
100 kb windows, one on each side of the cluster of ClinGen rows the region
belongs to, the nearest that stay clear of all of that, within 3 Mb; a side
with none (a chromosome end) is left out. Every row: depth pooled over the
targets against depth pooled over the controls, scale 2, a whole number
within 0.3 (round:0.3), at least 8x of control depth (4x on chrY), the two
flanks within 1.5x of each other; cross_check cnvpytor:cnvs (step 18).

Usage:
  scripts/ci/build-cn-intervals.py [--catalog bin/cn_loci.tsv] [--cache DIR] [--check]
  scripts/ci/build-cn-intervals.py --list   the downloads, one "NAME URL SHA256" per line

Rewrites the catalog's rows from the marker line on (MARKER), or with --check
exits 1 when the catalog is not what the downloads give. The downloads (about
150 KB) are checked against their sha256 and kept in --cache when given. The
e2e case tests/e2e/cn-intervals-1-catalog.sh runs --check.
"""
import argparse
import gzip
import hashlib
import os
import re
import sys
import tempfile
import time
import urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MARKER = "# --- interval rows: written by scripts/ci/build-cn-intervals.py; do not edit by hand ---"
GIAB = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all"
MASKS = [
    ("segdups.bed.gz", f"{GIAB}/SegmentalDuplications/GRCh38_segdups.bed.gz",
     "09ac029a7d0d681a9f0ed6581a6bc10fc708a834bc3f91c46a9bef504f6055c5"),
    ("gaps.bed.gz", f"{GIAB}/OtherDifficult/GRCh38_gaps_slop15kb.bed.gz",
     "b5fdbfc7e5b9cf4fa356f23841fd7b8a74fa239fd8cbc5b5c3463c9e7bdc0902"),
    ("satellites.bed.gz", f"{GIAB}/LowComplexity/GRCh38_satellites_slop5.bed.gz",
     "10207aa2dac0e1189ea553a374abbdcb4de25a81378d88dc4ba03aea8718b5a2"),
]
REGION_RGB, MASK_RGB = "255,111,0", "0,0,0"
CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
# GRCh38 chromosome lengths (the no-ALT analysis set's .fai).
LENGTHS = {"chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555, "chr5": 181538259,
           "chr6": 170805979, "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, "chr10": 133797422,
           "chr11": 135086622, "chr12": 133275309, "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
           "chr16": 90338345, "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
           "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415}
PIECE, FLANK, GAP, SEARCH, STEP, MIN_PIECE = 500_000, 100_000, 50_000, 3_000_000, 10_000, 500

# Rows the BED lacks: locus, chrom, [(name, 0-based start, end)] targets, [controls] or None (flanks),
# target_depth, whether segmental duplications are cut out of the targets, pmids. The label of the
# first target is the row's name in step 35's table.
Y_CONTROLS = [("ZFY", 2935200, 2982579, "ZFY, single-copy control on Yp"),
              ("NLGN4Y", 14522572, 14845654, "NLGN4Y, single-copy control on Yq")]
HAND_ROWS = [
    ("1p36_terminal", "chr1", [("1p36 terminal region (includes GABRD), ClinGen ISCA-37434", 898702, 6229913)],
     None, "pooled", True, "."),
    ("Xq28_MECP2", "chrX", [("Xq28 region (includes MECP2), ClinGen ISCA-46304: the MECP2 duplication",
                             154008528, 154110279)], None, "pooled", True, "."),
    ("GJB6_D13S1830", "chr13", [("del(GJB6-D13S1830), ClinVar 5546", 20223037, 20531806)], None, "pooled", True,
     "11807148"),
    ("CFTR_dele2_3", "chr7", [("CFTRdele2,3 (21 kb), ClinVar 66105", 117498312, 117519392)], None, "pooled", True,
     "."),
    ("HBB_exon3_del", "chr11", [("HBB deletion of exon 3 (621 bp), ClinVar 3769040", 5225254, 5225875)], None,
     "pooled", True, "."),
    ("HBB_619bp_del", "chr11", [("HBB 619 bp deletion (South Asian), HbVar 979", 5225388, 5226007)], None,
     "pooled", True, "287080"),
    ("HBD_HBB_Sicilian_del", "chr11", [("Sicilian (delta-beta)0 deletion (13.4 kb, HBD to HBB 3'), HbVar 1035",
                                        5219877, 5233280)], None, "pooled", True, "2307466"),
    ("AZFa", "chrY", [("AZFa (USP9Y and DDX3Y)", 12537649, 12860852), ("DDX3Y", 12904107, 12920478)],
     Y_CONTROLS, "pooled", False, "37674303"),
    ("AZFb", "chrY", [("AZFb (KDM5D and EIF1AY)", 19702297, 19745148), ("EIF1AY", 20575693, 20593154)],
     Y_CONTROLS, "pooled", False, "37674303"),
    ("AZFc_DAZ", "chrY", [("AZFc (DAZ1 to DAZ4 added up; gr/gr and b2/b3 leave two)", 23129354, 23199010),
                          ("DAZ2", 23219446, 23291356), ("DAZ3", 24763068, 24813492), ("DAZ4", 24833842, 24907040)],
     Y_CONTROLS, "sum", False, "37674303"),
]
COLUMNS = ["locus", "kind", "chrom", "targets", "controls", "details", "target_depth", "control_depth",
           "scale", "offset", "bedcov", "mapq0_min", "mapq0_max", "min_control_depth", "max_control_ratio",
           "call", "interval", "readings", "input", "owner", "cross_check", "pmids"]


class BuildError(ValueError):
    pass


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fetch(name, url, want, folder):
    path = os.path.join(folder, name)
    if not (os.path.exists(path) and sha256(path) == want):
        for i in range(3):
            try:
                with urllib.request.urlopen(url, timeout=120) as r, open(path + ".part", "wb") as f:
                    f.write(r.read())
                os.replace(path + ".part", path)
                break
            except OSError as e:
                # A missing or forbidden file stays so: no point waiting.
                if i == 2 or getattr(e, "code", None) in (403, 404):
                    raise BuildError(f"cannot download {url}: {e}")
                time.sleep(30)
    got = sha256(path)
    if got != want:
        raise BuildError(f"{url}: sha256 {got}, pinned {want}")
    return path


def versions_env(key):
    with open(os.path.join(REPO, "versions.env")) as f:
        for line in f:
            m = re.match(rf'^{key}="([^"]*)"', line)
            if m:
                return m.group(1)
    raise BuildError(f"no {key} in versions.env")


def norm_chrom(c):
    """chrx -> chrX (the ClinGen BED's spelling); anything not a GRCh38 chromosome stops the build."""
    c2 = re.sub(r"^chr([xy])$", lambda m: "chr" + m.group(1).upper(), c)
    if c2 not in CHROMS:
        raise BuildError(f"chromosome {c!r} is not a GRCh38 chromosome name")
    return c2


def read_clingen(lines):
    """(regions [(chrom, start, end, name)], masks [(chrom, start, end)]) of the BED's lines."""
    regions, masks = [], []
    for line in lines:
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9:
            raise BuildError(f"ClinGen BED: a row with {len(c)} columns: {line.strip()!r}")
        chrom, s, e, name, rgb = norm_chrom(c[0]), int(c[1]), int(c[2]), c[3].strip(), c[8].strip()
        if e <= s:
            raise BuildError(f"ClinGen BED: {name}: empty span")
        if rgb == REGION_RGB:
            regions.append((chrom, s, e, name))
        elif rgb == MASK_RGB:
            masks.append((chrom, s, e))
        else:
            raise BuildError(f"ClinGen BED: {name}: itemRgb {rgb!r} is neither a region nor a breakpoint row")
    if not regions:
        raise BuildError("ClinGen BED: no region row")
    return regions, masks


def read_bed(path):
    out = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.split("\t")
            out.append((c[0], int(c[1]), int(c[2])))
    return out


def read_parascopy(path):
    out = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            _, reg, _ = line.split()
            c, span = reg.split(":")
            s, e = span.split("-")
            out.append((c, int(s) - 1, int(e)))
    return out


def merge(spans):
    out = []
    for s, e in sorted(spans):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def subtract(s, e, cuts):
    """[(s, e)] of s-e outside the merged spans cuts."""
    out, pos = [], s
    for cs, ce in cuts:
        if ce <= pos or cs >= e:
            continue
        if cs > pos:
            out.append((pos, cs))
        pos = max(pos, ce)
    if pos < e:
        out.append((pos, e))
    return out


def pieces(spans):
    out = []
    for s, e in spans:
        if e - s < MIN_PIECE:
            continue
        n = -(-(e - s) // PIECE)
        step = -(-(e - s) // n)
        for i in range(n):
            out.append((s + i * step, min(e, s + (i + 1) * step)))
    return out


def clear(s, e, cuts):
    return all(ce <= s or cs >= e for cs, ce in cuts)


def flank(chrom, edge, side, cuts):
    """The nearest FLANK window GAP beyond edge on side (-1 left, +1 right), clear of cuts, within SEARCH."""
    for k in range(0, SEARCH // STEP + 1):
        if side < 0:
            e = edge - GAP - k * STEP
            s = e - FLANK
        else:
            s = edge + GAP + k * STEP
            e = s + FLANK
        if s < 0 or e > LENGTHS[chrom]:
            return None
        if clear(s, e, cuts):
            return (s, e)
    return None


def locus_id(name):
    band = name.split()[0]
    groups = [g for g in re.findall(r"\(([^()]*)\)", name) if not g.strip().startswith("includes")]
    words = re.sub(r"[^A-Za-z0-9.-]+", "_", " ".join(groups)).strip("_")
    return f"{band}_{words}" if words else band


def label(text):
    return re.sub(r"[;\t\n]", ",", text).strip()


def build(clingen_lines, segdups, gaps, sats, parascopy):
    regions, masks = read_clingen(clingen_lines)
    by_chrom = {}
    for name, items in (("mask", masks), ("segdup", segdups), ("gap", gaps), ("sat", sats), ("para", parascopy)):
        for c, s, e in items:
            by_chrom.setdefault(c, {}).setdefault(name, []).append((s, e))
    # Clusters: ClinGen rows (regions and breakpoints) that touch or overlap, per chromosome.
    clusters = {}
    for c in CHROMS:
        spans = [(s, e + 1) for cc, s, e, _ in regions if cc == c] + [(s, e + 1) for cc, s, e in masks if cc == c]
        clusters[c] = merge(spans)

    def cuts(c, with_segdups=True):
        d = by_chrom.get(c, {})
        keys = ("mask", "segdup", "gap", "sat", "para") if with_segdups else ("mask", "gap", "para")
        return merge([x for k in keys for x in d.get(k, [])])

    def flank_cuts(c):
        # A control never sits on a target: the ClinGen regions and the HAND_ROWS targets are cut out too.
        d = by_chrom.get(c, {})
        return merge([x for k in ("mask", "segdup", "gap", "sat", "para") for x in d.get(k, [])]
                     + [(s, e) for cc, s, e, _ in regions if cc == c]
                     + [(s, e) for _, cc, tl, *_ in HAND_ROWS if cc == c for _, s, e in tl])

    rows, seen = [], {}

    def add(locus, chrom, targets, controls, target_depth, pmids):
        if not targets:
            raise BuildError(f"{locus}: no target stretch is left outside the masks")
        if locus in seen:
            seen[locus] += 1
            locus = f"{locus}_{seen[locus]}"
        else:
            seen[locus] = 1
        y = chrom == "chrY"
        tg = ";".join(f"{locus}_t{i + 1}:{s}-{e}" + (f":{lab}" if lab else "") for i, (s, e, lab) in enumerate(targets))
        ct = ";".join(f"{locus}_{n}:{s}-{e}:{lab}" for n, s, e, lab in controls)
        rows.append({"locus": locus, "kind": "interval", "chrom": chrom, "targets": tg, "controls": ct,
                     "details": ".", "target_depth": target_depth, "control_depth": "pooled" if y else "flanks",
                     "scale": "2", "offset": "0", "bedcov": ".", "mapq0_min": ".", "mapq0_max": ".",
                     "min_control_depth": "4" if y else "8", "max_control_ratio": "." if y else "1.5",
                     "call": "round:0.3", "interval": ".", "readings": ".", "input": "bedcov", "owner": "cn-depth",
                     "cross_check": "cnvpytor:cnvs", "pmids": pmids})

    def flanks(chrom, s, e):
        cs, ce = next(((a, b) for a, b in clusters[chrom] if a <= s and e <= b), (s, e))
        fc = flank_cuts(chrom)
        out = []
        for side, edge, n in ((-1, cs, "L"), (1, ce, "R")):
            w = flank(chrom, edge, side, fc)
            if w:
                out.append((n, w[0], w[1], "flank " + ("left" if side < 0 else "right") + " of the region"))
        if not out:
            raise BuildError(f"{chrom}:{s}-{e}: no flank clear of the masks within {SEARCH // 1_000_000} Mb")
        return out

    for chrom, s, e, name in sorted(regions, key=lambda r: (sorted(CHROMS, key=_order).index(r[0]), r[1], r[2])):
        segs = pieces(subtract(s, e, cuts(chrom)))
        targets = [(a, b, label(name) if i == 0 else "") for i, (a, b) in enumerate(segs)]
        add(locus_id(name), chrom, targets, flanks(chrom, s, e), "pooled", ".")
    for locus, chrom, tlist, controls, td, with_sd, pmids in HAND_ROWS:
        c = cuts(chrom, with_sd)
        targets = []
        for name, s, e in tlist:
            segs = pieces(subtract(s, e, c)) if td == "pooled" else [(s, e)]
            if td == "sum" and not clear(s, e, merge(by_chrom.get(chrom, {}).get("para", []))):
                raise BuildError(f"{locus}: {name} overlaps a Parascopy model")
            targets += [(a, b, label(name) if i == 0 else "") for i, (a, b) in enumerate(segs)]
        if controls is None:
            ctl = flanks(chrom, tlist[0][1], tlist[-1][2])
        else:
            ctl = [(f"c{i + 1}", s, e, lab) for i, (_, s, e, lab) in enumerate(controls)]
        add(locus, chrom, targets, ctl, td, pmids)
    return rows


def _order(c):
    k = c[3:]
    return (0, int(k)) if k.isdigit() else (1, k)


def render(rows):
    return [MARKER + "\n"] + ["\t".join(r[c] for c in COLUMNS) + "\n" for r in rows]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--catalog", default=os.path.join(REPO, "bin", "cn_loci.tsv"))
    ap.add_argument("--cache")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args(argv)
    if a.list and a.check:
        ap.error("--list and --check do not go together")
    if a.list:
        try:
            files = [("clingen_recurrent_cnv.bed", versions_env("CLINGEN_RECURRENT_CNV_URL"),
                      versions_env("CLINGEN_RECURRENT_CNV_SHA256"))] + MASKS
        except BuildError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 2
        for n, u, h in files:
            print(n, u, h)
        return 0
    try:
        with tempfile.TemporaryDirectory() as tmp:
            folder = a.cache or tmp
            os.makedirs(folder, exist_ok=True)
            bed = fetch("clingen_recurrent_cnv.bed", versions_env("CLINGEN_RECURRENT_CNV_URL"),
                        versions_env("CLINGEN_RECURRENT_CNV_SHA256"), folder)
            sd, gp, st = (read_bed(fetch(n, u, h, folder)) for n, u, h in MASKS)
            with open(bed, errors="replace") as f:
                rows = build(f.readlines(), sd, gp, st,
                             read_parascopy(os.path.join(REPO, "assets", "parascopy_models_v1.7.txt")))
    except (OSError, BuildError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    with open(a.catalog) as f:
        lines = f.readlines()
    keep = lines[:lines.index(MARKER + "\n")] if MARKER + "\n" in lines else lines
    new = keep + render(rows)
    if a.check:
        if new != lines:
            old = set(lines) - set(new)
            print(f"{a.catalog} is not what the pinned downloads give ({len(old)} line(s) differ); "
                  "rerun scripts/ci/build-cn-intervals.py", file=sys.stderr)
            return 1
        print(f"{a.catalog}: the {len(rows)} interval rows are what the pinned downloads give")
        return 0
    with open(a.catalog + ".tmp", "w") as f:
        f.writelines(new)
    os.replace(a.catalog + ".tmp", a.catalog)
    print(f"{a.catalog}: {len(rows)} interval rows written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
