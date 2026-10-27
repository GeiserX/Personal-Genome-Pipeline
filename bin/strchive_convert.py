#!/usr/bin/env python3
"""strchive_convert.py: one repeat catalog for ExpansionHunter and Stranger from STRchive.

STRchive (Hiatt et al. 2025, PMID 40140942) publishes the disease repeats as
an ExpansionHunter catalog and, apart, a JSON of every locus with its
benign, intermediate and pathogenic sizes and its motifs. Stranger needs
both in one file: ExpansionHunter's fields plus NormalMax and
PathologicMin. ExpansionHunter reads only the fields it knows and Stranger
only its own, so one file serves both, and steps 09 and 9b cannot disagree
on the loci. This writes that file: one entry per ExpansionHunter locus.

  strchive_convert.py --eh STRchive-disease-loci-<v>.hg38.expansionhunter.json \\
                      --loci STRchive-loci-<v>.json --release <v> --out catalog.json \\
                      [--ranges str_ranges.tsv]

Each entry:
  LocusId        the gene in capitals (HTT, FMR1, C9ORF72, RFC1), the names
                 ExpansionHunter's bundled catalog uses, so the report finds
                 its key loci; STRchive's id where one gene has two loci
                 (EIEE1_ARX and PRTS_ARX, the three HOXA13 tracts)
  LocusStructure, ReferenceRegion, VariantType
                 STRchive's, unchanged
  VariantId      for a locus of several repeats: the pathogenic one (the
                 STRchive pathogenic_repeat) is named LocusId, which
                 ExpansionHunter writes as REPID and Stranger looks up; the
                 others LocusId_<motif>, which Stranger leaves without a status
  NormalMax, PathologicMin
                 Stranger's two thresholds: normal up to NormalMax,
                 pre_mutation above it and below PathologicMin, full_mutation
                 from PathologicMin. NormalMax is STRchive's benign_max and
                 PathologicMin its pathogenic_min. Where STRchive gives no
                 benign_max (the FAME loci, SCA31, XDP), NormalMax is the
                 reference allele's unit count, so a longer allele is never
                 called normal. A locus gets neither, and so no status, when
                 STRchive gives no pathogenic size, when its pathogenic size
                 is below its benign one (a contraction), or when the
                 reference repeat is longer than SHORT_READ_LIMIT bases: no
                 short read spans even the reference allele there, so
                 ExpansionHunter cannot size it. NoThresholds says why.
  PathologicRegion
                 the pathogenic repeat's region, for a locus of several
  Disease, InheritanceMode, Source, SourceId, SourceDisplay
                 STRchive's disease id, inheritance and locus id; Stranger
                 copies them into the VCF
  BenignMotifs, PathogenicMotifs
                 for a locus where the motif decides (RFC1, the FAME loci,
                 SCA31, SCA37): the motifs bin/rfc1_motif.py counts in the
                 reads. Benign: STRchive's benign and reference motifs;
                 pathogenic: its pathogenic ones. Only motifs of the
                 pathogenic motif's length, 2 to 6 bases. A locus is one of
                 these when its pathogenic motif is not a rotation of its
                 reference motif.

--ranges also writes the size ranges as a table, one row per entry: locus
(the LocusId), source_id (STRchive's id), normal_max and pathogenic_min
(the entry's NormalMax and PathologicMin, empty where it has none),
intermediate_min and intermediate_max (STRchive's own, as published, empty
where it gives none) and no_thresholds (the entry's NoThresholds).
tests/test_str_ranges.py checks it against the catalog.

One correction to STRchive's motif classes: at RFC1 STRchive lists AAAGG as
pathogenic, and its own notes say why: AAAGG expansions of 600 or more units
are seen in CANVAS, while 100 to 380 units are seen in unaffected people
(PMID 37450567). Short reads cannot size such an allele, so the read count
keeps AAAGG benign, as Cortese et al. 2019 (PMID 30926972) classed it.
STRchive's AAAGGG is six bases long and is not counted at a five-base locus.

Exits 2, writing nothing (neither file), when an ExpansionHunter locus has no STRchive
record, two loci get one name, or the pathogenic repeat of a locus cannot be
told.
"""
import argparse
import json
import os
import re
import sys

SHORT_READ_LIMIT = 300  # bases: twice a 150-base read
MIN_PERIOD, MAX_PERIOD = 2, 6  # a motif read in the reads is a short tandem repeat unit, not a VNTR's
SIZE_DEPENDENT = {"CANVAS_RFC1": {"AAAGG"}}  # listed pathogenic, benign at the sizes short reads see
REGION = re.compile(r"^(chr[^:]+):(\d+)-(\d+)$")
COMP = str.maketrans("ACGTN", "TGCAN")


class Bad(Exception):
    pass


def as_list(x):
    return x if isinstance(x, list) else [x]


def region(text):
    m = REGION.match(text)
    if not m:
        raise Bad(f"not a region: {text!r}")
    return m.group(1), int(m.group(2)), int(m.group(3))


def forms(m):
    """Every rotation of a motif on both strands."""
    out = set()
    for x in (m, m.translate(COMP)[::-1]):
        out.update(x[i:] + x[:i] for i in range(len(x)))
    return out


def same_motif(a, b):
    return a.upper() in forms(b.upper())


def motifs_of(structure):
    return re.findall(r"\(([^)]*)\)\*", structure)


def names(eh, loci):
    """STRchive id -> LocusId: the gene in capitals, or the STRchive id where one gene has two loci."""
    genes = {}
    for e in eh:
        genes.setdefault(loci[e["LocusId"]]["gene"].upper(), []).append(e["LocusId"])
    out = {}
    for gene, ids in genes.items():
        for i in ids:
            out[i] = gene if len(ids) == 1 else i
    return out


def pathogenic_index(e, rec):
    """Which of a locus's repeats is the pathogenic one."""
    regions = as_list(e["ReferenceRegion"])
    if len(regions) == 1:
        return 0
    mot = motifs_of(e["LocusStructure"])
    if len(mot) != len(regions):
        raise Bad(f"{e['LocusId']}: {len(mot)} motifs in {e['LocusStructure']} for {len(regions)} regions")
    want = [s["motif"] for s in rec.get("locus_structure") or [] if s.get("type") == "pathogenic_repeat"]
    hits = [i for i, m in enumerate(mot) if any(same_motif(m, w) for w in want)]

    def overlap(i):
        _, a, b = region(regions[i])
        return max(0, min(b, rec["stop_hg38"]) - max(a, rec["start_hg38"]))
    pool = hits or list(range(len(regions)))
    best = max(overlap(i) for i in pool)
    top = [i for i in pool if overlap(i) == best]
    if len(top) != 1:
        raise Bad(f"{e['LocusId']}: cannot tell its pathogenic repeat (candidates {top})")
    return top[0]


def motif_classes(sid, rec):
    """(benign, pathogenic) motif lists for a locus where the motif decides, else None."""
    path = [m.upper() for m in rec.get("pathogenic_motif_reference_orientation") or []]
    ref = [m.upper() for m in rec.get("reference_motif_reference_orientation") or []]
    ben = [m.upper() for m in rec.get("benign_motif_reference_orientation") or []]
    if not path or any("N" in m for m in path + ref):
        return None
    size = len(path[0])
    if not MIN_PERIOD <= size <= MAX_PERIOD or any(len(m) != size for m in path):
        return None
    moved = SIZE_DEPENDENT.get(sid, set())
    pathogenic = [m for m in path if m not in moved]
    benign = []
    for m in ref + ben + sorted(moved):
        if len(m) == size and "N" not in m and not any(same_motif(m, p) for p in pathogenic) \
                and not any(same_motif(m, b) for b in benign):
            benign.append(m)
    if all(any(same_motif(p, r) for r in ref) for p in pathogenic):
        return None  # the pathogenic motif is the reference one: the size decides, not the motif
    return benign, pathogenic


def convert(eh, loci, release):
    by = {r["id"]: r for r in loci}
    missing = [e["LocusId"] for e in eh if e["LocusId"] not in by]
    if missing:
        raise Bad(f"ExpansionHunter loci without a STRchive record: {', '.join(missing)}")
    ids = names(eh, by)
    out, seen = [], set()
    for e in eh:
        sid = e["LocusId"]
        rec = by[sid]
        lid = ids[sid]
        regions = as_list(e["ReferenceRegion"])
        entry = {"LocusId": lid, "LocusStructure": e["LocusStructure"],
                 "ReferenceRegion": e["ReferenceRegion"], "VariantType": e["VariantType"]}
        pi = pathogenic_index(e, rec)
        if len(regions) > 1:
            vids = []
            for i, m in enumerate(motifs_of(e["LocusStructure"])):
                v = lid if i == pi else f"{lid}_{m}"
                vids.append(v if v not in vids else f"{v}_{i + 1}")
            entry["VariantId"] = vids
        for v in as_list(entry.get("VariantId", lid)):
            if v in seen:
                raise Bad(f"two variants are named {v}")
            seen.add(v)
        _, start, end = region(regions[pi])
        size = len(motifs_of(e["LocusStructure"])[pi]) or 1
        span = rec["stop_hg38"] - rec["start_hg38"]
        bmax, pmin = rec.get("benign_max"), rec.get("pathogenic_min")
        if span > SHORT_READ_LIMIT:
            entry["NoThresholds"] = (f"the reference repeat is {span} bases, longer than {SHORT_READ_LIMIT}: "
                                     "short reads cannot size it")
        elif pmin is None:
            entry["NoThresholds"] = "STRchive gives no pathogenic size"
        elif bmax is not None and pmin <= bmax:
            entry["NoThresholds"] = (f"pathogenic from {pmin} units, below the benign maximum {bmax}: "
                                     "a contraction, which Stranger does not grade")
        else:
            entry["NormalMax"] = bmax if bmax is not None else (end - start) // size
            entry["PathologicMin"] = pmin
            if entry["NormalMax"] >= pmin:
                raise Bad(f"{sid}: normal maximum {entry['NormalMax']} is not below pathogenic {pmin}")
        if len(regions) > 1:
            entry["PathologicRegion"] = regions[pi]
        entry["Disease"] = rec.get("disease_id") or sid.split("_")[0]
        if rec.get("inheritance"):
            entry["InheritanceMode"] = "/".join(rec["inheritance"])
        entry.update(Source="STRchive", SourceId=sid, SourceDisplay=f"STRchive-{release}")
        cls = motif_classes(sid, rec)
        if cls:
            entry["BenignMotifs"], entry["PathogenicMotifs"] = cls
        out.append(entry)
    return out


RANGE_COLUMNS = ["locus", "source_id", "normal_max", "intermediate_min", "intermediate_max", "pathogenic_min",
                 "no_thresholds"]


def ranges(cat, loci, release):
    """The --ranges table as text: one row per catalog entry."""
    by = {r["id"]: r for r in loci}

    def cell(x):
        return "" if x is None else str(x)
    lines = [f"# STRchive {release} repeat ranges, one row per ExpansionHunter locus (bin/strchive_convert.py)",
             "\t".join(RANGE_COLUMNS)]
    for e in cat:
        r = by[e["SourceId"]]
        lines.append("\t".join([e["LocusId"], e["SourceId"], cell(e.get("NormalMax")), cell(r.get("intermediate_min")),
                                cell(r.get("intermediate_max")), cell(e.get("PathologicMin")),
                                e.get("NoThresholds", "")]))
    return "\n".join(lines) + "\n"


def write(path, text):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        f.write(text)
    os.replace(tmp, path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--eh", required=True, help="STRchive's ExpansionHunter catalog (JSON)")
    ap.add_argument("--loci", required=True, help="STRchive's loci file (JSON)")
    ap.add_argument("--release", required=True, help="the STRchive release, written into SourceDisplay")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ranges", help="also write the size ranges of each locus to this TSV")
    a = ap.parse_args(argv)
    try:
        with open(a.eh) as f:
            eh = json.load(f)
        with open(a.loci) as f:
            loci = json.load(f)
        cat = convert(eh, loci, a.release)
        table = ranges(cat, loci, a.release) if a.ranges else None
    except (OSError, ValueError, KeyError, Bad) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    if table is not None:
        write(a.ranges, table)
    write(a.out, json.dumps(cat, indent=1) + "\n")
    graded = sum("NormalMax" in e for e in cat)
    print(f"STRchive {a.release}: {len(cat)} loci, {graded} with Stranger thresholds, "
          f"{sum('PathogenicMotifs' in e for e in cat)} read by motif: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
