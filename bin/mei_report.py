#!/usr/bin/env python3
"""mei_report.py: the Alu, L1 and SVA insertions MEGAnE calls, and the known disease insertions.

  mei_report.py --sample S --vcf S_megane.vcf --sites assets/mei_sites.tsv
                --depth S_mei_depth.tsv --out-calls S_mei.tsv --out-sites S_mei_sites.tsv

scripts/50-mei.sh and the MEI process run MEGAnE (Kojima et al. 2023, Nat Genet,
PMID 37169872) on the sample's BAM. Its VCF has one record per insertion it
found, ALT <INS:ME>, the repeat class in INFO SVTYPE (SINE/Alu, LINE/L1,
Retroposon/SVA, and older classes too), MEPRED=PASS when the inserted sequence
was identified as that element, the subfamily in INFO MEI (or MEI_lbp and
MEI_rbp, one per breakpoint), and the copy number in FORMAT CN (0, 1 or 2).

--out-calls: one row per insertion of Alu, L1 or SVA with MEPRED=PASS and at
least one copy. Columns: chrom, pos, end (MEGAnE's right breakpoint), element,
subfamily, copies, filter (MEGAnE's FILTER: PASS, or its flags such as LC, low
confidence), confidence (CN_conf: H or L, how sure the copy number is). Other
classes (DNA transposons, old LTRs) are left out: none of them is active in
human genomes, so a call of one is read as an artefact.

--out-sites: every row of the catalog (assets/mei_sites.tsv) with three
columns added:
  status      seen        a listed insertion of the same element inside the
                          row's window (start to end)
              not seen    none, with at least MIN_DEPTH reads of depth over
                          the window
              not callable  the BAM has no reads of that contig, or the mean
                          depth over the window is under MIN_DEPTH: too few
                          reads to tell
  mean_depth  the mean depth over the window (samtools bedcov -c: summed depth
              / window length), "." when the contig is not in the BAM
  call        the insertion seen: position, subfamily, copies, filter; else "."
--depth is samtools bedcov -c over the catalog windows on the BAM's contigs:
chrom, start (0-based), end, gene, summed depth, reads. MIN_DEPTH is this
pipeline's choice, not a published value. "Not seen" is not proof of absence:
an SVA is about 2 kb of repeat, and short reads find its edges, not its body.

Used by scripts/50-mei.sh and the MEI_REPORT process (bin/ is on the task
PATH); tests/test_mei_report.py checks it. Standard library only.
"""
import argparse
import csv
import sys

MIN_DEPTH = 10  # mean depth over a catalog window before "not seen" is said
ELEMENTS = {"SINE/Alu": "Alu", "LINE/L1": "L1", "Retroposon/SVA": "SVA"}
CALL_HEADER = ["chrom", "pos", "end", "element", "subfamily", "copies", "filter", "confidence"]
SITE_COLUMNS = ["gene", "chrom", "start", "end", "site", "element", "inheritance", "condition", "clinvar", "pmid"]


def die(msg):
    sys.exit(f"mei_report.py: {msg}")


def info_dict(text):
    return dict(x.split("=", 1) if "=" in x else (x, "") for x in text.split(";") if x)


def subfamily(info):
    """The first subfamily MEGAnE names: from MEI, else from either breakpoint."""
    for key in ("MEI", "MEI_lbp", "MEI_rbp"):
        name = info.get(key, "").split("|")[0].split(",")[0]
        if name and name not in ("NA", "pA"):
            return name
    return "."


def read_calls(path):
    calls, header = [], False
    with open(path) as f:
        for n, line in enumerate(f, 1):
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = True
                continue
            if not line.strip():
                continue
            if not header:
                die(f"{path}: a record before the #CHROM line (line {n}): not a VCF")
            f_ = line.rstrip("\n").split("\t")
            if len(f_) < 10:
                die(f"{path} line {n}: {len(f_)} columns, a VCF with one sample has 10")
            info = info_dict(f_[7])
            if "SVTYPE" not in info or "MEPRED" not in info:
                die(f"{path} line {n}: no SVTYPE or MEPRED in INFO: not MEGAnE's VCF")
            element = ELEMENTS.get(info["SVTYPE"])
            fmt = dict(zip(f_[8].split(":"), f_[9].split(":")))
            cn = fmt.get("CN", ".")
            if element is None or info["MEPRED"] != "PASS" or not cn.isdigit() or int(cn) < 1:
                continue
            calls.append({"chrom": f_[0], "pos": int(f_[1]), "end": info.get("MEI_rpos", f_[1]),
                          "element": element, "subfamily": subfamily(info), "copies": cn,
                          "filter": f_[6], "confidence": info.get("CN_conf", ".")})
    if not header:
        die(f"{path}: no #CHROM line: not a VCF")
    return calls


def read_sites(path):
    rows = [r for r in open(path) if r.strip() and not r.startswith("#")]
    reader = list(csv.DictReader(rows, delimiter="\t"))
    if not reader or any(c not in reader[0] for c in SITE_COLUMNS):
        die(f"{path}: no rows, or a column of {', '.join(SITE_COLUMNS)} is missing")
    return reader


def read_depth(path):
    """samtools bedcov -c over the 4-column BED: chrom, start, end, gene, sum, reads."""
    depth = {}
    with open(path) as f:
        for n, line in enumerate(f, 1):
            if not line.strip():
                continue
            f_ = line.rstrip("\n").split("\t")
            if len(f_) < 6 or not f_[1].isdigit() or not f_[2].isdigit():
                die(f"{path} line {n}: not samtools bedcov -c output of a 4-column BED")
            length = int(f_[2]) - int(f_[1])
            depth[(f_[3], f_[0])] = float(f_[4]) / length if length > 0 else 0.0
    return depth


def chrom_rank(chrom):
    """chr1 to chr22, then chrX, chrY and the rest, as a karyotype is read."""
    c = chrom[3:] if chrom.startswith("chr") else chrom
    return int(c) if c.isdigit() else {"X": 23, "Y": 24}.get(c, 25)


def site_status(site, calls, depth):
    key = (site["gene"], site["chrom"])
    if key not in depth:
        return "not callable", ".", "."
    mean = depth[key]
    lo, hi = int(site["start"]), int(site["end"])
    hits = [c for c in calls if c["chrom"] == site["chrom"] and c["element"] == site["element"] and lo <= c["pos"] <= hi]
    if hits:
        c = hits[0]
        call = f"{c['chrom']}:{c['pos']} {c['subfamily']} {c['copies']} cop{'y' if c['copies'] == '1' else 'ies'} {c['filter']}"
        return "seen", f"{mean:.1f}", call
    if mean < MIN_DEPTH:
        return "not callable", f"{mean:.1f}", "."
    return "not seen", f"{mean:.1f}", "."


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sample", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--depth", required=True)
    ap.add_argument("--out-calls", required=True)
    ap.add_argument("--out-sites", required=True)
    a = ap.parse_args(argv)
    calls = read_calls(a.vcf)
    sites = read_sites(a.sites)
    depth = read_depth(a.depth)
    calls.sort(key=lambda c: (chrom_rank(c["chrom"]), c["chrom"], c["pos"]))
    with open(a.out_calls, "w") as out:
        out.write("\t".join(CALL_HEADER) + "\n")
        for c in calls:
            out.write("\t".join(str(c[k]) for k in CALL_HEADER) + "\n")
    with open(a.out_sites, "w") as out:
        out.write("\t".join(SITE_COLUMNS + ["status", "mean_depth", "call"]) + "\n")
        for s in sites:
            out.write("\t".join([s[k] for k in SITE_COLUMNS] + list(site_status(s, calls, depth))) + "\n")
    by = {}
    for c in calls:
        by[c["element"]] = by.get(c["element"], 0) + 1
    print(f"{a.sample}: {len(calls)} insertions ({', '.join(f'{k} {v}' for k, v in sorted(by.items())) or 'none'})")


if __name__ == "__main__":
    main()
