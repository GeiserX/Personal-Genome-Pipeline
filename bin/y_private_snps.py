#!/usr/bin/env python3
"""y_private_snps.py: the chrY SNVs of a sample that no public tree names (step 37).

A haplogroup call places a paternal line on the deepest branch the public
trees already have. The variants that would place it further down, next to
the men closest to it, are the ones on no public list yet. Step 37 lists
them here, from the sample's own calls, so the file can go with the chrY BAM
to a service that compares it with other men.

What goes in the list, and why each gate is there:
  * chrY only, outside the pseudoautosomal regions (the analysis set masks
    PAR on chrY, and a PAR call is really chrX)
  * one REF base and one ALT base, no other ALT allele at the site: a
    private SNP is a single-base change, and a multiallelic record at one
    chrY position is a mapping artefact on a haploid chromosome
  * the genotype carries the ALT allele and no REF allele (`1`, `1/1`,
    `1|1`): chrY is haploid outside PAR, so a heterozygous call there is an
    artefact, not a variant of this line
  * FILTER PASS (or `.`), depth at least --min-depth (10, Yleaf's -r
    default), genotype quality at least --min-gq (20), and at least
    --min-fraction (0.9, Yleaf's -b default) of the reads carrying ALT
  * depth no more than --max-depth-factor (2.0) times the median depth of
    the sample's own chrY calls: where reads from more than one place in the
    genome pile up, depth runs high and the call is not single-copy. Reads
    are what this pipeline has, so the repetitive part of the Y is found in
    the sample's own depth, not from a region list of another study. 0 turns
    the gate off.
  * not on a public SNP list: every position in the two GRCh38 SNP tables
    the current-tree folder installs (`positions_ftdna_hg38.txt`, made for
    FamilyTreeDNA's tree, and `positions_yfull_hg38.txt`, the one Yleaf 4
    uses for its YFull tree) is removed. That subtraction is what makes the
    list private: without it the file would be mostly the markers the
    haplogroup call is already built on.

Usage:
  y_private_snps.py --vcf VCF --tree-dir DIR --sample S --out-prefix P
                    [--min-depth N] [--min-gq N] [--min-fraction F]
                    [--max-depth-factor F]

Output: P_y_private_snps.tsv, one row per SNV (Pos, Ref, Alt, Depth,
Alt_fraction, GQ, Genotype), and a count on stdout.

A private SNP is unconfirmed until another man is found who carries it. Some
of these rows are sequencing or mapping errors, and some are mutations that
only this sample carries, which no tree can place. The number itself is not
a result: it rises with coverage and falls with every public tree release
that names one of them. Standard library only.
"""
import gzip
import sys

PAR1_END = 2781479     # GRCh38 chrY PAR1 is 10,001-2,781,479
PAR2_START = 56887903  # GRCh38 chrY PAR2 is 56,887,903-57,217,415
MIN_DEPTH = 10         # Yleaf's -r default
MIN_GQ = 20
MIN_FRACTION = 0.9     # Yleaf's -b default (90%)
MAX_DEPTH_FACTOR = 2.0
TABLES = ("positions_ftdna_hg38.txt", "positions_yfull_hg38.txt")
BASES = set("ACGT")
CHRY = ("chrY", "Y")


def public_positions(tree_dir):
    """Every chrY position the installed SNP tables name, whatever the tree
    uses them for. The tables are TSV with the position in column 4
    (bin/y_tree_place.py reads the same files)."""
    from pathlib import Path
    d = Path(tree_dir)
    pos = set()
    found = 0
    for t in TABLES:
        p = d / t
        if not p.is_file():
            continue
        found += 1
        with open(p) as f:
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 4 and c[3].isdigit():
                    pos.add(int(c[3]))
    if found < len(TABLES):
        sys.exit(f"ERROR: {tree_dir} has no {' and no '.join(TABLES)}: install the current Y tree "
                 "with scripts/setup.sh --yleaf-data <genome_dir>")
    if not pos:
        sys.exit(f"ERROR: the SNP tables in {tree_dir} name no position")
    return pos


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def fmt_field(fmt, sample, key):
    keys = fmt.split(":")
    vals = sample.split(":")
    if key not in keys:
        return None
    i = keys.index(key)
    return vals[i] if i < len(vals) else None


def carries_alt(gt):
    """True when the genotype is the ALT allele alone: haploid '1', or '1/1'
    or '1|1'. A call with a REF allele, a no-call or another allele is not
    this line's variant."""
    if gt is None:
        return False
    alleles = gt.replace("|", "/").split("/")
    return bool(alleles) and all(a == "1" for a in alleles)


def depth_and_fraction(fmt, sample, ref_alt_ad_only=True):
    """(depth, alt fraction) from AD when the record has it, else from DP and
    VAF. (None, None) when neither is usable."""
    ad = fmt_field(fmt, sample, "AD")
    if ad:
        parts = [p for p in ad.split(",")]
        if all(p.isdigit() for p in parts) and len(parts) >= 2:
            counts = [int(p) for p in parts]
            total = sum(counts)
            if total > 0:
                return total, counts[1] / total
    dp = fmt_field(fmt, sample, "DP")
    if dp and dp.isdigit() and int(dp) > 0:
        vaf = fmt_field(fmt, sample, "VAF")
        try:
            return int(dp), float(vaf)
        except (TypeError, ValueError):
            return int(dp), None
    return None, None


def candidates(vcf_path):
    """Every chrY SNV record of the one sample, before the gates: dicts with
    pos, ref, alt, depth, fraction, gq, gt, filt."""
    out = []
    with opener(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10 or c[0] not in CHRY or not c[1].isdigit():
                continue
            ref, alt = c[3].upper(), c[4].upper()
            if ref not in BASES or alt not in BASES:
                continue
            gq = fmt_field(c[8], c[9], "GQ")
            depth, fraction = depth_and_fraction(c[8], c[9])
            out.append({"pos": int(c[1]), "ref": ref, "alt": alt, "depth": depth,
                        "fraction": fraction, "gq": int(gq) if (gq or "").isdigit() else None,
                        "gt": fmt_field(c[8], c[9], "GT"), "filt": c[6]})
    return out


def median(values):
    v = sorted(values)
    n = len(v)
    if not n:
        return None
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2


def private_snps(vcf_path, tree_dir, min_depth=MIN_DEPTH, min_gq=MIN_GQ,
                 min_fraction=MIN_FRACTION, max_depth_factor=MAX_DEPTH_FACTOR):
    """(rows, counts): the private SNVs, and how many records each gate dropped."""
    public = public_positions(tree_dir)
    rows = []
    counts = dict.fromkeys(("records", "par", "filtered", "genotype", "depth", "quality",
                            "fraction", "excess_depth", "public"), 0)
    kept = []
    for r in candidates(vcf_path):
        counts["records"] += 1
        if r["pos"] <= PAR1_END or r["pos"] >= PAR2_START:
            counts["par"] += 1
            continue
        if r["filt"] not in ("PASS", ".", ""):
            counts["filtered"] += 1
            continue
        if not carries_alt(r["gt"]):
            counts["genotype"] += 1
            continue
        if r["depth"] is None or r["depth"] < min_depth:
            counts["depth"] += 1
            continue
        if r["gq"] is None or r["gq"] < min_gq:
            counts["quality"] += 1
            continue
        if r["fraction"] is None or r["fraction"] < min_fraction:
            counts["fraction"] += 1
            continue
        kept.append(r)
    # The single-copy gate needs the sample's own chrY depth, so it runs after
    # the rest: the median of the calls that passed them.
    med = median([r["depth"] for r in kept])
    limit = med * max_depth_factor if (med and max_depth_factor > 0) else None
    for r in kept:
        if limit is not None and r["depth"] > limit:
            counts["excess_depth"] += 1
            continue
        if r["pos"] in public:
            counts["public"] += 1
            continue
        rows.append(r)
    rows.sort(key=lambda r: r["pos"])
    counts["median_depth"] = med if med is not None else 0
    counts["private"] = len(rows)
    return rows, counts


def write(rows, counts, sample, prefix, min_depth, min_gq, min_fraction, max_depth_factor):
    with open(f"{prefix}_y_private_snps.tsv", "w") as o:
        o.write(f"# sample\t{sample}\n")
        o.write(f"# gates\tmin_depth={min_depth} min_gq={min_gq} min_fraction={min_fraction} "
                f"max_depth_factor={max_depth_factor} median_chrY_depth={counts['median_depth']}\n")
        o.write(f"# dropped\tPAR={counts['par']} not_PASS={counts['filtered']} genotype={counts['genotype']} "
                f"depth={counts['depth']} quality={counts['quality']} alt_fraction={counts['fraction']} "
                f"excess_depth={counts['excess_depth']} on_public_lists={counts['public']}\n")
        o.write("Pos\tRef\tAlt\tDepth\tAlt_fraction\tGQ\tGenotype\n")
        for r in rows:
            frac = "." if r["fraction"] is None else f"{r['fraction']:.3f}"
            o.write(f"{r['pos']}\t{r['ref']}\t{r['alt']}\t{r['depth']}\t{frac}\t{r['gq']}\t{r['gt']}\n")


def main(argv):
    import argparse
    ap = argparse.ArgumentParser(prog="y_private_snps.py", description=__doc__.split("\n\n")[0])
    for a in ("--vcf", "--tree-dir", "--sample", "--out-prefix"):
        ap.add_argument(a, required=True)
    ap.add_argument("--min-depth", type=int, default=MIN_DEPTH)
    ap.add_argument("--min-gq", type=int, default=MIN_GQ)
    ap.add_argument("--min-fraction", type=float, default=MIN_FRACTION)
    ap.add_argument("--max-depth-factor", type=float, default=MAX_DEPTH_FACTOR)
    a = ap.parse_args(argv)
    if a.min_depth < 1:
        sys.exit("ERROR: --min-depth must be 1 or more")
    if not 0 < a.min_fraction <= 1:
        sys.exit("ERROR: --min-fraction must be above 0 and at most 1")
    if a.max_depth_factor < 0:
        sys.exit("ERROR: --max-depth-factor must be 0 or more")
    rows, counts = private_snps(a.vcf, a.tree_dir, a.min_depth, a.min_gq, a.min_fraction,
                                a.max_depth_factor)
    write(rows, counts, a.sample, a.out_prefix, a.min_depth, a.min_gq, a.min_fraction,
          a.max_depth_factor)
    print(f"Private chrY SNVs: {len(rows)} past the public lists "
          f"({counts['records']} chrY SNV records read, {counts['public']} already named by a public list, "
          f"median chrY depth {counts['median_depth']})")
    print("  Each one is unconfirmed until another man is found who carries it; some are read errors "
          "(docs/37-y-haplogroup.md).")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
