#!/usr/bin/env python3
"""ancestry_fine_report.py: the Python half of step 26b (fine-grained ancestry).

Standard library only (it runs in the plain python image). Four commands:

  panel-map --meta TSV --bcf-samples LIST [--exclude IDS] [--paint-per-group N] --out samples.tsv
      The panel's sample table, written once by setup.sh --ancestry-fine-panel:
      each sample of the panel with its population, its continental and its
      regional group (GROUPS below, '-' when it is in none), whether it is a
      reference individual (unrelated and in a group) and whether it is one of
      the reference individuals of the painting (at most N per continental
      group, taken in turn from each population).

  order --pop FILE --q FILE
      Prints the group codes in the column order of an ADMIXTURE .Q (or .P)
      file written by `admixture --supervised`: each column is the group whose
      reference individuals have their largest share in it.

  sample-bed --bim ref.bim --vcf genotypes.vcf.gz --sample S --level LEVEL:K:P [...] --out-dir DIR
      The sample as a one-person PLINK .bed/.bim/.fam on the SNVs of the
      reference .bim, in its order and allele coding, and beside it the rows
      of each level's ADMIXTURE frequencies (.P) for the same SNVs, the .P.in
      `admixture -P` projects the sample on. ADMIXTURE stops on a SNV without
      a genotype, so a site the VCF has no genotype for, or a genotype with an
      allele the .bim does not have, is left out of both.

  report --sample S --model DIR --admixture DIR --paint DIR --panel NAME --out DIR
      The step's outputs from ADMIXTURE's .Q and .Q_se of each level and
      RFMix's .msp.tsv of each chromosome: <S>_ancestry_fine.tsv (the
      fractions with their bootstrap spread, and the share of the painted
      genome per group), <S>_painting.bed, <S>_painting.png and the plain
      text block <S>_ancestry_fine.txt.

A group name is ours: a set of reference populations we put together and
named, not a nationality. docs/26b-ancestry-fine.md lists them.
"""
import argparse
import csv
import gzip
import math
import os
import struct
import sys
import zlib

# population of the panel -> (continental group, regional group). The panel
# is gnomAD's HGDP + 1000 Genomes callset (population labels of its v2
# metadata). Admixed populations (between the continental groups) are in no
# group, so they are never a reference: ACB, ASW, CLM, MXL, PEL, PUR, Hazara
# and Uygur.
GROUPS = {
    # Africa
    "YRI": ("AFR", "WAF"), "ESN": ("AFR", "WAF"), "GWD": ("AFR", "WAF"), "MSL": ("AFR", "WAF"),
    "Yoruba": ("AFR", "WAF"), "Mandenka": ("AFR", "WAF"),
    "LWK": ("AFR", "EAF"), "BantuKenya": ("AFR", "EAF"),
    "Biaka": ("AFR", "CAF"), "Mbuti": ("AFR", "CAF"),
    "San": ("AFR", "SAF"), "BantuSouthAfrica": ("AFR", "SAF"),
    # Europe, the Middle East and North Africa
    "CEU": ("WEA", "NWE"), "GBR": ("WEA", "NWE"), "Orcadian": ("WEA", "NWE"), "French": ("WEA", "NWE"),
    "FIN": ("WEA", "NEE"), "Russian": ("WEA", "NEE"),
    "IBS": ("WEA", "IBE"), "Basque": ("WEA", "IBE"),
    "TSI": ("WEA", "ITA"), "Tuscan": ("WEA", "ITA"), "BergamoItalian": ("WEA", "ITA"),
    "Sardinian": ("WEA", "SAR"),
    "Adygei": ("WEA", "CAU"),
    "Druze": ("WEA", "LEV"), "Palestinian": ("WEA", "LEV"),
    "Bedouin": ("WEA", "BED"),
    "Mozabite": ("WEA", "NAF"),
    # Central and South Asia
    "GIH": ("CSA", "CSA"), "PJL": ("CSA", "CSA"), "BEB": ("CSA", "CSA"), "STU": ("CSA", "CSA"),
    "ITU": ("CSA", "CSA"), "Balochi": ("CSA", "CSA"), "Brahui": ("CSA", "CSA"), "Makrani": ("CSA", "CSA"),
    "Pathan": ("CSA", "CSA"), "Sindhi": ("CSA", "CSA"), "Burusho": ("CSA", "CSA"), "Kalash": ("CSA", "CSA"),
    # East Asia and Siberia
    "CHB": ("EAS", "EAS"), "CHS": ("EAS", "EAS"), "CDX": ("EAS", "EAS"), "JPT": ("EAS", "EAS"),
    "KHV": ("EAS", "EAS"), "Han": ("EAS", "EAS"), "NorthernHan": ("EAS", "EAS"), "Japanese": ("EAS", "EAS"),
    "Dai": ("EAS", "EAS"), "Daur": ("EAS", "EAS"), "Hezhen": ("EAS", "EAS"), "Lahu": ("EAS", "EAS"),
    "Miao": ("EAS", "EAS"), "Mongolian": ("EAS", "EAS"), "Naxi": ("EAS", "EAS"), "Oroqen": ("EAS", "EAS"),
    "She": ("EAS", "EAS"), "Tu": ("EAS", "EAS"), "Tujia": ("EAS", "EAS"), "Xibo": ("EAS", "EAS"),
    "Yi": ("EAS", "EAS"), "Cambodian": ("EAS", "EAS"),
    "Yakut": ("EAS", "SIB"),
    # The Americas (HGDP's Native American populations)
    "Karitiana": ("AMR", "AMR"), "Maya": ("AMR", "AMR"), "Pima": ("AMR", "AMR"), "Surui": ("AMR", "AMR"),
    "Colombian": ("AMR", "AMR"),
    # Oceania
    "Bougainville": ("OCE", "OCE"), "PapuanHighlands": ("OCE", "OCE"), "PapuanSepik": ("OCE", "OCE"),
}

LABELS = {
    "continent": {
        "AFR": "Africa",
        "WEA": "West Eurasia (Europe, Middle East, North Africa)",
        "CSA": "Central and South Asia",
        "EAS": "East Asia",
        "AMR": "Americas (Native American)",
        "OCE": "Oceania",
    },
    "region": {
        "WAF": "West Africa", "EAF": "East Africa", "CAF": "Central African rainforest", "SAF": "Southern Africa",
        "NWE": "Northwestern Europe", "NEE": "Northeastern Europe", "IBE": "Iberia", "ITA": "Italy",
        "SAR": "Sardinia", "CAU": "Caucasus", "LEV": "Levant", "BED": "Bedouin", "NAF": "North Africa",
        "CSA": "Central and South Asia", "EAS": "East Asia", "SIB": "Siberia",
        "AMR": "Americas (Native American)", "OCE": "Oceania",
    },
}
LEVELS = ("continent", "region")

# Painting colours of the continental groups, and grey where nothing is painted.
COLOURS = {"AFR": (230, 159, 0), "WEA": (0, 114, 178), "CSA": (0, 158, 115), "EAS": (204, 121, 167),
           "AMR": (213, 94, 0), "OCE": (86, 180, 233)}
UNPAINTED = (225, 225, 225)

# GRCh38 lengths of chr1 to chr22, for the painting picture.
CHROM_LEN = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
             138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
             83257441, 80373285, 58617616, 64444167, 46709983, 50818468]

# The meta TSV columns read by panel-map.
META_SAMPLE, META_POP, META_RELATED = "s", "population", "sample_filters.all_samples_related"


def die(msg):
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


def label(level, code):
    return LABELS.get(level, {}).get(code, code)


# --- panel-map ---------------------------------------------------------------------

def panel_map(meta_rows, bcf_samples, exclude=(), paint_per_group=100):
    """Rows (sample, population, continent, region, reference, paint) for each
    panel sample, in the panel's order."""
    meta = {r[META_SAMPLE]: r for r in meta_rows}
    missing = [s for s in bcf_samples if s not in meta]
    if missing:
        raise ValueError(f"{len(missing)} panel samples have no metadata row, e.g. {missing[0]}")
    rows = []
    for s in bcf_samples:
        if s in exclude:
            continue
        r = meta[s]
        pop = r[META_POP]
        cont, reg = GROUPS.get(pop, ("-", "-"))
        ref = cont != "-" and r[META_RELATED] == "false"
        rows.append([s, pop, cont, reg, "1" if ref else "0", "0"])
    # Painting references: in turn from each population of a continental
    # group, by sample id, until the group has paint_per_group.
    by_group = {}
    for row in rows:
        if row[4] == "1":
            by_group.setdefault(row[2], {}).setdefault(row[1], []).append(row)
    for pops in by_group.values():
        queues = [sorted(v, key=lambda x: x[0]) for _, v in sorted(pops.items())]
        taken = 0
        while taken < paint_per_group and any(queues):
            for q in queues:
                if q and taken < paint_per_group:
                    q.pop(0)[5] = "1"
                    taken += 1
    return rows


def panel_map_main(a):
    with open(a.meta, newline="") as f:
        meta_rows = list(csv.DictReader(f, delimiter="\t"))
    for col in (META_SAMPLE, META_POP, META_RELATED):
        if meta_rows and col not in meta_rows[0]:
            die(f"{a.meta} has no column {col}")
    with open(a.bcf_samples) as f:
        samples = [line.strip() for line in f if line.strip()]
    exclude = {x for x in (a.exclude or "").replace(",", " ").split() if x}
    unknown = sorted(exclude - set(samples))
    if unknown:
        die(f"--exclude names samples the panel does not have: {', '.join(unknown)}")
    try:
        rows = panel_map(meta_rows, samples, exclude, a.paint_per_group)
    except ValueError as e:
        die(str(e))
    write_rows(a.out, ["sample", "population", "continent", "region", "reference", "paint"], rows)
    n_ref = sum(r[4] == "1" for r in rows)
    print(f"{len(rows)} panel samples, {n_ref} reference individuals, {sum(r[5] == '1' for r in rows)} painting references")


def write_rows(path, header, rows):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    os.replace(tmp, path)


# --- order -------------------------------------------------------------------------

def read_matrix(path):
    with open(path) as f:
        return [[float(x) for x in line.split()] for line in f if line.strip()]


def column_order(pops, q):
    """The group of each column of a supervised ADMIXTURE .Q: the group whose
    reference individuals have their largest mean share in it."""
    if len(pops) != len(q):
        raise ValueError(f"{len(pops)} .pop lines but {len(q)} .Q rows")
    k = len(q[0]) if q else 0
    groups = sorted({p for p in pops if p not in ("-", "")})
    if len(groups) != k:
        raise ValueError(f"{len(groups)} groups in the .pop file but K = {k}")
    order = []
    for j in range(k):
        means = {}
        for g in groups:
            vals = [row[j] for p, row in zip(pops, q) if p == g]
            means[g] = sum(vals) / len(vals)
        best = max(means, key=means.get)
        if means[best] < 0.5:
            raise ValueError(f"column {j + 1} has no group with a mean share above 0.5")
        order.append(best)
    if len(set(order)) != k:
        raise ValueError(f"two columns map to the same group: {order}")
    return order


def order_main(a):
    with open(a.pop) as f:
        pops = [line.strip() for line in f]
    try:
        print(",".join(column_order(pops, read_matrix(a.q))))
    except ValueError as e:
        die(f"{a.q}: {e}")


# --- sample-bed --------------------------------------------------------------------

def norm_chrom(c):
    """plink2 writes chr22 as 22 in a .bim; the VCFs say chr22."""
    return c[3:] if c.startswith("chr") else c


def read_bim(path):
    out = []
    with open(path) as f:
        for line in f:
            c = line.split()
            if len(c) >= 6:
                out.append((norm_chrom(c[0]), int(c[3]), c[4], c[5], line))   # chrom, pos, A1, A2, the line
    return out


def gt_alleles(gt):
    parts = gt.replace("|", "/").split("/")
    if len(parts) != 2 or "." in parts:
        return None
    return parts


def read_vcf_genotypes(path):
    """{(chrom, pos): (ref, [alts], [allele indices])} of the one sample of a VCF."""
    out = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            fmt = c[8].split(":")
            if "GT" not in fmt:
                continue
            idx = gt_alleles(c[9].split(":")[fmt.index("GT")])
            if idx is None:
                continue
            alts = [] if c[4] == "." else c[4].split(",")
            out[(norm_chrom(c[0]), int(c[1]))] = (c[3], alts, idx)
    return out


def bed_code(rec, a1, a2):
    """The PLINK .bed code of one genotype for a .bim line with alleles A1, A2:
    0 = homozygous A1, 2 = heterozygous, 3 = homozygous A2, None = no call
    (no genotype, or an allele the .bim does not have)."""
    if rec is None:
        return None
    ref, alts, idx = rec
    alleles = [ref] + alts
    try:
        called = [alleles[int(i)] for i in idx]
    except (ValueError, IndexError):
        return None
    n1 = sum(x == a1 for x in called)
    n2 = sum(x == a2 for x in called)
    if n1 + n2 != 2:
        return None
    return {2: 0, 1: 2, 0: 3}[n1]


def sample_bed(bim, genotypes):
    """The called SNVs of the reference .bim: (row indices, .bed bytes).
    ADMIXTURE stops on a SNV with no genotype at all, which for one person
    is every missing site, so those rows are left out (and out of P)."""
    keep, out = [], bytearray(b"\x6c\x1b\x01")
    for i, (chrom, pos, a1, a2, _) in enumerate(bim):
        code = bed_code(genotypes.get((chrom, pos)), a1, a2)
        if code is not None:
            keep.append(i)
            out.append(code)
    return keep, bytes(out)


def sample_bed_main(a):
    bim = read_bim(a.bim)
    if not bim:
        die(f"{a.bim} has no SNV")
    keep, data = sample_bed(bim, read_vcf_genotypes(a.vcf))
    if len(keep) < a.min_sites:
        die(f"only {len(keep)} of the {len(bim)} reference SNVs have a genotype in {a.vcf} (at least {a.min_sites} needed)")
    os.makedirs(a.out_dir, exist_ok=True)
    for spec in a.level:
        level, k, p_path = spec.split(":", 2)
        base = os.path.join(a.out_dir, f"{a.sample}_{level}")
        with open(p_path) as f:
            p_rows = f.readlines()
        if len(p_rows) != len(bim):
            die(f"{p_path} has {len(p_rows)} rows but {a.bim} {len(bim)} SNVs")
        with open(base + ".bed", "wb") as f:
            f.write(data)
        with open(base + ".bim", "w") as f:
            f.writelines(bim[i][4] for i in keep)
        with open(base + ".fam", "w") as f:
            f.write(f"{a.sample} {a.sample} 0 0 0 -9\n")
        with open(f"{base}.{k}.P.in", "w") as f:
            f.writelines(p_rows[i] for i in keep)
    with open(os.path.join(a.out_dir, f"{a.sample}_sites.txt"), "w") as f:
        f.write(f"{len(keep)} of {len(bim)} reference SNVs genotyped\n")
    print(f"{len(keep)} of {len(bim)} reference SNVs genotyped")


# --- report ------------------------------------------------------------------------

def read_levels(model):
    """[(level, K, [codes in .P column order])] from model/levels.tsv."""
    out = []
    with open(os.path.join(model, "levels.tsv")) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 3 and not c[0].startswith("#"):
                out.append((c[0], int(c[1]), c[2].split(",")))
    if not out:
        die(f"{model}/levels.tsv lists no level")
    return out


def fractions(level, k, codes, adm_dir, sample):
    q = read_matrix(os.path.join(adm_dir, f"{sample}_{level}.{k}.Q"))
    if len(q) != 1 or len(q[0]) != k:
        die(f"{adm_dir}/{sample}_{level}.{k}.Q is not one row of {k} values")
    se_path = os.path.join(adm_dir, f"{sample}_{level}.{k}.Q_se")
    se = read_matrix(se_path)[0] if os.path.isfile(se_path) else [None] * k
    rows = []
    for code, f, s in zip(codes, q[0], se):
        lo = hi = None
        if s is not None:
            lo, hi = max(0.0, f - 1.96 * s), min(1.0, f + 1.96 * s)
        rows.append({"level": level, "code": code, "group": label(level, code), "fraction": f,
                     "se": s, "low": lo, "high": hi})
    rows.sort(key=lambda r: -r["fraction"])
    return rows


def read_msp(path):
    """RFMix's .msp.tsv: (codes {int: name}, [(chrom, start, end, hap1, hap2)])."""
    codes, segs = {}, []
    with open(path) as f:
        for line in f:
            if line.startswith("#Subpopulation order/codes:"):
                for item in line.split(":", 1)[1].split():
                    name, _, num = item.rpartition("=")
                    codes[int(num)] = name
            elif line.startswith("#"):
                continue
            else:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 8:
                    segs.append((c[0], int(c[1]), int(c[2]), int(c[6]), int(c[7])))
    if not codes:
        raise ValueError(f"{path} has no '#Subpopulation order/codes' line")
    return codes, segs


def painting(paint_dir):
    """BED rows (chrom, start, end, haplotype, code): RFMix's windows, the
    neighbouring windows of one haplotype with the same group merged."""
    rows = []
    for name in sorted(os.listdir(paint_dir)):
        if not name.endswith(".msp.tsv"):
            continue
        codes, segs = read_msp(os.path.join(paint_dir, name))
        for hap in (1, 2):
            cur = None
            for chrom, start, end, h1, h2 in segs:
                code = codes[h1 if hap == 1 else h2]
                if cur and cur[0] == chrom and cur[4] == code and start <= cur[2] + 1:
                    cur[2] = end
                else:
                    if cur:
                        rows.append(cur)
                    cur = [chrom, start - 1, end, hap, code]
            if cur:
                rows.append(cur)
    rows.sort(key=lambda r: (chrom_key(r[0]), r[1], r[3]))
    return rows


def chrom_key(c):
    c = c[3:] if c.startswith("chr") else c
    return (0, int(c)) if c.isdigit() else (1, c)


def painting_shares(rows):
    total = sum(r[2] - r[1] for r in rows)
    by = {}
    for r in rows:
        by[r[4]] = by.get(r[4], 0) + r[2] - r[1]
    return total, {k: v / total for k, v in by.items()} if total else {}


def png_bytes(width, height, pixels):
    """A plain RGB PNG of rows of (r, g, b) tuples."""
    raw = b"".join(b"\x00" + bytes(v for px in row for v in px) for row in pixels)

    def chunk(kind, data):
        return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data) & 0xffffffff)

    return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b""))


def painting_png(rows, width=1000, bar=6, gap=6):
    """Chromosomes 1 to 22 as two bars each (the two haplotypes), drawn to
    scale, each painted segment in its group's colour, the rest grey."""
    scale = (width - 20) / CHROM_LEN[0]
    height = gap + 22 * (2 * bar + 1 + gap)
    img = [[(255, 255, 255)] * width for _ in range(height)]
    for i, n in enumerate(CHROM_LEN):
        y0 = gap + i * (2 * bar + 1 + gap)
        for hap in (0, 1):
            y = y0 + hap * (bar + 1)
            for yy in range(y, y + bar):
                for x in range(10, 10 + int(n * scale)):
                    img[yy][x] = UNPAINTED
    for chrom, start, end, hap, code in rows:
        key = chrom_key(chrom)
        if key[0] != 0 or not 1 <= key[1] <= 22:
            continue
        i = key[1] - 1
        y = gap + i * (2 * bar + 1 + gap) + (hap - 1) * (bar + 1)
        x0 = 10 + int(start * scale)
        x1 = max(x0 + 1, 10 + int(end * scale))
        colour = COLOURS.get(code, (0, 0, 0))
        for yy in range(y, y + bar):
            for x in range(x0, min(x1, width - 10)):
                img[yy][x] = colour
    return png_bytes(width, height, img)


def pct(x):
    return "." if x is None else f"{100 * x:.0f}%"


def text_block(sample, panel, fr, shares, painted_bp, n_sites):
    out = [f"Fine-grained ancestry of {sample} (step 26b)",
           f"Reference panel: {panel}; {n_sites}",
           ""]
    for level in LEVELS:
        rows = [r for r in fr if r["level"] == level]
        if not rows:
            continue
        out.append(f"{'Continental' if level == 'continent' else 'Regional'} groups (ADMIXTURE, projection on "
                   f"{len(rows)} reference groups), largest first:")
        for r in rows:
            if r["fraction"] < 0.01 and level == "region":
                continue
            spread = f" (bootstrap range {pct(r['low'])} to {pct(r['high'])})" if r["low"] is not None else ""
            out.append(f"  {r['group']:<50} {pct(r['fraction']):>5}{spread}")
        out.append("")
    if shares:
        out.append(f"Chromosome painting (RFMix, continental groups), share of the "
                   f"{painted_bp / 1e6:.1f} Mb painted on both haplotypes:")
        for code, s in sorted(shares.items(), key=lambda kv: -kv[1]):
            out.append(f"  {label('continent', code):<50} {pct(s):>5}")
        out.append("")
    out += [
        "What this means: each name is a group of reference populations we put together and named,",
        "not a nationality. The fractions say which groups your genome is most similar to, with the",
        "bootstrap range of each. Continental groups are well separated; neighbouring regions are not:",
        "Iberian region-level splits are weak for every method, and Greek and Southern Italian ancestry",
        "overlap. Read a regional fraction together with its range. See docs/26b-ancestry-fine.md.",
    ]
    return "\n".join(out) + "\n"


def report_main(a):
    levels = read_levels(a.model)
    fr = []
    for level, k, codes in levels:
        if level not in LABELS:
            die(f"{a.model}/levels.tsv: unknown level {level}")
        fr += fractions(level, k, codes, a.admixture, a.sample)
    rows = painting(a.paint) if a.paint and os.path.isdir(a.paint) else []
    painted_bp, shares = painting_shares(rows)
    n_sites = "."
    stats = os.path.join(a.admixture, f"{a.sample}_sites.txt")
    if os.path.isfile(stats):
        with open(stats) as f:
            n_sites = f.read().strip()
    os.makedirs(a.out, exist_ok=True)
    base = os.path.join(a.out, a.sample)
    fmt = lambda x: "." if x is None else f"{x:.4f}"  # noqa: E731
    table = [[r["level"], r["code"], r["group"], fmt(r["fraction"]), fmt(r["se"]), fmt(r["low"]), fmt(r["high"])]
             for r in fr]
    table += [["painting", code, label("continent", code), fmt(s), ".", ".", "."]
              for code, s in sorted(shares.items(), key=lambda kv: -kv[1])]
    write_rows(base + "_ancestry_fine.tsv", ["level", "code", "group", "fraction", "se", "low", "high"], table)
    write_rows(base + "_painting.bed", ["#chrom", "start", "end", "haplotype", "group"], rows)
    with open(base + "_painting.png.tmp", "wb") as f:
        f.write(painting_png(rows))
    os.replace(base + "_painting.png.tmp", base + "_painting.png")
    with open(base + "_ancestry_fine.txt.tmp", "w") as f:
        f.write(text_block(a.sample, a.panel, fr, shares, painted_bp, n_sites))
    os.replace(base + "_ancestry_fine.txt.tmp", base + "_ancestry_fine.txt")
    sys.stdout.write(text_block(a.sample, a.panel, fr, shares, painted_bp, n_sites))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("panel-map")
    p.add_argument("--meta", required=True)
    p.add_argument("--bcf-samples", required=True)
    p.add_argument("--exclude", default="")
    p.add_argument("--paint-per-group", type=int, default=100)
    p.add_argument("--out", required=True)
    p = sub.add_parser("order")
    p.add_argument("--pop", required=True)
    p.add_argument("--q", required=True)
    p = sub.add_parser("sample-bed")
    p.add_argument("--bim", required=True)
    p.add_argument("--vcf", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--level", action="append", required=True, help="LEVEL:K:P_FILE, once per level")
    p.add_argument("--min-sites", type=int, default=100)
    p.add_argument("--out-dir", required=True)
    p = sub.add_parser("report")
    p.add_argument("--sample", required=True)
    p.add_argument("--model", required=True)
    p.add_argument("--admixture", required=True)
    p.add_argument("--paint", default="")
    p.add_argument("--panel", default="")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"panel-map": panel_map_main, "order": order_main, "sample-bed": sample_bed_main,
     "report": report_main}[a.cmd](a)


if __name__ == "__main__":
    main()
