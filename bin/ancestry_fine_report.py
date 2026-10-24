#!/usr/bin/env python3
"""ancestry_fine_report.py: the Python half of step 26b (fine-grained ancestry).

Standard library only (it runs in the plain python image). Six commands:

  panel-map --meta TSV --bcf-samples LIST [--exclude IDS] [--paint-per-group N] --out samples.tsv
      The panel's sample table, written once by setup.sh --ancestry-fine-panel:
      each person of the panel with population, continental and regional group
      (GROUPS below, '-' when in none), whether it is a reference individual
      (unrelated, in a group, and its regional group has MIN_REFERENCES
      unrelated people; a smaller group is left out of the model) and whether
      it is one of the painting references (at most N per continental group,
      taken in turn from each population).

  check-map --map samples.tsv
      Stops on a sample table given with ANCESTRY_FINE_SAMPLE_MAP that has a
      regional group with 1 to MIN_REFERENCES - 1 reference people, a
      continental group with reference people but fewer than 2 painting
      references, a painting reference who is no reference, or a reference
      person without both groups: RFMix fails on
      a group of one, and a smaller group stands for a handful of people.

  order --pop FILE --q FILE
      The group codes in the column order of the .Q (and .P) file
      `admixture --supervised` wrote: each column is the group whose reference
      individuals have their largest share in it.

  sample-bed --bim ref.bim --vcf genotypes.vcf.gz --sample S --level LEVEL:K:P [...]
             [--bootstrap B] --out-dir DIR
      ADMIXTURE's input for one person: a .bed/.bim/.fam of the reference
      .bim's SNVs the sample has a genotype for, in the .bim's order and allele
      coding, and each level's frequencies (.P) for the same SNVs (the .P.in
      `admixture -P` reads). ADMIXTURE stops on a SNV without a genotype, so a
      missing site, or a genotype with an allele the .bim does not have, is
      left out. Also run_admixture.sh, which the ADMIXTURE image runs (it has
      sh and awk, no python): the projection at each level, then B bootstrap
      replicates of it. ADMIXTURE's own -B crashes with -P, so a replicate is
      a draw, with replacement, of 100 blocks of neighbouring SNVs, written as
      a PLINK .ped (allele 2 is the .bim's A1) and projected the same way.

  chroms --counts counts.tsv --maps DIR [--min-sites N] [--min-cm X]
      The chromosomes to phase and paint: those with N genotyped SNVs or more
      (default 100) spanning X cM or more (default 5). RFMix stops with a
      segmentation fault on a span of a few cM; on a whole genome every
      autosome spans more than 40.

  report --sample S --model DIR --admixture DIR --paint DIR --panel NAME --out DIR
      The step's outputs from ADMIXTURE's .Q and bootstrap replicates of each
      level and RFMix's .msp.tsv of each chromosome: <S>_ancestry_fine.tsv
      (each group's fraction with its bootstrap standard error and 2.5 to 97.5
      percentile range, each group's share of the painted genome, a left_out
      row for each group the model leaves out and a not_painted row for each
      chromosome RFMix failed on, from the paint folder's failed.txt),
      <S>_painting.bed, <S>_painting.png and the plain text <S>_ancestry_fine.txt.

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
    "Bedouin": ("WEA", "BDN"),
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
        "SAR": "Sardinia", "CAU": "Caucasus", "LEV": "Levant", "BDN": "Bedouin", "NAF": "North Africa",
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

# The fewest unrelated reference people a group needs to be modelled. RFMix
# fits its CRF weight on people it simulates from floor(n / 10) + 1 of each
# group's references, and leaves those parents out of the reference while it
# does (select_parents in its gensamples.cpp): a group of one person then has
# no reference haplotype at all, RFMix asks for a zero-size block
# ("mm::allocate() ... requesting zero or negative size (0)") and can
# segfault. Two people is its hard floor. We ask for 10 because a group is
# meant to stand for a population: below 10 people (20 haplotypes) the
# supervised ADMIXTURE frequencies and the forest's reference are a handful
# of people's genotypes. A regional group below it leaves the model, its
# continental group too when none of its regional groups is left, and the
# report names them.
MIN_REFERENCES = 10

# The meta TSV columns read by panel-map.
META_SAMPLE, META_POP, META_RELATED = "s", "population", "sample_filters.all_samples_related"


def die(msg):
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


def label(level, code):
    return LABELS.get(level, {}).get(code, code)


# --- panel-map ---------------------------------------------------------------------

def panel_map(meta_rows, bcf_samples, exclude=(), paint_per_group=100, min_references=MIN_REFERENCES):
    """Rows (sample, population, continent, region, reference, paint) for each
    panel sample, in the panel's order. A person is a reference when unrelated,
    in a group, and its regional group has min_references unrelated people.
    Every continental group is the union of its regional groups, so one that
    keeps any regional group keeps min_references people too, and one that
    keeps none leaves the model with them (left_out_groups names both kinds
    in the report)."""
    # The painting references of a group are its RFMix reference, so the cap
    # must not undo the minimum: a cap of 1 would give every group one person.
    least = max(2, min_references)
    if paint_per_group < least:
        raise ValueError(f"--paint-per-group is {paint_per_group}; each group needs at least {least} painting references")
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
    n = {}
    for row in rows:
        if row[4] == "1":
            n[row[3]] = n.get(row[3], 0) + 1
    for row in rows:
        if row[4] == "1" and n[row[3]] < min_references:
            row[4] = "0"
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


def map_problems(rows, min_references=MIN_REFERENCES):
    """What makes a sample table (rows of six columns, no header) unfit for
    the model: a regional group with some but fewer than min_references
    reference people, a continental group with reference people but fewer
    than 2 painting references, a painting reference who is no reference, a
    reference person without both groups."""
    refs, conts, paint, stray, nogroup = {}, set(), {}, [], []
    for r in rows:
        if r[4] == "1" and (r[2] in ("-", "") or r[3] in ("-", "")):
            nogroup.append(r[0])
        elif r[4] == "1":
            refs[r[3]] = refs.get(r[3], 0) + 1
            conts.add(r[2])
            if r[5] == "1":
                paint[r[2]] = paint.get(r[2], 0) + 1
        elif r[5] == "1":
            stray.append(r[0])
    return ([f"regional group {g} has {n} of the {min_references} reference people it needs"
             for g, n in sorted(refs.items()) if n < min_references]
            + [f"continental group {g} has {paint.get(g, 0)} of the 2 painting references RFMix needs"
               for g in sorted(conts) if paint.get(g, 0) < 2]
            + ([f"painting references who are no reference people: {len(stray)}, e.g. {stray[0]}"] if stray else [])
            + ([f"reference people without a continental and a regional group: {len(nogroup)}, e.g. {nogroup[0]}"]
               if nogroup else []))


def check_map_main(a):
    with open(a.map) as f:
        rows = [line.rstrip("\n").split("\t") for line in f][1:]
    bad = [r for r in rows if len(r) != 6]
    if bad:
        die(f"{a.map}: a row without six columns, e.g. {bad[0]}")
    bad = [r for r in rows if r[4] not in ("0", "1") or r[5] not in ("0", "1")]
    if bad:
        die(f"{a.map}: reference and paint must be 0 or 1, e.g. {bad[0]}")
    problems = map_problems(rows)
    if problems:
        die(f"{a.map}: " + "; ".join(problems))


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
    codes = {level: {r[col] for r in rows if r[4] == "1"} for level, col in (("continent", 2), ("region", 3))}
    gone = left_out_groups(rows, codes)
    if gone:
        print(f"Left out of the model, fewer than {MIN_REFERENCES} unrelated reference people: {left_out_text(gone)}")


def left_out_groups(rows, codes):
    """[(level, code)] of the groups samples.tsv rows have people of but the
    model (codes: {level: group codes it has}) leaves out."""
    out = []
    for level, col in (("continent", 2), ("region", 3)):
        if level in codes:
            out += [(level, c) for c in sorted({r[col] for r in rows if r[col] not in ("-", "")})
                    if c not in codes[level]]
    return out


def left_out_group_names(gone):
    """Each left-out group as one name: 'Oceania (continental and regional
    group)' when both levels leave out a group of the same name."""
    kinds = {}
    for level, code in gone:
        kinds.setdefault(label(level, code), []).append("continental" if level == "continent" else "regional")
    return [f"{name} ({' and '.join(k)} group)" for name, k in kinds.items()]


def left_out_text(gone):
    return ", ".join(left_out_group_names(gone))


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


# Allele pair of a .bed code in a PLINK "12" .ped: allele 2 is the .bim's A1.
# (ADMIXTURE's projection of the .ped equals the .bed's only this way round.)
PED_ALLELES = {0: "2 2", 2: "1 2", 3: "1 1"}

# run_admixture.sh: the projection at each level, then the bootstrap. awk
# draws nb blocks of bs SNVs among the nblk blocks that hold SNVs (with bs
# rounded up, nblk can be below nb); srand(r) makes replicate r the same on
# every run.
BOOT_SH = r"""#!/bin/sh
# Written by bin/ancestry_fine_report.py sample-bed. Usage: sh run_admixture.sh THREADS
set -eu
t=${1:-1}
while read -r level k; do
  admixture -P -j"$t" "@S@_${level}.bed" "$k" > "@S@_${level}.log" 2>&1 \
    || { tail -n 20 "@S@_${level}.log"; exit 1; }
  : > "@S@_${level}.boot.Q"
  r=1
  while [ "$r" -le @B@ ]; do
    awk -v r="$r" -v nb=@NB@ -v k="$k" '
      FNR == 1 { f++ }
      f == 1 { g[n++] = $0; next }
      f == 2 { m[n2++] = $0; next }
      f == 3 { p[n3++] = $0; next }
      END {
        srand(r); bs = int((n + nb - 1) / nb); nblk = int((n + bs - 1) / bs); printf "S S 0 0 0 -9" > "rep.ped"
        for (b = 0; b < nb; b++) {
          s = int(rand() * nblk) * bs
          for (i = s; i < s + bs && i < n; i++) { printf " %s", g[i] > "rep.ped"; print m[i] > "rep.map"; print p[i] > ("rep." k ".P.in") }
        }
        printf "\n" > "rep.ped"
      }' "@S@.geno" "@S@.map" "@S@_${level}.${k}.P.in"
    admixture -P "rep.ped" "$k" > rep.log 2>&1 || { tail -n 20 rep.log; exit 1; }
    cat "rep.${k}.Q" >> "@S@_${level}.boot.Q"
    rm -f rep.ped rep.map "rep.${k}.P.in" "rep.${k}.P" "rep.${k}.Q"
    r=$((r + 1))
  done
done < levels.txt
"""


def sample_bed_main(a):
    bim = read_bim(a.bim)
    if not bim:
        die(f"{a.bim} has no SNV")
    keep, data = sample_bed(bim, read_vcf_genotypes(a.vcf))
    if len(keep) < a.min_sites:
        die(f"only {len(keep)} of the {len(bim)} reference SNVs have a genotype in {a.vcf} (at least {a.min_sites} needed)")
    os.makedirs(a.out_dir, exist_ok=True)
    out = lambda name: os.path.join(a.out_dir, name)  # noqa: E731
    levels = []
    for spec in a.level:
        level, k, p_path = spec.split(":", 2)
        base = out(f"{a.sample}_{level}")
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
        levels.append(f"{level} {k}\n")
    codes = data[3:]
    with open(out(f"{a.sample}.geno"), "w") as f:
        f.writelines(PED_ALLELES[c] + "\n" for c in codes)
    with open(out(f"{a.sample}.map"), "w") as f:
        f.writelines(" ".join(bim[i][4].split()[:4]) + "\n" for i in keep)
    with open(out("levels.txt"), "w") as f:
        f.writelines(levels)
    nb = max(1, min(100, len(keep) // 5))
    with open(out("run_admixture.sh"), "w") as f:
        f.write(BOOT_SH.replace("@S@", a.sample).replace("@B@", str(a.bootstrap)).replace("@NB@", str(nb)))
    with open(out(f"{a.sample}_sites.txt"), "w") as f:
        f.write(f"{len(keep)} of the {len(bim)} SNVs ADMIXTURE uses genotyped\n")
    print(f"{len(keep)} of the {len(bim)} SNVs ADMIXTURE uses genotyped")


# --- chroms ------------------------------------------------------------------------

def read_map(path):
    """[(position, cM)] of a PLINK .map (chromosome, id, cM, position)."""
    out = []
    with open(path) as f:
        for line in f:
            c = line.split()
            if len(c) >= 4:
                out.append((int(c[3]), float(c[2])))
    out.sort()
    return out


def cm_at(gmap, pos):
    """The genetic position of pos, linear between the map's markers."""
    if not gmap:
        return 0.0
    if pos <= gmap[0][0]:
        return gmap[0][1]
    if pos >= gmap[-1][0]:
        return gmap[-1][1]
    lo, hi = 0, len(gmap) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if gmap[mid][0] <= pos:
            lo = mid
        else:
            hi = mid
    (p0, c0), (p1, c1) = gmap[lo], gmap[hi]
    return c0 + (c1 - c0) * (pos - p0) / (p1 - p0) if p1 > p0 else c0


def paint_chroms(counts, maps_dir, min_sites=100, min_cm=5.0):
    """Chromosomes of counts.tsv (chrom, SNVs, first, last) to phase and paint."""
    out = []
    for chrom, n, first, last in counts:
        path = os.path.join(maps_dir, f"{chrom}.map")
        if n < min_sites or not os.path.isfile(path):
            continue
        gmap = read_map(path)
        if cm_at(gmap, last) - cm_at(gmap, first) >= min_cm:
            out.append(chrom)
    return sorted(out, key=chrom_key)


def read_counts(path):
    out = []
    with open(path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 4 and c[0].startswith("chr"):
                out.append((c[0], int(c[1]), int(c[2]), int(c[3])))
    return out


def chroms_main(a):
    print(" ".join(paint_chroms(read_counts(a.counts), a.maps, a.min_sites, a.min_cm)))


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


def percentile(xs, q):
    xs = sorted(xs)
    if not xs:
        return None
    i = q * (len(xs) - 1)
    lo = int(math.floor(i))
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (i - lo)


def fractions(level, k, codes, adm_dir, sample):
    """Rows of one level: each group's fraction, and from the bootstrap
    replicates its standard error and 2.5 to 97.5 percentile range."""
    q = read_matrix(os.path.join(adm_dir, f"{sample}_{level}.{k}.Q"))
    if len(q) != 1 or len(q[0]) != k:
        die(f"{adm_dir}/{sample}_{level}.{k}.Q is not one row of {k} values")
    boot_path = os.path.join(adm_dir, f"{sample}_{level}.boot.Q")
    boot = read_matrix(boot_path) if os.path.isfile(boot_path) else []
    if any(len(r) != k for r in boot):
        die(f"{boot_path} has a row that is not {k} values")
    rows = []
    for j, (code, f) in enumerate(zip(codes, q[0])):
        col = [r[j] for r in boot]
        se = lo = hi = None
        if len(col) >= 2:
            mean = sum(col) / len(col)
            se = math.sqrt(sum((x - mean) ** 2 for x in col) / (len(col) - 1))
            lo, hi = percentile(col, 0.025), percentile(col, 0.975)
        rows.append({"level": level, "code": code, "group": label(level, code), "fraction": f,
                     "se": se, "low": lo, "high": hi})
    rows.sort(key=lambda r: -r["fraction"])
    return rows, len(boot)


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


LIMITS = [
    "What this means: each name is a group of reference populations we put together and named,",
    "not a nationality. The fractions say which groups your genome is most similar to; the range",
    "after each is the 2.5 to 97.5 percentile of the bootstrap. Continental groups are well apart;",
    "neighbouring regions are not: splits within Iberia are weak for every method, and Greek and",
    "Southern Italian ancestry overlap (the panel has no Greek population at all). Read a regional",
    "fraction with its range, never as one exact percentage. See docs/26b-ancestry-fine.md.",
]


def text_block(sample, panel, fr, n_boot, shares, painted_bp, n_sites, painted_chroms, gone=(), failed=()):
    out = [f"Fine-grained ancestry of {sample} (step 26b)", f"Reference panel: {panel}; {n_sites}", ""]
    if gone:
        out += [f"Left out of the model: {left_out_text(gone)}. The panel has fewer than {MIN_REFERENCES} unrelated",
                "reference people for it, too few to stand for a population, so the model has no such group and",
                "any ancestry from it is spread over the groups below.", ""]
    for level in LEVELS:
        rows = [r for r in fr if r["level"] == level]
        if not rows:
            continue
        name = "Continental" if level == "continent" else "Regional"
        out.append(f"{name} groups (ADMIXTURE projection on {len(rows)} groups, {n_boot[level]} bootstrap replicates):")
        for r in rows:
            if level == "region" and r["fraction"] < 0.01 and (r["high"] or 0) < 0.05:
                continue
            spread = f"  ({pct(r['low'])} to {pct(r['high'])})" if r["low"] is not None else ""
            out.append(f"  {r['group']:<50} {pct(r['fraction']):>5}{spread}")
        out.append("")
    if shares:
        out.append(f"Chromosome painting (RFMix, continental groups) of {', '.join(painted_chroms)}: "
                   f"share of the {painted_bp / 1e6:.1f} Mb painted on the two haplotypes")
        for code, sh in sorted(shares.items(), key=lambda kv: -kv[1]):
            out.append(f"  {label('continent', code):<50} {pct(sh):>5}")
    elif not failed:
        out.append("Chromosome painting: none. No chromosome had enough genotyped SNVs spanning 5 cM.")
    else:
        out.append("Chromosome painting: none.")
    if failed:
        out.append(f"  Not painted: RFMix failed on {', '.join(failed)} (the step's log has its last lines).")
    out.append("")
    return "\n".join(out + LIMITS) + "\n"


def report_main(a):
    levels = read_levels(a.model)
    fr, n_boot = [], {}
    for level, k, codes in levels:
        if level not in LABELS:
            die(f"{a.model}/levels.tsv: unknown level {level}")
        rows, n_boot[level] = fractions(level, k, codes, a.admixture, a.sample)
        fr += rows
    rows = painting(a.paint) if a.paint and os.path.isdir(a.paint) else []
    painted_bp, shares = painting_shares(rows)
    painted_chroms = sorted({r[0] for r in rows}, key=chrom_key)
    failed = []
    failed_path = os.path.join(a.paint, "failed.txt") if a.paint else ""
    if failed_path and os.path.isfile(failed_path):
        with open(failed_path) as f:
            failed = sorted({line.strip() for line in f if line.strip()}, key=chrom_key)
    gone = []
    samples_path = os.path.join(os.path.dirname(os.path.normpath(a.model)), "samples.tsv")
    if os.path.isfile(samples_path):
        with open(samples_path) as f:
            samples = [line.rstrip("\n").split("\t") for line in f][1:]
        gone = left_out_groups([r for r in samples if len(r) >= 4], {level: set(codes) for level, _, codes in levels})
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
    table += [["painting", code, label("continent", code), fmt(sh), ".", ".", "."]
              for code, sh in sorted(shares.items(), key=lambda kv: -kv[1])]
    table += [["left_out", code, left_out_text([(level, code)]), ".", ".", ".", "."] for level, code in gone]
    table += [["not_painted", c, "RFMix failed", ".", ".", ".", "."] for c in failed]
    write_rows(base + "_ancestry_fine.tsv", ["level", "code", "group", "fraction", "se", "low", "high"], table)
    write_rows(base + "_painting.bed", ["#chrom", "start", "end", "haplotype", "group"], rows)
    with open(base + "_painting.png.tmp", "wb") as f:
        f.write(painting_png(rows))
    os.replace(base + "_painting.png.tmp", base + "_painting.png")
    text = text_block(a.sample, a.panel, fr, n_boot, shares, painted_bp, n_sites, painted_chroms, gone, failed)
    with open(base + "_ancestry_fine.txt.tmp", "w") as f:
        f.write(text)
    os.replace(base + "_ancestry_fine.txt.tmp", base + "_ancestry_fine.txt")
    sys.stdout.write(text)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("panel-map")
    p.add_argument("--meta", required=True)
    p.add_argument("--bcf-samples", required=True)
    p.add_argument("--exclude", default="")
    p.add_argument("--paint-per-group", type=int, default=100)
    p.add_argument("--out", required=True)
    p = sub.add_parser("check-map")
    p.add_argument("--map", required=True)
    p = sub.add_parser("order")
    p.add_argument("--pop", required=True)
    p.add_argument("--q", required=True)
    p = sub.add_parser("sample-bed")
    p.add_argument("--bim", required=True)
    p.add_argument("--vcf", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--level", action="append", required=True, help="LEVEL:K:P_FILE, once per level")
    p.add_argument("--bootstrap", type=int, default=100)
    p.add_argument("--min-sites", type=int, default=100)
    p.add_argument("--out-dir", required=True)
    p = sub.add_parser("chroms")
    p.add_argument("--counts", required=True)
    p.add_argument("--maps", required=True)
    p.add_argument("--min-sites", type=int, default=100)
    p.add_argument("--min-cm", type=float, default=5.0)
    p = sub.add_parser("report")
    p.add_argument("--sample", required=True)
    p.add_argument("--model", required=True)
    p.add_argument("--admixture", required=True)
    p.add_argument("--paint", default="")
    p.add_argument("--panel", default="")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"panel-map": panel_map_main, "check-map": check_map_main, "order": order_main, "sample-bed": sample_bed_main,
     "chroms": chroms_main, "report": report_main}[a.cmd](a)


if __name__ == "__main__":
    main()
