#!/usr/bin/env python3
"""Step 44: where a genome falls among ancient genomes on a PCA of present-day
West Eurasians (smartpca, EIGENSOFT), with the ancient samples of the Allen
Ancient DNA Resource (AADR) projected onto the same axes.

Subcommands, in the order they run:

  select       (setup.sh) read the AADR .anno and .ind files and choose the
               people: the present-day West Eurasians of MODERN_GROUPS build
               the axes, and the ancient West Eurasians with enough sites
               covered are projected. Writes the relabelled .ind (everyone
               else is Ignore), the population list and ancients.tsv.
  lift         (setup.sh) move the 1240K sites of the AADR .snp file from
               GRCh37 to GRCh38 with UCSC's chain file, and keep a site only
               when one of its two alleles is the GRCh38 reference base.
  person-geno  (step 44) the person's genotype at each site, from bcftools'
               output on the gVCF, appended as one more column to the
               present-day panel's EIGENSTRAT genotype file.
  align        (setup.sh, step 44) check that the axes of two smartpca runs
               on the same present-day panel agree, up to the sign of each.
  report       (step 44) the nearest ancient populations, the summary and
               the plot.

Only the Python standard library is used. Every number the report shows comes
from smartpca's .evec files and the AADR annotation.
"""
import argparse
import bisect
import collections
import csv
import gzip
import html
import math
import os
import statistics
import sys

# Present-day groups of the AADR 1240K file (Group ID, as written in the
# .anno file of v66.p1) that build the axes: Europe, the Caucasus and the Near
# East, from the 1000 Genomes Project, HGDP and SGDP shotgun genomes. North
# African groups (Mozabite) and groups with recent East Eurasian ancestry
# (Saami, Bashkir, Tatar, Yakut) are left out, so the axes describe West
# Eurasian variation only.
MODERN_GROUPS = (
    "Abkhasian", "Adygei", "Albanian", "Armenian", "Basque", "BedouinA", "BedouinB",
    "Bulgarian", "CEU", "Chechen", "Czech", "Druze", "English", "Estonian", "FIN",
    "Finnish", "French", "GBR", "Georgian", "Greek", "Greek_Crete", "Hungarian", "IBS",
    "Icelandic", "Iranian", "Italian_North", "Jew_Iraqi", "Jordanian", "Lezgin",
    "Norwegian", "Orcadian", "Ossetian", "Palestinian", "Polish", "Russian", "Sardinian",
    "Scottish", "Spanish", "Spanish_North", "TSI", "Turkish",
)
MODERN_DATA_TYPE = "Shotgun.diploid"
# Ancient samples projected: dated, passing the AADR's own assessment, with
# this many of the 1240K autosomal sites covered, inside this box.
MIN_ANCIENT_SNPS = 100_000
WEST_EURASIA = {"lat": (25.0, 75.0), "lon": (-25.0, 75.0)}
PASSING = ("PASS", "PROVISIONAL_PASS", "MERGE_PASS")
ANCIENT_LABEL = "Ancient"
PERSON_LABEL = "Person"
# PCs the distance uses, each weighted by the square root of its eigenvalue.
DISTANCE_PCS = 4
TOP = 20
# Two runs on the same panel give the same axes; any PC whose present-day
# coordinates correlate less than this between the two is a panel mismatch.
MIN_AXIS_CORRELATION = 0.999
CAVEAT = ("Closeness on a PCA is genetic similarity, not descent: a genome near an ancient group "
          "shares ancestry with it, it does not descend from those people.")
QPADM_NOTE = "qpAdm models of ancestry proportions are not part of this step."

# Header cells of the .anno file this script reads, matched by their start
# (the AADR writes long explanations after the name).
ANNO_COLUMNS = {
    "id": "Genetic ID",
    "individual": "Individual ID",
    "date_bp": "Date mean in BP",
    "date_sd": "Date standard deviation in BP",
    "full_date": "Full Date",
    "group": "Group ID",
    "locality": "Locality",
    "country": "Political Entity",
    "lat": "Latitude",
    "lon": "Longitude",
    "data_type": "Data type",
    "snps": "SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)",
    "assessment": "ASSESSMENT",
}
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


class Fail(Exception):
    pass


def opener(path):
    return gzip.open(path, "rt", errors="replace") if path.endswith(".gz") else open(path, errors="replace")


# --- select ---------------------------------------------------------------------

def anno_columns(header):
    cells = [c.strip().strip('"').strip() for c in header]
    out = {}
    for key, name in ANNO_COLUMNS.items():
        # ASSESSMENT must not match ASSESSMENT WARNINGS: the exact cell first.
        hits = [i for i, c in enumerate(cells) if c == name] or [i for i, c in enumerate(cells) if c.startswith(name)]
        if key == "assessment":
            hits = [i for i in hits if cells[i] == name]
        if not hits:
            raise Fail(f"the .anno file has no column '{name}'")
        out[key] = hits[0]
    return out


def to_float(s):
    try:
        v = float(s)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def passing(assessment):
    return assessment.strip().upper().startswith(PASSING)


def read_anno(path):
    with opener(path) as f:
        rows = list(csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE))
    if not rows:
        raise Fail(f"{path} is empty")
    col = anno_columns(rows[0])
    out = []
    for r in rows[1:]:
        if len(r) <= max(col.values()):
            continue
        out.append({k: r[i].strip() for k, i in col.items()})
    return out


def choose(anno, min_snps=MIN_ANCIENT_SNPS):
    """Return {genetic id: ('modern', group) | ('ancient', row)}: the
    present-day reference people and the ancient people to project, one
    genetic id per individual (the one with the most sites)."""
    modern_groups = set(MODERN_GROUPS)
    best = {}
    for r in anno:
        if not passing(r["assessment"]):
            continue
        bp = to_float(r["date_bp"])
        snps = to_float(r["snps"]) or 0
        if bp == 0 and r["group"] in modern_groups and r["data_type"] == MODERN_DATA_TYPE:
            kind = "modern"
        elif bp is not None and bp > 0 and snps >= min_snps:
            lat, lon = to_float(r["lat"]), to_float(r["lon"])
            if lat is None or lon is None:
                continue
            if not (WEST_EURASIA["lat"][0] <= lat <= WEST_EURASIA["lat"][1]
                    and WEST_EURASIA["lon"][0] <= lon <= WEST_EURASIA["lon"][1]):
                continue
            kind = "ancient"
        else:
            continue
        key = (kind, r["individual"] or r["id"])
        if key not in best or snps > best[key][0]:
            best[key] = (snps, r)
    out = {}
    for (kind, _), (_, r) in best.items():
        out[r["id"]] = (kind, r["group"] if kind == "modern" else r)
    return out


def cmd_select(a):
    anno = read_anno(a.anno)
    chosen = choose(anno, a.min_snps)
    n_mod = n_anc = 0
    groups = set()
    with open(a.ind) as f, open(a.out_ind + ".tmp", "w") as o:
        for line in f:
            c = line.split()
            if len(c) < 3:
                continue
            kind, val = chosen.get(c[0], (None, None))
            if kind == "modern":
                label = val
                n_mod += 1
                groups.add(val)
            elif kind == "ancient":
                label = ANCIENT_LABEL
                n_anc += 1
            else:
                label = "Ignore"
            o.write(f"{c[0]}\t{c[1]}\t{label}\n")
    missing = sorted(set(MODERN_GROUPS) - groups)
    if n_mod < 50 or n_anc < 1:
        os.remove(a.out_ind + ".tmp")
        raise Fail(f"the .ind file matches {n_mod} present-day reference people and {n_anc} ancient ones; "
                   "is it the .ind file of the same release as the .anno file?")
    os.replace(a.out_ind + ".tmp", a.out_ind)
    with open(a.poplist, "w") as o:
        for g in sorted(groups):
            o.write(g + "\n")
    with open(a.ancients, "w") as o:
        o.write("id\tgroup\tdate_bp\tdate_sd\tfull_date\tlocality\tcountry\tlat\tlon\tsnps\n")
        for gid, (kind, r) in sorted(chosen.items()):
            if kind != "ancient":
                continue
            o.write("\t".join(clean(r[k]) for k in ("id", "group", "date_bp", "date_sd", "full_date",
                                                    "locality", "country", "lat", "lon", "snps")) + "\n")
    print(f"Present-day reference people: {n_mod} in {len(groups)} groups"
          + (f" (not in this file: {', '.join(missing)})" if 0 < len(missing) <= 10 else
             f" ({len(missing)} of the {len(MODERN_GROUPS)} groups are not in this file)" if missing else ""))
    print(f"Ancient people projected: {n_anc} (at least {a.min_snps:,} sites, West Eurasia)")
    return 0


def clean(s):
    return " ".join(str(s).replace("\t", " ").split()) or "."


# --- lift -----------------------------------------------------------------------

def read_chain(path):
    """UCSC chain file -> {source chrom: sorted [(start, end, chain no, target
    chrom, target start, target strand, target size)]}, 0-based half-open
    ungapped blocks."""
    blocks = collections.defaultdict(list)
    n = -1
    with opener(path) as f:
        cur = None
        for line in f:
            c = line.split()
            if not c:
                cur = None
                continue
            if c[0] == "chain":
                n += 1
                t_name, t_strand, t_start = c[2], c[4], int(c[5])
                q_name, q_size, q_strand, q_start = c[7], int(c[8]), c[9], int(c[10])
                if t_strand != "+":
                    raise Fail(f"{path}: chain {c[12] if len(c) > 12 else n} has source strand {t_strand}")
                cur = [t_name, t_start, q_name, q_start, q_strand, q_size, n]
                continue
            if cur is None:
                raise Fail(f"{path}: a block line outside a chain: {line.strip()}")
            size = int(c[0])
            t_name, t_pos, q_name, q_pos, q_strand, q_size, k = cur
            blocks[t_name].append((t_pos, t_pos + size, k, q_name, q_pos, q_strand, q_size))
            if len(c) >= 3:
                cur[1] = t_pos + size + int(c[1])
                cur[3] = q_pos + size + int(c[2])
            else:
                cur = None
    for v in blocks.values():
        v.sort()
    return blocks


def index_blocks(blocks):
    """read_chain's blocks -> {chrom: (blocks, their starts, longest block)}."""
    return {k: (v, [x[0] for x in v], max(x[1] - x[0] for x in v)) for k, v in blocks.items()}


def lift_point(index, chrom, pos1):
    """1-based GRCh37 position -> (chrom, 1-based position, strand) on the
    target, or None when no block, or blocks of more than one chain, cover it."""
    if chrom not in index:
        return None
    b, starts, longest = index[chrom]
    p = pos1 - 1
    hits = []
    # Blocks of different chains can overlap: every block that starts at or
    # before p, back to the longest block's length, may still reach it.
    j = bisect.bisect_right(starts, p) - 1
    while j >= 0 and b[j][0] > p - longest:
        if b[j][1] > p:
            hits.append(b[j])
        j -= 1
    if len({h[2] for h in hits}) != 1:
        return None
    s, _, _, q_name, q_pos, q_strand, q_size = hits[0]
    q = q_pos + (p - s)
    if q_strand == "-":
        q = q_size - 1 - q
    return q_name, q + 1, q_strand


class Fasta:
    """Random access to a FASTA with a samtools .fai index."""

    def __init__(self, path):
        self.f = open(path, "rb")
        self.idx = {}
        with open(path + ".fai") as fai:
            for line in fai:
                c = line.split("\t")
                self.idx[c[0]] = (int(c[1]), int(c[2]), int(c[3]), int(c[4]))

    def base(self, chrom, pos1):
        if chrom not in self.idx:
            return None
        length, offset, line_bases, line_width = self.idx[chrom]
        if not 1 <= pos1 <= length:
            return None
        p = pos1 - 1
        self.f.seek(offset + (p // line_bases) * line_width + p % line_bases)
        return self.f.read(1).decode("ascii", "replace").upper()


def eigen_chrom(c):
    """The chromosome column of an EIGENSOFT .snp line -> 1..22, or None."""
    try:
        n = int(c)
    except ValueError:
        return None
    return n if 1 <= n <= 22 else None


def cmd_lift(a):
    blocks = index_blocks(read_chain(a.chain))
    fasta = Fasta(a.fasta)
    kept = []
    dropped = collections.Counter()
    n = 0
    with opener(a.snp) as f:
        for line in f:
            c = line.split()
            if len(c) < 6:
                continue
            n += 1
            chrom = eigen_chrom(c[1])
            a1, a2 = c[4].upper(), c[5].upper()
            if chrom is None:
                dropped["not an autosome"] += 1
                continue
            if a1 not in COMPLEMENT or a2 not in COMPLEMENT or a1 == a2:
                dropped["alleles not two different bases"] += 1
                continue
            hit = lift_point(blocks, f"chr{chrom}", int(float(c[3])))
            if hit is None:
                dropped["no single place on GRCh38"] += 1
                continue
            q_chrom, q_pos, strand = hit
            if q_chrom != f"chr{chrom}":
                dropped["moves to another chromosome"] += 1
                continue
            if strand == "-":
                a1, a2 = COMPLEMENT[a1], COMPLEMENT[a2]
            ref = fasta.base(q_chrom, q_pos)
            if ref not in (a1, a2):
                dropped["neither allele is the GRCh38 base"] += 1
                continue
            kept.append((chrom, q_pos, c[0], a1, a2))
    kept.sort()
    with open(a.out + ".tmp", "w") as o:
        o.write("snp\tchrom\tpos\tallele1\tallele2\n")
        for chrom, pos, sid, a1, a2 in kept:
            o.write(f"{sid}\tchr{chrom}\t{pos}\t{a1}\t{a2}\n")
    os.replace(a.out + ".tmp", a.out)
    print(f"Sites of the .snp file: {n:,}; on GRCh38: {len(kept):,}")
    for why, k in sorted(dropped.items()):
        print(f"  dropped, {why}: {k:,}")
    if not kept:
        raise Fail("no site could be placed on GRCh38: check the chain file and the reference")
    return 0


# --- person-geno ------------------------------------------------------------------

def read_sites(path):
    """sites_grch38.tsv -> {(chrom, pos): [(snp, allele1, allele2)]}."""
    out = collections.defaultdict(list)
    with open(path) as f:
        next(f, None)
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 5:
                out[(c[1], int(c[2]))].append((c[0], c[3], c[4]))
    return out


def dosage(gt, a1, a2):
    """Copies of allele1 in a bcftools %TGT genotype, or None (missing, a
    haploid or half call, or an allele that is neither of the two)."""
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2 or any(x not in (a1, a2) for x in alleles):
        return None
    return sum(1 for x in alleles if x == a1)


def person_dosages(sites, genotype_lines):
    """bcftools query lines (CHROM POS REF ALT GT) -> ({snp: dosage}, counts)."""
    out = {}
    counts = collections.Counter()
    for line in genotype_lines:
        c = line.rstrip("\n").split("\t")
        if len(c) < 5:
            continue
        here = sites.get((c[0], int(c[1])))
        if not here:
            continue
        for sid, a1, a2 in here:
            d = dosage(c[4].upper(), a1, a2)
            if d is None:
                counts["other_allele"] += 1
            else:
                out[sid] = d
                counts["called"] += 1
    return out, counts


def cmd_person_geno(a):
    sites = read_sites(a.sites)
    with open(a.genotypes) as f:
        dos, counts = person_dosages(sites, f)
    n = used = 0
    with open(a.snp) as snp, open(a.moderns_geno) as geno, open(a.out + ".tmp", "w") as o:
        width = None
        for sline in snp:
            c = sline.split()
            if not c:
                continue
            g = geno.readline().rstrip("\n")
            if not g:
                raise Fail(f"{a.moderns_geno} has fewer lines than {a.snp}")
            if width is None:
                width = len(g)
            elif len(g) != width:
                raise Fail(f"{a.moderns_geno} line {n + 1} has {len(g)} genotypes, not {width}")
            d = dos.get(c[0])
            if d is not None:
                used += 1
            o.write(g + ("9" if d is None else str(d)) + "\n")
            n += 1
        if geno.readline():
            raise Fail(f"{a.moderns_geno} has more lines than {a.snp}")
    os.replace(a.out + ".tmp", a.out)
    with open(a.stats, "w") as o:
        o.write(f"panel_sites\t{n}\nsites_genotyped\t{counts['called'] + counts['other_allele']}\n"
                f"other_allele\t{counts['other_allele']}\nsites_used\t{used}\n")
    print(f"Sites of the panel: {n:,}; genotyped from the gVCF: {counts['called'] + counts['other_allele']:,}; "
          f"used: {used:,} (an allele that is neither of a site's two: {counts['other_allele']:,})")
    if used == 0:
        raise Fail("none of the panel's sites was genotyped: is the gVCF on GRCh38 with chr names?")
    return 0


# --- evec, align ------------------------------------------------------------------

def read_evec(path):
    """smartpca .evec -> (eigenvalues, {id: (label, [pc...])})."""
    evals, rows = [], {}
    with open(path) as f:
        for line in f:
            c = line.split()
            if not c:
                continue
            if c[0].startswith("#"):
                if c[0] == "#eigvals:":
                    evals = [float(x) for x in c[1:]]
                continue
            rows[c[0]] = (c[-1], [float(x) for x in c[1:-1]])
    if not evals or not rows:
        raise Fail(f"{path} is not a smartpca .evec file (no #eigvals: line or no rows)")
    return evals, rows


def correlation(x, y):
    mx, my = statistics.fmean(x), statistics.fmean(y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sx = math.sqrt(sum((a - mx) ** 2 for a in x))
    sy = math.sqrt(sum((b - my) ** 2 for b in y))
    return sxy / (sx * sy) if sx > 0 and sy > 0 else 0.0


def axis_signs(panel_rows, run_rows, n_pcs):
    """The sign that turns each of the run's PCs into the panel's, from the
    present-day people both runs used to build the axes. A correlation below
    MIN_AXIS_CORRELATION on any PC stops: the two are not the same axes."""
    shared = sorted(k for k, (lab, _) in run_rows.items()
                    if k in panel_rows and lab not in (ANCIENT_LABEL, PERSON_LABEL)
                    and panel_rows[k][0] == lab)
    if len(shared) < 10:
        raise Fail(f"only {len(shared)} present-day people are in both runs; the panel and this run do not match")
    signs = []
    for k in range(n_pcs):
        r = correlation([panel_rows[i][1][k] for i in shared], [run_rows[i][1][k] for i in shared])
        if abs(r) < MIN_AXIS_CORRELATION:
            raise Fail(f"PC{k + 1} of this run correlates {r:.4f} with the panel's over {len(shared)} present-day "
                       "people: the axes differ. Install the panel again (setup.sh --ancient-panel).")
        signs.append(1.0 if r > 0 else -1.0)
    return signs


def cmd_align(a):
    _, panel = read_evec(a.panel_evec)
    _, run = read_evec(a.run_evec)
    signs = axis_signs(panel, run, a.pcs)
    print("Axes agree with the panel's: " + " ".join(f"PC{i + 1}{'+' if s > 0 else '-'}" for i, s in enumerate(signs)))
    return 0


# --- report -----------------------------------------------------------------------

def year_label(bp):
    """A date in years before 1950 -> 'N BCE' or 'N CE', the AADR's own
    reading (its 2600-2000 BCE has a mean of 4250 BP): BCE = BP - 1950 and
    CE = 1950 - BP, with 1950 BP read as 1 BCE."""
    y = 1950 - int(round(bp))
    return f"{y} CE" if y >= 1 else f"{max(-y, 1)} BCE"


def range_label(old_bp, young_bp):
    """The span of a population's dates, oldest first, e.g. '5300-4900 BCE',
    '200-450 CE', or '300 BCE-200 CE'."""
    a, b = year_label(old_bp), year_label(young_bp)
    if a == b:
        return a
    (na, ea), (nb, eb) = a.split(" "), b.split(" ")
    return f"{na}-{nb} {ea}" if ea == eb else f"{a}-{b}"


def read_ancients(path):
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE):
            out[r["id"]] = r
    return out


def weights(evals, n):
    return [math.sqrt(max(v, 0.0)) for v in evals[:n]]


def distance(p, q, w):
    return math.sqrt(sum((wi * (a - b)) ** 2 for a, b, wi in zip(p, q, w)))


def most_common(values):
    vals = [v for v in values if v and v != "."]
    if not vals:
        return "."
    c = collections.Counter(vals).most_common()
    top = [v for v, k in c if k == c[0][1]]
    return sorted(top)[0]


def nearest_populations(person, panel_rows, ancients, evals, n_pcs=DISTANCE_PCS, top=TOP):
    """Rank the ancient groups by the weighted distance from the person to the
    group's mean position on the first n_pcs PCs."""
    w = weights(evals, n_pcs)
    groups = collections.defaultdict(list)
    for gid, r in ancients.items():
        if gid in panel_rows and panel_rows[gid][0] == ANCIENT_LABEL:
            groups[r["group"]].append((panel_rows[gid][1], r))
    rows = []
    for g, members in groups.items():
        k = len(members[0][0])
        centre = [statistics.fmean(m[0][i] for m in members) for i in range(min(k, n_pcs))]
        dates = [to_float(m[1]["date_bp"]) for m in members]
        dates = [d for d in dates if d is not None]
        rows.append({
            "group": g,
            "n": len(members),
            "distance": distance(person[:n_pcs], centre, w),
            "pc1": centre[0],
            "pc2": centre[1] if len(centre) > 1 else 0.0,
            "date_bp": statistics.median(dates) if dates else None,
            "dates": range_label(max(dates), min(dates)) if dates else ".",
            "country": most_common(m[1]["country"] for m in members),
            "locality": most_common(m[1]["locality"] for m in members),
        })
    rows.sort(key=lambda r: (r["distance"], r["group"]))
    for i, r in enumerate(rows, 1):
        r["rank"] = i
    return rows[:top], len(rows)


AGE_BINS = (  # (oldest BP, colour, legend), oldest first; dark-background palette
    (math.inf, "#5aa9e6", "before 8000 BCE"),
    (9950, "#7fc97f", "8000-4000 BCE"),
    (5950, "#f0a35e", "4000-1000 BCE"),
    (2950, "#e377c2", "after 1000 BCE"),
)


def age_colour(bp):
    colour = AGE_BINS[0][1]
    for limit, c, _ in AGE_BINS:
        if bp is not None and bp < limit:
            colour = c
    return colour


def svg_plot(person, panel_rows, ancients, sample, nearest, width=900, height=720):
    """The PC1/PC2 plane: present-day people grey, ancient people coloured by
    age, the person a large ring. Dark background, so it reads in a dark
    report and on its own."""
    pts = [v[1][:2] for v in panel_rows.values()] + [person[:2]]
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    pad = 0.04
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)
    dx, dy = (x1 - x0) or 1.0, (y1 - y0) or 1.0
    left, right, top, bottom = 60, 20, 50, 110
    pw, ph = width - left - right, height - top - bottom

    def X(v):
        return left + ((v - x0) / dx * (1 - 2 * pad) + pad) * pw

    def Y(v):
        return top + (1 - ((v - y0) / dy * (1 - 2 * pad) + pad)) * ph

    e = html.escape
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
           f'viewBox="0 0 {width} {height}" font-family="sans-serif" font-size="13">',
           '<style>.m{fill:#9a9a9a;fill-opacity:.55}.a{fill-opacity:.7}text{fill:#e6e6e6}</style>',
           f'<rect width="{width}" height="{height}" fill="#15171c"/>',
           f'<text x="{left}" y="28" font-size="16">{e(sample)} among ancient West Eurasians '
           '(AADR, projected onto present-day PCs)</text>',
           f'<rect x="{left}" y="{top}" width="{pw}" height="{ph}" fill="none" stroke="#444"/>']
    for label, pcs in panel_rows.values():
        if label != ANCIENT_LABEL:
            out.append(f'<circle class="m" cx="{X(pcs[0]):.1f}" cy="{Y(pcs[1]):.1f}" r="1.6"/>')
    for gid, (label, pcs) in panel_rows.items():
        if label == ANCIENT_LABEL and gid in ancients:
            c = age_colour(to_float(ancients[gid]["date_bp"]))
            out.append(f'<circle class="a" cx="{X(pcs[0]):.1f}" cy="{Y(pcs[1]):.1f}" r="2" fill="{c}"/>')
    for r in nearest[:5]:
        out.append(f'<text x="{X(r["pc1"]) + 6:.1f}" y="{Y(r["pc2"]) - 4:.1f}" font-size="11">{e(r["group"])}</text>')
    px, py = X(person[0]), Y(person[1])
    out.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="9" fill="none" stroke="#ffffff" stroke-width="3"/>')
    out.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="3" fill="#ffd400"/>')
    out.append(f'<text x="{px + 12:.1f}" y="{py + 4:.1f}" font-weight="bold">{e(sample)}</text>')
    out.append(f'<text x="{left + pw / 2:.0f}" y="{top + ph + 22}" text-anchor="middle">PC1</text>')
    out.append(f'<text x="18" y="{top + ph / 2:.0f}" transform="rotate(-90 18 {top + ph / 2:.0f})" '
               'text-anchor="middle">PC2</text>')
    lx, ly = left, top + ph + 48
    out.append(f'<circle class="m" cx="{lx + 5}" cy="{ly - 4}" r="4"/><text x="{lx + 14}" y="{ly}">present day</text>')
    lx += 120
    for _, c, name in AGE_BINS:
        out.append(f'<circle cx="{lx + 5}" cy="{ly - 4}" r="4" fill="{c}"/><text x="{lx + 14}" y="{ly}">{e(name)}</text>')
        lx += 150
    out.append(f'<text x="{left}" y="{ly + 26}" font-size="12">{e(CAVEAT)}</text>')
    out.append("</svg>")
    return "\n".join(out) + "\n"


def cmd_report(a):
    evals, panel = read_evec(a.panel_evec)
    _, run = read_evec(a.run_evec)
    if a.sample not in run or run[a.sample][0] != PERSON_LABEL:
        raise Fail(f"{a.run_evec} has no row for {a.sample} labelled {PERSON_LABEL}")
    n_pcs = min(DISTANCE_PCS, len(evals), len(run[a.sample][1]))
    signs = axis_signs(panel, run, max(2, n_pcs))
    person = [s * v for s, v in zip(signs, run[a.sample][1])]
    ancients = read_ancients(a.ancients)
    nearest, n_groups = nearest_populations(person, panel, ancients, evals, n_pcs)
    stats = {}
    if a.stats:
        with open(a.stats) as f:
            for line in f:
                k, _, v = line.rstrip("\n").partition("\t")
                stats[k] = v
    n_anc = sum(1 for lab, _ in panel.values() if lab == ANCIENT_LABEL)
    n_mod = sum(1 for lab, _ in panel.values() if lab not in (ANCIENT_LABEL, PERSON_LABEL))
    with open(a.out_tsv + ".tmp", "w") as o:
        o.write("rank\tgroup\tn\tdistance\tpc1\tpc2\tdate_bp\tdates\tcountry\tlocality\n")
        for r in nearest:
            o.write(f"{r['rank']}\t{clean(r['group'])}\t{r['n']}\t{r['distance']:.5f}\t{r['pc1']:.5f}\t{r['pc2']:.5f}\t"
                    f"{'.' if r['date_bp'] is None else int(round(r['date_bp']))}\t{r['dates']}\t"
                    f"{clean(r['country'])}\t{clean(r['locality'])}\n")
    os.replace(a.out_tsv + ".tmp", a.out_tsv)
    info = [
        ("release", a.release),
        ("sites_used", stats.get("sites_used", ".")),
        ("panel_sites", stats.get("panel_sites", ".")),
        ("present_day_people", n_mod),
        ("ancient_people", n_anc),
        ("ancient_groups", n_groups),
        ("pcs_in_distance", n_pcs),
        ("pc1", f"{person[0]:.5f}"),
        ("pc2", f"{person[1]:.5f}"),
    ]
    with open(a.out_info + ".tmp", "w") as o:
        o.write("key\tvalue\n")
        for k, v in info:
            o.write(f"{k}\t{v}\n")
    os.replace(a.out_info + ".tmp", a.out_info)
    with open(a.out_svg + ".tmp", "w") as o:
        o.write(svg_plot(person, panel, ancients, a.sample, nearest))
    os.replace(a.out_svg + ".tmp", a.out_svg)
    lines = [f"Ancient populations nearest to {a.sample} (AADR {a.release}, smartpca with lsqproject)",
             f"Sites used: {stats.get('sites_used', '.')} of {stats.get('panel_sites', '.')}; axes from {n_mod} "
             f"present-day West Eurasians; {n_anc} ancient people in {n_groups} groups projected.",
             f"Distance: the first {n_pcs} PCs, each weighted by the square root of its eigenvalue.", ""]
    for r in nearest[:10]:
        lines.append(f"  {r['rank']:>2}. {r['group']:<45} {r['dates']:<22} {r['country']:<20} n={r['n']:<4} "
                     f"d={r['distance']:.4f}")
    lines += ["", CAVEAT, QPADM_NOTE]
    with open(a.out_txt, "w") as o:
        o.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select")
    s.add_argument("--anno", required=True)
    s.add_argument("--ind", required=True)
    s.add_argument("--out-ind", required=True)
    s.add_argument("--poplist", required=True)
    s.add_argument("--ancients", required=True)
    s.add_argument("--min-snps", type=int, default=MIN_ANCIENT_SNPS)
    s = sub.add_parser("lift")
    s.add_argument("--snp", required=True)
    s.add_argument("--chain", required=True)
    s.add_argument("--fasta", required=True)
    s.add_argument("--out", required=True)
    s = sub.add_parser("person-geno")
    s.add_argument("--sites", required=True)
    s.add_argument("--genotypes", required=True)
    s.add_argument("--snp", required=True)
    s.add_argument("--moderns-geno", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--stats", required=True)
    s = sub.add_parser("align")
    s.add_argument("--panel-evec", required=True)
    s.add_argument("--run-evec", required=True)
    s.add_argument("--pcs", type=int, default=DISTANCE_PCS)
    s = sub.add_parser("report")
    s.add_argument("--panel-evec", required=True)
    s.add_argument("--run-evec", required=True)
    s.add_argument("--ancients", required=True)
    s.add_argument("--sample", required=True)
    s.add_argument("--release", required=True)
    s.add_argument("--stats")
    s.add_argument("--out-tsv", required=True)
    s.add_argument("--out-info", required=True)
    s.add_argument("--out-txt", required=True)
    s.add_argument("--out-svg", required=True)
    a = p.parse_args(argv)
    try:
        return {"select": cmd_select, "lift": cmd_lift, "person-geno": cmd_person_geno,
                "align": cmd_align, "report": cmd_report}[a.cmd](a)
    except Fail as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
