#!/usr/bin/env python3
"""couple_report.py: the couple carrier screen of step 39 (scripts/39-couple.sh)
and of the Nextflow COUPLE_CARRIER process (workflows/family.nf). Both run the
same command. Standard library only, so it runs in the plain python image.

  couple_report.py screen --pedigree P.tsv --genes assets/carrier_genes.tsv
                          --clinvar ID=HITS.vcf ... [--delly ID=SV.vcf.gz ...]
                          [--eh ID=EH.vcf ...] [--smn ID=SMN.tsv ...]
                          [--exons EXONS.bed] [--declared-sex ID=SEX ...] --out-dir D
      For every two samples the pedigree marks as partners and that both have
      a ClinVar screen (step 06's <id>_clinvar_hits.vcf), collects each
      partner's carrier findings in the genes of the carrier list: ClinVar
      pathogenic and likely pathogenic hits, Delly deletions that remove an
      exon (with --exons), SMN1 copy number (Parascopy, step 35) and the FMR1
      and FXN repeats (ExpansionHunter, step 09). It writes, per couple A and
      B (sorted by name):
        D/couple_A_B.tsv          one row per gene where either partner has a
                                  finding: shared (both carry, autosomal
                                  recessive), x_linked, or one_partner
        D/couple_A_B.md           the same in words, what was not assessed,
                                  and the residual risk
        D/A_couple_carrier.tsv    a copy of the table for each partner's own
        D/B_couple_carrier.tsv    report (bin/collect_summary.py reads it)
      Exits 1 when no couple could be screened, 2 on an input problem.

  couple_report.py genes --table assets/carrier_genes.tsv --orphanet-genes en_product6.xml
                         --orphanet-inheritance en_product9_ages.xml
                         --gtf gencode.v<N>.basic.annotation.gtf.gz --par-bed assets/par_grch38.bed --out NEW.tsv
      Rebuilds the carrier list: the ACMG and ACOG rows of --table are kept
      as they are written, every other row is made again from Orphanet (genes
      of disorders with autosomal recessive or X-linked inheritance), and
      every row gets its GRCh38 span from the GENCODE gene lines.
"""
import argparse
import gzip
import os
import re
import sys
import xml.etree.ElementTree as ET
from bisect import bisect_left, bisect_right

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import collect_summary  # noqa: E402  ClinVar hit rows and stars, as both reports read them
import family_report  # noqa: E402  the pedigree rules of step 38
import paraphase_report  # noqa: E402  step 35's SMN1 rule

GENE_COLUMNS = ("gene", "chrom", "start", "end", "inheritance", "lists", "condition")
TABLE_COLUMNS = ("gene", "inheritance", "lists", "condition", "result", "offspring_risk",
                 "partner_a", "findings_a", "partner_b", "findings_b", "explanation")
# Repeat loci of the carrier list ExpansionHunter's catalog sizes: the
# smallest repeat count of an allele counted as a carrier finding, and its name.
# FMR1: 55 to 200 CGG is a premutation, more than 200 a full mutation.
# FXN: 66 GAA or more is in the range that causes Friedreich ataxia.
REPEATS = {"FMR1": (55, "CGG"), "FXN": (66, "GAA")}
FMR1_FULL = 201
SVTYPE = re.compile(r"(?:^|;)SVTYPE=([^;]+)")
END = re.compile(r"(?:^|;)END=(\d+)")
# A Delly deletion longer than this needs read depth to agree: its FORMAT RDCN
# (Delly's read-depth copy number) must be below the copies a normal region
# has. 1,000 bp is Delly's own line between small and large SVs in its read
# counts (src/delly.h, indelsize). Short-read Delly can mark PASS a long
# deletion whose RDCN shows no loss; a short one rests on split reads alone.
DEPTH_MIN_SIZE = 1000


def die(lines, code=2):
    family_report.die(lines, code)


def kv_args(values, what):
    out = {}
    for v in values or []:
        k, sep, p = v.partition("=")
        if not sep or not k or not p:
            die(["--%s takes ID=PATH, got '%s'" % (what, v)])
        if k in out:
            die(["--%s names '%s' twice" % (what, k)])
        out[k] = p
    return out


# --- the carrier list --------------------------------------------------------

def read_genes(path):
    """The carrier list: {gene: row}. A list without rows stops the screen:
    an empty list would report no shared gene for every couple."""
    genes, header = {}, None
    with open(path, encoding="utf-8") as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            vals = line.split("\t")
            if header is None:
                header = vals
                if tuple(header) != GENE_COLUMNS:
                    die(["%s: the header must be %s" % (path, "\t".join(GENE_COLUMNS))])
                continue
            if len(vals) != len(header):
                die(["%s line %d: %d fields, the header has %d" % (path, n, len(vals), len(header))])
            r = dict(zip(header, vals))
            if r["inheritance"] not in ("AR", "XL") or not r["start"].isdigit() or not r["end"].isdigit():
                die(["%s line %d: inheritance must be AR or XL, start and end whole numbers" % (path, n)])
            r["start"], r["end"] = int(r["start"]), int(r["end"])
            genes[r["gene"]] = r
    if not genes:
        die(["%s holds no gene: the screen would find nothing. Restore assets/carrier_genes.tsv." % path])
    return genes


class GeneIndex:
    """Genes by chromosome, for the genes a deletion overlaps."""

    def __init__(self, genes):
        self.by = {}
        for g in genes.values():
            self.by.setdefault(g["chrom"], []).append((g["start"], g["end"], g["gene"]))
        for c in self.by:
            self.by[c].sort()
        self.longest = {c: max(e - s for s, e, _ in v) for c, v in self.by.items()}

    def overlapping(self, chrom, s, e):
        """Genes whose 1-based span [start, end] overlaps the 1-based [s, e]."""
        rows = self.by.get(chrom, [])
        lo = bisect_left(rows, (s - self.longest.get(chrom, 0) - 1,))
        hi = bisect_right(rows, (e, float("inf")))
        return [(gs, ge, g) for gs, ge, g in rows[lo:hi] if gs <= e and ge >= s]


def read_exons(path):
    """GENCODE exons as a BED (0-based start, end), merged per chromosome."""
    by = {}
    with open(path) as f:
        for line in f:
            c = line.split("\t")
            if len(c) < 3 or line.startswith(("#", "track")):
                continue
            by.setdefault(c[0], []).append((int(c[1]), int(c[2])))
    merged = {}
    for chrom, iv in by.items():
        iv.sort()
        m = []
        for s, e in iv:
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        merged[chrom] = ([s for s, _ in m], [e for _, e in m])
    return merged


def exons_hit(exons, chrom, s0, e0):
    """Merged exons overlapping the 0-based half-open [s0, e0)."""
    starts, ends = exons.get(chrom, ([], []))
    i = bisect_right(ends, s0)
    n = 0
    while i < len(starts) and starts[i] < e0:
        n += 1
        i += 1
    return n


# --- findings ----------------------------------------------------------------

def finding(gene, text, copies):
    """copies: how many of the person's copies the finding is on (1 or 2)."""
    return {"gene": gene, "text": text, "copies": copies}


def zygosity(gt, chrom, sex, par):
    """het, hom, hemizygous, or None (no alt allele)."""
    a = [x for x in re.split(r"[/|]", gt) if x not in (".", "")]
    alt = [x for x in a if x != "0"]
    if not alt:
        return None
    if chrom in ("chrX", "chrY") and sex == "male" and not par:
        return "hemizygous"
    if len(a) == 1:
        return "hemizygous" if chrom in ("chrX", "chrY") else "hom"
    return "hom" if len(alt) == len(a) else "het"


def stars_text(n):
    return "%d star%s" % (n, "" if n == 1 else "s")


def clinvar_findings(path, genes, sex, par):
    """Step 06's hits in genes of the list, and the count of all hits."""
    out, total = [], 0
    for r in collect_summary.vcf_records(path):
        if len(r) < 8:
            continue
        total += 1
        h = collect_summary.clinvar_row(r)
        gt = r[9].split(":")[0] if len(r) > 9 else ""
        z = zygosity(gt, h["chrom"], sex, family_report.in_par(par, h["chrom"], h["pos"]))
        if z is None:
            continue
        for g in h["gene"].split(","):
            if g in genes:
                ident = "" if r[2] in (".", "") else ", ClinVar %s" % r[2]
                text = "%s:%d %s>%s %s, %s, %s%s" % (h["chrom"], h["pos"], h["ref"], h["alt"], z,
                                                     h["significance"], stars_text(h["stars"]), ident)
                out.append(finding(g, text, 1 if z == "het" else 2))
    return out, total


def delly_findings(path, index, exons, sex, par):
    """PASS Delly deletions that remove at least one exon of a gene of the list,
    and the count of long ones set aside because read depth shows no loss."""
    out, no_loss = [], 0
    for r in collect_summary.vcf_records(path):
        if len(r) < 10 or r[6] != "PASS":
            continue
        m = SVTYPE.search(r[7])
        e = END.search(r[7])
        if not m or m.group(1) != "DEL" or not e:
            continue
        fmt, vals = r[8].split(":"), r[9].split(":")
        if "FT" in fmt and fmt.index("FT") < len(vals) and vals[fmt.index("FT")] not in ("PASS", "."):
            continue
        gt = vals[0] if fmt and fmt[0] == "GT" else ""
        chrom, pos, end = r[0], int(r[1]), int(e.group(1))
        in_par = family_report.in_par(par, chrom, pos)
        z = zygosity(gt, chrom, sex, in_par)
        if z is None:
            continue
        size = end - pos
        # Read depth must show the loss for a long deletion: RDCN below the
        # copies a normal region has (1 on a male's chrX or chrY outside the
        # PARs, 2 elsewhere). A missing or negative RDCN (no reads in the
        # flanks) is no evidence of a loss. A VCF without the RDCN field (not
        # Delly 2.x) is read as before.
        if size > DEPTH_MIN_SIZE and "RDCN" in fmt:
            i = fmt.index("RDCN")
            cn = vals[i] if i < len(vals) else "."
            normal = 1 if chrom in ("chrX", "chrY") and sex == "male" and not in_par else 2
            if not (cn.isdigit() and int(cn) < normal):
                if any(exons_hit(exons, chrom, max(pos, gs - 1), min(end, ge))
                       for gs, ge, _ in index.overlapping(chrom, pos + 1, end)):
                    no_loss += 1
                continue
        # The deleted bases are POS+1 to END (POS is the base before them).
        for gs, ge, g in index.overlapping(chrom, pos + 1, end):
            n = exons_hit(exons, chrom, max(pos, gs - 1), min(end, ge))
            if not n:
                continue
            text = "deletion %s:%d-%d (%s) %s, %d exon%s of the gene (Delly)" % (
                chrom, pos + 1, end, ("%.1f kb" % (size / 1000.0)) if size >= 1000 else "%d bp" % size,
                z, n, "" if n == 1 else "s")
            out.append(finding(g, text, 1 if z == "het" else 2))
    return out, no_loss


def eh_findings(path, genes):
    """FMR1 and FXN alleles in their expanded range (ExpansionHunter)."""
    out = []
    for r in collect_summary.vcf_records(path):
        if len(r) < 10:
            continue
        info = collect_summary.info_map(r[7])
        locus = info.get("VARID") or info.get("REPID") or ""
        if locus not in REPEATS or locus not in genes:
            continue
        fmt, vals = r[8].split(":"), r[9].split(":")
        if "REPCN" not in fmt or fmt.index("REPCN") >= len(vals):
            continue
        cn = vals[fmt.index("REPCN")]
        sizes = [int(x) for x in cn.split("/") if x.isdigit()]
        least, unit = REPEATS[locus]
        big = [x for x in sizes if x >= least]
        if not big:
            continue
        if locus == "FMR1":
            what = "full mutation (more than 200)" if max(big) >= FMR1_FULL else "premutation (55 to 200)"
        else:
            what = "expanded (66 or more)"
        text = "%s %s repeat %s, %s allele%s %s (ExpansionHunter estimate)" % (
            locus, unit, cn, "two" if len(big) > 1 else "one", "s" if len(big) > 1 else "", what)
        out.append(finding(locus, text, len(big)))
    return out


def smn_finding(path, genes):
    """SMN1 copy number (step 35): ([findings], text of what was read)."""
    p = paraphase_report.parascopy_smn1(path)
    if p["cn"] is None or "SMN1" not in genes:
        return [], "SMN1 copy number %s" % p["text"]
    if p["cn"] >= 2:
        return [], "SMN1 %d copies" % p["cn"]
    text = "SMN1 %d cop%s (Parascopy: %s)" % (p["cn"], "y" if p["cn"] == 1 else "ies", p["text"])
    return [finding("SMN1", text, 1 if p["cn"] == 1 else 2)], "SMN1 %d cop%s" % (p["cn"], "y" if p["cn"] == 1 else "ies")


# --- the screen --------------------------------------------------------------

def odds(p):
    return {0.25: "1 in 4", 0.5: "1 in 2", 1.0: "every"}.get(p, "%.2f" % p)


def ar_row(g, a, b, fa, fb):
    ta = 1.0 if any(x["copies"] == 2 for x in fa) else 0.5
    tb = 1.0 if any(x["copies"] == 2 for x in fb) else 0.5
    p = ta * tb
    risk = "%s pregnancy" % odds(p) if p == 1.0 else "%s per pregnancy" % odds(p)
    more = [n for n, f in ((a, fa), (b, fb)) if len(f) > 1 and all(x["copies"] == 1 for x in f)]
    why = ("Both partners carry a change in %s, linked to %s (autosomal recessive). A child who inherits "
           "the change of each parent can have the condition: %s." % (g["gene"], g["condition"], risk))
    if more:
        why += (" %s %s two changes in the gene; when they sit on different copies, every child gets one, "
                "and the chance is higher." % (" and ".join(more), "has" if len(more) == 1 else "have"))
    return "shared", risk, why


def xl_lines(g, who, sex, f):
    two = any(x["copies"] == 2 for x in f)
    fmr1 = g["gene"] == "FMR1" and any("premutation" in x["text"] for x in f)
    if sex == "female":
        if two:
            risk = "%s: every son, every daughter a carrier" % who
            why = ("%s carries the change on both X chromosomes: every son inherits it and can have %s; "
                   "every daughter is a carrier." % (who, g["condition"]))
        else:
            risk = "%s: sons 1 in 2" % who
            why = ("%s carries a change in %s on one X chromosome (X-linked: %s). Each son has a 1 in 2 "
                   "chance of inheriting it and can have the condition; each daughter has a 1 in 2 chance of "
                   "being a carrier." % (who, g["gene"], g["condition"]))
        if fmr1:
            why += " A premutation can grow into a full mutation when a mother passes it on."
    elif sex == "male":
        risk = "%s: daughters carriers" % who
        why = ("%s has a change in %s on his X chromosome (X-linked: %s). Every daughter inherits it and "
               "is a carrier; a son gets his Y chromosome, not this X." % (who, g["gene"], g["condition"]))
    else:
        risk = "%s: sex unknown" % who
        why = ("%s has a change in %s (X-linked: %s). The pedigree gives no sex for %s, so the chance for "
               "sons and daughters cannot be said." % (who, g["gene"], g["condition"], who))
    return risk, why


def screen_couple(a, b, sexes, found, genes):
    """Rows of the couple table, shared genes first, then X-linked, then the rest."""
    rows = []
    for name in sorted(set(x["gene"] for x in found[a] + found[b])):
        g = genes[name]
        fa = [x for x in found[a] if x["gene"] == name]
        fb = [x for x in found[b] if x["gene"] == name]
        if g["inheritance"] == "XL":
            parts = [xl_lines(g, w, sexes.get(w), f) for w, f in ((a, fa), (b, fb)) if f]
            result, risk, why = "x_linked", "; ".join(p[0] for p in parts), " ".join(p[1] for p in parts)
        elif fa and fb:
            result, risk, why = ar_row(g, a, b, fa, fb)
        else:
            who = a if fa else b
            other = b if fa else a
            result, risk = "one_partner", "none shared"
            why = ("Only %s carries a change in %s (%s, autosomal recessive). A child can be a carrier like "
                   "%s, not affected, unless %s carries a change this screen does not see."
                   % (who, name, g["condition"], who, other))
        rows.append({"gene": name, "inheritance": g["inheritance"], "lists": g["lists"],
                     "condition": g["condition"], "result": result, "offspring_risk": risk,
                     "partner_a": a, "findings_a": "; ".join(x["text"] for x in fa) or ".",
                     "partner_b": b, "findings_b": "; ".join(x["text"] for x in fb) or ".",
                     "explanation": why})
    order = {"shared": 0, "x_linked": 1, "one_partner": 2}
    rows.sort(key=lambda r: (order[r["result"]], r["gene"]))
    return rows


def cell(v):
    return v.replace("\t", " ").replace("\n", " ")


def write_table(path, rows):
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(TABLE_COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(cell(str(r[k])) for k in TABLE_COLUMNS) + "\n")


RESIDUAL = [
    "Changes ClinVar does not list as pathogenic or likely pathogenic: a new loss-of-function change, or one of "
    "uncertain significance, is not counted.",
    "Genes outside the carrier list (assets/carrier_genes.tsv).",
    "SMN1 silent carriers: two copies on one chromosome and none on the other read as two copies.",
    "Changes in genes with a near-identical copy (HBA1 and HBA2, CYP21A2, GBA1, SMN1): short reads place them "
    "poorly, and a gene conversion is not seen. A deletion there counts only when Delly calls it.",
    "Repeat expansions other than FMR1 and FXN, and the exact size of a long FMR1 or FXN allele: short reads do "
    "not span it.",
    "Deletions Delly does not call, deletions longer than %d bp that read depth (Delly's RDCN) does not confirm, "
    "duplications, inversions and insertions." % DEPTH_MIN_SIZE,
    "Low-level mosaicism, and whether two changes of one person in one gene sit on the same copy.",
]
COUNSEL = ("This is a screen, not a diagnosis, and not medical advice. Before acting on a result, or on its "
           "absence, see a genetic counsellor: a clinical carrier test confirms a finding, and the counsellor "
           "can say what a negative screen leaves out for your ancestry.")


def write_markdown(path, a, b, sexes, rows, assessed, gaps, ngenes):
    def sx(s):
        return sexes.get(s) or "sex unknown"
    out = ["# Couple carrier screen: %s and %s" % (a, b), "",
           "%s (%s) and %s (%s), partners in the pedigree. The screen covers %d genes of autosomal recessive "
           "and X-linked conditions (assets/carrier_genes.tsv: the ACMG and ACOG carrier screening genes and the "
           "recessive and X-linked genes of Orphanet)." % (a, sx(a), b, sx(b), ngenes), "", COUNSEL, ""]

    def section(title, keep, none):
        out.extend(["## " + title, ""])
        sel = [r for r in rows if r["result"] == keep]
        if not sel:
            out.extend([none, ""])
            return
        for r in sel:
            out.append("- **%s** (%s; %s): %s" % (r["gene"], r["condition"], r["lists"], r["offspring_risk"]))
            for p, k in ((r["partner_a"], "findings_a"), (r["partner_b"], "findings_b")):
                if r[k] != ".":
                    out.append("  - %s: %s" % (p, r[k]))
            out.append("  - " + r["explanation"])
        out.append("")

    section("Shared risk: both partners carry a change in the same gene", "shared",
            "No gene of the list where both partners carry a change.")
    section("X-linked findings", "x_linked", "No change in an X-linked gene of the list.")
    section("Genes where only one partner carries a change", "one_partner", "None.")
    out.extend(["## What was assessed", "", "| Partner | ClinVar hits (in listed genes / all) | Deletions | SMN1 | FMR1 and FXN repeats |",
                "|---|---|---|---|---|"])
    for s in (a, b):
        x = assessed[s]
        out.append("| %s | %s | %s | %s | %s |" % (s, x["clinvar"], x["delly"], x["smn"], x["eh"]))
    out.extend(["", "## What this screen does not see", ""])
    out.extend("- " + g for g in gaps + RESIDUAL)
    out.extend(["", "A screen that finds no shared gene lowers the chance of a child with one of these conditions; "
                "it does not remove it.", ""])
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out))


def cmd_screen(a):
    genes = read_genes(a.genes)
    index = GeneIndex(genes)
    rows, errors = family_report.read_pedigree(a.pedigree)
    declared = {}
    for v in a.declared_sex or []:
        k, _, s = v.partition("=")
        if s.strip().lower() in family_report.SEX and family_report.SEX[s.strip().lower()]:
            declared[k] = family_report.SEX[s.strip().lower()]
    members, more = family_report.build_members(rows, declared, a.pedigree)
    if errors or more:
        die(errors + more)
    clinvar = kv_args(a.clinvar, "clinvar")
    delly = kv_args(a.delly, "delly")
    eh = kv_args(a.eh, "eh")
    smn = kv_args(a.smn, "smn")
    for what, d in (("clinvar", clinvar), ("delly", delly), ("eh", eh), ("smn", smn)):
        for k, p in d.items():
            if k not in members:
                die(["--%s %s: '%s' is not in the pedigree" % (what, k, k)])
            if not os.path.isfile(p):
                die(["--%s %s: no file %s" % (what, k, p)])
    exons = read_exons(a.exons) if a.exons else None
    par = family_report.read_par(a.par_bed)
    sexes = {s: m["sex"] for s, m in members.items()}

    couples = sorted({tuple(sorted((s, m["partner"]))) for s, m in members.items() if m["partner"]})
    if not couples:
        die(["%s marks no two people as partners (the partner column): there is no couple to screen." % a.pedigree], 1)
    os.makedirs(a.out_dir, exist_ok=True)
    done = 0
    for x, y in couples:
        lack = [s for s in (x, y) if s not in clinvar]
        if lack:
            print("%s and %s: not screened, %s %s no ClinVar screen (step 06's <sample>_clinvar_hits.vcf)"
                  % (x, y, " and ".join(lack), "has" if len(lack) == 1 else "have"))
            continue
        found, assessed, gaps = {}, {}, []
        for s in (x, y):
            f, total = clinvar_findings(clinvar[s], genes, sexes[s], par)
            info = {"clinvar": "%d / %d" % (len(f), total)}
            if s not in delly:
                info["delly"] = "not assessed (no Delly calls, step 19)"
                gaps.append("Deletions of %s: no Delly calls (step 19)." % s)
            elif exons is None:
                info["delly"] = "not assessed (no exon list)"
                gaps.append("Deletions of %s: Delly ran, but without the GENCODE exon list a deletion cannot be "
                            "placed on exons." % s)
            else:
                d, no_loss = delly_findings(delly[s], index, exons, sexes[s], par)
                info["delly"] = "%d over exons of listed genes" % len(d)
                if no_loss:
                    info["delly"] += "; %d longer than %d bp set aside, read depth shows no loss" % (
                        no_loss, DEPTH_MIN_SIZE)
                f += d
            if s in smn:
                m, txt = smn_finding(smn[s], genes)
                info["smn"] = txt
                f += m
                if not m and ("not resolved" in txt or "no row" in txt):
                    gaps.append("SMN1 copy number of %s: %s." % (s, txt))
            else:
                info["smn"] = "not assessed (no Parascopy, step 35)"
                gaps.append("SMN1 copy number of %s: no Parascopy (step 35)." % s)
            if s in eh:
                r = eh_findings(eh[s], genes)
                info["eh"] = "%d expanded" % len(r)
                f += r
            else:
                info["eh"] = "not assessed (no ExpansionHunter, step 09)"
                gaps.append("FMR1 and FXN repeats of %s: no ExpansionHunter (step 09)." % s)
            found[s], assessed[s] = f, info
        rows = screen_couple(x, y, sexes, found, genes)
        base = os.path.join(a.out_dir, "couple_%s_%s" % (x, y))
        write_table(base + ".tsv", rows)
        for s in (x, y):
            write_table(os.path.join(a.out_dir, "%s_couple_carrier.tsv" % s), rows)
        write_markdown(base + ".md", x, y, sexes, rows, assessed, gaps, len(genes))
        n = {k: sum(1 for r in rows if r["result"] == k) for k in ("shared", "x_linked", "one_partner")}
        print("%s and %s: %d shared gene%s, %d X-linked, %d with one partner only (%s.tsv)" % (
            x, y, n["shared"], "" if n["shared"] == 1 else "s", n["x_linked"], n["one_partner"], base))
        for r in rows:
            if r["result"] == "shared":
                print("  shared: %s (%s): %s" % (r["gene"], r["condition"], r["offspring_risk"]))
        done += 1
    if not done:
        die(["no couple of the pedigree has a ClinVar screen for both partners: run step 06 for each"], 1)


# --- the carrier list, rebuilt ---------------------------------------------------

ORPHANET_KEEP = {"Disease-causing germline mutation(s) in",
                 "Disease-causing germline mutation(s) (loss of function) in",
                 "Disease-causing germline mutation(s) (gain of function) in"}
ORPHANET_MODE = {"Autosomal recessive": "AR", "X-linked recessive": "XL", "X-linked dominant": "XL"}
PRIMARY = ["chr%s" % c for c in list(range(1, 23)) + ["X", "Y"]]


def gencode_genes(path):
    """{gene_name: (chrom, start, end)} of GENCODE's gene lines (1-based), and
    the release. A name on several chromosomes keeps the first of chr1 to
    chrY in that order (a PAR gene: chrX)."""
    found, release = {}, "?"
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                m = re.search(r"version (\d+)", line)
                if m:
                    release = m.group(1)
                continue
            c = line.split("\t", 9)
            if len(c) < 9 or c[2] != "gene" or c[0] not in PRIMARY:
                continue
            m = re.search(r'gene_name "([^"]+)"', c[8])
            if m:
                found.setdefault(m.group(1), []).append((c[0], int(c[3]), int(c[4])))
    out = {}
    for g, locs in found.items():
        locs.sort(key=lambda x: (PRIMARY.index(x[0]), x[1]))
        out[g] = locs[0]
    return out, release


def cmd_genes(a):
    head, panel = [], []
    with open(a.table, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            panel.append(line.rstrip("\r\n").split("\t"))
    cols, panel = panel[0], [dict(zip(panel[0], r)) for r in panel[1:] if r and r[0]]
    if tuple(cols) != GENE_COLUMNS:
        die(["%s: the header must be %s" % (a.table, "\t".join(GENE_COLUMNS))])
    panel = [r for r in panel if "ACMG" in r["lists"] or "ACOG" in r["lists"]]
    gtf, release = gencode_genes(a.gtf)
    upper = {g.upper(): g for g in gtf}
    par = family_report.read_par(a.par_bed)

    def place(g):
        return gtf.get(g) or gtf.get(upper.get(g.upper(), ""))

    missing = [r["gene"] for r in panel if not place(r["gene"])]
    if missing:
        die(["GENCODE has no gene line for %s; give the current HGNC symbol" % ", ".join(missing)])

    inh = {}
    p9 = ET.parse(a.orphanet_inheritance).getroot()
    for d in p9.iter("Disorder"):
        inh[d.findtext("OrphaCode")] = {ORPHANET_MODE[t.findtext("Name")] for t in d.iter("TypeOfInheritance")
                                        if t.findtext("Name") in ORPHANET_MODE}
    p6 = ET.parse(a.orphanet_genes).getroot()
    disorders = {}
    for d in p6.iter("Disorder"):
        modes = inh.get(d.findtext("OrphaCode"), set())
        if not modes:
            continue
        for x in d.iter("DisorderGeneAssociation"):
            if (x.findtext("DisorderGeneAssociationType/Name") not in ORPHANET_KEEP
                    or x.findtext("DisorderGeneAssociationStatus/Name") != "Assessed"):
                continue
            disorders.setdefault(x.findtext("Gene/Symbol"), []).append((d.findtext("Name"), modes))

    rows = {}
    for r in panel:
        c, s, e = place(r["gene"])
        rows[r["gene"]] = dict(r, chrom=c, start=str(s), end=str(e))
    for sym, ds in disorders.items():
        loc = place(sym)
        if not loc:
            continue
        gene = sym if sym in gtf else upper[sym.upper()]
        if gene in rows:
            continue
        c, s, e = loc
        # The chromosome decides how a gene is passed on: X outside the PARs is
        # X-linked, an autosome or a PAR is autosomal.
        mode = "XL" if c == "chrX" and not family_report.in_par(par, c, s) else "AR"
        names = sorted({n for n, m in ds if mode in m})
        if not names or c == "chrY":
            continue
        cond = names[0] if len(names) == 1 else "%s; and %d other%s" % (names[0], len(names) - 1,
                                                                        "" if len(names) == 2 else "s")
        rows[gene] = {"gene": gene, "chrom": c, "start": str(s), "end": str(e), "inheritance": mode,
                      "lists": "Orphanet", "condition": cond.replace("\t", " ")}
    date6 = p6.get("date", "?").split()[0]
    head = [
        "# carrier_genes.tsv: the genes of the couple carrier screen (step 39, bin/couple_report.py).",
        "# One row per gene. chrom, start and end: the gene's GRCh38 span (1-based, inclusive) from the",
        "# gene lines of GENCODE release %s (basic annotation). inheritance: AR (autosomal, or in a" % release,
        "# pseudoautosomal region) or XL (chrX outside the PARs). lists: where the gene comes from:",
        "#   ACMG      Gregg et al., Genet Med 2021;23(10):1793-1806 (doi:10.1038/s41436-021-01203-z),",
        "#             tier 3: the 97 autosomal recessive genes of its tables 1 to 5 and the 16 X-linked",
        "#             genes of table 6; condition names from those tables",
        "#   ACOG      the conditions of ACOG Committee Opinions 690 and 691 (2017): cystic fibrosis, spinal",
        "#             muscular atrophy, the hemoglobinopathies, fragile X and the Ashkenazi Jewish panel",
        "#   Orphanet  every other gene of an Orphanet disorder with autosomal recessive or X-linked",
        "#             inheritance (Orphadata en_product6 and en_product9_ages of %s, CC BY 4.0," % date6,
        "#             https://www.orphadata.com), with a disease-causing germline association that",
        "#             Orphanet assessed; condition: the first such disorder by name",
        "# Rebuilt with: couple_report.py genes (docs/39-couple.md says how).",
    ]
    with open(a.out, "w", encoding="utf-8") as f:
        f.write("\n".join(head) + "\n" + "\t".join(GENE_COLUMNS) + "\n")
        for g in sorted(rows, key=lambda g: (PRIMARY.index(rows[g]["chrom"]), int(rows[g]["start"]), g)):
            f.write("\t".join(rows[g][k] for k in GENE_COLUMNS) + "\n")
    print("%d genes: %d ACMG or ACOG, %d Orphanet" % (len(rows), len(panel), len(rows) - len(panel)))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("screen")
    s.add_argument("--pedigree", required=True)
    s.add_argument("--genes", required=True)
    s.add_argument("--clinvar", action="append")
    s.add_argument("--delly", action="append")
    s.add_argument("--eh", action="append")
    s.add_argument("--smn", action="append")
    s.add_argument("--exons")
    s.add_argument("--par-bed")
    s.add_argument("--declared-sex", action="append")
    s.add_argument("--out-dir", required=True)
    g = sub.add_parser("genes")
    g.add_argument("--table", required=True)
    g.add_argument("--orphanet-genes", required=True)
    g.add_argument("--orphanet-inheritance", required=True)
    g.add_argument("--gtf", required=True)
    g.add_argument("--par-bed", required=True)
    g.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "screen":
        cmd_screen(a)
    else:
        cmd_genes(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
