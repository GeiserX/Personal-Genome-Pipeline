#!/usr/bin/env python3
"""pgx_outside_calls.py: the PharmCAT outside-call file of one sample.

  pgx_outside_calls.py --calls OUT_CALLS.tsv --consensus OUT_CONSENSUS.tsv
                       [--hla T1K_genotype.tsv] [--pypgx PYPGX_summary.tsv]
                       [--cyrius CYRIUS.tsv] [--depth-check CYP2D6_depth_check.tsv]

writes the outside-call file PharmCAT reads with -po (gene<TAB>diplotype, one
line per gene, empty when nothing is passed on) and a consensus table (Gene,
Result, Outside_call, Reason, Evidence) that says, for HLA-A, HLA-B and
CYP2D6, what each caller said and what reached PharmCAT.

  HLA-A, HLA-B  T1K's two alleles (step 08), cut to two fields (*57:01), when
                both have a quality above 0; one allele reads as homozygous,
                as T1K reports a homozygous gene.
  CYP2D6        only when pypgx (step 32) and Cyrius (step 21) give the same
                diplotype and the depth check (bin/cyp2d6_depth_check.py)
                passed. One caller alone, a disagreement, a failed or missing
                depth check: 'indeterminate', and nothing reaches PharmCAT.
                PharmCAT itself calls no CYP2D6 from a VCF (3.4.0 reports it
                with callSource NONE), so it is not a second caller here.

A CYP2D6 call of pypgx alone is never passed on. The reports show it as
'probable (one caller: pypgx)' (probable_cyp2d6) when the depth check passed
and every star-defining variant pypgx used is a PASS call at or above the GQ
floor; step 36's table still says 'indeterminate'.

  pgx_outside_calls.py pypgx-summary --outbase DIR --bam-genes "G ..." --vcf-genes "G ..."
                       --vcf VCF --out SUMMARY.tsv [--not-callable GENE ...] [--gq-floor 20]

writes step 32's summary (scripts/32-pypgx.sh and the PYPGX module run it in
the pypgx image): per gene the diplotype, phenotype, pypgx's copy-number call
and, for each star allele called, the variants that define it (pypgx's core
variants) looked up in the VCF: how many it has, the lowest GQ among them,
the weakest one with its AD, their FILTER, and 'low-confidence: check the
reads' when one has a GQ below the floor. A gene in --not-callable (GSTT1 on
a reference without ALT contigs) is written as not callable on this
reference, not as FAILED.

Used by scripts/36-pgx-consensus.sh and the PGX_CONSENSUS module (bin/ is on
the task PATH); tests/test_pgx_outside_calls.py checks it. Standard library
only, except pypgx-summary, which imports pypgx and pysam.
"""
import argparse
import csv
import os
import re
import subprocess
import sys

CONSENSUS_HEADER = ["Gene", "Result", "Outside_call", "Reason", "Evidence"]
HLA_GENES = ("HLA-A", "HLA-B")
NO_CALL = {"", "none", "n/a", "na", "failed", "indeterminate", "not called", "unknown", "."}
ONE_CALLER = "one caller only: a second caller must agree (Cyrius: setup.sh --cyrius)"
PROBABLE = "probable (one caller: pypgx)"
# Step 32's summary: the status of a gene the reference cannot hold, the
# quality flag, and the GQ below which a defining variant raises it.
NOT_CALLABLE = "not callable on this reference"
LOW_GQ = "low-confidence: check the reads"
NO_GQ = "GQ not in the VCF: check the reads"
DEFAULT_GQ_FLOOR = 20
SUMMARY_HEADER = ["Gene", "Diplotype", "Phenotype", "CNV_call", "Source", "Status",
                  "Defining_variants", "Min_GQ", "Weakest_variant", "Defining_filter", "Quality"]


# --- HLA ---------------------------------------------------------------------

def two_fields(allele, gene):
    """'HLA-B*57:01:01:02' -> '*57:01'."""
    name = allele[len(gene):] if allele.startswith(gene) else allele
    if not name.startswith("*"):
        name = "*" + name
    return ":".join(name.split(":")[:2])


def read_t1k(path):
    """{gene: [(allele, quality)]} from T1K's genotype TSV: gene, number of
    distinct alleles, then allele, abundance, quality twice ('.', 0, -1 when
    there is no second allele)."""
    out = {}
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or fields[0] not in HLA_GENES:
                continue
            alleles = []
            for i in (2, 5):
                if i + 2 < len(fields) and fields[i] not in ("", "."):
                    try:
                        q = float(fields[i + 2])
                    except ValueError:
                        q = -1.0
                    alleles.append((fields[i], q))
            out[fields[0]] = alleles
    return out


def hla_row(gene, t1k):
    """(consensus row, outside-call diplotype or None) for one HLA gene."""
    if t1k is None:
        return [gene, "not typed", "no", "HLA typing (step 08) did not run", "-"], None
    alleles = t1k.get(gene) or []
    if not alleles:
        return [gene, "not typed", "no", "T1K reported no allele", "-"], None
    evidence = "; ".join(f"{a} (quality {q:g})" for a, q in alleles)
    if any(q <= 0 for _, q in alleles):
        return [gene, "indeterminate", "no",
                "T1K gives an allele quality 0 or below (T1K says to ignore it)", evidence], None
    names = [two_fields(a, gene) for a, _ in alleles]
    reason = "T1K (step 08)"
    if len(names) == 1:
        names *= 2
        reason += ": one allele, read as homozygous"
    dip = "/".join(names)
    return [gene, dip, "yes", reason, evidence], dip


# --- CYP2D6 ------------------------------------------------------------------

def is_call(diplotype):
    d = (diplotype or "").strip()
    return d.lower() not in NO_CALL and "/" in d and ";" not in d


def same_diplotype(a, b):
    """Equal up to the order of the two alleles."""
    def key(d):
        return sorted(x.strip() for x in d.split("/"))
    return key(a) == key(b)


def read_pypgx(path):
    """pypgx's CYP2D6 diplotype ('' when its summary has no CYP2D6 row)."""
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if (row.get("Gene") or "").strip() == "CYP2D6":
                return (row.get("Diplotype") or "").strip()
    return ""


def read_cyrius(path):
    """(genotype, filter) of the first sample row of Cyrius's TSV."""
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            return (row.get("Genotype") or "").strip(), (row.get("Filter") or "").strip()
    return "", ""


def read_depth_check(path):
    """(status, message) from bin/cyp2d6_depth_check.py's TSV."""
    kv = {}
    with open(path) as f:
        for line in f:
            k, _, v = line.rstrip("\n").partition("\t")
            kv[k] = v
    return kv.get("status", ""), kv.get("message", "")


def cyp2d6_row(pypgx_path, cyrius_path, depth_path):
    """(consensus row, outside-call diplotype or None) for CYP2D6."""
    pg = read_pypgx(pypgx_path) if pypgx_path else None
    cy = read_cyrius(cyrius_path) if cyrius_path else None
    dc = read_depth_check(depth_path) if depth_path else None

    pg_ok = pg is not None and is_call(pg)
    cy_ok = cy is not None and is_call(cy[0]) and cy[1] == "PASS"
    if pg is None:
        pg_txt = "not run"
    elif pg_ok:
        pg_txt = pg
    else:
        pg_txt = f"no call ({pg or 'no CYP2D6 row'})"
    if cy is None:
        cy_txt = "not run"
    elif cy_ok:
        cy_txt = cy[0]
    else:
        cy_txt = f"no call ({', '.join(x for x in cy if x) or 'no row'})"
    dc_txt = "not run" if dc is None else (dc[0] or "unreadable")
    evidence = f"pypgx: {pg_txt}; Cyrius: {cy_txt}; depth check: {dc_txt}"

    def indeterminate(reason):
        return ["CYP2D6", "indeterminate", "no", reason, evidence], None

    if dc is not None and dc[0] != "ok":
        return indeterminate(dc[1] or "the CYP2D6 depth check failed")
    if not pg_ok and not cy_ok:
        return indeterminate("no caller made a call")
    if not (pg_ok and cy_ok):
        return indeterminate(ONE_CALLER)
    if not same_diplotype(pg, cy[0]):
        return indeterminate("pypgx and Cyrius disagree")
    if dc is None:
        return indeterminate("the CYP2D6 depth check did not run")
    return ["CYP2D6", pg, "yes", "pypgx and Cyrius agree; the depth check passed", evidence], pg


def probable_cyp2d6(row, cyrius, depth_status):
    """{'diplotype', 'phenotype', 'label'} for pypgx's CYP2D6 call when it is
    the only usable call, else None. row: pypgx's CYP2D6 summary row (a
    dict); cyrius: (genotype, filter), or None when Cyrius did not run;
    depth_status: the depth check's status, or None without its file (step 32
    turns the row into Indeterminate when the check fails, so a call in the
    row passed it).

    Not probable: no call, a usable Cyrius call (then the two callers decide),
    a depth check that did not pass, a defining variant that is missing from
    the VCF or not PASS (a summary written before these columns has none), a
    star allele pypgx's table lacks, or a defining variant with a GQ below the
    floor or no GQ."""
    if not row or not is_call(row.get("Diplotype")):
        return None
    if cyrius is not None and is_call(cyrius[0]) and cyrius[1] == "PASS":
        return None
    if depth_status is not None and depth_status != "ok":
        return None
    if (row.get("Defining_filter") or "") not in ("PASS", "none") or row.get("Quality") not in ("ok", "N/A"):
        return None
    return {"diplotype": row["Diplotype"].strip(), "phenotype": (row.get("Phenotype") or "").strip(),
            "label": PROBABLE}


# --- step 32's summary ---------------------------------------------------------

def star_alleles(diplotype):
    """The star alleles of a pypgx diplotype, each once: '*2x2/*4' gives
    ['*2', '*4'] (a copy count is not an allele), '*36+*10/*4' gives
    ['*36', '*10', '*4']; a no-call gives []."""
    d = (diplotype or "").strip()
    if d.lower() in NO_CALL or "/" not in d:
        return []
    out = []
    for half in d.split("/"):
        for a in half.split("+"):
            a = re.sub(r"x(\d+|N)$", "", a.strip())
            if a and a not in out:
                out.append(a)
    return out


def defining_evidence(diplotype, variants_of, lookup, floor=DEFAULT_GQ_FLOOR):
    """[Defining_variants, Min_GQ, Weakest_variant, Defining_filter, Quality]
    for one call. variants_of(allele): the variants pypgx defines the allele
    by ('22-42128945-C-T'), [] for one defined by none (*1, a deletion), None
    for one pypgx's table lacks. lookup(variant): {'filter', 'gq', 'ad'} of
    the VCF record, or None when the VCF has none."""
    if not star_alleles(diplotype):
        return ["N/A", "N/A", "N/A", "N/A", "N/A"]
    variants, unknown = [], []
    for a in star_alleles(diplotype):
        vs = variants_of(a)
        if vs is None:
            unknown.append(a)
            continue
        for v in vs:
            if v not in variants:
                variants.append(v)
    not_in_table = [f"not in pypgx's table: {', '.join(unknown)}"] if unknown else []
    if not variants:
        return (["unknown", "N/A", "N/A", not_in_table[0], "N/A"] if unknown
                else ["none", "N/A", "N/A", "none", "N/A"])
    found, missing, not_pass = [], [], []
    for v in variants:
        rec = lookup(v)
        name = "chr{}:{} {}>{}".format(*v.replace("chr", "", 1).split("-", 3))
        if rec is None:
            missing.append(name)
            continue
        gq = rec.get("gq")
        found.append((int(gq) if isinstance(gq, (int, float)) else None, name, rec.get("ad") or "."))
        if rec.get("filter") not in ("PASS", "."):
            not_pass.append(f"{name} {rec.get('filter')}")
    with_gq = [f for f in found if isinstance(f[0], int)]
    weakest = min(with_gq) if with_gq else None
    filt = "PASS" if not (missing or not_pass or unknown) else "; ".join(
        ([f"not PASS: {', '.join(not_pass)}"] if not_pass else [])
        + ([f"not in the VCF: {', '.join(missing)}"] if missing else []) + not_in_table)
    low = weakest is not None and weakest[0] < floor
    quality = LOW_GQ if low else NO_GQ if len(with_gq) < len(found) else "ok"
    return [f"{len(found)}/{len(variants)} in the VCF",
            str(weakest[0]) if weakest else "N/A",
            f"{weakest[1]} GQ {weakest[0]} AD {weakest[2]}" if weakest else "N/A",
            filt, quality]


def vcf_lookup(path):
    """lookup(variant) for defining_evidence over a bgzipped, indexed VCF
    (pysam), on its first sample."""
    import pysam
    vf = pysam.VariantFile(path)
    contigs = set(vf.header.contigs)

    def lookup(v):
        chrom, pos, ref, alt = v.split("-", 3)
        chrom = chrom.replace("chr", "", 1)
        c = f"chr{chrom}" if f"chr{chrom}" in contigs else chrom if chrom in contigs else None
        if c is None:
            return None
        pos = int(pos)
        for rec in vf.fetch(c, pos - 1, pos):
            if rec.pos != pos or rec.ref != ref or alt not in (rec.alts or ()):
                continue
            smp = rec.samples[0] if len(rec.samples) else {}
            gq, ad = smp.get("GQ"), smp.get("AD")
            return {"filter": ";".join(rec.filter.keys()) or ".",
                    "gq": int(gq) if isinstance(gq, (int, float)) else None,
                    "ad": ",".join("." if x is None else str(x) for x in ad) if ad else "."}
        return None
    return lookup


def pypgx_summary_main(argv):
    p = argparse.ArgumentParser(prog="pgx_outside_calls.py pypgx-summary",
                                description="Write step 32's summary from pypgx's results (pypgx image).")
    p.add_argument("--outbase", required=True, help="the folder with <gene>/results.zip")
    p.add_argument("--bam-genes", required=True)
    p.add_argument("--vcf-genes", required=True)
    p.add_argument("--vcf", required=True, help="the VCF pypgx read (bgzipped, indexed)")
    p.add_argument("--out", required=True)
    p.add_argument("--not-callable", action="append", default=[], help="a gene this reference cannot hold")
    p.add_argument("--gq-floor", type=int, default=DEFAULT_GQ_FLOOR)
    a = p.parse_args(argv)
    bam_genes, vcf_genes = a.bam_genes.split(), a.vcf_genes.split()
    import pypgx
    # pypgx writes Indeterminate for a gene it has no phenotype rule for (no
    # PhenotypeMethod in its own gene table): say that instead, so the label
    # is not read as an undecided call.
    try:
        no_rule = {g for g in bam_genes + vcf_genes if not pypgx.has_phenotype(g)}
    except Exception as e:
        print(f"WARNING: could not read pypgx gene table ({e}); phenotypes left as pypgx wrote them", file=sys.stderr)
        no_rule = set()
    try:
        lookup = vcf_lookup(a.vcf)
    except Exception as e:
        print(f"WARNING: cannot read {a.vcf} ({e}); the defining variants are not looked up", file=sys.stderr)
        lookup = None
    rows = []
    for gene in bam_genes + vcf_genes:
        source = "BAM" if gene in bam_genes else "VCF"
        if gene in a.not_callable:
            rows.append([gene, "N/A", "N/A", "N/A", source, NOT_CALLABLE] + ["N/A"] * 5)
            continue
        results_zip = os.path.join(a.outbase, gene, "results.zip")
        if not os.path.isfile(results_zip):
            rows.append([gene, "FAILED", "N/A", "N/A", source, "failed"] + ["N/A"] * 5)
            continue
        diplotype = phenotype = cnv = "N/A"
        # pypgx print-data results.zip prints a TSV whose first column is the
        # sample and whose named columns include Genotype, Phenotype and CNV
        # (the copy-number call pypgx made from read depth, e.g. Normal or
        # WholeDel1).
        out = subprocess.run(["pypgx", "print-data", results_zip], capture_output=True, text=True)
        if out.returncode != 0:
            print(f"WARNING: pypgx print-data failed for {gene}: {out.stderr.strip()}", file=sys.stderr)
        else:
            lines = out.stdout.rstrip().split("\n")
            if len(lines) >= 2:
                col = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
                diplotype = col.get("Genotype") or "N/A"
                phenotype = col.get("Phenotype") or "N/A"
                if gene in no_rule:
                    phenotype = "no pypgx phenotype rule"
                # Copy number is called from read depth, so it is meaningful only
                # for the BAM-based genes: pypgx's own value, never a guess.
                if source == "BAM":
                    cnv = col.get("CNV")
                    if cnv is None:
                        print(f"WARNING: {gene} results have no CNV column", file=sys.stderr)
                        cnv = "N/A"

        def variants_of(allele, gene=gene):
            try:
                return pypgx.list_variants(gene, alleles=[allele], mode="core", assembly="GRCh38")
            except Exception:
                return None
        evidence = (defining_evidence(diplotype, variants_of, lookup, a.gq_floor) if lookup
                    else ["not looked up"] + ["N/A"] * 4)
        rows.append([gene, diplotype, phenotype, cnv, source, "called" if is_call(diplotype) else "no call"] + evidence)
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(SUMMARY_HEADER)
        w.writerows(rows)
    print(f"Genes called: {sum(1 for r in rows if r[5] == 'called')}/{len(rows)}")
    low = [f"{r[0]} {r[1]} ({r[8]})" for r in rows if r[10] == LOW_GQ]
    if low:
        print(f"WARNING: star alleles resting on a variant below GQ {a.gq_floor}, check the reads: {'; '.join(low)}")
    return 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ["pypgx-summary"]:
        return pypgx_summary_main(argv[1:])
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--calls", required=True, help="outside-call file to write (PharmCAT -po)")
    p.add_argument("--consensus", required=True, help="consensus table to write")
    p.add_argument("--hla", help="T1K genotype TSV (step 08 / HLA_TYPING)")
    p.add_argument("--pypgx", help="pypgx summary TSV (step 32 / PYPGX)")
    p.add_argument("--cyrius", help="Cyrius TSV (step 21 / CYRIUS)")
    p.add_argument("--depth-check", help="CYP2D6 depth check TSV (bin/cyp2d6_depth_check.py)")
    args = p.parse_args(argv)

    try:
        t1k = read_t1k(args.hla) if args.hla else None
        rows, calls = [], []
        for gene in HLA_GENES:
            row, dip = hla_row(gene, t1k)
            rows.append(row)
            if dip:
                calls.append((gene, dip))
        row, dip = cyp2d6_row(args.pypgx, args.cyrius, args.depth_check)
        rows.append(row)
        if dip:
            calls.append(("CYP2D6", dip))
    except (OSError, UnicodeDecodeError, csv.Error) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    for c in calls:
        if not re.fullmatch(r"[^\t\n]+", c[1]):
            print(f"ERROR: refusing to write a malformed call {c!r}", file=sys.stderr)
            return 1
    with open(args.calls, "w") as f:
        f.writelines(f"{g}\t{d}\n" for g, d in calls)
    with open(args.consensus, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(CONSENSUS_HEADER)
        w.writerows(rows)
    for r in rows:
        print(f"{r[0]}: {r[1]} ({'passed to PharmCAT' if r[2] == 'yes' else 'not passed'}; {r[3]})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
