#!/usr/bin/env python3
"""collect_summary.py: read one sample's step outputs into one JSON summary.

Both reports (scripts/24-html-report.sh and scripts/generate-report.sh) are
rendered from this summary by bin/render_report.py, so they cannot disagree
about a count. Standard library only (it runs in the plain python image).

  collect_summary.py --sample S --sample-dir DIR [--out S/summary.json]

DIR is a bash output folder (GENOME_DIR/S) or a Nextflow one (outdir/S): both
publish names are read (vcf/ or roh/ for ROH, mosdepth/ or coverage/ for
depth, hla_t1k/ or hla/ for HLA, vcf/ or pharmcat/ for PharmCAT).

Each section has a state:
  ok          read from its file
  missing     no file: the step has not run
  failed      no file, and its step failed in the latest run-all.sh run
              (logs/run_status.tsv)
  stale       the file is older than the latest run-all.sh run and its step was
              not ok in that run (logs/run_status.tsv); the values are shown
              with the file's date
  unreadable  the file exists but could not be read (the reason is in note)

The run manifest (S/run_manifest.tsv, bin/write_manifest.sh) is copied in for
the report footer, with the ClinVar release date and the HLA database release.

  collect_summary.py sample-qc --sample S --somalier-samples F --selfsm F [...]

writes the sample identity and contamination table of step 33 (and of the
Nextflow SAMPLE_QC process): the sex somalier infers from the reads, the
declared sex, VerifyBamID2's FREEMIX against a warning threshold, and the
other samples somalier finds to be the same person. Both callers stop on the
sex_check value this writes, so the rule lives here once.

  collect_summary.py prs-format --scores DIR --labels assets/pgs_scores.tsv --out DIR --alleles F
  collect_summary.py prs-table --sample S --results DIR --sampleset NAME --scores DIR --input-kind gvcf --out F

are the two ends of step 25 (and of the Nextflow PRS_PREPARE and PRS_SUMMARY
processes) around pgsc_calc: prs-format writes each GRCh38-harmonised PGS
Catalog file as the custom GRCh38 scoring file pgsc_calc reads, labelled with
the catalog's trait, and lists every effect and other allele so the score
positions can be genotyped from the gVCF; prs-table reads pgsc_calc's match
summary and scores into <sample>_prs_summary.tsv, with the ancestry-adjusted
percentile and the reference group when pgsc_calc ran with an ancestry panel,
and then writes the principal components and population of step 26.
"""
import argparse
import csv
import glob
import gzip
import json
import os
import re
import sys
import zlib
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pgx_outside_calls  # noqa: E402  step 36's rule for a CYP2D6 call and an HLA allele

# 2: cyp2d6.consensus, Cyrius's Filter applied, HLA quality per allele,
# Stranger's status on the repeat loci, the CPIC unclassified count, and the
# section state 'failed'.
SCHEMA_VERSION = 2

# ClinVar review status -> stars. bin/clinvar_hits.awk holds the same table
# (tests/test_collect_summary.py checks that the two agree).
STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
}

# Steps run-all.sh runs and records in logs/run_status.tsv. A section whose
# step is not here (variant calling, step 03, which run-all reuses on purpose)
# is never marked stale.
RUN_ALL_STEPS = {"04", "05", "06", "07", "08", "09", "09b", "10", "11", "12", "13", "14",
                 "15", "16", "16b", "17", "18", "19", "20", "21", "22", "23", "25", "26",
                 "27", "29", "30", "31", "32", "37", "04b"}

# What this pipeline does not assess, whatever ran.
NOT_ASSESSED = [
    "Mosaic and low-fraction variants (only chrM heteroplasmy is reported)",
    "Methylation and imprinting",
    "Repeat expansions outside the ExpansionHunter catalog",
    "Phase: two variants in one gene are compound-het candidates, not confirmed",
    "Variants of uncertain significance (not interpreted)",
]

# Added to NOT_ASSESSED unless step 25 ran with an ancestry reference panel.
PRS_NOT_ADJUSTED = "Polygenic score percentiles: no ancestry reference panel was used, so the scores are raw sums"

EH_LOCI = ["HTT", "FMR1", "C9ORF72", "ATXN1", "DMPK"]

# ACMG SF v3.3 (Lee et al., Genet Med 2025;27(8)): the 84 genes in which the
# ACMG recommends reporting pathogenic and likely pathogenic variants found
# by chance, as ClinGen lists them (search.clinicalgenome.org/kb/genes/acmgsf).
# A new list version means a new set here and a new ACMG_SF_VERSION.
ACMG_SF_VERSION = "ACMG SF v3.3"
ACMG_SF_GENES = frozenset([
    "ABCD1", "ACTA2", "ACTC1", "ACVRL1", "APC", "APOB", "ATP7B", "BAG3", "BMPR1A", "BRCA1", "BRCA2",
    "BTD", "CACNA1S", "CALM1", "CALM2", "CALM3", "CASQ2", "COL3A1", "CYP27A1", "DES", "DSC2",
    "DSG2", "DSP", "ENG", "FBN1", "FLNC", "GAA", "GLA", "HFE", "HNF1A", "KCNH2", "KCNQ1", "LDLR",
    "LMNA", "MAX", "MEN1", "MLH1", "MSH2", "MSH6", "MUTYH", "MYBPC3", "MYH11", "MYH7", "MYL2",
    "MYL3", "NF2", "OTC", "PALB2", "PCSK9", "PKP2", "PLN", "PMS2", "PRKAG2", "PTEN", "RB1", "RBM20",
    "RET", "RPE65", "RYR1", "RYR2", "SCN5A", "SDHAF2", "SDHB", "SDHC", "SDHD", "SMAD3", "SMAD4",
    "STK11", "TGFBR1", "TGFBR2", "TMEM127", "TMEM43", "TNNC1", "TNNI3", "TNNT2", "TP53", "TPM1",
    "TRDN", "TSC1", "TSC2", "TTN", "TTR", "VHL", "WT1",
])


class Unreadable(Exception):
    pass


def iso(ts):
    return datetime.fromtimestamp(ts, tz=timezone.utc).strftime("%Y-%m-%d %H:%M UTC")


def open_text(path):
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", errors="replace")
    return open(path, errors="replace")


def first_existing(base, rels):
    for rel in rels:
        for p in sorted(glob.glob(os.path.join(base, rel), recursive=True)):
            if os.path.isfile(p):
                return p
    return None


def newest(base, rels):
    found = [p for rel in rels for p in glob.glob(os.path.join(base, rel), recursive=True) if os.path.isfile(p)]
    return max(found, key=os.path.getmtime) if found else None


def vcf_records(path):
    """Yield the split data lines of a VCF (plain or gzip)."""
    try:
        with open_text(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                yield line.rstrip("\n").split("\t")
    except (OSError, EOFError, gzip.BadGzipFile, zlib.error) as e:
        raise Unreadable(f"{os.path.basename(path)}: {e}") from e


def info_map(field):
    out = {}
    for kv in field.split(";"):
        k, _, v = kv.partition("=")
        out[k] = v
    return out


def read_tsv(path):
    with open_text(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


# --- sections -------------------------------------------------------------------

def sec_variants(d, s):
    p = first_existing(d, [f"vcf/{s}.vcf.gz"])
    if not p:
        return None, {}
    total = passed = snps = indels = 0
    for r in vcf_records(p):
        if len(r) < 7:
            continue
        total += 1
        if r[6] == "PASS":
            passed += 1
        alts = [a for a in r[4].split(",") if a not in (".", "*", "<NON_REF>")]
        if any(len(a) == len(r[3]) == 1 for a in alts):
            snps += 1
        if any(len(a) != len(r[3]) and not a.startswith("<") for a in alts):
            indels += 1
    return p, {"total": total, "pass": passed, "snps": snps, "indels": indels}


def clinvar_row(r):
    info = info_map(r[7]) if len(r) > 7 else {}
    rev = info.get("CLNREVSTAT", "")
    genes = [g.split(":")[0] for g in info.get("GENEINFO", "").split("|") if g]
    gt = r[9].split(":")[0] if len(r) > 9 else ""
    zyg = {"0/1": "het", "1/0": "het", "0|1": "het", "1|0": "het", "1/1": "hom", "1|1": "hom"}.get(gt, gt)
    return {
        "stars": STARS.get(rev, 0),
        "chrom": r[0], "pos": int(r[1]), "ref": r[3], "alt": r[4],
        "genotype": zyg,
        "gene": ",".join(genes) or ".",
        "significance": info.get("CLNSIG", "").replace("_", " ") or ".",
        "review_status": rev.replace("_", " ") or ".",
    }


def sec_clinvar(d, s):
    p = first_existing(d, [f"clinvar/{s}_clinvar_hits.vcf"])
    if not p:
        return None, {}
    rows = [clinvar_row(r) for r in vcf_records(p) if len(r) >= 8]
    rows.sort(key=lambda x: (-x["stars"], x["chrom"], x["pos"]))
    by = {str(k): sum(1 for x in rows if x["stars"] == k) for k in (4, 3, 2, 1, 0)}
    return p, {"count": len(rows), "by_stars": by, "hits": rows}


def pharmcat_json(d, s):
    return newest(d, [f"pharmcat/{s}.report.json", "pharmcat/*.report.json",
                      f"vcf/{s}.report.json", "vcf/*.report.json"])


def sec_pharmcat(d, s):
    p = pharmcat_json(d, s)
    if not p:
        return None, {}
    try:
        with open(p) as f:
            data = json.load(f)
    except (OSError, ValueError) as e:
        raise Unreadable(f"{os.path.basename(p)}: {e}")
    html = os.path.isfile(p[: -len(".json")] + ".html")
    return p, {"version": str(data.get("pharmcatVersion", "")) if isinstance(data, dict) else "",
               "html_report": html}


def sec_cpic(d, s):
    p = first_existing(d, [f"cpic/{s}_phenotypes.tsv"])
    if not p:
        return None, {}
    genes = []
    for row in read_tsv(p):
        genes.append({"gene": row.get("Gene", ""), "diplotype": row.get("Diplotype", ""),
                      "phenotype": row.get("Phenotype", ""), "status": row.get("Status", "")})
    rec = os.path.join(os.path.dirname(p), f"{s}_cpic_recommendations.txt")
    warnings, failed = [], False
    if os.path.isfile(rec):
        with open(rec, errors="replace") as f:
            for line in f:
                if "PARSING FAILED" in line:
                    failed = True
                if line.strip().startswith("WARNING:"):
                    warnings.append(line.strip()[len("WARNING:"):].strip())
    return p, {"genes": genes,
               "non_normal": sum(1 for g in genes if g["status"] == "non-normal"),
               "ambiguous": sum(1 for g in genes if g["status"] == "ambiguous"),
               "not_called": sum(1 for g in genes if g["status"] == "not called"),
               "unclassified": sum(1 for g in genes if g["status"] == "unclassified"),
               "parse_failed": failed or not genes,
               "warnings": warnings}


def sec_pypgx(d, s):
    p = first_existing(d, [f"pypgx/{s}_pypgx_summary.tsv"])
    if not p:
        return None, {}
    rows = read_tsv(p)
    # A no-call by step 36's list (Indeterminate, FAILED, N/A...), in any case
    called = [r for r in rows if (r.get("Diplotype") or "").strip().lower() not in pgx_outside_calls.NO_CALL]
    cyp = next((r.get("Diplotype", "") for r in rows if r.get("Gene") == "CYP2D6"), "")
    out = {"genes_total": len(rows), "genes_called": len(called), "cyp2d6": cyp or "not called"}
    comp = os.path.join(os.path.dirname(p), f"{s}_pharmcat_comparison.tsv")
    if os.path.isfile(comp):
        crow = read_tsv(comp)
        out["comparison"] = {
            "conflicts": sum(1 for r in crow if r.get("Match") == "No"),
            "one_tool_only": sum(1 for r in crow if r.get("Match") in ("pypgx only", "PharmCAT only")),
            "concordant": sum(1 for r in crow if r.get("Match") == "Yes"),
        }
    return p, out


def sec_cyrius(d, s):
    p = first_existing(d, [f"cyrius/{s}_cyp2d6.tsv", "cyrius/*_cyp2d6.tsv"])
    if not p:
        return None, {}
    rows = read_tsv(p)
    r = rows[0] if rows else {}
    return p, {"genotype": r.get("Genotype", "") or "none", "filter": r.get("Filter", "")}


def sec_hla(d, s):
    p = first_existing(d, [f"hla_t1k/{s}_hla_genotype.tsv", f"hla/{s}_hla_genotype.tsv", "hla/*_hla_genotype.tsv"])
    if not p:
        return None, {}
    loci = []
    with open(p, errors="replace") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3 or c[0].startswith("#"):
                continue
            # T1K: allele, abundance, quality, twice; '.', 0, -1 when there is no
            # second allele. Each quality stays with its allele.
            alleles, quals = [], []
            for i in (2, 5):
                a = c[i] if len(c) > i else ""
                if a and a != ".":
                    alleles.append(a)
                    quals.append(c[i + 2] if len(c) > i + 2 else "")
            # Step 36's rule (bin/pgx_outside_calls.py): an allele quality of 0 or
            # below withholds the whole gene; only HLA-A and HLA-B reach PharmCAT.
            low = any(hla_quality(q) <= 0 for q in quals)
            loci.append({"gene": c[0], "alleles": alleles, "quality": quals, "low_confidence": low,
                         "withheld_from_pharmcat": low and c[0] in pgx_outside_calls.HLA_GENES})
    return p, {"loci": loci}


def hla_quality(q):
    """T1K's quality as a number; one that is not a number counts as -1, as in step 36."""
    try:
        return float(q)
    except ValueError:
        return -1.0


def sec_prs(d, s):
    p = first_existing(d, [f"prs/{s}_prs_summary.tsv"])
    if not p:
        return None, {}
    rows = []
    for r in read_tsv(p):
        pct = r.get("Percentile", "") or ""
        rows.append({"condition": r.get("Condition", ""), "pgs_id": r.get("PGS_ID", ""),
                     "score": r.get("Score_SUM", ""), "matched": r.get("Variants_Matched", ""),
                     "total": r.get("Variants_Total", ""),
                     "chrx": r.get("ChrX", ""),
                     "percentile": "" if pct == "NA" else pct,
                     "group": "" if r.get("Ancestry_Group", "NA") == "NA" else r.get("Ancestry_Group", "")})
    data = {"scores": rows,
            # Percentiles exist only when step 25 ran pgsc_calc with an ancestry panel.
            "adjusted": any(r["percentile"] for r in rows)}
    a = first_existing(d, [f"ancestry/{s}_ancestry.tsv", f"prs/{s}_ancestry.tsv"])
    if a:
        kv = {r.get("key"): r.get("value", "") for r in read_tsv(a) if r.get("key")}
        data["ancestry"] = {"population": kv.get("population", ""), "panel": kv.get("reference_panel", ""),
                            "low_confidence": kv.get("population_low_confidence", "")}
    return p, data


# --- PRS with pgsc_calc (steps 25 and 26, and the PRS processes) ------------------
# The PGS Catalog's pgsc_calc scores the sample. Step 25 and the Nextflow
# PRS_PREPARE and PRS_SUMMARY processes run the two commands below, so the
# score files pgsc_calc reads and the table the reports read are made once.

PRS_COLUMNS = ["Condition", "PGS_ID", "Score_SUM", "Variants_Matched", "Variants_Total",
               "Matched_Pct", "Percentile", "Ancestry_Group", "ChrX", "Input"]
# prs-format writes it beside the scores: per score, its chrX rows and
# whether they were kept (scored) or left out (no sex given).
CHRX_FILE = "chrx_rows.tsv"
# Columns of a scoring file that make a score non-additive: pgsc_calc scores
# them, but this pipeline's custom file keeps effect_weight only.
NON_ADDITIVE = ("dosage_0_weight", "dosage_1_weight", "dosage_2_weight")
AUTOSOMES = frozenset(str(i) for i in range(1, 23))


def pgs_id_of(path):
    """PGS000018 for PGS000018.txt.gz, PGS000018_hmPOS_GRCh38.txt.gz and the like."""
    name = os.path.basename(path)
    for suffix in (".gz", ".txt", ".tsv"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name.split("_hmPOS_")[0]


def read_labels(path):
    """ID -> trait label of assets/pgs_scores.tsv (columns pgs_id, trait_reported)."""
    out = {}
    if not path:
        return out
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader((line for line in f if not line.startswith("#")), delimiter="\t"):
            if r.get("pgs_id"):
                out[r["pgs_id"].strip()] = (r.get("trait_reported") or "").strip()
    return out


def single_allele(a):
    a = (a or "").strip().upper()
    return a if a and "/" not in a and a != "." else ""


def format_pgs(path, out_dir, label=None, keep_x=False):
    """Write the GRCh38-harmonised scoring file PATH as a custom GRCh38 file
    pgsc_calc reads (chr_name and chr_position from hm_chr and hm_pos), and
    return (pgs_id, rows written, rows dropped off the autosomes,
    [(chrom, pos, allele), ...], chrX rows) with the effect and other alleles
    of every row written. Chromosomes 1 to 22 are kept, and chrX with keep_x
    (the sample's sex is known): pgsc_calc converts the sample with plink2,
    which stops on a chrX record when no sex is given. chrY and MT rows are
    always dropped. Raises Unreadable for a file that is not harmonised to
    GRCh38, lacks a column, is not additive or has no row."""
    header, cols, rows, alleles = {}, None, [], []
    with open_text(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                k, _, v = line[1:].partition("=")
                header[k.strip()] = v.strip()
                continue
            if cols is None:
                cols = line.split("\t")
                continue
            if line:
                rows.append(line.split("\t"))
    name = os.path.basename(path)
    if header.get("HmPOS_build") != "GRCh38":
        raise Unreadable(f"{name} has #HmPOS_build='{header.get('HmPOS_build', '')}', expected GRCh38 "
                         "(use the PGS Catalog's Harmonized/<id>_hmPOS_GRCh38 file)")
    c = {k: i for i, k in enumerate(cols or [])}
    missing = [k for k in ("hm_chr", "hm_pos", "effect_allele", "effect_weight") if k not in c]
    if missing:
        raise Unreadable(f"{name} has no {', '.join(missing)} column")
    if any(k in c for k in NON_ADDITIVE):
        raise Unreadable(f"{name} has dosage weights (a non-additive score); only additive scores are supported")

    def get(r, k):
        return r[c[k]].strip() if k in c and c[k] < len(r) else ""

    for k in ("is_dominant", "is_recessive"):
        if k in c and any(get(r, k) == "True" for r in rows):
            raise Unreadable(f"{name} has {k}=True rows (a non-additive score); only additive scores are supported")
    pid = header.get("pgs_id") or pgs_id_of(path)
    trait = (label or header.get("trait_reported") or pid).replace("=", "-").replace("\n", " ")
    out = os.path.join(out_dir, f"{pid}.txt.gz")
    n = off = nx = 0
    with gzip.open(out + ".tmp", "wt", encoding="utf-8") as w:
        w.write(f"#pgs_id={pid}\n#pgs_name={pid}\n#trait_reported={trait}\n#genome_build=GRCh38\n")
        w.write("chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n")
        for r in rows:
            chrom, pos, ea, ew = get(r, "hm_chr"), get(r, "hm_pos"), get(r, "effect_allele").upper(), get(r, "effect_weight")
            if not (chrom and pos and ea and ew):
                continue   # a row the catalog could not place on GRCh38
            chrom = chrom[3:] if chrom.startswith("chr") else chrom
            if chrom == "X":
                nx += 1
            if chrom not in AUTOSOMES and not (keep_x and chrom == "X"):
                off += 1   # chrY, MT, and chrX unless the sex is known: plink2 refuses chrX without it
                continue
            oa = single_allele(get(r, "other_allele")) or single_allele(get(r, "hm_inferOtherAllele"))
            w.write(f"{chrom}\t{pos}\t{ea}\t{oa}\t{ew}\n")
            n += 1
            alleles.append((f"chr{chrom}", pos, ea))
            if oa:
                alleles.append((f"chr{chrom}", pos, oa))
    if n == 0:
        os.remove(out + ".tmp")
        raise Unreadable(f"{name} has no autosomal row with hm_chr, hm_pos, effect_allele and effect_weight")
    os.replace(out + ".tmp", out)
    return pid, n, off, alleles, nx


def prs_format_main(argv):
    ap = argparse.ArgumentParser(prog="collect_summary.py prs-format",
                                 description="Write the scoring files as the custom GRCh38 files pgsc_calc reads.")
    ap.add_argument("--scores", required=True, help="folder of GRCh38-harmonised PGS Catalog scoring files")
    ap.add_argument("--ids", default="", help="comma-separated PGS ids to use (default: every file in --scores)")
    ap.add_argument("--labels", help="assets/pgs_scores.tsv: the trait label of each id")
    ap.add_argument("--out", required=True, help="folder for the files pgsc_calc reads")
    ap.add_argument("--alleles", required=True,
                    help="written: CHROM POS ALLELE of every effect and other allele, sorted and unique")
    ap.add_argument("--keep-x", action="store_true",
                    help="keep the chrX rows (the sample's sex is known and reaches pgsc_calc's plink2); "
                         "without it they are dropped with chrY and MT")
    a = ap.parse_args(argv)
    labels = read_labels(a.labels)
    files = sorted(glob.glob(os.path.join(a.scores, "*.txt.gz")) + glob.glob(os.path.join(a.scores, "*.txt")))
    if a.ids:
        want = [x.strip() for x in a.ids.split(",") if x.strip()]
        by_id = {pgs_id_of(p): p for p in files}
        missing = [x for x in want if x not in by_id]
        if missing:
            print(f"ERROR: no scoring file for {', '.join(missing)} in {a.scores}", file=sys.stderr)
            return 1
        files = [by_id[x] for x in want]
    if not files:
        print(f"ERROR: no scoring file (*.txt.gz or *.txt) in {a.scores}", file=sys.stderr)
        return 1
    os.makedirs(a.out, exist_ok=True)
    allele_set, chrx = set(), []
    try:
        for p in files:
            pid, n, off, al, nx = format_pgs(p, a.out, labels.get(pgs_id_of(p)), a.keep_x)
            allele_set.update(al)
            chrx.append((pid, nx))
            kept = f", {nx} on chrX scored" if a.keep_x and nx else ""
            dropped = f", {off} off the autosomes dropped" if off else ""
            print(f"  {pid}: {n} GRCh38 rows{kept}{dropped} ({labels.get(pid) or 'label from the file'})")
    except (Unreadable, OSError, EOFError, zlib.error) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    with open(os.path.join(a.out, CHRX_FILE), "w") as w:
        w.write("pgs_id\tchrx_rows\tchrx\n")
        for pid, nx in chrx:
            w.write(f"{pid}\t{nx}\t{'scored' if a.keep_x else 'left out'}\n")

    def key(t):
        return (t[0], int(t[1]) if t[1].isdigit() else 0, t[2])
    with open(a.alleles + ".tmp", "w") as w:
        for t in sorted(allele_set, key=key):
            w.write("\t".join(t) + "\n")
    os.replace(a.alleles + ".tmp", a.alleles)
    return 0


def fmt_num(x):
    """A number as pgsc_calc wrote it, without a trailing .0: 62.0 -> 62."""
    try:
        return f"{float(x):.10g}"
    except (TypeError, ValueError):
        return "NA"


def read_gz_tsv(path):
    with open_text(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


# pgscatalog-match's log line for a score under the minimum overlap.
BELOW_RE = re.compile(r"Score (\S+) fails minimum matching threshold \(([0-9.]+)% variants match\)")


def prs_table(results, sampleset, scores_dir, input_kind, zero_matches=False, below_log=None):
    """Rows of the summary table (PRS_COLUMNS) and the ancestry rows (key,
    value), from pgsc_calc's output folder RESULTS for SAMPLESET. Totals and
    labels come from the files prs-format wrote in SCORES_DIR. With
    zero_matches (pgsc_calc stopped because no score variant is in the
    sample's genotypes) every score is reported unmatched. With below_log
    (pgsc_calc stopped because every score matched under its minimum
    overlap, so it published no match summary) no score has a sum, the
    matched count is unknown (NA), and the match rate is read from that log
    when pgscatalog-match printed it."""
    totals, labels = {}, {}
    for p in sorted(glob.glob(os.path.join(scores_dir, "*.txt.gz"))):
        pid, n = pgs_id_of(p), 0
        with open_text(p) as f:
            for line in f:
                if line.startswith("#trait_reported="):
                    labels[pid] = line.rstrip("\n").split("=", 1)[1]
                elif not line.startswith("#"):
                    n += 1
        totals[pid] = n - 1   # minus the column header
    if not totals:
        raise Unreadable(f"no formatted scoring file in {scores_dir}")
    # "8 scored", "8 left out" or "0"; empty for scores formatted before the file existed
    chrx = {}
    if os.path.isfile(os.path.join(scores_dir, CHRX_FILE)):
        for r in read_tsv(os.path.join(scores_dir, CHRX_FILE)):
            n = r.get("chrx_rows", "")
            chrx[r.get("pgs_id", "")] = f"{n} {r.get('chrx', '')}" if n not in ("", "0") else n
    matched, sums, pct, group = {}, {}, {}, {}
    pops, rate = {}, {}
    if below_log:
        with open(below_log, errors="replace") as f:
            for m in BELOW_RE.finditer(f.read()):
                rate[pgs_id_of(m.group(1))] = f"{float(m.group(2)):.1f}"
    elif not zero_matches:
        summ = os.path.join(results, sampleset, "match", f"{sampleset}_summary.csv")
        if not os.path.isfile(summ):
            raise Unreadable(f"pgsc_calc wrote no match summary ({summ})")
        with open(summ, newline="") as f:
            for r in csv.DictReader(f):
                pid = pgs_id_of(r.get("accession", ""))
                if r.get("match_status") == "matched":
                    matched[pid] = matched.get(pid, 0) + int(float(r.get("count") or 0))
        score_dir = os.path.join(results, sampleset, "score")
        adjusted = os.path.join(score_dir, f"{sampleset}_pgs.txt.gz")
        plain = os.path.join(score_dir, "aggregated_scores.txt.gz")
        src = adjusted if os.path.isfile(adjusted) else plain
        if not os.path.isfile(src):
            raise Unreadable(f"pgsc_calc wrote no scores ({plain})")
        for r in read_gz_tsv(src):
            if r.get("sampleset") != sampleset:
                continue   # the reference panel's own samples
            pid = pgs_id_of(r.get("PGS", ""))
            sums[pid] = fmt_num(r.get("SUM"))
            if r.get("percentile_MostSimilarPop") not in (None, ""):
                v = fmt_num(r["percentile_MostSimilarPop"])
                pct[pid] = "NA" if v == "NA" else f"{float(v):.1f}"
        popsim = os.path.join(score_dir, f"{sampleset}_popsimilarity.txt.gz")
        if os.path.isfile(popsim):
            target = [r for r in read_gz_tsv(popsim)
                      if r.get("sampleset") == sampleset and r.get("REFERENCE", "False") != "True"]
            if len(target) != 1:
                raise Unreadable(f"{os.path.basename(popsim)} has {len(target)} rows for {sampleset}, expected 1")
            pops = target[0]
            for pid in sums:
                group[pid] = pops.get("MostSimilarPop") or "NA"
    rows = []
    for pid in sorted(totals):
        m, t = matched.get(pid, 0), totals[pid]
        if below_log:
            rows.append([labels.get(pid, pid), pid, "NA", "NA", str(t), rate.get(pid, "NA"), "NA", "NA",
                         chrx.get(pid, ""), input_kind])
            continue
        rows.append([labels.get(pid, pid), pid, sums.get(pid, "NA"), str(m), str(t),
                     f"{100 * m / t:.1f}" if t else "0.0", pct.get(pid, "NA"), group.get(pid, "NA"),
                     chrx.get(pid, ""), input_kind])
    ancestry = []
    if pops:
        pcs = sorted((k for k in pops if k.startswith("PC") and k[2:].isdigit()), key=lambda k: int(k[2:]))
        ancestry.append(("population", pops.get("MostSimilarPop", "")))
        ancestry.append(("population_low_confidence", pops.get("MostSimilarPop_LowConfidence", "")))
        for k in sorted(k for k in pops if k.startswith("RF_P_")):
            ancestry.append((f"probability_{k[5:]}", fmt_num(pops[k])))
        for k in pcs:
            ancestry.append((k, fmt_num(pops[k])))
    return rows, ancestry


def prs_table_main(argv):
    ap = argparse.ArgumentParser(prog="collect_summary.py prs-table",
                                 description="Write step 25's summary (and step 26's ancestry) table from pgsc_calc's output.")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--results", required=True, help="pgsc_calc's --outdir")
    ap.add_argument("--sampleset", required=True, help="the sampleset name given to pgsc_calc")
    ap.add_argument("--scores", required=True, help="the folder prs-format wrote")
    ap.add_argument("--input-kind", required=True, choices=["gvcf", "vcf"],
                    help="gvcf: the score positions were genotyped from the gVCF; vcf: variant sites only")
    ap.add_argument("--panel", default="", help="name of the ancestry reference panel, when one was used")
    ap.add_argument("--zero-matches", action="store_true",
                    help="pgsc_calc stopped because no score variant is in the sample's genotypes")
    ap.add_argument("--below-threshold", metavar="LOG",
                    help="pgsc_calc stopped because every score matched under its minimum overlap; LOG is its console output")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ancestry-out", help="written when pgsc_calc ran with the panel: PCs and population")
    a = ap.parse_args(argv)
    try:
        rows, ancestry = prs_table(a.results, a.sampleset, a.scores, a.input_kind, a.zero_matches, a.below_threshold)
    except (Unreadable, OSError, ValueError, KeyError, EOFError, csv.Error, zlib.error) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    with open(a.out + ".tmp", "w", encoding="utf-8") as f:
        f.write("\t".join(PRS_COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    os.replace(a.out + ".tmp", a.out)
    if a.ancestry_out and ancestry:
        with open(a.ancestry_out + ".tmp", "w") as f:
            f.write("key\tvalue\n")
            for k, v in [("sample", a.sample), ("reference_panel", a.panel)] + ancestry:
                f.write(f"{k}\t{v}\n")
        os.replace(a.ancestry_out + ".tmp", a.ancestry_out)
    for r in rows:
        extra = f", percentile {r[6]} among {r[7]}" if r[6] != "NA" else ""
        if r[3] == "NA":
            print(f"  {r[0]} ({r[1]}): no sum, under pgsc_calc's minimum overlap ({r[5]}% of {r[4]} variants matched)")
        else:
            print(f"  {r[0]} ({r[1]}): sum {r[2]}, {r[3]} of {r[4]} variants matched ({r[5]}%){extra}")
    return 0


def count_vcf(p):
    total = passed = 0
    for r in vcf_records(p):
        total += 1
        if len(r) > 6 and r[6] == "PASS":
            passed += 1
    return total, passed


def sec_manta(d, s):
    p = first_existing(d, ["manta/results/variants/diploidSV.vcf.gz", "manta/*diploidSV.vcf.gz",
                           "manta/**/diploidSV.vcf.gz"])
    if not p:
        return None, {}
    t, ps = count_vcf(p)
    return p, {"total": t, "pass": ps}


def sec_delly(d, s):
    p = first_existing(d, [f"delly/{s}_sv.vcf.gz"])
    if not p:
        return None, {}
    t, ps = count_vcf(p)
    return p, {"total": t, "pass": ps}


def sec_cnvpytor(d, s):
    p = first_existing(d, [f"cnvpytor/{s}_cnvs.txt"])
    if not p:
        return None, {}
    total = dels = dups = sig = 0
    with open(p, errors="replace") as f:
        for line in f:
            c = line.split()
            if not c:
                continue
            total += 1
            dels += c[0] == "deletion"
            dups += c[0] == "duplication"
            try:
                sig += float(c[4]) < 0.01
            except (IndexError, ValueError):
                pass
    return p, {"total": total, "deletions": dels, "duplications": dups, "significant": sig}


def sec_sv_consensus(d, s):
    p = first_existing(d, [f"sv_merged/{s}_sv_consensus.vcf.gz"])
    if not p:
        return None, {}
    t, _ = count_vcf(p)
    return p, {"consensus": t}


def sec_expansions(d, s):
    eh = first_existing(d, [f"expansion_hunter/{s}_eh.vcf", f"expansion_hunter/{s}.vcf",
                            "expansion_hunter/*_eh.vcf"])
    st = first_existing(d, [f"expansion_hunter/{s}_eh_stranger.vcf", "expansion_hunter/*_eh_stranger.vcf"])
    # Stranger (step 9b) copies ExpansionHunter's records and adds STR_STATUS.
    # Its file is read unless it is older than the ExpansionHunter file it
    # would have been made from; then the report says to rerun step 9b.
    outdated = bool(st and eh and os.path.getmtime(st) < os.path.getmtime(eh))
    if outdated:
        st = None
    p = st or eh
    if not p:
        return None, {}
    loci, tested, flagged, no_status = {}, 0, [], 0
    for r in vcf_records(p):
        if len(r) < 10:
            continue
        tested += 1
        info = info_map(r[7])
        rep = info.get("REPID") or info.get("VARID") or ""
        fmt = r[8].split(":")
        val = r[9].split(":")
        repcn = val[fmt.index("REPCN")] if "REPCN" in fmt and fmt.index("REPCN") < len(val) else ""
        if rep in EH_LOCI and rep not in loci:
            loci[rep] = repcn or "."
        if st:
            # One value per record, the most severe of normal, pre_mutation and
            # full_mutation; a locus outside Stranger's catalog has none.
            status = [x for x in info.get("STR_STATUS", "").split(",") if x and x != "."]
            if not status:
                no_status += 1
            elif any(x != "normal" for x in status):
                flagged.append({"locus": rep or f"{r[0]}:{r[1]}", "repeat_count": repcn or ".",
                                "status": ",".join(status)})
    out = {"records": tested, "key_loci": [{"locus": k, "repeat_count": loci.get(k, "not in output")}
                                           for k in EH_LOCI],
           "stranger": bool(st), "stranger_outdated": outdated}
    if st:
        out.update(flagged=flagged, no_status=no_status)
    return p, out


def sec_telomere(d, s):
    p = first_existing(d, [f"telomere/{s}/{s}/{s}_summary.tsv", f"telomere/{s}/{s}_summary.tsv",
                           "telomere/**/*_summary.tsv"])
    if not p:
        return None, {}
    with open(p, errors="replace") as f:
        lines = [l.rstrip("\n").split("\t") for l in f if l.strip()]
    if len(lines) < 2:
        return p, {"tel_content": "."}
    head, row = lines[0], lines[1]
    # By header name; older TelomereHunter versions put it in column 11.
    i = head.index("tel_content") if "tel_content" in head else 10
    return p, {"tel_content": row[i] if i < len(row) else "."}


def sec_roh(d, s):
    p = first_existing(d, [f"vcf/{s}_roh.txt", f"roh/{s}_roh.txt", "roh/*_roh.txt"])
    if not p:
        return None, {}
    total = largest = 0.0
    segs, big = 0, []
    with open(p, errors="replace") as f:
        for line in f:
            c = line.split()
            if len(c) < 6 or c[0] != "RG":
                continue
            try:
                length = float(c[5])
            except ValueError:
                continue
            segs += 1
            total += length
            largest = max(largest, length)
            if length > 5e6 and c[2] not in ("chrX", "chrY", "X", "Y"):
                big.append({"region": f"{c[2]}:{c[3]}-{c[4]}", "mb": round(length / 1e6, 1)})
    return p, {"segments": segs, "total_mb": round(total / 1e6, 1), "largest_mb": round(largest / 1e6, 1),
               "autosomal_over_5mb": big}


def sec_haplogroup(d, s):
    p = first_existing(d, [f"mito/{s}_haplogroup.txt"])
    if not p:
        return None, {}
    rows = read_tsv(p)
    r = rows[0] if rows else {}
    hg = (r.get("Haplogroup") or (list(r.values())[1] if len(r) > 1 else "") or "").strip('"')
    out = {"haplogroup": hg or "."}
    # haplocheck (step 12 with step 20's Mutect2 calls): a second haplogroup
    # in the allele fractions is another person's DNA in the sample.
    hc = first_existing(d, [f"mito/{s}_haplocheck.txt"])
    if hc:
        rows = [{k.strip('"'): (v or "").strip('"') for k, v in row.items() if k} for row in read_tsv(hc)]
        if not rows or "Contamination Status" not in rows[0]:
            raise Unreadable(f"{os.path.basename(hc)} has no Contamination Status column")
        out["contamination_status"] = rows[0]["Contamination Status"] or "."
        out["contamination_level"] = rows[0].get("Contamination Level") or "."
    return p, out


def sec_y_haplogroup(d, s):
    p = first_existing(d, [f"y_haplogroup/{s}_y_haplogroup.txt"])
    if not p:
        return None, {}
    rows = read_tsv(p)
    if not rows or "Hg" not in rows[0]:
        raise Unreadable(f"{os.path.basename(p)} has no Hg column or no sample row")
    r = rows[0]
    hg = (r.get("Hg") or "").strip()
    return p, {"haplogroup": hg if hg not in ("", "NA") else "insufficient markers",
               "valid_markers": r.get("Valid_markers") or ".", "qc_score": r.get("QC-score") or "."}


def sec_mito(d, s):
    p = first_existing(d, [f"mito/{s}_chrM_filtered.vcf.gz"])
    if not p:
        return None, {}
    passed = het = 0
    for r in vcf_records(p):
        if len(r) < 10 or r[6] != "PASS":
            continue
        passed += 1
        fmt = r[8].split(":")
        val = r[9].split(":")
        if "AF" not in fmt or fmt.index("AF") >= len(val):
            continue
        try:
            af = max(float(x) for x in val[fmt.index("AF")].split(",") if x not in ("", "."))
        except ValueError:
            continue
        # Heteroplasmic: 5% to 95%. Below 5% NUMT reads and noise dominate.
        if 0.05 <= af < 0.95:
            het += 1
    return p, {"pass": passed, "heteroplasmic": het, "heteroplasmy_floor": 0.05}


def sec_cpsr(d, s):
    html = first_existing(d, [f"cpsr/{s}.cpsr.grch38.html", "cpsr/*.cpsr.*.html"])
    tsv = first_existing(d, [f"cpsr/{s}.cpsr.grch38.classification.tsv.gz", "cpsr/*.cpsr.*.classification.tsv.gz"])
    p = tsv or html
    if not p:
        return None, {}
    out = {"html_report": bool(html)}
    if tsv:
        try:
            with gzip.open(tsv, "rt", errors="replace") as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
                cols = rows[0].keys() if rows else []
        except (OSError, EOFError, gzip.BadGzipFile, zlib.error) as e:
            raise Unreadable(f"{os.path.basename(tsv)}: {e}") from e
        col = next((c for c in ("CLASSIFICATION", "FINAL_CLASSIFICATION", "CPSR_CLASSIFICATION") if c in cols), None)
        if col:
            counts = {}
            for r in rows:
                counts[r.get(col) or "."] = counts.get(r.get(col) or ".", 0) + 1
            out["classification_column"] = col
            out["classification"] = dict(sorted(counts.items(), key=lambda kv: -kv[1]))
        else:
            out["classification_column"] = None
    return p, out


def csq_high_impact_genes(path):
    """{(chrom, pos, ref, alt): (genotype, {gene: consequence})} of the records
    with a HIGH-impact VEP consequence (any transcript), from the CSQ field."""
    fmt, out = None, {}
    try:
        with open_text(path) as f:
            for line in f:
                if line.startswith("##INFO=<ID=CSQ"):
                    m = re.search(r"Format: ([^\"]+)", line)
                    fmt = m.group(1).split("|") if m else None
                    continue
                if line.startswith("#") or not fmt or "IMPACT" not in fmt or "SYMBOL" not in fmt:
                    continue
                r = line.rstrip("\n").split("\t")
                if len(r) < 8:
                    continue
                csq = info_map(r[7]).get("CSQ", "")
                genes = {}
                for tr in csq.split(","):
                    v = tr.split("|")
                    if len(v) != len(fmt) or v[fmt.index("IMPACT")] != "HIGH":
                        continue
                    cons = v[fmt.index("Consequence")] if "Consequence" in fmt else ""
                    genes.setdefault(v[fmt.index("SYMBOL")], cons.replace("&", ","))
                if genes:
                    gt = r[9].split(":")[0] if len(r) > 9 else ""
                    out[(r[0], int(r[1]), r[3], r[4])] = (gt, genes)
    except (OSError, EOFError, gzip.BadGzipFile, zlib.error) as e:
        raise Unreadable(f"{os.path.basename(path)}: {e}") from e
    return out


def acmg_sf_tier(d, s, clinical_vcf):
    """The variants in ACMG SF genes: the ClinVar P/LP hits of step 6, and the
    clinical filter's rare HIGH-impact records."""
    clinvar, high = [], []
    hits = first_existing(d, [f"clinvar/{s}_clinvar_hits.vcf"])
    if hits:
        for r in vcf_records(hits):
            if len(r) < 8:
                continue
            row = clinvar_row(r)
            on = [g for g in row["gene"].split(",") if g in ACMG_SF_GENES]
            if on:
                clinvar.append({"gene": ",".join(on), "variant": f"{row['chrom']}:{row['pos']} {row['ref']}>{row['alt']}",
                                "genotype": row["genotype"], "significance": row["significance"], "stars": row["stars"]})
    if clinical_vcf:
        for (c, pos, ref, alt), (gt, genes) in sorted(csq_high_impact_genes(clinical_vcf).items()):
            for g, cons in sorted(genes.items()):
                if g in ACMG_SF_GENES:
                    zyg = {"0/1": "het", "1/0": "het", "0|1": "het", "1|0": "het", "1/1": "hom", "1|1": "hom"}.get(gt, gt)
                    high.append({"gene": g, "variant": f"{c}:{pos} {ref}>{alt}", "genotype": zyg, "consequence": cons})
    return {"version": ACMG_SF_VERSION, "genes_on_list": len(ACMG_SF_GENES),
            "clinvar_hits": clinvar, "clinvar_checked": bool(hits),
            "high_impact": high, "high_impact_checked": bool(clinical_vcf)}


def sec_clinical(d, s):
    p = first_existing(d, [f"clinical/{s}_clinical_summary.tsv"])
    v = first_existing(d, [f"clinical/{s}_clinical.vcf.gz"])
    if not p:
        # The Nextflow report gets the clinical VCF, not the summary table:
        # the count only (one record per variant, as in the table).
        if not v:
            return None, {}
        return v, {"variants": sum(1 for _ in vcf_records(v)), "acmg_sf": acmg_sf_tier(d, s, v)}
    rows = read_tsv(p)
    by = {}
    for r in rows:
        by[r.get("IMPACT") or "."] = by.get(r.get("IMPACT") or ".", 0) + 1
    genes = sorted({r.get("GENE") for r in rows if r.get("GENE") not in (None, "", ".")})
    return p, {"variants": len(rows), "by_impact": by, "genes": len(genes), "acmg_sf": acmg_sf_tier(d, s, v)}


def sec_slivar(d, s):
    p = first_existing(d, [f"slivar/{s}_slivar_summary.tsv"])
    if not p:
        # The Nextflow report gets the prioritized VCF: the count only.
        v = first_existing(d, [f"slivar/{s}_prioritized.vcf.gz"])
        if not v:
            return None, {}
        return v, {"prioritized": sum(1 for _ in vcf_records(v))}
    rows = read_tsv(p)
    out = {"prioritized": len(rows)}
    ch = os.path.join(os.path.dirname(p), f"{s}_compound_hets.tsv")
    if os.path.isfile(ch):
        crow = read_tsv(ch)
        out["compound_het_variants"] = len(crow)
        out["compound_het_genes"] = len({r.get("GENE") for r in crow if r.get("GENE")})
    return p, out


def sec_coverage(d, s):
    p = first_existing(d, [f"mosdepth/{s}.mosdepth.summary.txt", f"coverage/{s}.mosdepth.summary.txt"])
    if not p:
        return None, {}
    mean = None
    for r in read_tsv(p):
        if r.get("chrom") == "total":
            mean = r.get("mean")
    return p, {"mean_depth": float(mean) if mean not in (None, "") else None}


def sec_sex_check(d, s):
    p = first_existing(d, ["indexcov/indexcov-indexcov.ped", "indexcov/*-indexcov.ped"])
    if not p:
        return None, {}
    with open(p, errors="replace") as f:
        lines = [l.rstrip("\n").split("\t") for l in f if l.strip()]
    head = [h.lstrip("#") for h in lines[0]] if lines else []
    if "sex" not in head or len(lines) < 2:
        raise Unreadable(f"{os.path.basename(p)} has no sex column or no sample row")
    row = dict(zip(head, lines[-1]))
    sex = {"1": "male", "2": "female"}.get(row.get("sex", ""), "unknown")
    return p, {"inferred_sex": sex, "cn_chrX": row.get("CNchrX", ""), "cn_chrY": row.get("CNchrY", "")}


# --- sample identity and contamination (step 33) ----------------------------------

FREEMIX_WARN = 0.03   # VerifyBamID2 FREEMIX above this is reported as possible contamination
SAME_PERSON = 0.9     # somalier relatedness at or above this: the same person (a duplicate or a swap)
SOMALIER_SEX = {"1": "male", "2": "female"}   # its samples.tsv sex column; -9 unknown, -2 X and Y disagree


def read_rows(path):
    """Rows of a tab-separated file whose header may start with '#'."""
    with open_text(path) as f:
        lines = [l.rstrip("\n").split("\t") for l in f if l.strip()]
    if not lines:
        return []
    head = [lines[0][0].lstrip("#")] + lines[0][1:]
    return [dict(zip(head, r)) for r in lines[1:]]


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if v == v else None   # NaN is no value


def somalier_x_sex(row):
    """male or female when somalier's own rule calls the sex from chrX, else
    None. The rule of relate --infer (relate.nim, add_parents_and_check_sex):
    more than 10 chrX sites with depth, heterozygous / homozygous-ALT sites
    below 0.05 for male or above 0.4 for female, and fewer than 6% of all
    sites with an allele balance outside 0.1 to 0.9."""
    n, het, hom = num(row.get("X_n")), num(row.get("X_het")), num(row.get("X_hom_alt"))
    # A missing column passes (older somalier, test stubs); a present value
    # that is not a number (nan from 0/0 when no autosomal site has a call)
    # fails, as somalier's own `< 0.06` does.
    raw_mid = row.get("p_middling_ab")
    mid = num(raw_mid)
    if n is None or het is None or hom is None or n <= 10 or (
            raw_mid is not None and (mid is None or mid >= 0.06)):
        return None
    if hom == 0:
        return "female" if het > 0 else None
    ratio = het / hom
    return "male" if ratio < 0.05 else "female" if ratio > 0.4 else None


def sample_qc_table(sample, samples_tsv, selfsm=None, pairs_tsv=None, declared_sex=None,
                    freemix_warn=FREEMIX_WARN, somalier_id=None, marker_check=""):
    """The step 33 verdict for one sample, as ordered (key, value) pairs."""
    sid = somalier_id or sample
    rows = read_rows(samples_tsv)
    row = next((r for r in rows if r.get("sample_id") == sid), None)
    if row is None:
        raise Unreadable(f"{os.path.basename(samples_tsv)} has no row for sample '{sid}' "
                         f"(it has: {', '.join(r.get('sample_id', '?') for r in rows) or 'none'})")
    inferred = SOMALIER_SEX.get(row.get("sex", ""), "unknown")
    declared = (declared_sex or "").lower() or None
    # With the declared sex as its pedigree (--ped), somalier starts its sex
    # column from that sex and overwrites it only when the reads tell, so the
    # pedigree's sex left standing is not a call: unknown, as without a
    # pedigree. -2 from a pedigree female with chrY reads is not one either
    # unless chrX looks female too.
    ped = {"male": "male", "female": "female"}.get(row.get("original_pedigree_sex", "").lower())
    x_says = somalier_x_sex(row) if ped else None
    if ped and inferred == ped and x_says != ped:
        inferred = "unknown"
    x_and_y = row.get("sex", "") == "-2" and (not ped or x_says == "female")
    # -2 first: it is worth saying even when no sex was declared
    if x_and_y:
        sex_check = "not_checked"
        why = ("chrX is heterozygous like a female sample but chrY has reads: "
               "a sex-chromosome aneuploidy or a mixed sample")
    elif not declared:
        sex_check, why = "not_checked", "no declared sex"
    elif inferred == "unknown":
        sex_check = "not_checked"
        why = (f"somalier could not tell the sex from {row.get('X_n', '0')} chrX sites "
               f"(it needs more than 10 with reads, and allele balances that look like one person)")
    elif inferred == declared:
        sex_check, why = "ok", f"declared {declared}, somalier infers {inferred}"
    else:
        sex_check, why = "mismatch", f"declared {declared}, somalier infers {inferred}"
    gt = num(row.get("gt_depth_mean"))
    ydp = num(row.get("Y_depth_mean"))
    xdp = num(row.get("X_depth_mean"))
    sites = sum(int(num(row.get(k)) or 0) for k in ("n_hom_ref", "n_het", "n_hom_alt"))
    out = [("sample", sample), ("somalier_id", sid), ("declared_sex", declared or ""),
           ("inferred_sex", inferred), ("sex_check", sex_check), ("sex_check_reason", why),
           ("sites_genotyped", str(sites)), ("depth_at_sites", row.get("gt_depth_mean", "")),
           ("x_sites", row.get("X_n", "")), ("x_het", row.get("X_het", "")),
           ("x_hom_alt", row.get("X_hom_alt", "")),
           ("x_depth_ratio", f"{xdp / gt:.2f}" if gt and xdp is not None else ""),
           ("y_sites", row.get("Y_n", "")),
           ("y_depth_ratio", f"{ydp / gt:.2f}" if gt and ydp is not None else "")]
    freemix, status, markers = None, "not_run", ""
    if selfsm:
        sm = read_rows(selfsm)
        freemix = num(sm[0].get("FREEMIX")) if sm else None
        if freemix is None:
            raise Unreadable(f"{os.path.basename(selfsm)} has no FREEMIX value")
        status = "warn" if freemix > freemix_warn else "ok"
        # #SNPS of the selfSM is the size of the panel, not the markers with reads
        markers = sm[0].get("#SNPS", "")
    # marker_check: passed, or skipped when fewer than 1,000 panel markers had
    # reads and VerifyBamID2 ran with --DisableSanityCheck.
    out += [("freemix", f"{freemix:.4f}" if freemix is not None else ""),
            ("freemix_warn_above", f"{freemix_warn:g}"), ("contamination", status),
            ("panel_markers", markers), ("verifybamid2_marker_check", marker_check)]
    same = []
    if pairs_tsv:
        for r in read_rows(pairs_tsv):
            a, b = r.get("sample_a"), r.get("sample_b")
            rel = num(r.get("relatedness"))
            if sid in (a, b) and a != b and rel is not None and rel >= SAME_PERSON:
                same.append(f"{b if a == sid else a} ({rel:.2f})")
    out.append(("same_person_as", ", ".join(same)))
    return out


def sec_sample_qc(d, s):
    p = first_existing(d, [f"qc/{s}_sample_qc.tsv"])
    if not p:
        return None, {}
    data = {}
    for r in read_tsv(p):
        if r.get("key"):
            data[r["key"]] = r.get("value", "")
    if "inferred_sex" not in data:
        raise Unreadable(f"{os.path.basename(p)} has no inferred_sex row")
    return p, data


SECTIONS = [
    # key, title, step, reader
    ("coverage", "Coverage (mosdepth)", "16b", sec_coverage),
    ("sex_check", "Sex check (indexcov)", "16", sec_sex_check),
    ("sample_qc", "Sample identity and contamination (somalier, VerifyBamID2)", "33", sec_sample_qc),
    ("variants", "Variant calling", "03", sec_variants),
    ("clinvar", "ClinVar screen", "06", sec_clinvar),
    ("pharmcat", "PharmCAT", "07", sec_pharmcat),
    ("cpic", "CPIC drug-gene recommendations", "27", sec_cpic),
    ("pypgx", "pypgx", "32", sec_pypgx),
    ("cyrius", "Cyrius (CYP2D6)", "21", sec_cyrius),
    ("hla", "HLA typing (T1K)", "08", sec_hla),
    ("prs", "Polygenic risk scores", "25", sec_prs),
    ("manta", "Structural variants (Manta)", "04", sec_manta),
    ("delly", "Structural variants (Delly)", "19", sec_delly),
    ("cnvpytor", "Copy number (CNVpytor)", "18", sec_cnvpytor),
    ("sv_consensus", "SV consensus merge", "22", sec_sv_consensus),
    ("expansions", "Repeat expansions (ExpansionHunter)", "09", sec_expansions),
    ("telomere", "Telomere content (TelomereHunter)", "10", sec_telomere),
    ("roh", "Runs of homozygosity", "11", sec_roh),
    ("haplogroup", "Mitochondrial haplogroup", "12", sec_haplogroup),
    ("y_haplogroup", "Y-chromosome haplogroup (Yleaf)", "37", sec_y_haplogroup),
    ("mito", "Mitochondrial variants", "20", sec_mito),
    ("cpsr", "Cancer predisposition (CPSR)", "17", sec_cpsr),
    ("clinical", "Clinical variant filter", "23", sec_clinical),
    ("slivar", "Variant prioritization (slivar)", "31", sec_slivar),
]


# --- run status and manifest -----------------------------------------------------

def read_run_status(d):
    p = os.path.join(d, "logs", "run_status.tsv")
    if not os.path.isfile(p):
        return None
    st = {"started_epoch": None, "started_utc": None, "declared_sex": None, "steps": {}}
    with open(p, errors="replace") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3 or c[0].startswith("#"):
                continue
            if c[0] == "meta" and c[1] == "started_epoch":
                try:
                    st["started_epoch"] = float(c[2])
                except ValueError:
                    pass
            elif c[0] == "meta":
                st[c[1]] = c[2]
            elif c[0] == "step":
                st["steps"][c[1]] = c[2]
    if st["started_epoch"] is not None:
        st["started_utc"] = iso(st["started_epoch"])
    return st


def read_manifest(d):
    p = os.path.join(d, "run_manifest.tsv")
    if not os.path.isfile(p):
        return None
    rows = []
    with open(p, errors="replace") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 3 and not c[0].startswith("#"):
                rows.append([c[0], c[1], "\t".join(c[2:])])
    return rows


def manifest_value(rows, section, key):
    for r in rows or []:
        if r[0] == section and r[1] == key:
            return r[2]
    return None


def read_consensus(d, s, sections):
    """Step 36's CYP2D6 row (pgx_consensus/<id>_pgx_consensus.tsv), or None
    when the table is missing or older than a caller's result it was made from."""
    p = first_existing(d, [f"pgx_consensus/{s}_pgx_consensus.tsv"])
    if not p:
        return None
    for k in ("pypgx", "cyrius"):
        src = sections[k]["source"]
        if src and os.path.getmtime(os.path.join(d, src)) > os.path.getmtime(p):
            return None
    row = next((r for r in read_tsv(p) if (r.get("Gene") or "").strip() == "CYP2D6"), None)
    if row is None:
        return None
    return {"result": row.get("Result") or "", "passed_to_pharmcat": row.get("Outside_call") == "yes",
            "reason": row.get("Reason") or "", "source": os.path.relpath(p, d)}


def cyp2d6_block(sample, sample_dir, sections):
    """CYP2D6 from each caller, side by side, judged by step 36's rule
    (bin/pgx_outside_calls.py): a no-call in any case (Indeterminate, None,
    '.'...) is not a call, and Cyrius counts only with Filter PASS. PharmCAT
    calls no CYP2D6 from a VCF: its column is what step 36 passed to it, so
    agreement is between pypgx and Cyrius. With step 36's table its verdict
    is kept as well."""
    pc = next((g for g in sections["cpic"]["data"].get("genes", []) if g["gene"] == "CYP2D6"), None)
    pg = sections["pypgx"]["data"].get("cyp2d6") if sections["pypgx"]["state"] == "ok" else None
    cy = sections["cyrius"]["data"] if sections["cyrius"]["state"] == "ok" else None
    # A caller whose step failed says so, not 'not run'.
    failed = {k: f"failed (step {sections[k]['step']})" for k in ("cpic", "pypgx", "cyrius")
              if sections[k]["state"] == "failed"}
    pg_ok = pg is not None and pgx_outside_calls.is_call(pg)
    cy_ok = cy is not None and pgx_outside_calls.is_call(cy.get("genotype")) and cy.get("filter") == "PASS"
    cy_txt = None
    if cy is not None:
        cy_txt = cy.get("genotype") or "none"
        if pgx_outside_calls.is_call(cy_txt) and not cy_ok:
            cy_txt += f" (not usable: Filter {cy.get('filter') or 'missing'})"
    calls = {
        "PharmCAT": (pc["diplotype"] if pc["status"] not in ("not called", "ambiguous") else pc["status"]) if pc
        else failed.get("cpic"),
        "pypgx": failed.get("pypgx", pg),
        "Cyrius": failed.get("cyrius", cy_txt),
    }
    agree = pgx_outside_calls.same_diplotype(pg, cy["genotype"]) if pg_ok and cy_ok else None
    try:
        consensus = read_consensus(sample_dir, sample, sections)
    except (OSError, csv.Error) as e:
        consensus = {"result": "unreadable", "passed_to_pharmcat": False, "reason": str(e), "source": None}
    return {"calls": calls, "agree": agree, "consensus": consensus}


def collect(sample, sample_dir, declared_sex=None):
    status = read_run_status(sample_dir)
    manifest = read_manifest(sample_dir)
    sections = {}
    for key, title, step, reader in SECTIONS:
        sec = {"title": title, "step": step, "state": "missing", "source": None,
               "file_date": None, "note": None, "data": {}}
        try:
            path, data = reader(sample_dir, sample)
        except Unreadable as e:
            sec["state"], sec["note"] = "unreadable", str(e)
            sections[key] = sec
            continue
        except (OSError, ValueError, KeyError, IndexError, EOFError, csv.Error, zlib.error) as e:
            sec["state"], sec["note"] = "unreadable", f"{type(e).__name__}: {e}"
            sections[key] = sec
            continue
        if not path and status and (status["steps"].get(step) or "").startswith("failed"):
            # run-all.sh records 'failed' for a step it ran that exited non-zero
            sec["state"] = "failed"
            sec["note"] = (f"step {step} failed in the run of {status['started_utc'] or 'unknown date'}; "
                           "its log is in logs/")
        if path:
            mtime = os.path.getmtime(path)
            sec.update(state="ok", source=os.path.relpath(path, sample_dir), file_date=iso(mtime), data=data)
            if status and status["started_epoch"] is not None and step in RUN_ALL_STEPS:
                result = status["steps"].get(step)
                if result != "ok" and mtime < status["started_epoch"]:
                    sec["state"] = "stale"
                    sec["note"] = (f"from an earlier run ({iso(mtime)}); step {step} was "
                                   f"{result or 'not run'} in the run of {status['started_utc']}")
        sections[key] = sec

    cyp2d6 = cyp2d6_block(sample, sample_dir, sections)

    clinvar_date = manifest_value(manifest, "data", "clinvar_file_date")
    summary = {
        "schema_version": SCHEMA_VERSION,
        "sample": sample,
        "generated_utc": iso(datetime.now(tz=timezone.utc).timestamp()),
        "run": {
            "status_file": status is not None,
            "started_utc": status["started_utc"] if status else None,
            "declared_sex": (declared_sex or (status or {}).get("declared_sex")
                             or manifest_value(manifest, "run", "declared_sex")
                             or sections["sample_qc"]["data"].get("declared_sex") or None),
        },
        "databases": {
            "clinvar_release": clinvar_date,
            "hla_database": manifest_value(manifest, "data", "hla_database"),
        },
        "manifest": manifest or [],
        "sections": sections,
        "cyp2d6": cyp2d6,
        "not_run": [sections[k]["title"] + (f" (step {sections[k]['step']} failed)" if sections[k]["state"] == "failed" else "")
                    for k, *_ in SECTIONS if sections[k]["state"] in ("missing", "failed")],
        "not_assessed": NOT_ASSESSED + ([] if sections["prs"]["data"].get("adjusted") else [PRS_NOT_ADJUSTED]),
    }
    return summary


def sample_qc_main(argv):
    ap = argparse.ArgumentParser(prog="collect_summary.py sample-qc",
                                 description="Write the step 33 table (key, value) for one sample.")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--somalier-samples", required=True, help="somalier relate's samples.tsv")
    ap.add_argument("--somalier-pairs", help="somalier relate's pairs.tsv (other samples of the run)")
    ap.add_argument("--somalier-id", help="the sample's name in somalier's files (the BAM's @RG SM in the bash step, the samplesheet id in Nextflow; default --sample)")
    ap.add_argument("--selfsm", help="VerifyBamID2's .selfSM")
    ap.add_argument("--declared-sex", choices=["male", "female", ""], default="")
    ap.add_argument("--freemix-warn", type=float, default=FREEMIX_WARN)
    ap.add_argument("--marker-check", choices=["passed", "skipped", ""], default="",
                    help="skipped: VerifyBamID2 ran with --DisableSanityCheck (fewer than 1,000 markers had reads)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    try:
        rows = sample_qc_table(a.sample, a.somalier_samples, a.selfsm, a.somalier_pairs, a.declared_sex,
                               a.freemix_warn, a.somalier_id, a.marker_check)
    except (Unreadable, OSError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    tmp = a.out + ".tmp"
    with open(tmp, "w") as f:
        f.write("key\tvalue\n")
        for k, v in rows:
            f.write(f"{k}\t{v}\n")
    os.replace(tmp, a.out)
    for k, v in rows:
        print(f"  {k}: {v}")
    return 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ["sample-qc"]:
        return sample_qc_main(argv[1:])
    if argv[:1] == ["prs-format"]:
        return prs_format_main(argv[1:])
    if argv[:1] == ["prs-table"]:
        return prs_table_main(argv[1:])
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sample", required=True)
    ap.add_argument("--sample-dir", required=True)
    ap.add_argument("--declared-sex", help="the samplesheet's sex, when no run status or manifest records it")
    ap.add_argument("--out", help="summary JSON (default: SAMPLE_DIR/summary.json)")
    a = ap.parse_args(argv)
    summary = collect(a.sample, a.sample_dir, a.declared_sex)
    out = a.out or os.path.join(a.sample_dir, "summary.json")
    tmp = out + ".tmp"
    with open(tmp, "w") as f:
        json.dump(summary, f, indent=1)
        f.write("\n")
    os.replace(tmp, out)
    print(f"Summary written: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
