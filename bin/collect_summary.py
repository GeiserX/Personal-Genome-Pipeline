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
"""
import argparse
import csv
import glob
import gzip
import json
import os
import sys
import zlib
from datetime import datetime, timezone

SCHEMA_VERSION = 1

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
                 "27", "29", "30", "31", "32", "04b"}

# What this pipeline does not assess, whatever ran.
NOT_ASSESSED = [
    "Mosaic and low-fraction variants (only chrM heteroplasmy is reported)",
    "Methylation and imprinting",
    "Repeat expansions outside the ExpansionHunter catalog",
    "Phase: two variants in one gene are compound-het candidates, not confirmed",
    "Polygenic scores against an ancestry-matched reference (raw scores only)",
    "Variants of uncertain significance (not interpreted)",
]

EH_LOCI = ["HTT", "FMR1", "C9ORF72", "ATXN1", "DMPK"]


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
               "parse_failed": failed or not genes,
               "warnings": warnings}


def sec_pypgx(d, s):
    p = first_existing(d, [f"pypgx/{s}_pypgx_summary.tsv"])
    if not p:
        return None, {}
    rows = read_tsv(p)
    called = [r for r in rows if (r.get("Diplotype") or "") not in ("", "FAILED", "N/A")]
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
            alleles = [a for a in (c[2] if len(c) > 2 else "", c[5] if len(c) > 5 else "") if a and a != "."]
            quals = [q for q in (c[4] if len(c) > 4 else "", c[7] if len(c) > 7 else "") if q and q != "."]
            loci.append({"gene": c[0], "alleles": alleles, "quality": quals})
    return p, {"loci": loci}


def sec_prs(d, s):
    p = first_existing(d, [f"prs/{s}_prs_summary.tsv"])
    if not p:
        return None, {}
    rows = []
    for r in read_tsv(p):
        rows.append({"condition": r.get("Condition", ""), "pgs_id": r.get("PGS_ID", ""),
                     "score": r.get("Score_SUM", ""), "matched": r.get("Variants_Matched", ""),
                     "total": r.get("Variants_Total", "")})
    return p, {"scores": rows}


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
    p = first_existing(d, [f"expansion_hunter/{s}_eh.vcf", f"expansion_hunter/{s}.vcf",
                           "expansion_hunter/*_eh.vcf"])
    if not p:
        return None, {}
    loci, tested = {}, 0
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
    return p, {"records": tested, "key_loci": [{"locus": k, "repeat_count": loci.get(k, "not in output")}
                                                for k in EH_LOCI]}


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
    return p, {"haplogroup": hg or "."}


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


def sec_clinical(d, s):
    p = first_existing(d, [f"clinical/{s}_clinical_summary.tsv"])
    if not p:
        # The Nextflow report gets the clinical VCF, not the summary table:
        # the count only (one record per variant, as in the table).
        v = first_existing(d, [f"clinical/{s}_clinical.vcf.gz"])
        if not v:
            return None, {}
        return v, {"variants": sum(1 for _ in vcf_records(v))}
    rows = read_tsv(p)
    by = {}
    for r in rows:
        by[r.get("IMPACT") or "."] = by.get(r.get("IMPACT") or ".", 0) + 1
    genes = sorted({r.get("GENE") for r in rows if r.get("GENE") not in (None, "", ".")})
    return p, {"variants": len(rows), "by_impact": by, "genes": len(genes)}


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
    if not declared:
        sex_check, why = "not_checked", "no declared sex"
    elif row.get("sex", "") == "-2":
        sex_check = "not_checked"
        why = ("chrX is heterozygous like a female sample but chrY has reads: "
               "a sex-chromosome aneuploidy or a mixed sample")
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

    # CYP2D6 from the three callers, side by side.
    pc = next((g for g in sections["cpic"]["data"].get("genes", []) if g["gene"] == "CYP2D6"), None)
    calls = {
        "PharmCAT": (pc["diplotype"] if pc["status"] not in ("not called", "ambiguous") else pc["status"]) if pc else None,
        "pypgx": sections["pypgx"]["data"].get("cyp2d6") if sections["pypgx"]["state"] == "ok" else None,
        "Cyrius": sections["cyrius"]["data"].get("genotype") if sections["cyrius"]["state"] == "ok" else None,
    }
    called = {k: v for k, v in calls.items()
              if v and v not in ("not called", "ambiguous", "none", "None/None", "FAILED", "N/A")}
    norm = {k: "/".join(sorted(v.split("/"))) for k, v in called.items()}
    cyp2d6 = {"calls": calls,
              "agree": (len(set(norm.values())) == 1) if len(norm) >= 2 else None}

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
        "not_run": [sections[k]["title"] for k, *_ in SECTIONS if sections[k]["state"] == "missing"],
        "not_assessed": NOT_ASSESSED,
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
