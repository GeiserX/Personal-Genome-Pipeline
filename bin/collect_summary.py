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
                                [--keep-x] [--keep-ambiguous --reference FASTA]
  collect_summary.py prs-table --sample S --results DIR --sampleset NAME --scores DIR --input-kind gvcf --out F
                               [--sex male|female]

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
import base64
import csv
import glob
import gzip
import json
import os
import re
import struct
import sys
import zlib
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pgx_outside_calls  # noqa: E402  step 36's rule for a CYP2D6 call and an HLA allele
import ancestry_fine_report  # noqa: E402  step 26b's group names and painting colours
import paraphase_report  # noqa: E402  step 35b's tables and its SMN1 rule beside Parascopy
import clip_check  # noqa: E402  CPSR's variant ids and class names, as step 17's read check reads them
import annotsv_postfilter  # noqa: E402  step 05's class 4-5 checks and the step 22 support join

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

# The call floor of step 6 (scripts/06-clinvar-screen.sh writes a floor other
# than this default into the hits VCF as ##PGP_CALL_FLOOR=GQ=..,DP=..): a P/LP
# hit below either is a low-confidence call, listed but counted apart. A hits
# file without the line used this default, or was written before the floor.
CALL_FLOOR = {"GQ": 30, "DP": 10}
# gnomAD allele frequency above which a P/LP hit is marked common.
COMMON_AF = 0.01
# The genotype check of the clinical filter (bin/csq_rank.awk, step 23 and the
# CLINICAL_FILTER module): step 6's GQ and DP floor, and the alt read fraction
# (AD) a heterozygote must lie within and a homozygous alt call must reach.
# bin/csq_rank.awk writes the floors it used into the clinical VCF as
# ##PGP_CALL_FLOOR=GQ=..,DP=..,HET_AB=lo-hi,HOM_AB=..; these are its defaults.
AB_FLOOR = {"HET_AB": (0.2, 0.8), "HOM_AB": 0.8}
IMPACT_RANK = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}

# Steps run-all.sh runs and records in logs/run_status.tsv. A section whose
# step is not here (variant calling, step 03, which run-all reuses on purpose)
# is never marked stale.
RUN_ALL_STEPS = {"04", "05", "06", "07", "08", "09", "09b", "10", "11", "12", "13", "14",
                 "15", "16", "16b", "17", "18", "19", "20", "21", "22", "23", "25", "26", "26b",
                 "27", "29", "30", "31", "32", "35", "35b", "37", "40", "41", "42", "47", "04b"}

# What this pipeline does not assess, whatever ran.
NOT_ASSESSED = [
    "Mosaic and low-fraction variants (only chrM heteroplasmy is reported)",
    "Methylation and imprinting",
    "Repeat expansions outside the ExpansionHunter catalog",
    "Phase: two variants in one gene are compound-het candidates, not confirmed",
    "Variants of uncertain significance (not interpreted)",
]

# Added to NOT_ASSESSED unless step 35's depth count (bin/hba_copy_number.py,
# short reads) or step 35b (Paraphase, PacBio HiFi, both genes) made a call.
HBA_NOT_ASSESSED = ("HBA1/HBA2 copy number (alpha-thalassaemia deletions): neither step 35's read-depth count "
                    "(short reads) nor step 35b (Paraphase, PacBio HiFi) made a call")

# Added to NOT_ASSESSED unless step 25 ran with an ancestry reference panel.
PRS_NOT_ADJUSTED = "Polygenic score percentiles: no ancestry reference panel was used, so the scores are raw sums"

EH_LOCI = ["HTT", "FMR1", "C9ORF72", "ATXN1", "DMPK"]

# ACMG SF v3.3 (Lee et al., Genet Med 2025;27(8):101454): the 84 genes in which
# the ACMG recommends reporting pathogenic and likely pathogenic variants found
# by chance, and the reporting rule of each gene (Table 1), are the rows of
# assets/acmg_sf.tsv. A new list version means new rows there and a new
# ACMG_SF_VERSION. The sample's carrier list reads assets/carrier_genes.tsv.
# Both are read from ASSETS_DIR: the repository's assets/ folder beside bin/,
# or the folder render_report.py --assets names (the containers of step 24 and
# of HTML_REPORT mount bin/ only).
ACMG_SF_VERSION = "ACMG SF v3.3"
ASSETS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "assets")
SF_REPORT_IF = ("any", "biallelic", "hom_variant", "hemizygous", "het_female_too")
SEVERITY_KINDS = ("hypomorph", "pseudodeficiency")


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


def clinvar_variants(v, sep):
    """'chr1:10:C:T|chr1:15:G:A' -> 'chr1:10 C>T<sep>chr1:15 G>A'."""
    out = []
    for x in v.split("|"):
        f = (x.split(":") + ["", "", "", ""])[:4]
        out.append(f"{f[0]}:{f[1]} {f[2]}>{f[3]}")
    return sep.join(out)


def clinvar_combined_note(info, gt=""):
    """Step 6's note on a hit that is a combined allele of two neighbouring
    records, or part of one, from the tags bin/clinvar_join_adjacent.awk
    writes; bin/clinvar_hits.awk words it the same way. "." when none.
    On an unphased het call (0/1) the genotype cannot say the two parts share
    a haplotype, so the note says they may not."""
    a = gt.split("/")
    phase = ("; unphased, so the parts may sit on different haplotypes"
             if len(a) == 2 and all(x.isdigit() for x in a) and a[0] != a[1] else "")
    if info.get("PGP_JOINED_FROM"):
        return "combined allele of " + clinvar_variants(info["PGP_JOINED_FROM"], " + ") + phase
    if info.get("PGP_PART_OF"):
        return "part of a combined allele: " + clinvar_variants(info["PGP_PART_OF"], ", ") + phase
    if info.get("PGP_ONE_ALLELE"):
        return ("possibly one allele with " + clinvar_variants(info["PGP_ONE_ALLELE"], ", ")
                + ": check ClinVar for a combined record")
    return "."


def call_floor(path):
    """Step 6's call floor from its hits VCF header, CALL_FLOOR without the line."""
    floor = dict(CALL_FLOOR)
    try:
        with open_text(path) as f:
            for line in f:
                if not line.startswith("##"):
                    break
                if line.startswith("##PGP_CALL_FLOOR="):
                    for k, v in re.findall(r"(GQ|DP)=(\d+)", line):
                        floor[k] = int(v)
    except (OSError, EOFError, gzip.BadGzipFile, zlib.error) as e:
        raise Unreadable(f"{os.path.basename(path)}: {e}") from e
    return floor


def low_confidence(dp, gq, floor):
    """True when the call's GQ or DP is known and below the floor."""
    return (gq is not None and gq < floor["GQ"]) or (dp is not None and dp < floor["DP"])


def clinvar_row(r, floor=None):
    info = info_map(r[7]) if len(r) > 7 else {}
    rev = info.get("CLNREVSTAT", "")
    genes = [g.split(":")[0] for g in info.get("GENEINFO", "").split("|") if g]
    gt = r[9].split(":")[0] if len(r) > 9 else ""
    zyg = {"0/1": "het", "1/0": "het", "0|1": "het", "1|0": "het", "1/1": "hom", "1|1": "hom"}.get(gt, gt)
    fmt = dict(zip(r[8].split(":"), r[9].split(":"))) if len(r) > 9 else {}
    dp, gq = (int(fmt[k]) if fmt.get(k, "").isdigit() else None for k in ("DP", "GQ"))
    return {
        "stars": STARS.get(rev, 0),
        "chrom": r[0], "pos": int(r[1]), "ref": r[3], "alt": r[4],
        "genotype": zyg,
        "gene": ",".join(genes) or ".",
        "significance": info.get("CLNSIG", "").replace("_", " ") or ".",
        "review_status": rev.replace("_", " ") or ".",
        "combined_allele": clinvar_combined_note(info, gt),
        "dp": dp, "gq": gq,
        "low_confidence": low_confidence(dp, gq, floor or CALL_FLOOR),
    }


def clinvar_cdna_checks(d, s, p=None):
    """Step 06's cDNA read check (bin/clip_check.py), keyed by (chrom, pos,
    ref, alt), or the table at p. None when the file is absent: the Nextflow
    pipeline does not write it."""
    p = p or first_existing(d, [f"clinvar/{s}_clinvar_cdna_check.tsv"])
    if not p:
        return None
    try:
        rows = read_tsv(p)
    except (OSError, EOFError, gzip.BadGzipFile, zlib.error) as e:
        raise Unreadable(f"{os.path.basename(p)}: {e}") from e

    def num(v):
        return int(v) if (v or "").isdigit() else None
    out = {}
    for r in rows:
        if not (r.get("pos") or "").isdigit():
            continue
        out[(r.get("chrom"), int(r["pos"]), r.get("ref"), r.get("alt"))] = {
            "check": r.get("cdna_check") or "not_checked", "note": r.get("cdna_note") or "",
            "alt_reads": num(r.get("alt_reads")), "alt_clipped": num(r.get("alt_clipped")),
            "ref_reads": num(r.get("ref_reads")), "ref_clipped": num(r.get("ref_clipped"))}
    return out


def tabix_linear_index(path):
    """{contig: linear index} of a bgzipped VCF's .tbi: for each 16 kb window,
    the virtual file offset of the first record that can overlap it (standard
    library only, as the rest of this script)."""
    with gzip.open(path + ".tbi", "rb") as f:
        b = f.read()
    if b[:4] != b"TBI\x01":
        raise ValueError(f"{os.path.basename(path)}.tbi is not a tabix index")
    n_ref = struct.unpack_from("<i", b, 4)[0]
    l_nm = struct.unpack_from("<i", b, 32)[0]
    names = [x.decode() for x in b[36:36 + l_nm].split(b"\0")[:n_ref]]
    off, out = 36 + l_nm, {}
    for name in names:
        n_bin = struct.unpack_from("<i", b, off)[0]
        off += 4
        for _ in range(n_bin):
            n_chunk = struct.unpack_from("<i", b, off + 4)[0]
            off += 8 + 16 * n_chunk
        n_intv = struct.unpack_from("<i", b, off)[0]
        off += 4
        out[name] = struct.unpack_from(f"<{n_intv}Q", b, off)
        off += 8 * n_intv
    return out


def tabix_lines(path, index, chrom, pos):
    """The data lines of a bgzipped VCF that start at chrom:pos, read from the
    offset its tabix linear index gives for pos's window. A window with no
    offset recorded (0, or all ones) reads from the start of the file and
    skips the other contigs, which is slower but still right."""
    ioff = index.get(chrom)
    if not ioff:
        return []
    voff = ioff[min((pos - 1) >> 14, len(ioff) - 1)]
    if voff == 0xFFFFFFFFFFFFFFFF:
        voff = 0
    out, seen = [], False
    with open(path, "rb") as raw:
        raw.seek(voff >> 16)
        g = gzip.GzipFile(fileobj=raw)
        g.read(voff & 0xFFFF)
        for line in g:
            if line.startswith(b"#"):
                continue
            f = line.split(b"\t", 2)
            if len(f) < 3 or not f[1].isdigit():
                continue
            if f[0].decode() != chrom:
                if seen:
                    break
                continue
            seen = True
            if int(f[1]) > pos:
                break
            if int(f[1]) == pos:
                out.append(line.decode(errors="replace").rstrip("\n").split("\t"))
    return out


def vep_gnomad_af(vep, hits):
    """gnomAD exome and genome AF of each hit from step 13's VEP output, as
    {(chrom, pos, ref, alt): (exome, genome)}; a hit with no record is left
    out. The fields are VEP's gnomADe_AF and gnomADg_AF (--everything)."""
    fmt = None
    with open_text(vep) as f:
        for line in f:
            if not line.startswith("#"):
                break
            if line.startswith("##INFO=<ID=CSQ"):
                m = re.search(r"Format: ([^\"]+)", line)
                fmt = m.group(1).split("|") if m else None
    if not fmt or "Allele" not in fmt:
        return {}
    cols = [fmt.index(c) if c in fmt else None for c in ("gnomADe_AF", "gnomADg_AF")]

    def num(v):
        vals = []
        for x in (v or "").split("&"):
            try:
                vals.append(float(x))
            except ValueError:
                pass
        return max(vals) if vals else None
    out = {}
    keys = sorted({(h["chrom"], h["pos"], h["ref"], h["alt"]) for h in hits})
    index = tabix_linear_index(vep) if keys else {}
    for key in keys:
        chrom, pos, ref, alt = key
        for r in tabix_lines(vep, index, chrom, pos):
            alts = r[4].split(",") if len(r) > 7 else []
            if r[3] != ref or alt not in alts:
                continue
            # VEP's Allele drops the base an indel's alleles share.
            snv = len(ref) == 1 and all(len(a) == 1 for a in alts)
            allele = alt if snv or any(a[:1] != ref[:1] for a in alts) else (alt[1:] or "-")
            got = [None, None]
            for tr in info_map(r[7]).get("CSQ", "").split(","):
                v = tr.split("|")
                if len(v) != len(fmt) or (len(alts) > 1 and v[fmt.index("Allele")] != allele):
                    continue
                for i, c in enumerate(cols):
                    if c is not None and got[i] is None:
                        got[i] = num(v[c])
            out[key] = tuple(got)
    return out


def clinvar_frequencies(d, s, rows):
    """Put the gnomAD AF from step 13's VEP output on each ClinVar hit, and
    say where it came from or why there is none."""
    vep = first_existing(d, [f"vep/{s}_vep.vcf.gz"])
    note = None
    if not vep:
        note = "no VEP output (step 13) in this folder"
    elif not os.path.exists(vep + ".tbi"):
        note = f"{os.path.relpath(vep, d)} has no .tbi index"
    found = {}
    if note is None:
        try:
            found = vep_gnomad_af(vep, rows)
        except (OSError, EOFError, ValueError, struct.error, gzip.BadGzipFile, zlib.error) as e:
            note = f"{os.path.relpath(vep, d)} could not be read: {e}"
    for x in rows:
        e, g = found.get((x["chrom"], x["pos"], x["ref"], x["alt"]), (None, None))
        known = [v for v in (e, g) if v is not None]
        x["gnomad_exome_af"], x["gnomad_genome_af"] = e, g
        x["af"] = max(known) if known else None
        x["common"] = x["af"] is not None and x["af"] > COMMON_AF
    return {"source": os.path.relpath(vep, d) if note is None else None, "note": note, "common_af": COMMON_AF}


def sec_clinvar(d, s):
    p = first_existing(d, [f"clinvar/{s}_clinvar_hits.vcf"])
    if not p:
        return None, {}
    floor = call_floor(p)
    rows = [clinvar_row(r, floor) for r in vcf_records(p) if len(r) >= 8]
    rows.sort(key=lambda x: (-x["stars"], x["chrom"], x["pos"]))
    # The headline count and the stars leave the low-confidence calls out;
    # they stay in the hits, marked, and are counted on their own.
    sure = [x for x in rows if not x["low_confidence"]]
    by = {str(k): sum(1 for x in sure if x["stars"] == k) for k in (4, 3, 2, 1, 0)}
    data = {"count": len(sure), "by_stars": by, "hits": rows, "low_confidence": len(rows) - len(sure),
            "call_floor": floor, "frequency": clinvar_frequencies(d, s, rows)}
    # The cDNA read check, on each hit when step 06 wrote it.
    cdna = clinvar_cdna_checks(d, s)
    if cdna is not None:
        for x in rows:
            x["cdna"] = cdna.get((x["chrom"], x["pos"], x["ref"], x["alt"])) or {
                "check": "not_checked", "note": "not in the read-check table", "alt_reads": None,
                "alt_clipped": None, "ref_reads": None, "ref_clipped": None}
        data["cdna_flagged"] = sum(1 for x in rows if x["cdna"]["check"] == "flagged")
    return p, data


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


# The rows of the CYP2D6 depth check (bin/cyp2d6_depth_check.py) the reports
# show: its verdict and the MAPQ shares beside it.
DEPTH_CHECK_KEYS = ("status", "gene_mapq1_share", "gene_mapq20_share", "body_mapq1_share", "body_mapq20_share")


def sec_pypgx(d, s):
    p = first_existing(d, [f"pypgx/{s}_pypgx_summary.tsv"])
    if not p:
        return None, {}
    rows = read_tsv(p)
    # A no-call by step 36's list (Indeterminate, FAILED, N/A...), in any case
    called = [r for r in rows if (r.get("Diplotype") or "").strip().lower() not in pgx_outside_calls.NO_CALL]
    cyp_row = next((r for r in rows if r.get("Gene") == "CYP2D6"), None)
    cyp = (cyp_row or {}).get("Diplotype", "")
    out = {"genes_total": len(rows), "genes_called": len(called), "cyp2d6": cyp or "not called"}
    # pypgx's copy-number call from read depth, BAM-based genes only; 'N/A'
    # (a VCF gene, GSTT1 without its contig, CYP2D6 after a failed depth check)
    # is none.
    out["cnv_call"] = {r["Gene"]: r["CNV_call"].strip() for r in rows
                       if (r.get("Source") or "").strip().upper() == "BAM"
                       and (r.get("CNV_call") or "").strip() not in ("", ".")
                       and not r["CNV_call"].strip().startswith("N/A")}
    out["not_callable"] = [r["Gene"] for r in rows
                           if (r.get("Status") or "").strip() == pgx_outside_calls.NOT_CALLABLE]
    out["low_gq"] = [{"gene": r["Gene"], "diplotype": (r.get("Diplotype") or "").strip(),
                      "min_gq": int(r["Min_GQ"]) if (r.get("Min_GQ") or "").strip().isdigit() else None,
                      "weakest": (r.get("Weakest_variant") or "").strip()}
                     for r in rows if (r.get("Quality") or "").strip() == pgx_outside_calls.LOW_GQ]
    if cyp_row is not None:
        out["cyp2d6_evidence"] = {"phenotype": (cyp_row.get("Phenotype") or "").strip(),
                                  "defining_filter": (cyp_row.get("Defining_filter") or "").strip(),
                                  "quality": (cyp_row.get("Quality") or "").strip()}
    dc = os.path.join(os.path.dirname(p), f"{s}_cyp2d6_depth_check.tsv")
    if os.path.isfile(dc):
        with open(dc, errors="replace") as f:
            kv = dict(line.rstrip("\n").split("\t", 1) for line in f if "\t" in line)
        out["depth_check"] = {k: kv[k] for k in DEPTH_CHECK_KEYS if k in kv}
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


# The Paraphase regions the report shows gene by gene, in this order; the
# others are counted.
PARALOG_REGIONS = ("smn1", "hba", "rccx", "GBA", "pms2", "strc", "ncf1", "OTOA")


def sec_paralogs(d, s):
    """Step 35b (Paraphase, PacBio HiFi): the genes of PARALOG_REGIONS, a
    count of the other regions, the variant rows per region, and SMN1 beside
    step 35's Parascopy call when that table is there too."""
    p = first_existing(d, [f"paralogs/{s}_paraphase_copy_number.tsv"])
    if not p:
        return None, {}
    rows = read_tsv(p)
    if rows and not {"region", "gene", "copy_number", "status"} <= set(rows[0]):
        raise Unreadable(f"{os.path.basename(p)} lacks the region, gene, copy_number or status column")
    genes = [{"region": r["region"], "gene": r["gene"], "copy_number": r["copy_number"],
              "status": r["status"], "detail": r.get("detail") or "."}
             for r in rows if r["region"] in PARALOG_REGIONS]
    genes.sort(key=lambda g: PARALOG_REGIONS.index(g["region"]))
    others = [r for r in rows if r["region"] not in PARALOG_REGIONS]
    variants = {}
    vp = os.path.join(os.path.dirname(p), f"{s}_paraphase_variants.tsv")
    if os.path.isfile(vp):
        for r in read_tsv(vp):
            variants[r.get("region", "")] = variants.get(r.get("region", ""), 0) + 1
    smn1 = next((g["copy_number"] for g in genes if g["gene"] == "SMN1"), None)
    smn1 = int(smn1) if smn1 and smn1.isdigit() else None
    parascopy = None
    pp = os.path.join(os.path.dirname(p), f"{s}_smn_copy_number.tsv")
    if os.path.isfile(pp):
        parascopy = paraphase_report.parascopy_smn1(pp)
    return p, {"genes": genes,
               "other_regions": len(others),
               "other_regions_called": sum(1 for r in others if r["status"] == "called"),
               "variants": variants,
               "smn1": {"paraphase": smn1, "parascopy": parascopy,
                        "agree": paraphase_report.compare_smn1(smn1, parascopy)}}


# SMN2 exon 7, c.840 (the base that tells SMN2 from SMN1), GRCh38 1-based;
# SMN1's is paraphase_report.SMN1_EXON7.
SMN2_EXON7 = ("chr5", 70076526)


def parascopy_cn_at(rows, chrom, pos):
    """The copy number at one position from step 35's table: the
    paralog-specific value (psCN) of the row, or of its homologous region,
    that holds the position, with the row's aggregate value (agCN). As
    paraphase_report.parascopy_smn1: cn is None unless psCN's filter is PASS
    and its quality at least 20 (Parascopy's threshold)."""
    for r in rows:
        try:
            regions = [(r["chrom"], int(r["start"]) + 1, int(r["end"]))]
        except (KeyError, ValueError, TypeError):  # TypeError: a short row (None fields)
            continue
        for h in (r.get("homologous_regions") or "*").split(","):
            m = re.match(r"^([^:]+):(\d+)-(\d+)", h.strip())
            if m:
                regions.append((m.group(1), int(m.group(2)), int(m.group(3))))
        idx = next((i for i, (c, b, e) in enumerate(regions) if c == chrom and b <= pos <= e), None)
        if idx is None:
            continue
        cns = (r.get("psCN") or "").split(",")
        quals = (r.get("psCN_qual") or "").split(",")
        cn_txt = cns[idx] if idx < len(cns) else "?"
        q_txt = quals[idx] if idx < len(quals) else "*"
        filt = r.get("psCN_filter") or "."
        try:
            q = float(q_txt)
        except ValueError:
            q = None
        ok = filt == "PASS" and q is not None and q >= 20 and cn_txt.isdigit()
        return {"cn": int(cn_txt) if ok else None,
                "text": f"{cn_txt} (quality {q_txt}, filter {filt})" + ("" if ok else ", not resolved"),
                "agcn": r.get("agCN") or ".", "agcn_qual": r.get("agCN_qual") or ".",
                "agcn_filter": r.get("agCN_filter") or "."}
    return {"cn": None, "text": "no row covers exon 7", "agcn": ".", "agcn_qual": ".", "agcn_filter": "."}


def smn_reading(smn1):
    """What an SMN1 copy number means for SMA, in one line."""
    if smn1 is None:
        return "SMN1 not resolved (quality under 20 or filter not PASS): no carrier reading"
    if smn1 == 0:
        return "no SMN1 copy: the pattern of spinal muscular atrophy (SMA) itself"
    if smn1 == 1:
        return "one SMN1 copy: the pattern of an SMA carrier"
    if smn1 == 2:
        return ("two SMN1 copies: not a carrier by copy number; a 2+0 carrier (both copies on one chromosome) "
                "cannot be ruled out")
    return f"{smn1} SMN1 copies: not a carrier by copy number"


def sec_parascopy(d, s):
    """Step 35 (Parascopy, short reads): SMN1 and SMN2 copy number at exon 7,
    the aggregate copy number of the pair there, and the SMA reading."""
    p = first_existing(d, [f"paralogs/{s}_smn_copy_number.tsv"])
    if not p:
        return None, {}
    rows = read_tsv(p)
    if rows and not {"chrom", "start", "end", "agCN", "psCN", "psCN_qual", "psCN_filter"} <= set(rows[0]):
        raise Unreadable(f"{os.path.basename(p)} lacks the chrom, start, end, agCN, psCN, psCN_qual or "
                         "psCN_filter column")
    smn1 = parascopy_cn_at(rows, *paraphase_report.SMN1_EXON7)
    smn2 = parascopy_cn_at(rows, *SMN2_EXON7)
    hp = os.path.join(os.path.dirname(p), f"{s}_hba_copy_number.tsv")
    return p, {"smn1": smn1, "smn2": {"cn": smn2["cn"], "text": smn2["text"]},
               "agcn": smn1["agcn"], "agcn_qual": smn1["agcn_qual"], "agcn_filter": smn1["agcn_filter"],
               "rows": len(rows),
               # HBA1/HBA2 from depth (bin/hba_copy_number.py); None when
               # step 35 wrote no HBA table (a run from before the check)
               "hba": read_hba(hp) if os.path.isfile(hp) else None,
               "reading": ("no row of step 35's table covers SMN1 exon 7: no carrier reading"
                           if smn1["text"] == "no row covers exon 7" else smn_reading(smn1["cn"]))}


def read_hba(path):
    """Step 35's HBA table (metric, value) as a dict; a table without a
    status row reads as not assessed, with that reason."""
    out = {}
    with open(path, errors="replace") as f:
        for line in f:
            c = line.rstrip("\n").split("\t", 1)
            if len(c) == 2 and c[0] != "metric":
                out[c[0]] = c[1]
    if out.get("status") not in ("called", "not assessed"):
        out = {"status": "not assessed", "message": f"not assessed: {os.path.basename(path)} has no status row"}
    return out


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
            loci.append(hla_locus_extras({"gene": c[0], "alleles": alleles, "quality": quals, "low_confidence": low,
                                          "withheld_from_pharmcat": low and c[0] in pgx_outside_calls.HLA_GENES},
                                         c[1] if len(c) > 1 else ""))
    out = {"loci": loci, "derived": hla_derived(loci)}
    ot, kir = hla_optitype(d, s), kir_genotypes(d, s)
    if ot is not None:
        out["optitype"] = ot
    if kir is not None:
        out["kir"] = kir
    return p, out


def first_comment(path):
    """The first '#' line of a file step 08 wrote instead of a result, without the '#'."""
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                return line.lstrip("#").strip()
    return ""


def hla_optitype(d, s):
    """Step 08's second HLA caller (HLA_QC=true): per class I gene the row of
    its concordance table (pgx_outside_calls.py hla-concordance), or, when
    there is no table, why OptiType's result file says it did not run. None
    when OptiType was not asked for."""
    conc = first_existing(d, [f"hla_optitype/{s}_hla_concordance.tsv"])
    res = first_existing(d, [f"hla_optitype/{s}_optitype_result.tsv"])
    if not conc and not res:
        return None
    if not conc:
        return {"genes": [], "reason": first_comment(res) or "no concordance table: step 08 did not compare the callers"}
    genes = []
    with open(conc, errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if (r.get("Gene") or "").strip():
                genes.append({"gene": r["Gene"].strip(), "t1k": (r.get("T1K") or "").strip(),
                              "optitype": (r.get("OptiType") or "").strip(), "agree": (r.get("Agree") or "").strip(),
                              "note": (r.get("Note") or "").strip()})
    return {"genes": genes, "reason": "" if genes else "the concordance table has no rows"}


def kir_genotypes(d, s):
    """Step 08's KIR pass (KIR=true): the genes T1K gave an allele, with each
    allele's quality, and the genes it gave none; or the line step 08 wrote
    for a BAM with too few KIR reads. None when KIR was not typed."""
    p = first_existing(d, [f"kir_t1k/{s}_kir_genotype.tsv", f"kir/{s}_kir_genotype.tsv"])
    if not p:
        return None
    typed, not_called = [], []
    with open(p, errors="replace") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3 or c[0].startswith("#"):
                continue
            alleles, quals = [], []
            for i in (2, 5):
                if len(c) > i and c[i] not in ("", "."):
                    alleles.append(c[i])
                    quals.append(c[i + 2] if len(c) > i + 2 else "")
            if alleles:
                typed.append({"gene": c[0], "alleles": alleles, "quality": quals})
            else:
                not_called.append(c[0])
    reason = "" if typed or not_called else (first_comment(p) or "the KIR genotype file has no gene rows")
    return {"typed": typed, "not_called": not_called, "reason": reason}


def hla_quality(q):
    """T1K's quality as a number; one that is not a number counts as -1, as in step 36."""
    try:
        return float(q)
    except ValueError:
        return -1.0


# --- HLA: candidate caps, the HFE label, the DQ haplotypes and the allele table ---
# Rules in docs/08-hla-typing.md "What the report derives". The allele table is
# bin/hla_alleles.tsv (in bin/ because the report containers mount only bin/).

HLA_ALLELES_TSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hla_alleles.tsv")
HLA_MAX_CANDIDATES = 3   # a quality-0 allele with more candidates is "no confident call"
HLA_HFE_LABEL = ("T1K sequence allele, not the HFE C282Y/H63D genotype: that genotype comes from the "
                 "variant calls (ClinVar section)")


def hla_name(x):
    """'HLA-DQB1*02:01:01' or 'DQB1*02:01:01' -> 'DQB1*02:01:01'; 'HLA-DQB1' -> 'DQB1'."""
    return x.strip()[4:] if x.strip().startswith("HLA-") else x.strip()


def hla_candidates(allele):
    """T1K's allele field: one allele, or equally likely candidates separated by commas."""
    return [a for a in (hla_name(x) for x in allele.split(",")) if a]


def hla_locus_extras(locus, n_diff):
    """Cap a quality-0 candidate list to HLA_MAX_CANDIDATES (the full list stays
    in the TSV), count the candidates, mark one homozygous call and label the
    HFE row."""
    counts, kept = [], []
    for a, q in zip(locus["alleles"], locus["quality"] + [""] * len(locus["alleles"])):
        cands = hla_candidates(a)
        counts.append(len(cands))
        kept.append(",".join(cands[:HLA_MAX_CANDIDATES]) if hla_quality(q) <= 0 and len(cands) > HLA_MAX_CANDIDATES
                    else a)
    locus["alleles"] = kept
    locus["candidates"] = counts
    # T1K writes one allele and '. 0 -1' for a homozygous gene (or a second allele it did not see).
    locus["one_allele_type"] = n_diff.strip() == "1" and len(kept) == 1
    if hla_name(locus["gene"]) == "HFE":
        locus["label"] = HLA_HFE_LABEL
    return locus


def hla_match(allele, group):
    """'yes' when every candidate of the allele is in the group (B*57:01 holds
    B*57:01:01:02), 'no' when none is, 'mixed' when only some are."""
    hits = [c == group or c.startswith(group + ":") for c in hla_candidates(allele)]
    return "yes" if hits and all(hits) else ("mixed" if any(hits) else "no")


def hla_read_table(path=HLA_ALLELES_TSV):
    """bin/hla_alleles.tsv as a list of dicts (tests/test_hla_report.py checks its schema)."""
    rows, head = [], None
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if head is None:
                head = c
                continue
            rows.append(dict(zip(head, c)))
    return rows


def hla_derived(loci, table=None):
    """The DQ haplotypes, the DQB1*02 dose and each table allele present or absent."""
    by = {}
    for l in loci:
        by.setdefault(hla_name(l["gene"]), l)

    def usable(gene):
        l = by.get(gene)
        if not l or not l["alleles"]:
            return None, f"HLA-{gene} not typed"
        if l["low_confidence"]:
            return None, f"HLA-{gene} low confidence"
        return l, ""

    def carries(l, group):
        """'yes', 'no' or 'mixed' over both alleles of a locus."""
        m = [hla_match(a, group) for a in l["alleles"]]
        return "yes" if "yes" in m else ("mixed" if "mixed" in m else "no")

    def both(a, b):
        return "present" if a == b == "yes" else ("absent" if "no" in (a, b) else "uncertain")

    dqa, ra = usable("DQA1")
    dqb, rb = usable("DQB1")
    coeliac = {"status": "not_determined", "reason": "; ".join(r for r in (ra, rb) if r),
               "dq2_5": None, "dq2_5_phase": None, "dq2_2": None, "dq8": None, "dqb1_02_dose": None}
    if dqb:
        m = [hla_match(a, "DQB1*02") for a in dqb["alleles"]]
        if "mixed" not in m:
            n = m.count("yes")
            coeliac["dqb1_02_dose"] = 2 * n if dqb["one_allele_type"] else n
    if dqa and dqb:
        coeliac.update(status="derived", reason="")
        coeliac["dq2_5"] = both(carries(dqa, "DQA1*05"), carries(dqb, "DQB1*02"))
        coeliac["dq2_2"] = both(carries(dqa, "DQA1*02:01"), carries(dqb, "DQB1*02:02"))
        coeliac["dq8"] = both(carries(dqa, "DQA1*03"), carries(dqb, "DQB1*03:02"))
        if coeliac["dq2_5"] == "present":
            drb1, why = usable("DRB1")
            dr3 = carries(drb1, "DRB1*03") if drb1 else None
            coeliac["dq2_5_phase"] = ("likely cis (DRB1*03 present)" if dr3 == "yes"
                                      else "possible trans (no DRB1*03)" if dr3 == "no"
                                      else f"phase unknown ({why or 'DRB1*03 uncertain'})")
    alleles = []
    for r in (hla_read_table() if table is None else table):
        gene = r["allele"].split("*")[0]
        l, why = usable(gene)
        state = carries(l, r["allele"]) if l else None
        alleles.append({"allele": r["allele"], "kind": r["kind"], "association": r["association"],
                        "pmids": r["pmids"], "citation": r["citation"],
                        "state": {"yes": "present", "no": "absent", "mixed": "uncertain", None: "not_determined"}[state],
                        "reason": why})
    return {"coeliac": coeliac, "alleles": alleles}


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
                     # "12 dropped", "12 kept", "12 dropped, strand check failed", "0"; empty in older tables
                     "ambiguous": r.get("Ambiguous", ""),
                     # the sex a sex-specific score applies to, and the check against the declared sex
                     "applies_to": r.get("Applies_To", ""), "sex_check": r.get("Sex_Check", ""),
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

# Input stays last: readers take it as the last column.
PRS_COLUMNS = ["Condition", "PGS_ID", "Score_SUM", "Variants_Matched", "Variants_Total",
               "Matched_Pct", "Percentile", "Ancestry_Group", "ChrX", "Ambiguous", "Applies_To",
               "Sex_Check", "Input"]
# prs-format writes it beside the scores: per score, its chrX rows and
# whether they were kept (scored) or left out (no sex given).
CHRX_FILE = "chrx_rows.tsv"
# prs-format writes it beside the scores too: per score, the sex it applies
# to (applies_to of assets/pgs_scores.tsv) and its strand-ambiguous rows (an
# A/T or C/G pair) with what happened to them.
SCORE_NOTES_FILE = "score_notes.tsv"
SEXES = ("female", "male")
# Sex_Check of a sex-specific score: the declared sex is the one it applies
# to, another one (no percentile then), or none was declared.
SEX_APPLIES = "applies"
SEX_NOT_APPLICABLE = "not applicable to the declared sex"
SEX_NOT_DECLARED = "no sex declared"
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
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


def read_applies_to(path):
    """ID -> the sex a score applies to, from the optional applies_to column of
    assets/pgs_scores.tsv: 'female', 'male', or '' for a score of both sexes.
    Raises Unreadable for any other value."""
    out = {}
    if not path:
        return out
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader((line for line in f if not line.startswith("#")), delimiter="\t"):
            pid = (r.get("pgs_id") or "").strip()
            sex = (r.get("applies_to") or "").strip().lower()
            if not pid:
                continue
            if sex and sex not in SEXES:
                raise Unreadable(f"{os.path.basename(path)}: applies_to of {pid} is '{sex}', expected female, male or empty")
            out[pid] = sex
    return out


class Fasta:
    """Bases of an uncompressed FASTA with its .fai, read by seeking."""

    def __init__(self, path):
        self.idx = {}
        try:
            with open(path + ".fai") as f:
                for line in f:
                    p = line.rstrip("\n").split("\t")
                    if len(p) >= 5:
                        self.idx[p[0]] = tuple(int(x) for x in p[1:5])
            self.f = open(path, "rb")
        except (OSError, ValueError) as e:
            raise Unreadable(f"cannot read the reference {path} and its .fai: {e}")
        if self.f.read(2) == b"\x1f\x8b":
            raise Unreadable(f"{path} is compressed; the strand check reads an uncompressed FASTA")

    def base(self, chrom, pos):
        """The base at 1-based POS of CHROM (with or without chr), upper case;
        '' when the contig is not in the reference or POS is outside it."""
        c = self.idx.get(chrom) or self.idx.get("chr" + chrom)
        if not c or not 1 <= pos <= c[0]:
            return ""
        length, offset, bases, width = c
        i = pos - 1
        self.f.seek(offset + (i // bases) * width + i % bases)
        return self.f.read(1).decode("ascii", "replace").upper()


def ambiguous_pair(ea, oa):
    """True for a strand-ambiguous SNV: its two alleles are complements (A/T, C/G)."""
    return len(ea) == 1 and len(oa) == 1 and COMPLEMENT.get(ea) == oa


def strand_check(rows, ref):
    """Whether a score's alleles are given on the forward strand, so its
    strand-ambiguous rows can be kept: (passed, rows that agree with the
    reference, rows that agree only after a strand flip). ROWS are (chrom,
    pos, effect, other) of the rows written. A row counts when both alleles
    are single bases and not an ambiguous pair: it agrees when the reference
    base is one of its alleles, and is flipped when the base is the
    complement of one of them instead. A score passes with at least one
    agreeing row and no flipped one; a row whose contig is not in the
    reference is not counted."""
    agree = flipped = 0
    for chrom, pos, ea, oa in rows:
        if len(ea) != 1 or len(oa) != 1 or ambiguous_pair(ea, oa) or not pos.isdigit():
            continue
        b = ref.base(chrom, int(pos))
        if b not in COMPLEMENT:
            continue
        if b in (ea, oa):
            agree += 1
        elif b in (COMPLEMENT.get(ea), COMPLEMENT.get(oa)):
            flipped += 1
    return agree > 0 and flipped == 0, agree, flipped


def single_allele(a):
    a = (a or "").strip().upper()
    return a if a and "/" not in a and a != "." else ""


def format_pgs(path, out_dir, label=None, keep_x=False, ref=None):
    """Write the GRCh38-harmonised scoring file PATH as a custom GRCh38 file
    pgsc_calc reads (chr_name and chr_position from hm_chr and hm_pos), and
    return (pgs_id, rows written, rows dropped off the autosomes,
    [(chrom, pos, allele), ...], chrX rows, ambiguous) with the effect and
    other alleles of every row written. Chromosomes 1 to 22 are kept, and chrX
    with keep_x (the sample's sex is known): pgsc_calc converts the sample
    with plink2, which stops on a chrX record when no sex is given. chrY and
    MT rows are always dropped. AMBIGUOUS is (strand-ambiguous rows, what
    happens to them, agreeing rows, flipped rows): without REF they are
    written and pgsc_calc drops them ("dropped"); with REF (keep_ambiguous)
    they are written and kept when the score passes strand_check ("kept"),
    and left out of the file when it does not ("dropped, strand check
    failed"). Raises Unreadable for a file that is not harmonised to GRCh38,
    lacks a column, is not additive or has no row."""
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
    counts = {"off": 0, "x": 0}

    def placed():
        """(chrom, pos, effect, other, weight) of each row pgsc_calc gets; one
        pass over ROWS each call, so no second copy of a large score is held."""
        counts["off"] = counts["x"] = 0
        for r in rows:
            chrom, pos, ea, ew = get(r, "hm_chr"), get(r, "hm_pos"), get(r, "effect_allele").upper(), get(r, "effect_weight")
            if not (chrom and pos and ea and ew):
                continue   # a row the catalog could not place on GRCh38
            chrom = chrom[3:] if chrom.startswith("chr") else chrom
            if chrom == "X":
                counts["x"] += 1
            if chrom not in AUTOSOMES and not (keep_x and chrom == "X"):
                counts["off"] += 1   # chrY, MT, and chrX unless the sex is known: plink2 refuses chrX without it
                continue
            yield chrom, pos, ea, single_allele(get(r, "other_allele")) or single_allele(get(r, "hm_inferOtherAllele")), ew

    agree = flipped = 0
    amb, drop_amb = "dropped", False
    if ref is not None:
        passed, agree, flipped = strand_check(((c, p, ea, oa) for c, p, ea, oa, _ in placed()), ref)
        if passed:
            amb = "kept"
        else:
            # pgsc_calc's keep_ambiguous applies to every score of the run,
            # so a score that fails loses its ambiguous rows here.
            amb, drop_amb = "dropped, strand check failed", True
    n = n_amb = 0
    with gzip.open(out + ".tmp", "wt", encoding="utf-8") as w:
        w.write(f"#pgs_id={pid}\n#pgs_name={pid}\n#trait_reported={trait}\n#genome_build=GRCh38\n")
        w.write("chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n")
        for chrom, pos, ea, oa, ew in placed():
            if ambiguous_pair(ea, oa):
                n_amb += 1
                if drop_amb:
                    continue
            w.write(f"{chrom}\t{pos}\t{ea}\t{oa}\t{ew}\n")
            n += 1
            alleles.append((f"chr{chrom}", pos, ea))
            if oa:
                alleles.append((f"chr{chrom}", pos, oa))
    if n == 0:
        os.remove(out + ".tmp")
        raise Unreadable(f"{name} has no autosomal row with hm_chr, hm_pos, effect_allele and effect_weight")
    os.replace(out + ".tmp", out)
    return pid, n, counts["off"], alleles, counts["x"], (n_amb, amb, agree, flipped)


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
    ap.add_argument("--keep-ambiguous", action="store_true",
                    help="pgsc_calc runs with --keep_ambiguous: keep each score's strand-ambiguous rows (A/T, C/G) "
                         "when its other rows agree with --reference on the forward strand, and leave them out of "
                         "a score that fails that check")
    ap.add_argument("--reference", help="the uncompressed GRCh38 FASTA with its .fai, for --keep-ambiguous")
    a = ap.parse_args(argv)
    if a.keep_ambiguous and not a.reference:
        print("ERROR: --keep-ambiguous needs --reference, the FASTA the strand check reads", file=sys.stderr)
        return 1
    try:
        applies = read_applies_to(a.labels)
        ref = Fasta(a.reference) if a.keep_ambiguous else None
    except (Unreadable, OSError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
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
    allele_set, chrx, notes = set(), [], []
    try:
        for p in files:
            pid, n, off, al, nx, (n_amb, amb, agree, flipped) = format_pgs(
                p, a.out, labels.get(pgs_id_of(p)), a.keep_x, ref)
            allele_set.update(al)
            chrx.append((pid, nx))
            notes.append((pid, applies.get(pid, ""), n_amb, amb))
            kept = f", {nx} on chrX scored" if a.keep_x and nx else ""
            dropped = f", {off} off the autosomes dropped" if off else ""
            ambig = f", {n_amb} strand-ambiguous {amb}" if n_amb else ""
            if ref is not None:
                ambig += f" (strand check: {agree} rows agree with the reference, {flipped} flipped)"
            sex = f", {applies[pid]} only" if applies.get(pid) else ""
            print(f"  {pid}: {n} GRCh38 rows{kept}{dropped}{ambig}{sex} ({labels.get(pid) or 'label from the file'})")
    except (Unreadable, OSError, EOFError, zlib.error) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    with open(os.path.join(a.out, CHRX_FILE), "w") as w:
        w.write("pgs_id\tchrx_rows\tchrx\n")
        for pid, nx in chrx:
            w.write(f"{pid}\t{nx}\t{'scored' if a.keep_x else 'left out'}\n")
    with open(os.path.join(a.out, SCORE_NOTES_FILE), "w") as w:
        w.write("pgs_id\tapplies_to\tambiguous_rows\tambiguous\n")
        for pid, sex, n_amb, amb in notes:
            w.write(f"{pid}\t{sex}\t{n_amb}\t{amb}\n")

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


def chrx_matched(results, sampleset):
    """PGS id -> chrX rows pgsc_calc matched, from its match log; None when
    there is no log. In a run with an ancestry panel pgsc_calc keeps only
    the variants it found in the panel too (the rest are 'excluded'), so a
    chrX row the panel lacks is not scored even with the sex."""
    logs = sorted(glob.glob(os.path.join(results, sampleset, "match", "*_log.csv.gz")))
    if not logs:
        return None
    out = {}
    for p in logs:
        with open_text(p) as f:
            head = next(csv.reader([f.readline()]), [])
            c = {k: i for i, k in enumerate(head)}
            if not {"accession", "chr_name", "match_status"} <= set(c):
                return None
            for line in f:
                # Most lines are not on chrX; split only those that may be.
                if ",X," not in line and ",23," not in line:
                    continue
                r = next(csv.reader([line]))
                if len(r) > max(c.values()) and r[c["chr_name"]] in ("X", "23") and r[c["match_status"]] == "matched":
                    pid = pgs_id_of(r[c["accession"]])
                    out[pid] = out.get(pid, 0) + 1
    return out


def prs_table(results, sampleset, scores_dir, input_kind, zero_matches=False, below_log=None, sex="", panel=False):
    """Rows of the summary table (PRS_COLUMNS) and the ancestry rows (key,
    value), from pgsc_calc's output folder RESULTS for SAMPLESET. Totals and
    labels come from the files prs-format wrote in SCORES_DIR. With
    zero_matches (pgsc_calc stopped because no score variant is in the
    sample's genotypes) every score is reported unmatched. With below_log
    (pgsc_calc stopped because every score matched under its minimum
    overlap, so it published no match summary) no score has a sum, the
    matched count is unknown (NA), and the match rate is read from that log
    when pgscatalog-match printed it. SEX is the declared sex ('' for none):
    a score that applies to the other sex keeps its sum and loses its
    percentile. PANEL: pgsc_calc ran with an ancestry panel, so a chrX row
    counts as scored only when its match log matched it."""
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
    chrx, chrx_scored = {}, set()
    if os.path.isfile(os.path.join(scores_dir, CHRX_FILE)):
        for r in read_tsv(os.path.join(scores_dir, CHRX_FILE)):
            n = r.get("chrx_rows", "")
            chrx[r.get("pgs_id", "")] = f"{n} {r.get('chrx', '')}" if n not in ("", "0") else n
            if n not in ("", "0") and r.get("chrx") == "scored":
                chrx_scored.add(r.get("pgs_id", ""))
    # "" for a score of both sexes; empty for scores formatted before the file existed
    applies, ambig = {}, {}
    if os.path.isfile(os.path.join(scores_dir, SCORE_NOTES_FILE)):
        for r in read_tsv(os.path.join(scores_dir, SCORE_NOTES_FILE)):
            pid, n = r.get("pgs_id", ""), r.get("ambiguous_rows", "")
            applies[pid] = r.get("applies_to", "")
            ambig[pid] = f"{n} {r.get('ambiguous', '')}" if n not in ("", "0") else n
    sex = (sex or "").lower()
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
        if panel and chrx_scored:
            # chrX rows given to pgsc_calc that its matches did not keep:
            # "3 left out (not matched in the panel run)" for all of them,
            # "1 scored, 2 left out (not matched in the panel run)" for some
            xm = chrx_matched(results, sampleset)
            if xm is not None:
                for pid in chrx_scored:
                    n_x, got = int(chrx[pid].split()[0]), xm.get(pid, 0)
                    if not got:
                        chrx[pid] = f"{n_x} left out (not matched in the panel run)"
                    elif got < n_x:
                        chrx[pid] = f"{got} scored, {n_x - got} left out (not matched in the panel run)"
    rows = []
    for pid in sorted(totals):
        m, t = matched.get(pid, 0), totals[pid]
        a = applies.get(pid, "")
        check = "" if not a else (SEX_NOT_DECLARED if not sex else (SEX_APPLIES if sex == a else SEX_NOT_APPLICABLE))
        extra = [chrx.get(pid, ""), ambig.get(pid, ""), a, check, input_kind]
        if below_log:
            rows.append([labels.get(pid, pid), pid, "NA", "NA", str(t), rate.get(pid, "NA"), "NA", "NA"] + extra)
            continue
        # A score of the other sex keeps its raw sum and gets no percentile:
        # a percentile for an organ the person does not have reads as a risk.
        p = "NA" if check == SEX_NOT_APPLICABLE else pct.get(pid, "NA")
        rows.append([labels.get(pid, pid), pid, sums.get(pid, "NA"), str(m), str(t),
                     f"{100 * m / t:.1f}" if t else "0.0", p, group.get(pid, "NA")] + extra)
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
    ap.add_argument("--sex", default="", choices=["", "male", "female"],
                    help="the declared sex: a score that applies to the other one gets no percentile")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ancestry-out", help="written when pgsc_calc ran with the panel: PCs and population")
    a = ap.parse_args(argv)
    try:
        rows, ancestry = prs_table(a.results, a.sampleset, a.scores, a.input_kind, a.zero_matches, a.below_threshold,
                                   a.sex, bool(a.panel))
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
        if r[11] == SEX_NOT_APPLICABLE:
            extra = f", {r[10]}-specific: {SEX_NOT_APPLICABLE}, no percentile"
        elif r[11] == SEX_NOT_DECLARED:
            extra += f", {r[10]}-specific score: {SEX_NOT_DECLARED}"
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


def sec_ancestry_fine(d, s):
    """Step 26b: <s>_ancestry_fine.tsv (level, code, group, fraction, se, low,
    high; levels continent, region and painting, plus left_out for a group the
    model leaves out and not_painted for a chromosome RFMix failed on) and
    the painting picture."""
    p = first_existing(d, [f"ancestry_fine/{s}_ancestry_fine.tsv"])
    if not p:
        return None, {}
    rows = read_tsv(p)
    need = {"level", "code", "group", "fraction", "se", "low", "high"}
    if not rows or not need <= set(rows[0]):
        raise Unreadable(f"{os.path.basename(p)} has no rows or lacks a column of {sorted(need)}")
    out = {"continent": [], "region": [], "painting": [], "left_out": [], "not_painted": [],
           "min_references": ancestry_fine_report.MIN_REFERENCES}
    for r in rows:
        if r["level"] in out:
            out[r["level"]].append({"code": r["code"], "group": r["group"], "fraction": num(r["fraction"]),
                                    "low": num(r["low"]), "high": num(r["high"])})
    # The painting picture travels in the summary, so the HTML report shows it
    # from the summary alone (render_report.py --summary).
    png = first_existing(d, [f"ancestry_fine/{s}_painting.png"])
    out["painting_png"] = None
    out["legend"] = [{"code": c, "group": ancestry_fine_report.label("continent", c),
                      "colour": "#%02x%02x%02x" % rgb} for c, rgb in ancestry_fine_report.COLOURS.items()]
    if png and os.path.getsize(png) < 2_000_000:
        with open(png, "rb") as f:
            out["painting_png"] = base64.b64encode(f.read()).decode()
    return p, out


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


def gencode_gene(d, gene):
    """(chrom, start, end, release) of GENCODE's gene line for gene, 1-based,
    from the gene lines setup.sh installs beside the sample folders
    (reference/gencode.v<N>.basic.genes.gtf, the newest release there), or
    None when there are none or the gene is not in them."""
    files = glob.glob(os.path.join(d, os.pardir, "reference", "gencode.v*.basic.genes.gtf"))
    rel = {f: int(m.group(1)) for f in files for m in [re.search(r"gencode\.v(\d+)\.basic", os.path.basename(f))] if m}
    if not rel:
        return None
    path = max(rel, key=rel.get)
    want = f'gene_name "{gene}";'
    with open(path, errors="replace") as f:
        for line in f:
            c = line.split("\t", 9)
            if len(c) >= 9 and c[2] == "gene" and want in c[8] and re.fullmatch(r"chr([0-9]+|X|Y)", c[0]):
                return c[0], int(c[3]), int(c[4]), str(rel[path])
    return None


# CNVpytor's read-depth level (1.0 is two copies) below which a deletion
# leaves no copy, and below which it leaves one.
UGT2B17_NO_COPY, UGT2B17_ONE_COPY = 0.25, 0.75


def ugt2b17_copy_number(calls, span):
    """UGT2B17 copy number from the CNVpytor call that covers most of its
    GENCODE span. calls: (type, chrom, start, end, level, e-value)."""
    chrom, gs, ge = span[:3]
    best = None
    for t, c, cs, ce, level, ev in calls:
        ov = min(ce, ge) - max(cs, gs) + 1
        if c == chrom and ov > 0 and (best is None or ov > best[0]):
            best = (ov, t, c, cs, ce, level, ev)
    out = {"span": f"{chrom}:{gs}-{ge}", "gencode_release": span[3]}
    if best is None:
        return dict(out, call="two copies", note="no CNVpytor call over the gene")
    ov, t, c, cs, ce, level, ev = best
    if t == "duplication":
        call = "more than two copies"
    elif level < UGT2B17_NO_COPY:
        call = "no copy (deleted on both chromosomes)"
    elif level < UGT2B17_ONE_COPY:
        call = "one copy"
    else:
        call = "two copies"
    note = (f"CNVpytor {t} {c}:{cs}-{ce}, depth level {level:.2f}, e-value {ev:.2g}, "
            f"covers {ov / (ge - gs + 1):.0%} of the gene")
    if ev >= 0.01:
        note += "; not significant (e-value 0.01 or above)"
    return dict(out, call=call, note=note, level=level, e_value=ev)


# A CNVpytor call is kept in the filtered count when fewer than half of its
# reads have mapping quality 0 (q0, column 9) and fewer than half of its bases
# are N (pN, column 10), it touches no centromere or heterochromatin band
# (MASKED_BANDS), and it is not on chrY of a declared female (collect() takes
# those out, once the declared sex is known).
CNV_MAX_Q0 = CNV_MAX_PN = 0.5


def sec_cnvpytor(d, s):
    p = first_existing(d, [f"cnvpytor/{s}_cnvs.txt"])
    if not p:
        return None, {}
    bands = masked_bands()
    total = dels = dups = sig = kept = kept_y = 0
    dropped = {"q0": 0, "gap": 0, "band": 0}
    calls = []
    with open(p, errors="replace") as f:
        for line in f:
            c = line.split()
            if not c:
                continue
            total += 1
            dels += c[0] == "deletion"
            dups += c[0] == "duplication"
            m = re.fullmatch(r"(\S+):(\d+)-(\d+)", c[1]) if len(c) > 4 else None
            try:
                if m:
                    calls.append((c[0], m.group(1), int(m.group(2)), int(m.group(3)), float(c[3]), float(c[4])))
            except ValueError:
                pass
            try:
                if float(c[4]) >= 0.01:
                    continue
            except (IndexError, ValueError):
                continue
            sig += 1
            q0, pn = (num(c[8]), num(c[9])) if len(c) > 9 else (None, None)
            if q0 is not None and q0 >= CNV_MAX_Q0:
                dropped["q0"] += 1
            elif pn is not None and pn >= CNV_MAX_PN:
                dropped["gap"] += 1
            elif m and masked_overlap(bands, m.group(1), int(m.group(2)), int(m.group(3))):
                dropped["band"] += 1
            else:
                kept += 1
                kept_y += bool(m) and m.group(1) in ("chrY", "Y")
    span = gencode_gene(d, "UGT2B17")
    ugt = (ugt2b17_copy_number(calls, span) if span else
           {"call": "not assessed", "note": "no GENCODE gene lines beside the results "
            "(reference/gencode.v<N>.basic.genes.gtf, which setup.sh installs) to place the gene"})
    return p, {"total": total, "deletions": dels, "duplications": dups, "significant": sig,
               "significant_filtered": kept, "filtered_chry": kept_y, "dropped": dropped, "ugt2b17": ugt}


def cnvpytor_sex_filter(data, declared_sex):
    """Take the chrY calls out of the filtered count of a declared female;
    with no declared sex they stay, and the data says so."""
    if "significant_filtered" not in data:
        return
    data["declared_sex"] = declared_sex or None
    if declared_sex == "female":
        data["dropped"]["chry_female"] = data["filtered_chry"]
        data["significant_filtered"] -= data["filtered_chry"]


def sec_sv_consensus(d, s):
    p = first_existing(d, [f"sv_merged/{s}_sv_consensus.vcf.gz"])
    if not p:
        return None, {}
    t, _ = count_vcf(p)
    return p, {"consensus": t}


def sec_annotsv(d, s):
    """Step 05's class 4 and 5 SVs (annotsv/<s>_sv_class4_5_checked.tsv,
    bin/annotsv_postfilter.py), each with how many callers support it from
    step 22's consensus (sv_merged/<s>_sv_consensus.vcf.gz), joined here
    because step 22 runs after step 05 (and beside ANNOTSV in Nextflow)."""
    p = first_existing(d, [f"annotsv/{s}_sv_class4_5_checked.tsv"])
    if not p:
        return None, {}
    rows = read_tsv(p) if os.path.getsize(p) else []
    need = {"SV_chrom", "SV_start", "SV_end", "SV_type", "ACMG_class", "Caller", "Exon_edge", "Population_AF", "Common"}
    if rows and not need <= set(rows[0]):
        raise Unreadable(f"{os.path.basename(p)} lacks a column of {sorted(need)}")
    cons = first_existing(d, [f"sv_merged/{s}_sv_consensus.vcf.gz"])
    cons_note = None
    try:
        support = annotsv_postfilter.consensus_support(rows, cons)
    except (OSError, EOFError, zlib.error, UnicodeDecodeError) as e:
        support = annotsv_postfilter.consensus_support(rows, None)
        cons_note = f"{os.path.basename(cons)} could not be read: {e}"
        cons = None
    out = []
    for r, sp in zip(rows, support):
        try:
            af = float(r["Population_AF"]) if r["Population_AF"] else None
        except ValueError:
            af = None
        genes = [g for g in (r.get("Gene_name") or "").split(";") if g]
        out.append({"sv": f"{r['SV_chrom']}:{r['SV_start']}-{r['SV_end']}", "type": r["SV_type"],
                    "length": r.get("SV_length") or "", "acmg_class": r["ACMG_class"], "genes": genes,
                    "caller": r["Caller"], "support": sp["Support"],
                    "support_callers": sp["Support_callers"].split(",") if sp["Support_callers"] else [],
                    "exon_edge": r["Exon_edge"], "exon_edge_detail": r.get("Exon_edge_detail") or "",
                    "population_af": af, "population_source": r.get("Population_source") or "",
                    "common": r["Common"] == "yes"})
    return p, {"count": len(out), "rows": out,
               "consensus": os.path.relpath(cons, d) if cons else None, "consensus_note": cons_note,
               "single_caller": sum(x["support"] == 1 for x in out),
               "exon_edge": sum(x["exon_edge"] == "yes" for x in out),
               "exon_edge_checked": not any(x["exon_edge"] == "not checked" for x in out),
               "common": sum(x["common"] for x in out)}


STR_RANGES_TSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "str_intermediate_ranges.tsv")
# How well ExpansionHunter sized an allele, from its FORMAT SO, ADSP and
# REPCI fields, worst first; a locus takes the worst of its alleles.
STR_SIZING = ["no_call", "flanking_only", "in_repeat", "low_confidence", "sized", "unknown"]
STR_MIN_SPANNING = 2


def str_ranges(path=STR_RANGES_TSV):
    """bin/str_intermediate_ranges.tsv by locus (tests/test_str_ranges.py checks it)."""
    rows, head = {}, None
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if head is None:
                head = c
                continue
            r = dict(zip(head, c))
            rows[r["locus"]] = r
    return rows


def str_int(x):
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


def str_allele_sizing(so, adsp, ci, cn, path_min):
    """One allele's sizing: sized, low_confidence (fewer than
    STR_MIN_SPANNING spanning reads, or a confidence interval that reaches
    the catalog's pathogenic minimum from below it), in_repeat (in-repeat
    reads only: longer than a read), flanking_only, no_call, or unknown
    (no SO field)."""
    if str_int(cn) is None:
        return "no_call"
    if not so or so == ".":
        return "unknown"
    if so == "FLANKING":
        return "flanking_only"
    if so == "INREPEAT":
        return "in_repeat"
    hi = str_int(ci.split("-")[-1]) if ci else None
    if (str_int(adsp) is not None and str_int(adsp) < STR_MIN_SPANNING) or (
            path_min is not None and hi is not None and str_int(cn) < path_min <= hi):
        return "low_confidence"
    return "sized"


def str_field(fmt, val, key):
    return val[fmt.index(key)] if key in fmt and fmt.index(key) < len(val) else ""


def str_split(x, n):
    """A per-allele FORMAT value as n parts; one value serves every allele."""
    parts = x.split("/") if x else []
    return parts if len(parts) == n else (parts * n if len(parts) == 1 else [""] * n)


def rfc1_motif_row(d, s, after):
    """Step 9b's RFC1 motif count (bin/rfc1_motif.py), unless it is older than
    the Stranger file it was made from."""
    p = first_existing(d, [f"expansion_hunter/{s}_rfc1_motif.tsv"])
    if not p or (after and os.path.getmtime(p) < os.path.getmtime(after)):
        return None
    rows = read_tsv(p)
    if not rows:
        return None
    r = rows[0]
    return {k: r.get(k, ".") for k in ("verdict", "note", "dominant_motif", "dominant_class",
                                       "dominant_share", "motif_shares", "repeat_reads")}


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
    ranges = str_ranges()
    rfc1 = rfc1_motif_row(d, s, st) if st else None
    loci, tested, flagged, no_status = [], 0, [], 0
    for r in vcf_records(p):
        if len(r) < 10:
            continue
        tested += 1
        info = info_map(r[7])
        rep = info.get("REPID") or info.get("VARID") or f"{r[0]}:{r[1]}"
        fmt = r[8].split(":")
        val = r[9].split(":")
        repcn = str_field(fmt, val, "REPCN")
        counts = repcn.split("/") if repcn and repcn != "." else []
        n = max(len(counts), 1)
        path_min = str_int(info.get("STR_PATHOLOGIC_MIN"))
        sizing = [str_allele_sizing(so, ad, ci, cn, path_min) for so, ad, ci, cn in
                  zip(str_split(str_field(fmt, val, "SO"), n), str_split(str_field(fmt, val, "ADSP"), n),
                      str_split(str_field(fmt, val, "REPCI"), n), counts or ["."])]
        status = [x for x in info.get("STR_STATUS", "").split(",") if x and x != "."] if st else []
        # Stranger writes '.' where its catalog has no value.
        item = {"locus": rep, "repeat_count": repcn or ".", "status": ",".join(status),
                "normal_max": info.get("STR_NORMAL_MAX", "").strip("."),
                "pathologic_min": info.get("STR_PATHOLOGIC_MIN", "").strip("."),
                "support": str_field(fmt, val, "SO"), "adsp": str_field(fmt, val, "ADSP"),
                "ci": str_field(fmt, val, "REPCI"), "sizing": min(sizing, key=STR_SIZING.index),
                "allele_sizing": "/".join(sizing), "zone": "", "zone_label": "", "zone_citation": ""}
        rng = ranges.get(rep)
        sized = [str_int(c) for c, z in zip(counts, sizing) if z in ("sized", "unknown") and str_int(c) is not None]
        if rng and any(int(rng["intermediate_min"]) <= c <= int(rng["intermediate_max"]) for c in sized):
            item.update(zone=f"{rng['intermediate_min']}-{rng['intermediate_max']}", zone_label=rng["label"],
                        zone_citation=rng["citation"])
        if rep == "RFC1" and st and any(x != "normal" for x in status):
            item["rfc1_motif"] = rfc1 or {"verdict": "not_checked",
                                          "note": "no motif count beside Stranger's output (step 9b runs it "
                                                  "with the sample's BAM): confirm the flag with repeat-primed PCR"}
        loci.append(item)
        if st:
            # One value per record, the most severe of normal, pre_mutation and
            # full_mutation; a locus outside Stranger's catalog has none.
            if not status:
                no_status += 1
            elif any(x != "normal" for x in status):
                f = {"locus": item["locus"], "repeat_count": item["repeat_count"], "status": item["status"],
                     "sizing": item["sizing"]}
                if "rfc1_motif" in item:
                    f["rfc1_motif"] = item["rfc1_motif"]
                flagged.append(f)
    by = {l["locus"]: l for l in loci}
    # Every locus, the five key loci first, then by name.
    loci.sort(key=lambda l: (EH_LOCI.index(l["locus"]) if l["locus"] in EH_LOCI else len(EH_LOCI), l["locus"]))
    out = {"records": tested, "key_loci": [{"locus": k, "repeat_count": by[k]["repeat_count"] if k in by
                                            else "not in output"} for k in EH_LOCI],
           "loci": loci, "stranger": bool(st), "stranger_outdated": outdated,
           "unsized": [l["locus"] for l in loci if l["sizing"] not in ("sized", "unknown")],
           "intermediate": [l["locus"] for l in loci if l["zone"]]}
    if st:
        out.update(flagged=flagged, no_status=no_status)
    return p, out


# TelomereHunter's estimate rests on its intratelomeric reads. Below this many
# per 30x of depth, about 30 times fewer than a typical 30x short-read genome
# has, the value is too thin to compare (a library or platform problem, or
# reads too short for its repeat threshold). The floor scales with the depth
# TelomereHunter itself counted (total_reads x read length / GENOME_BASES).
TELOMERE_FLOOR_PER_30X = 2000
GENOME_BASES = 3.1e9


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
    out = {"tel_content": row[i] if i < len(row) else "."}

    def col(*names):
        for n in names:
            if n in head and head.index(n) < len(row):
                return row[head.index(n)]
        return None
    # TelomereHunter 1.1.0 names it intratel_reads.
    intratel = num(col("intratel_reads", "intratelomeric_reads"))
    total = num(col("total_reads"))
    lengths = [num(x) for x in re.split(r"[,;| ]+", col("read_lengths", "read_length") or "")]
    lengths = [x for x in lengths if x]
    if intratel is None:
        return p, dict(out, read_check="no intratelomeric read count in the summary")
    out["intratelomeric_reads"] = int(intratel)
    if not total or not lengths:
        return p, dict(out, read_check="depth unknown (no total_reads or read_lengths in the summary): not checked")
    depth = total * max(lengths) / GENOME_BASES
    floor = TELOMERE_FLOOR_PER_30X * depth / 30
    out.update(depth=round(depth, 1), read_floor=int(round(floor)), low_reads=intratel < floor)
    return p, out


# GRCh38's centromeres and heterochromatin: the acen, gvar and stalk rows of
# UCSC's cytoBand.hg38.txt (the file setup.sh installs as
# reference/cytoBand.hg38.txt, sha256 e514e48b... in scripts/lib/common.sh),
# neighbouring rows merged, as chrom:start-end with a 0-based start and an
# exclusive end. Reads there map badly, so a run of homozygosity or a copy
# number call over them is mostly an artefact. Regenerate with
#   gzip -cd cytoBand.txt.gz | awk '$1 ~ /^chr([0-9]+|X|Y)$/ && $5 ~ /^(acen|gvar|stalk)$/' \
#     | sort -k1,1V -k2,2n | awk '{ if ($1 == c && $2 == e) e = $3; else { if (c != "") printf "%s:%d-%d ", c, s, e;
#       c = $1; s = $2; e = $3 } } END { print c ":" s "-" e }'
# scripts/11-roh-analysis.sh and modules/local/roh/main.nf carry the same
# line (tests/test_collect_summary.py checks that the three agree).
MASKED_BANDS = ("chr1:121700000-143200000 chr2:91800000-96000000 chr3:87800000-98600000 chr4:48200000-51800000 "
                "chr5:46100000-51400000 chr6:58500000-62600000 chr7:58100000-62100000 chr8:43200000-47200000 "
                "chr9:42200000-61500000 chr10:38000000-41600000 chr11:51000000-55800000 chr12:33200000-37800000 "
                "chr13:0-18900000 chr14:0-18200000 chr15:0-20500000 chr16:35300000-47000000 "
                "chr17:22700000-27400000 chr18:15400000-21500000 chr19:19900000-31900000 chr20:25700000-30400000 "
                "chr21:0-13000000 chr22:0-17400000 chrX:58100000-63800000 chrY:10300000-10600000 "
                "chrY:26600000-57227415")


def masked_bands(text=MASKED_BANDS):
    """{chrom: [(start, end), ...]} of MASKED_BANDS, 0-based, end exclusive."""
    out = {}
    for x in text.split():
        c, _, se = x.partition(":")
        a, _, b = se.partition("-")
        out.setdefault(c, []).append((int(a), int(b)))
    return out


def masked_overlap(bands, chrom, start, end):
    """Bases of the 1-based, inclusive span start-end inside the masked bands."""
    return sum(max(0, min(end, b) - max(start - 1, a)) for a, b in bands.get(chrom, ()))


def sec_roh(d, s):
    """Step 11's bcftools roh output. Segments on chrX and chrY are left out
    (a man's single X is one long run), segments on any other contig that is
    not chr1 to chr22 (chrM, chrEBV, chrUn_*) are not counted, and so is the
    part of a segment that lies in a centromere or heterochromatin band
    (MASKED_BANDS), where few markers map. Total, Largest and the list of segments over 5 Mb are of the
    autosomal length that is left."""
    p = first_existing(d, [f"vcf/{s}_roh.txt", f"roh/{s}_roh.txt", "roh/*_roh.txt"])
    if not p:
        return None, {}
    bands = masked_bands()
    total = largest = masked = 0.0
    segs, sex_segs, big = 0, 0, []
    with open(p, errors="replace") as f:
        for line in f:
            c = line.split()
            if len(c) < 6 or c[0] != "RG":
                continue
            try:
                start, end, length = int(c[3]), int(c[4]), float(c[5])
            except ValueError:
                continue
            if c[2] in ("chrX", "chrY", "X", "Y"):
                sex_segs += 1
                continue
            if not (c[2].startswith("chr") and c[2][3:] in AUTOSOMES):
                continue
            inside = min(length, masked_overlap(bands, c[2], start, end))
            masked += inside
            length -= inside
            if length <= 0:
                continue
            segs += 1
            total += length
            largest = max(largest, length)
            if length > 5e6:
                row = {"region": f"{c[2]}:{c[3]}-{c[4]}", "mb": round(length / 1e6, 1)}
                if inside:
                    row["masked_mb"] = round(inside / 1e6, 1)
                big.append(row)
    return p, {"segments": segs, "total_mb": round(total / 1e6, 1), "largest_mb": round(largest / 1e6, 1),
               "autosomal_over_5mb": big, "masked_mb": round(masked / 1e6, 1), "sex_chromosome_segments": sex_segs}


# Pages that show more of a called branch, free and with no upload. The URL
# patterns were checked live on 2026-10-25 (docs/37-y-haplogroup.md,
# "Explore your branch further"): FamilyTreeDNA Discover answers its own
# HaplogroupSearch for R-M269 and H5a7, and YFull's tree pages answer 200 for
# a branch it has and 404 for one it does not.
EXPLORE_DOCS = "https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/docs/"


def haplogroup_explore(kind, hg):
    """[{name, url}] for a called branch: kind 'y' (Yleaf, e.g. R-M269 or
    R-M269*(xR-L23)) or 'mt' (haplogrep3, e.g. H5a7). Empty when there is no call."""
    from urllib.parse import quote
    hg = (hg or "").strip()
    if kind == "y":
        hg = hg.split("*", 1)[0]   # Yleaf's paragroup mark: the branch is the name before '*'
    if not hg or hg in (".", "NA", "insufficient markers"):
        return []
    q = quote(hg, safe="")
    if kind == "y":
        return [{"name": "FamilyTreeDNA Discover (free, no upload)", "url": f"https://discover.familytreedna.com/y-dna/{q}"},
                {"name": "YFull YTree", "url": f"https://www.yfull.com/tree/{q}/"},
                {"name": "More routes, and the files to share", "url": EXPLORE_DOCS + "37-y-haplogroup.md#explore-your-branch-further"}]
    return [{"name": "FamilyTreeDNA Discover (free, no upload)", "url": f"https://discover.familytreedna.com/mtdna/{q}"},
            {"name": "YFull MTree", "url": f"https://www.yfull.com/mtree/{q}/"},
            {"name": "More routes, and the files to share", "url": EXPLORE_DOCS + "12-mito-haplogroup.md#explore-your-branch-further"}]


def sec_haplogroup(d, s):
    p = first_existing(d, [f"mito/{s}_haplogroup.txt"])
    if not p:
        return None, {}
    rows = read_tsv(p)
    r = rows[0] if rows else {}
    hg = (r.get("Haplogroup") or (list(r.values())[1] if len(r) > 1 else "") or "").strip('"')
    out = {"haplogroup": hg or ".", "explore": haplogroup_explore("mt", hg)}
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
               "valid_markers": r.get("Valid_markers") or ".", "qc_score": r.get("QC-score") or ".",
               "explore": haplogroup_explore("y", hg)}


def sec_y_str(d, s):
    """Step 45's table (y_str/<id>_y_str.tsv, bin/ystr_report.py): one row per
    FamilyTreeDNA marker with its count and confidence, and the counts per
    confidence."""
    p = first_existing(d, [f"y_str/{s}_y_str.tsv"])
    if not p:
        return None, {}
    rows = read_tsv(p)
    need = ("order", "panel", "marker", "count", "confidence")
    if not rows or any(k not in rows[0] for k in need):
        raise Unreadable(f"{os.path.basename(p)} has no rows or lacks a column of {', '.join(need)}")
    keep = [{k: r.get(k) or "" for k in need} for r in rows]
    counts = {c: sum(1 for r in keep if r["confidence"] == c) for c in ("high", "low", "no call", "not callable")}
    if sum(counts.values()) != len(keep):
        raise Unreadable(f"{os.path.basename(p)} has a confidence other than high, low, no call or not callable")
    return p, {"markers": keep, **{k.replace(" ", "_"): v for k, v in counts.items()}}


def sec_consumer_export(d, s):
    """Step 40's summary table (key, value): how many chip sites got a genotype."""
    p = first_existing(d, [f"consumer_export/{s}_consumer_export.tsv"])
    if not p:
        return None, {}
    data = {r["key"]: r.get("value", "") for r in read_tsv(p) if r.get("key")}
    for k in ("sites", "called", "no_call"):
        if not (data.get(k) or "").isdigit():
            raise Unreadable(f"{os.path.basename(p)} has no whole-number '{k}' row")
    reasons = {}
    for k, v in data.items():
        if k.startswith("no_call: "):
            if not (v or "").isdigit():
                raise Unreadable(f"{os.path.basename(p)} has a '{k}' row without a whole-number count")
            reasons[k[len("no_call: "):]] = int(v)
    return p, {"sites": int(data["sites"]), "called": int(data["called"]), "no_call": int(data["no_call"]),
               "from_reference_blocks": int(data.get("from_reference_blocks") or 0),
               "template_sites": data.get("template_sites") or ".",
               "file_23andme": data.get("file_23andme") or ".",
               "file_ancestrydna": data.get("file_ancestrydna") or ".",
               "no_call_reasons": reasons}


def sec_bloodtype(d, s):
    """Step 41's table: one row per blood group system, predicted from the genome."""
    p = first_existing(d, [f"bloodtype/{s}_bloodtype.tsv"])
    if not p:
        return None, {}
    rows = read_tsv(p)
    if not rows or not {"system", "phenotype", "confidence"} <= set(rows[0]):
        raise Unreadable(f"{os.path.basename(p)} has no system, phenotype and confidence columns, or no row")
    systems = [{"system": r["system"], "phenotype": r["phenotype"] or ".", "genotype": r.get("genotype") or ".",
                "confidence": r["confidence"] or "."} for r in rows]
    by = {x["system"]: x["phenotype"] for x in systems}
    return p, {"systems": systems, "abo": by.get("ABO", "."), "rhd": by.get("RhD", ".")}


# The two LPA risk variants of assets/traits.tsv (rsid -> its alt, the allele
# that goes with a high Lp(a)); tests/test_lpa_kiv2.py checks they agree.
LPA_RISK_ALLELES = {"rs10455872": "G", "rs3798220": "C"}


def sec_lpa(d, s):
    """Step 46's KIV-2 count (lpa/<id>_lpa_kiv2.tsv, bin/lpa_kiv2.py), with
    step 42's two LPA risk variants when the traits table is there, and
    whether either points to a genetically high Lp(a)."""
    p = first_existing(d, [f"lpa/{s}_lpa_kiv2.tsv"])
    if not p:
        return None, {}
    kv = {r.get("metric"): r.get("value") for r in read_tsv(p)}
    if kv.get("status") not in ("estimated", "too little depth", "not assessed"):
        raise Unreadable(f"{os.path.basename(p)} has no status row (estimated, too little depth or not assessed)")
    data = {k: kv.get(k) or "." for k in ("status", "kiv2_copies", "ci95_low", "ci95_high", "class",
                                          "kiv2_depth", "flank_depth", "mapq0_share", "message", "advice")}
    variants = []
    t = first_existing(d, [f"traits/{s}_traits.tsv"])
    for r in (read_tsv(t) if t else []):
        alt = LPA_RISK_ALLELES.get(r.get("rsid"))
        if not alt:
            continue
        called = r.get("status") == "called"
        n = re.split(r"[/|]", r.get("genotype") or "").count(alt) if called else None
        variants.append({"rsid": r["rsid"], "genotype": r.get("genotype") or ".", "status": r.get("status") or ".",
                         "risk_alleles": n})
    risk = sum(v["risk_alleles"] for v in variants if v["risk_alleles"] is not None)
    reasons = ([f"{risk} risk allele{'s' if risk != 1 else ''} at rs10455872 and rs3798220"] if risk else []) \
        + (["a low KIV-2 repeat count"] if data["class"] == "low" else [])
    read = data["status"] == "estimated" or any(v["risk_alleles"] is not None for v in variants)
    data.update(variants=variants, risk_alleles=risk, flag="high" if reasons else ".",
                flag_text=("Genetic signs of a high Lp(a): " + " and ".join(reasons) + ". Measure Lp(a) in blood."
                           if reasons else
                           "No genetic sign of a high Lp(a) from these markers. They explain only part of it: "
                           "the KIV-2 count explained 21% to 27% of the variation in a Danish population study."
                           if read else
                           "Not determined: neither the KIV-2 count nor the two risk variants could be read. "
                           "Measure Lp(a) in blood."))
    return p, data


def sec_viral(d, s):
    """Step 47's table (viral/<id>_viral.tsv, bin/viral_report.py): one row per
    viral genome; the HHV-6 line and the genomes with reads (phiX included)."""
    p = first_existing(d, [f"viral/{s}_viral.tsv"])
    if not p:
        return None, {}
    rows = read_tsv(p)
    need = ("virus", "accession", "reads", "breadth", "relative_depth", "copies_per_cell", "call")
    if not rows or any(k not in rows[0] for k in need):
        raise Unreadable(f"{os.path.basename(p)} has no rows or lacks a column of {', '.join(need)}")
    keep = [{k: r[k] or "." for k in need + ("human_depth", "note")} for r in rows]
    ici = [r["virus"] for r in keep if r["call"] == "inherited integrated HHV-6 likely"]
    # No read on chr1 to chr22 leaves no human depth, so integration was never assessed.
    hhv6 = (f"inherited integrated {ici[0]} likely" if ici
            else "not assessed (no human depth)" if keep[0]["human_depth"] == "."
            else "no integrated HHV-6")
    return p, {"hhv6": hhv6,
               "human_depth": keep[0]["human_depth"],
               "found": [r for r in keep if r["reads"].isdigit() and int(r["reads"]) > 0], "genomes": len(keep)}


def sec_traits(d, s):
    """Step 42's traits table (traits/<id>_traits.tsv, bin/traits_report.py):
    every row with its genotype, status and reading, the counts per status,
    and the APOE line on its own (the reports fold it away)."""
    p = first_existing(d, [f"traits/{s}_traits.tsv"])
    if not p:
        return None, {}
    with open(p, encoding="utf-8", errors="replace") as f:
        rows = list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))
    need = ("rsid", "gene", "trait", "genotype", "status", "reading", "ancestry", "note", "pmids")
    if not rows or any(k not in rows[0] for k in need):
        raise Unreadable(f"{os.path.basename(p)} has no rows or lacks a column of {', '.join(need)}")
    # A row cut short (DictReader fills the missing fields with None) or with
    # extra fields (kept under the key None) is a damaged table, never a
    # called trait without its reading.
    bad = next((i for i, r in enumerate(rows, 2) if None in r or any(r.get(k) is None for k in need)), None)
    if bad:
        raise Unreadable(f"{os.path.basename(p)} line {bad} does not have one field per column")
    keep = [{k: r[k] for k in need} for r in rows]
    apoe = [r for r in keep if r["gene"] == "APOE"]
    counts = {st: sum(1 for r in keep if r["status"] == st and r["gene"] != "APOE")
              for st in ("called", "not_covered", "other_allele")}
    return p, {"rows": [r for r in keep if r["gene"] != "APOE"], "apoe": apoe, **counts}


def sec_carrier(d, s):
    """The sample's own carrier list (carrier_list) from step 06's ClinVar
    hits; its own section, so the ClinVar section reads the same whether or
    not the assets folder is there."""
    p = first_existing(d, [f"clinvar/{s}_clinvar_hits.vcf"])
    if not p:
        return None, {}
    return p, carrier_list_or_note(p)


def sec_couple(d, s):
    """Step 39's couple table (couple/<id>_couple_carrier.tsv, bin/couple_report.py):
    the partner, the shared and X-linked genes with the chance for a child, and
    how many genes only one partner carries a change in. A table without rows
    (no finding in either partner) names no partner."""
    p = first_existing(d, [f"couple/{s}_couple_carrier.tsv"])
    if not p:
        return None, {}
    with open(p, encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        rows = list(reader)
        cols = reader.fieldnames or []
    need = ("gene", "condition", "result", "offspring_risk", "partner_a", "findings_a", "partner_b",
            "findings_b", "explanation")
    if any(k not in cols for k in need):
        raise Unreadable(f"{os.path.basename(p)} lacks a column of {', '.join(need)}")
    bad = next((i for i, r in enumerate(rows, 2) if None in r or any(r.get(k) is None for k in need)
                or s not in (r["partner_a"], r["partner_b"]) or r["result"] not in ("shared", "x_linked", "one_partner")),
               None)
    if bad:
        raise Unreadable(f"{os.path.basename(p)} line {bad} is not a row of {s}'s couple table")

    def side(r):
        own = "a" if r["partner_a"] == s else "b"
        other = "b" if own == "a" else "a"
        return {"gene": r["gene"], "condition": r["condition"], "risk": r["offspring_risk"],
                "own": r["findings_" + own], "partner_findings": r["findings_" + other],
                "explanation": r["explanation"]}
    partner = next((r["partner_b"] if r["partner_a"] == s else r["partner_a"] for r in rows), None)
    return p, {"partner": partner,
               "shared": [side(r) for r in rows if r["result"] == "shared"],
               "x_linked": [side(r) for r in rows if r["result"] == "x_linked"],
               "one_partner": sum(1 for r in rows if r["result"] == "one_partner")}


# The D-loop poly-C tracts (rCRS m.303-315 and m.16184-16193). Length
# variants there, an extra or missing C, are found in most people, often as a
# mix of lengths that reads as heteroplasmy; they are counted apart and left
# out of the heteroplasmic count.
POLYC_TRACTS = ((303, 315), (16184, 16193))


def mito_span(pos, ref, alt):
    """The rCRS positions a chrM allele changes: an SNV or MNV its own bases,
    an indel the bases after its anchor base (Mutect2 writes the anchor)."""
    if len(ref) == len(alt):
        return pos, pos + len(ref) - 1
    return pos + 1, pos + max(len(ref) - 1, 1)


def in_polyc(pos, ref, alt):
    lo, hi = mito_span(pos, ref, alt)
    return any(lo <= b and hi >= a for a, b in POLYC_TRACTS)


def min_allele(pos, ref, alt):
    """(pos, ref, alt) without the bases ref and alt share at the end and then
    at the start, one base always kept: Mutect2's A>AC of a two-ALT record
    AC>A,ACC becomes the same key as ClinVar's and gnomAD's."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


MITO_ANNOT_COLS = ["chrom", "pos", "ref", "alt", "filter", "af", "clinvar_id", "clinvar_significance",
                   "clinvar_review", "clinvar_condition", "gnomad_af_hom", "gnomad_af_het", "gnomad_max_hl"]
# ClinVar classes: P/LP are counted; a record whose every class is benign or
# carries no classification is not listed unless the call is heteroplasmic.
CLINVAR_PLP = ("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")
CLINVAR_BENIGN = ("Benign", "Likely_benign", "Benign/Likely_benign", "not_provided")


def clinvar_terms(sig):
    return [t for t in re.split(r"[|,]", sig or "") if t and t != "."]


def vcf_file_date(path):
    with open_text(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            if line.startswith("##fileDate="):
                return line.strip().split("=", 1)[1]
    return None


def mito_annotate(vcf, clinvar=None, gnomad=None):
    """Step 20's annotation table: each ALT of each record of the chrM VCF,
    with its ClinVar record (CLNSIG, CLNREVSTAT, CLNDN of the same allele) and
    gnomAD v3.1's chrM sites (AF_hom, AF_het, max_hl of the same allele).
    Returns (header lines, rows); a file not given or not there is named in
    the header and its columns are '.'."""
    head, cv, gn = [], {}, {}
    if clinvar and os.path.isfile(clinvar) and os.path.isfile(clinvar + ".tbi"):
        head.append(f"# clinvar\t{os.path.basename(clinvar)}\t{vcf_file_date(clinvar) or 'fileDate unknown'}")
    else:
        head.append(f"# clinvar\tnot annotated\t{'no ' + clinvar + ' with its .tbi' if clinvar else 'no ClinVar file given'}")
        clinvar = None
    if gnomad and os.path.isfile(gnomad):
        head.append(f"# gnomad\t{os.path.basename(gnomad)}\t-")
        for r in vcf_records(gnomad):
            if len(r) < 8:
                continue
            i = info_map(r[7])
            for alt in r[4].split(","):
                gn[min_allele(int(r[1]), r[3], alt)] = (i.get("AF_hom", "."), i.get("AF_het", "."), i.get("max_hl", "."))
    else:
        head.append(f"# gnomad\tnot annotated\t{'no ' + gnomad if gnomad else 'no gnomAD chrM file given'}")
    recs = list(vcf_records(vcf))
    if clinvar:
        index = tabix_linear_index(clinvar)
        chrom = next((c for c in ("chrM", "MT") if c in index), None)
        wanted = {int(r[1]) for r in recs if len(r) >= 5} | {min_allele(int(r[1]), r[3], a)[0] for r in recs if len(r) >= 5 for a in r[4].split(",")}
        for pos in sorted(wanted) if chrom else ():
            # ClinVar's own record can start a base earlier or later than
            # Mutect2's for the same indel; the alleles decide.
            for q in (x for x in (pos - 1, pos, pos + 1) if x >= 1):
                for c in tabix_lines(clinvar, index, chrom, q):
                    i = info_map(c[7]) if len(c) > 7 else {}
                    for alt in c[4].split(","):
                        cv[min_allele(int(c[1]), c[3], alt)] = (c[2], i.get("CLNSIG", "."), i.get("CLNREVSTAT", "."),
                                                                i.get("CLNDN", "."))
    rows = []
    for r in recs:
        if len(r) < 8:
            continue
        fmt = r[8].split(":") if len(r) > 9 else []
        val = r[9].split(":") if len(r) > 9 else []
        afs = val[fmt.index("AF")].split(",") if "AF" in fmt and fmt.index("AF") < len(val) else []
        for k, alt in enumerate(r[4].split(",")):
            key = min_allele(int(r[1]), r[3], alt)
            rows.append([r[0], str(key[0]), key[1], key[2], r[6], afs[k] if k < len(afs) else ".",
                         *cv.get(key, (".", ".", ".", ".")), *gn.get(key, (".", ".", "."))])
    return head, rows


def mito_annotate_main(argv):
    ap = argparse.ArgumentParser(prog="collect_summary.py mito-annotate",
                                 description="Write step 20's chrM annotation table (ClinVar and gnomAD).")
    ap.add_argument("--vcf", required=True, help="the filtered chrM VCF (<sample>_chrM_filtered.vcf.gz)")
    ap.add_argument("--clinvar", default="", help="ClinVar's VCF with chr names, and its .tbi")
    ap.add_argument("--gnomad", default="", help="gnomAD v3.1's chrM sites VCF")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    try:
        head, rows = mito_annotate(a.vcf, a.clinvar or None, a.gnomad or None)
    except (Unreadable, OSError, ValueError, struct.error, zlib.error) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    tmp = a.out + ".tmp"
    with open(tmp, "w") as f:
        f.write("".join(h + "\n" for h in head))
        f.write("\t".join(MITO_ANNOT_COLS) + "\n")
        f.write("".join("\t".join(r) + "\n" for r in rows))
    os.replace(tmp, a.out)
    for h in head:
        c = h[2:].split("\t")
        print(f"  {c[0]}: {c[1]}" + (f" ({c[2]})" if c[2] != "-" else ""))
    n_cv = sum(r[7] != "." for r in rows)
    n_gn = sum(r[10] != "." for r in rows)
    print(f"  {len(rows)} chrM alleles: {n_cv} with a ClinVar record, {n_gn} in gnomAD")
    return 0


def mito_annotation(d, s, vcf):
    """Step 20's annotation table, as the report lists it, or a note on why
    there is none."""
    p = first_existing(d, [f"mito/{s}_chrM_annotated.tsv"])
    if not p:
        return {"note": "no annotation table (mito/<id>_chrM_annotated.tsv): step 20 ran before it annotated, "
                        "or the Nextflow MITO_VARIANTS_ANNOTATE task failed (see the Nextflow log)"}
    if os.path.getmtime(p) < os.path.getmtime(vcf):
        return {"note": f"{os.path.relpath(p, d)} is older than the chrM calls: rerun step 20"}
    sources, lines = {}, []
    with open(p, errors="replace") as f:
        for line in f:
            if line.startswith("# "):
                c = line[2:].rstrip("\n").split("\t")
                if len(c) >= 3:
                    sources[c[0]] = (f"not annotated: {c[2]}" if c[1] == "not annotated"
                                     else c[1] + (f" ({c[2]})" if c[2] != "-" else ""))
            elif line.strip():
                lines.append(line)
    rows = [r for r in csv.DictReader(lines, delimiter="\t") if r.get("filter") == "PASS"]
    return {"source": os.path.relpath(p, d), "clinvar": sources.get("clinvar"), "gnomad": sources.get("gnomad"),
            "rows": rows}


def sec_mito(d, s):
    p = first_existing(d, [f"mito/{s}_chrM_filtered.vcf.gz"])
    if not p:
        return None, {}
    passed = het = 0
    polyc = []
    for r in vcf_records(p):
        if len(r) < 10 or r[6] != "PASS":
            continue
        passed += 1
        fmt = r[8].split(":")
        val = r[9].split(":")
        if "AF" not in fmt or fmt.index("AF") >= len(val):
            continue
        try:
            afs = [float(x) for x in val[fmt.index("AF")].split(",") if x not in ("", ".")]
            af = max(afs)
        except ValueError:
            continue
        # Heteroplasmic: 5% to 95%. Below 5% NUMT reads and noise dominate.
        if not 0.05 <= af < 0.95:
            continue
        alts = r[4].split(",")
        if any(in_polyc(int(r[1]), r[3], a) for a in alts):
            polyc.append({"variant": f"{r[0]}:{r[1]} {r[3]}>{r[4]}", "af": round(af, 3)})
        else:
            het += 1
    ann = mito_annotation(d, s, p)
    listed, plp = [], 0
    for r in ann.pop("rows", []):
        terms = clinvar_terms(r.get("clinvar_significance"))
        is_plp = any(t in CLINVAR_PLP for t in terms)
        plp += is_plp
        try:
            het_row = 0.05 <= float(r.get("af")) < 0.95
        except (TypeError, ValueError):
            het_row = False
        if is_plp or het_row or (terms and not all(t in CLINVAR_BENIGN for t in terms)):
            listed.append({k: r.get(k) for k in MITO_ANNOT_COLS if k not in ("chrom", "filter")}
                          | {"polyc": in_polyc(int(r["pos"]), r["ref"], r["alt"])})
    if "source" in ann:
        ann.update(variants=listed, clinvar_plp=plp)
    return p, {"pass": passed, "heteroplasmic": het, "heteroplasmy_floor": 0.05, "polyc_heteroplasmic": polyc,
               "polyc_tracts": [f"m.{a}-{b}" for a, b in POLYC_TRACTS], "annotation": ann}


def cpsr_cdna_checks(d, s):
    """Step 17's cDNA read check of CPSR's P/LP rows (bin/clip_check.py
    --cpsr), keyed like step 06's; None when step 17 did not write it."""
    p = first_existing(d, [f"cpsr/{s}_cpsr_cdna_check.tsv"])
    return clinvar_cdna_checks(d, s, p) if p else None


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
            out.update(cpsr_compare(d, s, rows, cols, col))
        else:
            out["classification_column"] = None
    return p, out


def cpsr_compare(d, s, rows, cols, col):
    """CPSR's own class beside the final class and ClinVar's.

    The final class (CLASSIFICATION) is ClinVar's whenever ClinVar has a
    trusted record, so it hides CPSR's own call (CPSR_CLASSIFICATION). This
    lists the P/LP rows (final or own class), counts the final P/LP rows that
    CPSR's own rules class benign or likely benign, and lists every row where
    ClinVar and CPSR's own class disagree on P/LP. ClinVar's class is
    CLINVAR_CLASSIFICATION, or the final class where ASSERTION_AUTHORITY is
    ClinVar. A table with one class column has no disagreement to show.
    Low-confidence: DP_CONTROL below step 6's depth floor (CPSR writes no
    GQ). The cDNA read check of step 17 goes on each P/LP row."""
    own = "CPSR_CLASSIFICATION" if "CPSR_CLASSIFICATION" in cols and col != "CPSR_CLASSIFICATION" else None
    hits = first_existing(d, [f"clinvar/{s}_clinvar_hits.vcf"])
    floor = call_floor(hits) if hits else dict(CALL_FLOOR)
    cdna = cpsr_cdna_checks(d, s)
    norm = clip_check.norm_class
    plp, disagree, own_benign = [], [], 0
    for r in rows:
        final, mine = norm(r.get(col)), norm(r.get(own)) if own else ""
        authority = (r.get("ASSERTION_AUTHORITY") or "").strip()
        clinvar_raw = ((r.get("CLINVAR_CLASSIFICATION") or "").strip() if "CLINVAR_CLASSIFICATION" in cols
                       else (r.get(col) or "").strip() if authority.lower() == "clinvar" else "")
        clinvar = norm(clinvar_raw)
        if final in clip_check.PLP and mine in ("benign", "likely benign"):
            own_benign += 1
        is_plp = final in clip_check.PLP or mine in clip_check.PLP
        differ = bool(own and clinvar and mine) and (clinvar in clip_check.PLP) != (mine in clip_check.PLP)
        if not (is_plp or differ):
            continue
        v = clip_check.cpsr_variant(r)
        dp = (r.get("DP_CONTROL") or "").strip()
        dp = int(dp) if dp.isdigit() else None
        x = {"variant": f"{v[0]}:{v[1]} {v[2]}>{v[3]}" if v else (r.get("GENOMIC_CHANGE") or "."),
             "gene": r.get("SYMBOL") or ".", "final": r.get(col) or ".",
             "cpsr": (r.get(own) or ".") if own else None, "clinvar": clinvar_raw if clinvar else None,
             "authority": authority or None, "acmg_codes": (r.get("ACMG_CODE") or "").strip() or None,
             "dp": dp, "low_confidence": dp is not None and dp < floor["DP"]}
        if cdna is not None and is_plp:
            x["cdna"] = (cdna.get((v[0], v[1], v[2], v[3])) if v else None) or {
                "check": "not_checked", "note": "not in the read-check table", "alt_reads": None,
                "alt_clipped": None, "ref_reads": None, "ref_clipped": None}
        if is_plp:
            plp.append(x)
        if differ:
            disagree.append(x)
    out = {"plp": plp, "disagreements": disagree, "final_plp_own_benign": own_benign,
           "own_class_column": own, "low_confidence": sum(1 for x in plp if x["low_confidence"]),
           "call_floor_dp": floor["DP"]}
    if cdna is not None:
        out["cdna_flagged"] = sum(1 for x in plp if x["cdna"]["check"] == "flagged")
    return out


def csq_rank(fmt, csq):
    """The rule of bin/csq_rank.awk on one CSQ value: (class, entries, freq).
    class is MANE (entries of the MANE Select transcript), else CANONICAL,
    else ANY (every entry); freq is rare, common, unknown (no value) or none
    (the CSQ has no frequency field): MAX_AF, else gnomADe_AF and gnomADg_AF,
    the highest over the entries."""
    ix = {k: i for i, k in enumerate(fmt)}
    entries = [v for v in (t.split("|") for t in csq.split(",") if t) if len(v) == len(fmt)]

    def get(v, k):
        return v[ix[k]] if k in ix else ""
    mane = [v for v in entries if get(v, "MANE_SELECT")]
    canon = [v for v in entries if get(v, "CANONICAL") == "YES"]
    cls, chosen = ("MANE", mane) if mane else ("CANONICAL", canon) if canon else ("ANY", entries)
    fields = ["MAX_AF"] if "MAX_AF" in ix else [f for f in ("gnomADe_AF", "gnomADg_AF") if f in ix]
    afs = []
    for v in entries:
        for f in fields:
            for a in get(v, f).split("&"):
                try:
                    afs.append(float(a))
                except ValueError:
                    pass
    freq = "none" if not fields else "unknown" if not afs else "rare" if max(afs) < COMMON_AF else "common"
    return cls, [{"impact": get(v, "IMPACT"), "gene": get(v, "SYMBOL"),
                  "consequence": get(v, "Consequence").replace("&", ",")} for v in chosen], freq


def call_check(fmt, sample, floor):
    """bin/csq_rank.awk's genotype check of one sample: (weak, reasons, gq,
    dp, ad). Weak: GQ or DP below the floor, a heterozygote's alt read
    fraction outside HET_AB, a homozygous alt call's below HOM_AB."""
    f = dict(zip(fmt.split(":"), sample.split(":")))
    gq, dp = (int(f[k]) if f.get(k, "").isdigit() else None for k in ("GQ", "DP"))
    ad = f.get("AD", "")
    why = []
    if gq is not None and gq < floor["GQ"]:
        why.append(f"GQ {gq} < {floor['GQ']}")
    if dp is not None and dp < floor["DP"]:
        why.append(f"DP {dp} < {floor['DP']}")
    al = f.get("GT", "").replace("|", "/").split("/")
    depths = ad.split(",")
    if len(al) == 2 and all(a.isdigit() for a in al) and len(depths) > 1 and all(x.isdigit() for x in depths):
        tot = sum(int(x) for x in depths)
        if tot:
            ab = sum(int(x) for x in depths[1:]) / tot
            lo, hi = floor["HET_AB"]
            if (al[0] == "0") != (al[1] == "0") and not lo <= ab <= hi:
                why.append(f"alt reads {ab:.2f} of the depth on a heterozygote")
            elif al[0] == al[1] != "0" and ab < floor["HOM_AB"]:
                why.append(f"alt reads {ab:.2f} of the depth on a homozygous call")
    return bool(why), why, gq, dp, ad or None


def clinical_call_floor(path):
    """The floors bin/csq_rank.awk wrote into the clinical VCF, the defaults without the line."""
    floor = dict(CALL_FLOOR, **AB_FLOOR)
    try:
        with open_text(path) as f:
            for line in f:
                if not line.startswith("##"):
                    break
                if line.startswith("##PGP_CALL_FLOOR="):
                    for k, v in re.findall(r"\b(GQ|DP)=(\d+)", line):
                        floor[k] = int(v)
                    m = re.search(r"HET_AB=([0-9.]+)-([0-9.]+)", line)
                    if m:
                        floor["HET_AB"] = (float(m.group(1)), float(m.group(2)))
                    m = re.search(r"HOM_AB=([0-9.]+)", line)
                    if m:
                        floor["HOM_AB"] = float(m.group(1))
    except (OSError, EOFError, gzip.BadGzipFile, zlib.error) as e:
        raise Unreadable(f"{os.path.basename(path)}: {e}") from e
    return floor


def csq_high_impact_genes(path):
    """{(chrom, pos, ref, alt): record} of the records whose ranked CSQ entries
    (bin/csq_rank.awk's rule: the MANE Select transcript, else the canonical
    one, else any) hold a HIGH-impact consequence. record: genotype, genes
    ({gene: consequence}), transcript (MANE, CANONICAL or ANY), frequency
    (rare, common, unknown or none), weak (the call check failed), why, gq,
    dp and ad. A HIGH consequence on a transcript other than the MANE Select
    one does not count when the record has a MANE Select transcript."""
    floor = clinical_call_floor(path)
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
                cls, entries, freq = csq_rank(fmt, info_map(r[7]).get("CSQ", ""))
                genes = {}
                for e in entries:
                    if e["impact"] == "HIGH":
                        genes.setdefault(e["gene"], e["consequence"])
                if genes:
                    gt = r[9].split(":")[0] if len(r) > 9 else ""
                    weak, why, gq, dp, ad = call_check(r[8], r[9], floor) if len(r) > 9 else (False, [], None, None, None)
                    out[(r[0], int(r[1]), r[3], r[4])] = {
                        "genotype": gt, "genes": genes, "transcript": cls, "frequency": freq,
                        "weak": weak, "why": why, "gq": gq, "dp": dp, "ad": ad}
    except (OSError, EOFError, gzip.BadGzipFile, zlib.error) as e:
        raise Unreadable(f"{os.path.basename(path)}: {e}") from e
    return out


# --- the ACMG SF rules and the carrier list (assets/) --------------------------

_ASSETS = {}


def asset_rows(name):
    """The rows of assets/<name> (comment lines left out), read once."""
    p = os.path.join(ASSETS_DIR, name)
    if p in _ASSETS:
        return _ASSETS[p]
    if not os.path.isfile(p):
        raise Unreadable(f"{name} is not in {ASSETS_DIR}: the report needs the repository's assets/ folder "
                         "(render_report.py --assets)")
    with open(p, encoding="utf-8") as f:
        rows = list(csv.DictReader((x for x in f if not x.startswith("#")), delimiter="\t", quoting=csv.QUOTE_NONE))
    _ASSETS[p] = rows
    return rows


def allele_key(locus):
    """'chr6:26092913:G:A' -> ('chr6', 26092913, 'G', 'A')."""
    c, pos, ref, alt = locus.split(":")
    return (c, int(pos), ref, alt)


def parse_alleles(text, kinds=False):
    """The variants or severity cell of assets/acmg_sf.tsv and
    assets/carrier_genes.tsv: '.' or entries separated by ';', each
    '[kind ]chrom:pos:ref:alt name'. {key: {'kind', 'name'}}."""
    out = {}
    for e in (text or ".").split(";"):
        e = e.strip()
        if e in ("", "."):
            continue
        parts = e.split(" ", 2 if kinds else 1)
        kind = parts.pop(0) if kinds else None
        if kinds and kind not in SEVERITY_KINDS:
            raise ValueError(f"severity '{kind}' is not one of {', '.join(SEVERITY_KINDS)}: {e}")
        out[allele_key(parts[0])] = {"kind": kind, "name": parts[1] if len(parts) > 1 else parts[0]}
    return out


def acmg_sf_rules():
    """{gene: rule} of assets/acmg_sf.tsv."""
    rules = {}
    for n, r in enumerate(asset_rows("acmg_sf.tsv"), 1):
        try:
            g = r["gene"]
            if r["report_if"] not in SF_REPORT_IF or r["sv_cn"] not in ("yes", "no"):
                raise ValueError(f"report_if must be one of {', '.join(SF_REPORT_IF)}, sv_cn yes or no")
            if r["sv_cn"] == "yes" and r["report_if"] not in ("any", "het_female_too"):
                raise ValueError("sv_cn yes needs a rule that one variant meets (any, het_female_too)")
            regions = []
            for e in (r["sv_regions"] or ".").split(";"):
                e = e.strip()
                if e not in ("", "."):
                    loc, _, name = e.partition(" ")
                    c, _, span = loc.partition(":")
                    a, _, b = span.partition("-")
                    regions.append((c, int(a), int(b), name or loc))
            variants = parse_alleles(r["variants"])
            if r["report_if"] == "hom_variant" and not variants:
                raise ValueError("report_if hom_variant needs its variants")
            if g in rules:
                raise ValueError("the gene is on two rows")
            rules[g] = {"inheritance": r["inheritance"], "report_if": r["report_if"], "variants": variants,
                        "sv_cn": r["sv_cn"] == "yes", "sv_regions": regions,
                        "severity": parse_alleles(r["severity"], kinds=True),
                        "table_text": r["table_text"], "note": None if r["note"] in ("", ".") else r["note"]}
        except (KeyError, ValueError, TypeError) as e:
            raise Unreadable(f"acmg_sf.tsv row {n}: {e}") from e
    if not rules:
        raise Unreadable("acmg_sf.tsv has no gene")
    return rules


def carrier_genes():
    """{gene: row} of assets/carrier_genes.tsv, the severity cell parsed."""
    out = {}
    for n, r in enumerate(asset_rows("carrier_genes.tsv"), 1):
        try:
            out[r["gene"]] = {"inheritance": r["inheritance"], "lists": r["lists"], "condition": r["condition"],
                              "severity": parse_alleles(r.get("severity"), kinds=True)}
        except (KeyError, ValueError, TypeError) as e:
            raise Unreadable(f"carrier_genes.tsv row {n}: {e}") from e
    return out


def allele_copies(gt):
    """(alt copies, called copies) of a GT: '0/1' (1, 2), '1/1' (2, 2), '1' (1, 1)."""
    a = [x for x in re.split(r"[/|]", gt or "") if x not in (".", "")]
    return sum(1 for x in a if x != "0"), len(a)


def zygosity_word(copies, ploidy):
    if not copies:
        return "unknown"
    if ploidy == 1:
        return "hemizygous"
    return "hom" if copies >= 2 else "het"


def sf_verdicts(gene, rule, alleles):
    """Which of a gene's variants the ACMG rule reports. alleles: {key:
    (copies, ploidy)} of the gene's distinct small variants, structural ones
    keyed by their span. {key: None (reported) or the reason it is not}, and
    the partial-deficiency note of a reported gene with a hypomorph allele."""
    sev = rule["severity"]
    out, kept = {}, {}
    for k, cp in alleles.items():
        if sev.get(k, {}).get("kind") == "pseudodeficiency":
            out[k] = f"{sev[k]['name']}: a pseudodeficiency allele, not counted"
        else:
            kept[k] = cp
    n = sum(c for c, _ in kept.values())
    ri = rule["report_if"]
    if ri in ("any", "het_female_too"):
        ok = {k: True for k in kept}
    elif ri == "biallelic":
        ok = {k: n >= 2 for k in kept}
    elif ri == "hemizygous":
        hemi = any(p == 1 for _, p in kept.values())
        ok = {k: hemi or n >= 2 for k in kept}
    else:  # hom_variant
        ok = {k: k in rule["variants"] and c >= 2 and p == 2 for k, (c, p) in kept.items()}
    for k, good in ok.items():
        if good:
            out[k] = None
        elif ri == "hom_variant":
            out[k] = f"{gene} is reported only for: {rule['table_text']}"
        elif ri == "hemizygous":
            out[k] = f"one heterozygous variant; {gene} is reported for {rule['table_text']}"
        else:
            out[k] = f"one heterozygous variant; {gene} is reported with {rule['table_text']}"
    partial = None
    hypo = [sev[k]["name"] for k, why in out.items() if why is None and sev.get(k, {}).get("kind") == "hypomorph"]
    if hypo:
        partial = (f"partial, not profound: {'; '.join(hypo)} is a hypomorph allele, so the deficiency it causes "
                   "is partial")
    return out, partial


def sf_structural(d, s, rules):
    """Step 05's class 4 and 5 SVs (AnnotSV) in genes whose rule allows them
    (sv_cn yes), and the class 4 or 5 deletions over a gene's sv_regions (the
    3' end of EPCAM for MSH2). Common SVs are left out. (rows, checked, note)."""
    try:
        p, data = sec_annotsv(d, s)
    except (Unreadable, OSError, ValueError, KeyError) as e:
        return [], False, f"the AnnotSV table could not be read: {e}"
    if not p:
        return [], False, None
    out = []
    for x in data["rows"]:
        if x["common"] or str(x["acmg_class"]) not in ("4", "5"):
            continue
        chrom, _, span = x["sv"].partition(":")
        a, _, b = span.partition("-")
        chrom = chrom if chrom.startswith("chr") else "chr" + chrom
        hits = {}
        for g in x["genes"]:
            if g in rules and rules[g]["sv_cn"]:
                hits[g] = None
        for g, rule in rules.items():
            for c, rs, re_, name in rule["sv_regions"]:
                if (x["type"] == "DEL" and c == chrom and a.isdigit() and b.isdigit()
                        and int(a) <= re_ and int(b) >= rs and g not in hits):
                    hits[g] = name
        for g, via in sorted(hits.items()):
            out.append({"gene": g, "variant": f"{x['type']} {chrom}:{a}-{b}", "acmg_class": str(x["acmg_class"]),
                        "caller": x["caller"], "support": x["support"], "via": via,
                        "key": ("sv", chrom, a, b, x["type"])})
    return out, True, None


def acmg_sf_tier(d, s, clinical_vcf):
    """The variants that are ACMG secondary findings by the rule of their gene
    (assets/acmg_sf.tsv): the ClinVar P/LP hits of step 6, the clinical
    filter's rare HIGH-impact records, and step 05's class 4 and 5 SVs in the
    genes whose rule allows them. A variant in an SF gene that its rule does
    not report (one MUTYH variant, a heterozygous HFE p.C282Y) is listed
    under not_reported with the reason."""
    rules = acmg_sf_rules()
    clinvar, high = [], []
    hits = first_existing(d, [f"clinvar/{s}_clinvar_hits.vcf"])
    if hits:
        floor = call_floor(hits)
        for r in vcf_records(hits):
            if len(r) < 8:
                continue
            row = clinvar_row(r, floor)
            copies, ploidy = allele_copies(r[9].split(":")[0] if len(r) > 9 else "")
            for g in row["gene"].split(","):
                if g in rules:
                    clinvar.append({"gene": g, "variant": f"{row['chrom']}:{row['pos']} {row['ref']}>{row['alt']}",
                                    "genotype": row["genotype"], "significance": row["significance"],
                                    "stars": row["stars"], "low_confidence": row["low_confidence"],
                                    "key": (row["chrom"], row["pos"], row["ref"], row["alt"]),
                                    "copies": copies, "ploidy": ploidy})
    if clinical_vcf:
        # A common record is in the clinical VCF through the ClinVar tier,
        # not as a rare HIGH impact; one of unknown frequency is listed and
        # marked, never counted as rare; a weak call is listed and marked.
        for (c, pos, ref, alt), x in sorted(csq_high_impact_genes(clinical_vcf).items()):
            if x["frequency"] == "common":
                continue
            for g, cons in sorted(x["genes"].items()):
                if g in rules:
                    gt = x["genotype"]
                    copies, ploidy = allele_copies(gt)
                    zyg = {"0/1": "het", "1/0": "het", "0|1": "het", "1|0": "het", "1/1": "hom", "1|1": "hom"}.get(gt, gt)
                    high.append({"gene": g, "variant": f"{c}:{pos} {ref}>{alt}", "genotype": zyg, "consequence": cons,
                                 "transcript": x["transcript"], "frequency": x["frequency"],
                                 "frequency_unknown": x["frequency"] == "unknown",
                                 "weak_call": x["weak"], "call_note": "; ".join(x["why"]) or None,
                                 "gq": x["gq"], "dp": x["dp"], "ad": x["ad"],
                                 "key": (c, pos, ref, alt), "copies": copies, "ploidy": ploidy})
    svs, sv_checked, sv_note = sf_structural(d, s, rules)
    # One verdict per gene, over its distinct variants: a ClinVar hit that is
    # also a HIGH-impact record is one variant.
    alleles = {}
    for x in clinvar + high:
        alleles.setdefault(x["gene"], {}).setdefault(x["key"], (max(x["copies"], 1), x["ploidy"] or 2))
    for x in svs:
        alleles.setdefault(x["gene"], {}).setdefault(x["key"], (1, 2))
    verdict = {g: sf_verdicts(g, rules[g], a) for g, a in alleles.items()}
    reported = {"clinvar": [], "high": [], "sv": []}
    not_reported = []
    for kind, rows, source in (("clinvar", clinvar, "ClinVar"), ("high", high, "HIGH impact"), ("sv", svs, "AnnotSV")):
        for x in rows:
            why, partial = verdict[x["gene"]][0][x["key"]], verdict[x["gene"]][1]
            rule = rules[x["gene"]]
            for k in ("key", "copies", "ploidy"):
                x.pop(k, None)
            if why is None:
                x.update(rule=rule["table_text"], rule_note=rule["note"], severity=partial)
                reported[kind].append(x)
            else:
                not_reported.append({"gene": x["gene"], "variant": x["variant"], "genotype": x.get("genotype"),
                                     "source": source, "reason": why})
    high = reported["high"]
    sure = [x for x in high if not (x["weak_call"] or x["frequency_unknown"])]
    return {"version": ACMG_SF_VERSION, "genes_on_list": len(rules), "rules": "assets/acmg_sf.tsv",
            "clinvar_hits": reported["clinvar"], "clinvar_checked": bool(hits),
            "high_impact": high, "high_impact_checked": bool(clinical_vcf),
            "high_impact_rare": len(sure),
            "high_impact_weak": sum(1 for x in high if x["weak_call"]),
            "high_impact_frequency_unknown": sum(1 for x in high if x["frequency_unknown"]),
            "structural": reported["sv"], "structural_checked": sv_checked, "structural_note": sv_note,
            "not_reported": not_reported}


def carrier_list(p):
    """The sample's own carrier list: the ClinVar hits of step 6 (p) in the
    genes of assets/carrier_genes.tsv, one person, no partner needed. An
    allele the list marks hypomorph or pseudodeficiency is set aside, never a
    carrier finding on its own."""
    genes = carrier_genes()
    carriers, aside = [], []
    for r in vcf_records(p):
        if len(r) < 8:
            continue
        row = clinvar_row(r)
        copies, ploidy = allele_copies(r[9].split(":")[0] if len(r) > 9 else "")
        key = (row["chrom"], row["pos"], row["ref"], row["alt"])
        for g in row["gene"].split(","):
            if g not in genes:
                continue
            x = genes[g]
            item = {"gene": g, "variant": f"{row['chrom']}:{row['pos']} {row['ref']}>{row['alt']}",
                    "genotype": zygosity_word(copies, ploidy), "significance": row["significance"],
                    "stars": row["stars"], "inheritance": x["inheritance"], "condition": x["condition"]}
            sev = x["severity"].get(key)
            if sev:
                item.update(kind=sev["kind"], name=sev["name"])
                aside.append(item)
            else:
                carriers.append(item)
    return {"genes_on_list": len(genes), "list": "assets/carrier_genes.tsv",
            "carriers": carriers, "set_aside": aside}


def acmg_sf_or_note(d, s, clinical_vcf):
    """acmg_sf_tier, or the reason it could not run (a missing or broken
    assets/acmg_sf.tsv), which both reports print: the clinical counts stay."""
    try:
        return acmg_sf_tier(d, s, clinical_vcf)
    except Unreadable as e:
        return {"version": ACMG_SF_VERSION, "unavailable": str(e), "clinvar_hits": [], "high_impact": [],
                "structural": [], "not_reported": []}


def carrier_list_or_note(p):
    """carrier_list, or the reason it could not run, which both reports print:
    the ClinVar hits stay."""
    try:
        return carrier_list(p)
    except Unreadable as e:
        return {"unavailable": str(e)}


def sec_clinical(d, s):
    p = first_existing(d, [f"clinical/{s}_clinical_summary.tsv"])
    v = first_existing(d, [f"clinical/{s}_clinical.vcf.gz"])
    if not p:
        # The Nextflow report gets the clinical VCF, not the summary table:
        # the count only (one record per variant, as in the table).
        if not v:
            return None, {}
        return v, {"variants": sum(1 for _ in vcf_records(v)), "acmg_sf": acmg_sf_or_note(d, s, v)}
    rows = read_tsv(p)
    by = {}
    for r in rows:
        by[r.get("IMPACT") or "."] = by.get(r.get("IMPACT") or ".", 0) + 1
    genes = sorted({r.get("GENE") for r in rows if r.get("GENE") not in (None, "", ".")})
    out = {"variants": len(rows), "by_impact": by, "genes": len(genes), "acmg_sf": acmg_sf_or_note(d, s, v)}
    # Step 23's tier and call check columns; a table written before them has neither.
    if rows and "FREQ_TIER" in rows[0]:
        out["frequency_unknown"] = sum(1 for r in rows if r.get("FREQ_TIER") == "unknown")
    if rows and "CALL_CHECK" in rows[0]:
        out["weak_calls"] = sum(1 for r in rows if r.get("CALL_CHECK") == "weak")
    return p, out


def sec_slivar(d, s):
    p = first_existing(d, [f"slivar/{s}_slivar_summary.tsv"])
    if not p:
        # The Nextflow report gets the prioritized VCF: the count, plus the
        # phase when COMPHET_PHASE's pairs table sits beside it (HTML_REPORT
        # stages it there).
        v = first_existing(d, [f"slivar/{s}_prioritized.vcf.gz"])
        if not v:
            return None, {}
        out = {"prioritized": sum(1 for _ in vcf_records(v))}
        out.update(comphet_phase(os.path.join(os.path.dirname(v), f"{s}_compound_het_pairs.tsv")))
        return v, out
    rows = read_tsv(p)
    out = {"prioritized": len(rows)}
    # Step 31's frequency tier: prioritized for the compound-het search, not rare.
    if rows and "FREQ_TIER" in rows[0]:
        out["frequency_unknown"] = sum(1 for r in rows if r.get("FREQ_TIER") == "unknown")
    ch = os.path.join(os.path.dirname(p), f"{s}_compound_hets.tsv")
    if os.path.isfile(ch):
        crow = read_tsv(ch)
        out["compound_het_variants"] = len(crow)
        out["compound_het_genes"] = len({r.get("GENE") for r in crow if r.get("GENE")})
    out.update(comphet_phase(os.path.join(os.path.dirname(p), f"{s}_compound_het_pairs.tsv")))
    return p, out


# The trans pairs the report lists by name; the rest are counted.
COMPHET_TRANS_SHOWN = 20


def comphet_phase(path):
    """Step 31's pairs table (bin/comphet_phase.py): the pairs by phase from
    the reads, and the trans ones (both copies of the gene hit) by name. A run
    from before the table has none of these keys."""
    if not os.path.isfile(path):
        return {}
    rows = read_tsv(path)
    phase = {k: sum(1 for r in rows if r.get("PHASE") == k) for k in ("trans", "cis", "unknown")}

    def n(x):
        return int(x) if (x or "").isdigit() else 0
    trans = [{"gene": r.get("GENE", ""),
              "variants": f"{r.get('CHROM')}:{r.get('POS_1')} {r.get('REF_1')}>{r.get('ALT_1')} + "
                          f"{r.get('CHROM')}:{r.get('POS_2')} {r.get('REF_2')}>{r.get('ALT_2')}",
              "reads": f"{r.get('TRANS_READS')} of {n(r.get('TRANS_READS')) + n(r.get('CIS_READS'))}"}
             for r in rows if r.get("PHASE") == "trans"]
    return {"compound_het_pairs": len(rows), "compound_het_phase": phase,
            "compound_het_trans": trans[:COMPHET_TRANS_SHOWN]}


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


# --- mosaic loss of chrY or of one chrX (step 33, bin/mosaic_sex_loss.py) ----------

MOSAIC_CALLS = ("detected", "not_detected", "gain", "not_estimated")
MOSAIC_NUMBERS = ("depth_ratio", "depth_ratio_expected", "loss_depth", "loss_depth_low", "loss_depth_high",
                  "loss_af", "loss_af_low", "loss_af_high", "floor")


def sec_mosaic_sex_loss(d, s):
    """The mosaic_* rows of the step 33 table, without the prefix; none in a
    table written before the estimate existed."""
    p = first_existing(d, [f"qc/{s}_sample_qc.tsv"])
    if not p:
        return None, {}
    data = {r["key"][len("mosaic_"):]: r.get("value") or "" for r in read_tsv(p)
            if (r.get("key") or "").startswith("mosaic_")}
    if not data:
        return None, {}
    if data.get("call") not in MOSAIC_CALLS:
        raise Unreadable(f"{os.path.basename(p)} has no mosaic_call row ({', '.join(MOSAIC_CALLS)})")
    for k in MOSAIC_NUMBERS:
        v = num(data.get(k))
        data[k] = v if v is not None and abs(v) != float("inf") else None   # inf is no value either
    data["x_het_sites"] = int(data["x_het_sites"]) if (data.get("x_het_sites") or "").isdigit() else None
    return p, data


def sec_read_count(d, s):
    """Step 02's read count (and the ALIGN_MINIMAP2 task's): the reads that
    went into the aligner and the primary records of the BAM it wrote."""
    p = first_existing(d, [f"qc/{s}_read_count.tsv"])
    if not p:
        return None, {}
    data = {}
    for r in read_tsv(p):
        if r.get("key"):
            data[r["key"]] = r.get("value", "")
    for k in ("input_reads", "primary_records", "missing_reads"):
        if not re.fullmatch(r"-?[0-9]+", data.get(k, "")):
            raise Unreadable(f"{os.path.basename(p)} has no whole number for {k}")
        data[k] = int(data[k])
    if data.get("read_count_check") not in ("ok", "mismatch"):
        raise Unreadable(f"{os.path.basename(p)} has no read_count_check row (ok or mismatch)")
    return p, data


# Sections that are not a step of their own, or not one of one person (the
# couple screen needs two): no "not run" line without them.
NOT_A_STEP = {"read_count", "couple", "carrier"}

SECTIONS = [
    # key, title, step, reader
    ("coverage", "Coverage (mosdepth)", "16b", sec_coverage),
    ("sex_check", "Sex check (indexcov)", "16", sec_sex_check),
    ("sample_qc", "Sample identity and contamination (somalier, VerifyBamID2)", "33", sec_sample_qc),
    ("mosaic_sex_loss", "Mosaic loss of chromosome Y or X", "33", sec_mosaic_sex_loss),
    ("read_count", "Reads kept by alignment", "02", sec_read_count),
    ("variants", "Variant calling", "03", sec_variants),
    ("clinvar", "ClinVar screen", "06", sec_clinvar),
    ("couple", "Couple carrier screen", "39", sec_couple),
    ("carrier", "Carrier list", "06", sec_carrier),
    ("pharmcat", "PharmCAT", "07", sec_pharmcat),
    ("cpic", "CPIC drug-gene recommendations", "27", sec_cpic),
    ("pypgx", "pypgx", "32", sec_pypgx),
    ("cyrius", "Cyrius (CYP2D6)", "21", sec_cyrius),
    ("parascopy", "SMN1/SMN2 copy number (Parascopy, short reads)", "35", sec_parascopy),
    ("paralogs", "Paralog genes (Paraphase, PacBio HiFi)", "35b", sec_paralogs),
    ("hla", "HLA typing (T1K)", "08", sec_hla),
    ("prs", "Polygenic risk scores", "25", sec_prs),
    ("ancestry_fine", "Fine-grained ancestry (ADMIXTURE, RFMix)", "26b", sec_ancestry_fine),
    ("manta", "Structural variants (Manta)", "04", sec_manta),
    ("delly", "Structural variants (Delly)", "19", sec_delly),
    ("cnvpytor", "Copy number (CNVpytor)", "18", sec_cnvpytor),
    ("sv_consensus", "SV consensus merge", "22", sec_sv_consensus),
    ("annotsv", "Structural variants, class 4 and 5 (AnnotSV)", "05", sec_annotsv),
    ("expansions", "Repeat expansions (ExpansionHunter)", "09", sec_expansions),
    ("telomere", "Telomere content (TelomereHunter)", "10", sec_telomere),
    ("roh", "Runs of homozygosity", "11", sec_roh),
    ("haplogroup", "Mitochondrial haplogroup", "12", sec_haplogroup),
    ("y_haplogroup", "Y-chromosome haplogroup (Yleaf)", "37", sec_y_haplogroup),
    ("y_str", "Y-STR markers (HipSTR)", "45", sec_y_str),
    ("consumer_export", "Raw-data export (23andMe v5 layout)", "40", sec_consumer_export),
    ("bloodtype", "Blood groups (predicted)", "41", sec_bloodtype),
    ("lpa", "Lp(a): LPA KIV-2 repeat count", "46", sec_lpa),
    ("mito", "Mitochondrial variants", "20", sec_mito),
    ("traits", "Traits panel", "42", sec_traits),
    ("viral", "Viral DNA in the unmapped reads", "47", sec_viral),
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
    # pypgx's call alone, labelled, when step 36's rule allows it
    # (pgx_outside_calls.probable_cyp2d6); it never reaches PharmCAT.
    probable = None
    if pg_ok:
        pd = sections["pypgx"]["data"]
        ev = pd.get("cyp2d6_evidence") or {}
        probable = pgx_outside_calls.probable_cyp2d6(
            {"Diplotype": pg, "Phenotype": ev.get("phenotype"), "Defining_filter": ev.get("defining_filter"),
             "Quality": ev.get("quality")},
            (cy.get("genotype"), cy.get("filter")) if cy is not None else None,
            (pd.get("depth_check") or {}).get("status"))
    return {"calls": calls, "agree": agree, "consensus": consensus, "probable": probable}


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
    declared = (declared_sex or (status or {}).get("declared_sex") or manifest_value(manifest, "run", "declared_sex")
                or sections["sample_qc"]["data"].get("declared_sex") or None)
    cnvpytor_sex_filter(sections["cnvpytor"]["data"], declared)

    clinvar_date = manifest_value(manifest, "data", "clinvar_file_date")
    summary = {
        "schema_version": SCHEMA_VERSION,
        "sample": sample,
        "generated_utc": iso(datetime.now(tz=timezone.utc).timestamp()),
        "run": {
            "status_file": status is not None,
            "started_utc": status["started_utc"] if status else None,
            "declared_sex": declared,
        },
        "databases": {
            "clinvar_release": clinvar_date,
            "cpsr_clinvar_release": manifest_value(manifest, "data", "cpsr_clinvar_release"),
            "hla_database": manifest_value(manifest, "data", "hla_database"),
        },
        "manifest": manifest or [],
        "sections": sections,
        "cyp2d6": cyp2d6,
        # The read count is written by alignment, which a run from a BAM or
        # a VCF never does, so it is no step left undone.
        "not_run": [sections[k]["title"] + (f" (step {sections[k]['step']} failed)" if sections[k]["state"] == "failed" else "")
                    for k, *_ in SECTIONS if sections[k]["state"] in ("missing", "failed") and k not in NOT_A_STEP],
        "not_assessed": NOT_ASSESSED + ([] if sections["prs"]["data"].get("adjusted") else [PRS_NOT_ADJUSTED])
        # Step 35b on a region subset (--paraphase_regions) may leave hba out
        + ([] if {"HBA1", "HBA2"} <= {g["gene"] for g in sections["paralogs"]["data"].get("genes", [])
                                      if g["region"] == "hba" and g["status"] == "called"}
           or (sections["parascopy"]["data"].get("hba") or {}).get("status") == "called"
           else [HBA_NOT_ASSESSED]),
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
    if argv[:1] == ["mito-annotate"]:
        return mito_annotate_main(argv[1:])
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
