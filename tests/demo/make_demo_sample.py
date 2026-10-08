#!/usr/bin/env python3
"""make_demo_sample.py: write DEMO-001, the invented sample of the docs pictures.

  make_demo_sample.py --genome-dir DIR [--seed 1] [--records 4700000]
  make_demo_sample.py check --genome-dir DIR

The first form writes DIR/DEMO-001 with every file steps 24, 27 and 28 read,
as the steps that make them would have written them, so that every card of
the reports shows a value. Nothing is read from a real sample: the values come
from a random generator with a fixed seed and from a few hand-picked numbers
below. It also writes the three reference files the run manifest reads
(DIR/clinvar, DIR/t1k_idx, DIR/prs_scores), cut down to their header lines.

Borrowed, not rewritten:
  - the PharmCAT report starts from PharmCAT's own example report
    (tests/fixtures/pharmcat/pharmcat-docs-example.json, PharmCAT's invented
    example sample); HLA-A, HLA-B and CYP2D6 become the outside calls step 36
    passes on for DEMO-001, worked out by bin/pgx_outside_calls.py itself;
  - T1K's genotype table and the CYP2D6 depth check are the invented
    fixtures of tests/fixtures/pgx;
  - step 33's table is written by `collect_summary.py sample-qc`, and step 25's
    PRS and ancestry tables by `collect_summary.py prs-format` and `prs-table`,
    from made-up somalier, VerifyBamID2 and pgsc_calc files.

The second form reads DIR/DEMO-001/summary.json (step 24 writes it) and fails
when a section is not "ok", so a new report card that this script does not
feed yet is caught before a picture shows "Not run".

Standard library only. tests/demo/retake-screenshots.sh runs both forms.
"""
import argparse
import calendar
import copy
import gzip
import json
import math
import os
import random
import re
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "bin"))
import collect_summary  # noqa: E402
import pgx_outside_calls  # noqa: E402

SAMPLE = "DEMO-001"
SEX = "male"
PGX_FIXTURES = os.path.join(REPO, "tests", "fixtures", "pgx")
PHARMCAT_EXAMPLE = os.path.join(REPO, "tests", "fixtures", "pharmcat", "pharmcat-docs-example.json")

# GRCh38 primary chromosomes and their lengths.
CHROMS = [
    ("chr1", 248956422), ("chr2", 242193529), ("chr3", 198295559), ("chr4", 190214555),
    ("chr5", 181538259), ("chr6", 170805979), ("chr7", 159345973), ("chr8", 145138636),
    ("chr9", 138394717), ("chr10", 133797422), ("chr11", 135086622), ("chr12", 133275309),
    ("chr13", 114364328), ("chr14", 107043718), ("chr15", 101991189), ("chr16", 90338345),
    ("chr17", 83257441), ("chr18", 80373285), ("chr19", 58617616), ("chr20", 64444167),
    ("chr21", 46709983), ("chr22", 50818468), ("chrX", 156040895), ("chrY", 57227415),
]
AUTOSOMES = [c for c in CHROMS if c[0] not in ("chrX", "chrY")]
# Variants per base of each chromosome, relative to an autosome (one X, one Y).
DENSITY = {"chrX": 0.45, "chrY": 0.04}
MEAN_DEPTH = 31.3
# When the invented run started: 2026-10-01 09:00 UTC, the time of the demo's
# PharmCAT report too (logs/run_status.tsv, the report's "Latest run" line).
RUN_STARTED = calendar.timegm((2026, 10, 1, 9, 0, 0, 0, 0, 0))


def put(path, text, gz=False):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if gz:
        with gzip.open(path, "wt", compresslevel=6) as f:
            f.write(text)
    else:
        with open(path, "w") as f:
            f.write(text)


def vcf_header(extra=(), sample=True):
    lines = ["##fileformat=VCFv4.2"] + [f"##contig=<ID={c},length={n}>" for c, n in CHROMS] + list(extra)
    cols = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    lines.append(cols + ("\tFORMAT\t" + SAMPLE if sample else ""))
    return "\n".join(lines) + "\n"


def hla_release():
    """HLA_DB_RELEASE's default in scripts/lib/common.sh (3.65.0 -> '3.65.0')."""
    with open(os.path.join(REPO, "scripts", "lib", "common.sh")) as f:
        for line in f:
            m = re.match(r'HLA_DB_RELEASE=\$\{HLA_DB_RELEASE:-([0-9.]+)\}', line.strip())
            if m:
                return m.group(1)
    raise SystemExit("HLA_DB_RELEASE is not in scripts/lib/common.sh")


def image_tag(var):
    """The tag of an image line of versions.env (PHARMCAT_IMAGE -> 3.4.0)."""
    with open(os.path.join(REPO, "versions.env")) as f:
        for line in f:
            if line.startswith(var + "="):
                ref = line.split("=", 1)[1].split("#")[0].strip().strip('"')
                return ref.split("@")[0].rsplit(":", 1)[-1]
    raise SystemExit(f"{var} is not in versions.env")


# --- step 03: the variant calls ------------------------------------------------------

def variants(d, rng, n):
    """vcf/S.vcf.gz: N records at random positions, about 94% PASS, 85% SNVs."""
    weights = [length * DENSITY.get(c, 1.0) for c, length in CHROMS]
    total_w = sum(weights)
    path = os.path.join(d, "vcf", f"{SAMPLE}.vcf.gz")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    head = vcf_header(['##FILTER=<ID=RefCall,Description="Genotyping model thinks this site is reference.">',
                       '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                       '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">',
                       '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">'])
    bases = "ACGT"
    with gzip.open(path, "wt", compresslevel=1) as f:
        f.write(head)
        for (chrom, length), w in zip(CHROMS, weights):
            k = round(n * w / total_w)
            hap = chrom in ("chrX", "chrY") and SEX == "male"
            out = []
            for pos in sorted(rng.sample(range(10001, length - 10000), k)):
                ref = bases[rng.randrange(4)]
                r = rng.random()
                if r < 0.855:
                    alt = bases[(bases.index(ref) + 1 + rng.randrange(3)) % 4]
                elif r < 0.93:
                    alt = ref + "".join(bases[rng.randrange(4)] for _ in range(1 + rng.randrange(6)))
                else:
                    alt, ref = ref, ref + "".join(bases[rng.randrange(4)] for _ in range(1 + rng.randrange(6)))
                dp = max(3, int(rng.gauss(MEAN_DEPTH / (2 if hap else 1), 6)))
                if rng.random() < 0.942:
                    gt = "1/1" if hap or rng.random() < 0.38 else "0/1"
                    out.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{rng.randint(20, 66)}.{rng.randrange(10)}\tPASS\t.\t"
                               f"GT:GQ:DP\t{gt}:{rng.randint(20, 60)}:{dp}\n")
                else:
                    out.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{rng.randint(0, 9)}.{rng.randrange(10)}\tRefCall\t.\t"
                               f"GT:GQ:DP\t0/0:{rng.randint(5, 30)}:{dp}\n")
            f.write("".join(out))


# --- step 06: ClinVar screen -----------------------------------------------------------

def clinvar(d, genome_dir):
    """Three heterozygous hits in genes with invented names (the picture stops
    above the hits table, and no real gene is tied to an invented person)."""
    hits = [
        ("chr2", 47403191, "C", "T", "DEMOA:90001", "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts"),
        ("chr11", 108259012, "G", "A", "DEMOB:90002", "Likely_pathogenic", "criteria_provided,_single_submitter"),
        ("chr16", 2087914, "CT", "C", "DEMOC:90003", "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts"),
    ]
    body = "".join(f"{c}\t{p}\t{90000 + i}\t{r}\t{a}\t50\tPASS\tGENEINFO={g};CLNSIG={s};CLNREVSTAT={rev}\tGT\t0/1\n"
                   for i, (c, p, r, a, g, s, rev) in enumerate(hits))
    put(os.path.join(d, "clinvar", f"{SAMPLE}_clinvar_hits.vcf"), vcf_header() + body)
    # The ClinVar file itself: only the header line the run manifest reads.
    put(os.path.join(genome_dir, "clinvar", "clinvar.vcf.gz"),
        "##fileformat=VCFv4.1\n##fileDate=2026-09-28\n##source=ClinVar\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", gz=True)


# --- pharmacogenomics: steps 07, 08, 21, 32 (36 and 27 run for real) -------------------

def pgx(d, genome_dir, rng):
    """T1K, pypgx, Cyrius and the depth check, then the PharmCAT report that
    step 07 would write once step 36 passed their consensus on."""
    example = json.load(open(PHARMCAT_EXAMPLE))
    genes = example["genes"]

    def called(gene):
        g = genes[gene]["sourceDiplotypes"][0]
        return g["label"], "; ".join(g.get("phenotypes") or [])

    # T1K (step 08) and the CYP2D6 depth check: the invented fixtures.
    t1k = os.path.join(d, "hla_t1k", f"{SAMPLE}_hla_genotype.tsv")
    os.makedirs(os.path.dirname(t1k), exist_ok=True)
    shutil.copy(os.path.join(PGX_FIXTURES, "t1k_genotype.tsv"), t1k)
    depth = os.path.join(d, "pypgx", f"{SAMPLE}_cyp2d6_depth_check.tsv")
    os.makedirs(os.path.dirname(depth), exist_ok=True)
    shutil.copy(os.path.join(PGX_FIXTURES, "depth_check_ok.tsv"), depth)
    # The HLA database release the manifest prints: one header line of hla.dat,
    # naming the release setup.sh installs (HLA_DB_RELEASE in scripts/lib/common.sh).
    put(os.path.join(genome_dir, "t1k_idx", "hlaidx", "hla.dat"),
        f"ID   HLA00001; standard; DNA; HUM; 3503 BP.\nDT   15/07/2025 (Rel. {hla_release()}, Created)\n")

    # pypgx (step 32) and Cyrius (step 21) agree with PharmCAT's example calls,
    # CYP2D6 included, so the three CYP2D6 callers agree.
    cyp, cyp_phen = called("CYP2D6")
    rows = []
    for gene in ("CYP2B6", "CYP2C19", "CYP2C9", "CYP2D6", "CYP3A5", "NUDT15", "SLCO1B1", "TPMT", "UGT1A1"):
        dip, phen = called(gene)
        rows.append(f"{gene}\t{dip}\t{phen}\t{'Normal' if gene == 'CYP2D6' else 'N/A'}\t{'BAM' if gene == 'CYP2D6' else 'VCF'}\n")
    put(os.path.join(d, "pypgx", f"{SAMPLE}_pypgx_summary.tsv"), "Gene\tDiplotype\tPhenotype\tCNV_call\tSource\n" + "".join(rows))
    put(os.path.join(d, "cyrius", f"{SAMPLE}_cyp2d6.tsv"), f"Sample\tGenotype\tFilter\n{SAMPLE}\t{cyp}\tPASS\n")

    # What step 36 will pass to PharmCAT, from its own code.
    t1k_calls = pgx_outside_calls.read_t1k(t1k)
    outside = {}
    for gene in pgx_outside_calls.HLA_GENES:
        _, dip = pgx_outside_calls.hla_row(gene, t1k_calls)
        if dip:
            outside[gene] = dip
    _, dip = pgx_outside_calls.cyp2d6_row(os.path.join(d, "pypgx", f"{SAMPLE}_pypgx_summary.tsv"),
                                          os.path.join(d, "cyrius", f"{SAMPLE}_cyp2d6.tsv"), depth)
    if dip:
        outside["CYP2D6"] = dip

    # The PharmCAT report (step 07, run after step 36): PharmCAT's example,
    # with DEMO-001's outside calls. HLA phenotypes as PharmCAT writes them.
    report = copy.deepcopy(example)
    report.pop("_fixture_note", None)
    report["title"] = SAMPLE
    report["timestamp"] = "2026-10-01T09:00:00.000Z"
    report["pharmcatVersion"] = image_tag("PHARMCAT_IMAGE")
    hla_tests = {"HLA-A": ["*31:01"], "HLA-B": ["*15:02", "*57:01", "*58:01"]}
    changed = {}
    for gene, dip in outside.items():
        alleles = dip.split("/")
        if gene in hla_tests:
            phen = [f"{t} {'positive' if t in alleles else 'negative'}" for t in hla_tests[gene]]
        else:
            phen = [cyp_phen] if dip == cyp else ["Indeterminate"]
        if dip != called(gene)[0]:
            changed[gene] = dip
        report["genes"][gene]["callSource"] = "OUTSIDE"
        report["genes"][gene]["sourceDiplotypes"] = [{
            "gene": gene, "label": dip, "phenotypes": phen, "activityScore": None, "outsidePhenotype": False,
            "allele1": {"gene": gene, "name": alleles[0], "function": None},
            "allele2": {"gene": gene, "name": alleles[1], "function": None}}]
    # MT-RNR1 is an outside call in PharmCAT's example; this pipeline passes
    # none for it, so it is read from the VCF here, as the reference allele.
    mt = report["genes"]["MT-RNR1"]
    mt["callSource"] = "MATCHER"
    mt["sourceDiplotypes"] = [{"gene": "MT-RNR1", "label": "Reference",
                               "phenotypes": ["normal risk of aminoglycoside-induced hearing loss"],
                               "activityScore": None, "outsidePhenotype": False,
                               "allele1": {"gene": "MT-RNR1", "name": "Reference", "function": None}, "allele2": None}]
    changed["MT-RNR1"] = "Reference"
    prune_annotations(report, changed)
    put(os.path.join(d, "pharmcat", f"{SAMPLE}.report.json"), json.dumps(report, indent=1) + "\n")
    return outside


def prune_annotations(report, changed):
    """Drop the example's drug annotations written for a diplotype of a gene
    whose call changed: PharmCAT lists only those that fit the sample."""
    for src, by_drug in list(report.get("drugs", {}).items()):
        for name, rep in list(by_drug.items()):
            for gl in rep.get("guidelines") or []:
                keep = []
                for ann in gl.get("annotations") or []:
                    labels = {(dp.get("gene"), dp.get("label")) for gt in ann.get("genotypes") or []
                              for dp in (gt or {}).get("diplotypes") or [] if isinstance(dp, dict)}
                    stale = any(g in changed and lab != changed[g] for g, lab in labels)
                    stale = stale or any(g in changed and not any(x[0] == g for x in labels)
                                         for g in (ann.get("phenotypes") or {}))
                    if not stale:
                        keep.append(ann)
                gl["annotations"] = keep
            rep["guidelines"] = [gl for gl in rep.get("guidelines") or [] if gl.get("annotations")]
            if not rep["guidelines"]:
                del by_drug[name]
        if not by_drug:
            del report["drugs"][src]


# --- step 25: polygenic scores through pgsc_calc's files -------------------------------

def prs(d, genome_dir, rng):
    labels = collect_summary.read_labels(os.path.join(REPO, "assets", "pgs_scores.tsv"))
    harmonised = os.path.join(genome_dir, "prs_scores")
    for pid, trait in labels.items():
        n = rng.randint(250, 4000)
        rows = []
        for _ in range(n):
            c, length = AUTOSOMES[rng.randrange(len(AUTOSOMES))]
            ea, oa = rng.sample("ACGT", 2)
            rows.append(f"{c[3:]}\t{rng.randint(1, length)}\t{ea}\t{oa}\t{rng.gauss(0, 0.05):.5f}\t{c[3:]}\t{rng.randint(1, length)}\n")
        put(os.path.join(harmonised, f"{pid}.txt.gz"),
            f"#format_version=2.0\n#pgs_id={pid}\n#pgs_name={pid}\n#trait_reported={trait}\n#variants_number={n}\n"
            "#HmPOS_build=GRCh38\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\thm_chr\thm_pos\n"
            + "".join(rows), gz=True)
    work = os.path.join(d, "prs", "pgsc_calc")
    scores = os.path.join(work, "scores")
    if collect_summary.prs_format_main(["--scores", harmonised, "--labels", os.path.join(REPO, "assets", "pgs_scores.tsv"),
                                        "--out", scores, "--alleles", os.path.join(work, "score_alleles.tsv")]):
        raise SystemExit("prs-format failed")
    # pgsc_calc's output for sampleset "sample" (step 25's name) with the 1000 Genomes panel.
    res = os.path.join(work, "results", "sample")
    summ, score = ["dataset,accession,ambiguous,is_multiallelic,match_flipped,duplicate_best_match,"
                   "duplicate_ID,match_IDs,match_status,count,score_pass,match_rate"], []
    for pid in sorted(labels):
        with gzip.open(os.path.join(scores, f"{pid}.txt.gz"), "rt") as f:
            total = sum(1 for line in f if not line.startswith("#")) - 1
        m = round(total * rng.uniform(0.86, 0.99))
        summ.append(f"sample,{pid},false,false,false,false,false,NA,matched,{m},true,{m / total:.2f}")
        summ.append(f"sample,{pid},false,false,false,false,false,NA,unmatched,{total - m},true,{m / total:.2f}")
        score.append(f"sample\t{SAMPLE}\t{SAMPLE}\t{pid}\t{rng.gauss(0, 0.4):.4f}\t{rng.gauss(0, 1):.3f}\t"
                     f"{rng.uniform(18, 82):.2f}")
    put(os.path.join(res, "match", "sample_summary.csv"), "\n".join(summ) + "\n")
    put(os.path.join(res, "score", "sample_pgs.txt.gz"),
        "sampleset\tFID\tIID\tPGS\tSUM\tZ_MostSimilarPop\tpercentile_MostSimilarPop\n" + "\n".join(score) + "\n", gz=True)
    pcs = "\t".join(f"{rng.gauss(0, 0.01):.4f}" for _ in range(10))
    put(os.path.join(res, "score", "sample_popsimilarity.txt.gz"),
        "sampleset\tFID\tIID\t" + "\t".join(f"PC{i}" for i in range(1, 11))
        + "\tRF_P_AFR\tRF_P_AMR\tRF_P_EAS\tRF_P_EUR\tRF_P_SAS\tMostSimilarPop\tMostSimilarPop_LowConfidence\tREFERENCE\n"
        f"sample\t{SAMPLE}\t{SAMPLE}\t{pcs}\t0.00\t0.02\t0.00\t0.97\t0.01\tEUR\tFalse\tFalse\n", gz=True)
    os.makedirs(os.path.join(d, "ancestry"), exist_ok=True)
    if collect_summary.prs_table_main(["--sample", SAMPLE, "--results", os.path.join(work, "results"), "--sampleset", "sample",
                                       "--scores", scores, "--input-kind", "gvcf", "--panel", image_tag_panel(),
                                       "--ancestry-out", os.path.join(d, "ancestry", f"{SAMPLE}_ancestry.tsv"),
                                       "--out", os.path.join(d, "prs", f"{SAMPLE}_prs_summary.tsv")]):
        raise SystemExit("prs-table failed")


def image_tag_panel():
    with open(os.path.join(REPO, "versions.env")) as f:
        for line in f:
            if line.startswith("PGSC_PANEL="):
                return line.split("=", 1)[1].split("#")[0].strip().strip('"')
    return "pgsc_1000G_v1"


# --- steps 16, 16b, 33, 01b and samtools: QC ------------------------------------------

def coverage(d, rng):
    """mosdepth (step 16b) and indexcov (step 16)."""
    m = os.path.join(d, "mosdepth")
    summary = ["chrom\tlength\tbases\tmean\tmin\tmax"]
    dist = []
    tot_len = tot_bases = 0

    def cumulative(mean, sd, gap):
        # Share of bases with at least each depth: a normal curve, plus the
        # gaps of the reference (no reads) at depth 0.
        out = []
        for dep in range(0, int(mean * 2.6) + 1):
            p = (1 - gap) * 0.5 * math.erfc((dep - 0.5 - mean) / (sd * math.sqrt(2))) if dep else 1.0
            if p < 0.005 and dep:
                break
            out.append((dep, p))
        return out

    curves = {}
    for chrom, length in CHROMS:
        half = chrom in ("chrX", "chrY") and SEX == "male"
        mean = MEAN_DEPTH / (2 if half else 1) * rng.uniform(0.97, 1.03)
        gap = rng.uniform(0.01, 0.03) if chrom != "chrY" else 0.55
        bases = int(length * mean * (1 - gap))
        summary.append(f"{chrom}\t{length}\t{bases}\t{bases / length:.2f}\t0\t{int(mean * 9)}")
        tot_len += length
        tot_bases += bases
        curves[chrom] = cumulative(mean, mean * 0.24, gap)
    summary.append(f"total\t{tot_len}\t{tot_bases}\t{tot_bases / tot_len:.2f}\t0\t{int(MEAN_DEPTH * 9)}")
    for chrom, _ in CHROMS:
        dist += [f"{chrom}\t{dep}\t{p:.2f}" for dep, p in curves[chrom]]
    dist += [f"total\t{dep}\t{p:.2f}" for dep, p in (cumulative(MEAN_DEPTH, MEAN_DEPTH * 0.25, 0.02))]
    put(os.path.join(m, f"{SAMPLE}.mosdepth.summary.txt"), "\n".join(summary) + "\n")
    put(os.path.join(m, f"{SAMPLE}.mosdepth.global.dist.txt"), "\n".join(dist) + "\n")

    i = os.path.join(d, "indexcov")
    put(os.path.join(i, "indexcov-indexcov.ped"),
        "#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\tCNchrX\tCNchrY\tbins.out\tbins.lo\tbins.hi\t"
        "bins.in\tslope\tp.out\tPC1\tPC2\tPC3\tPC4\tPC5\n"
        f"unknown\t{SAMPLE}\t-9\t-9\t{1 if SEX == 'male' else 2}\t-9\t1.01\t0.98\t312\t205\t107\t174820\t0.993\t0.0018\t"
        "0\t0\t0\t0\t0\n")
    roc = ["#chrom\tcov\t" + SAMPLE]
    for chrom, _ in CHROMS:
        half = chrom in ("chrX", "chrY") and SEX == "male"
        centre = 0.5 if half else 1.0
        for k in range(0, 41):
            cov = k * 0.05
            frac = 0.5 * math.erfc((cov - centre) / (0.12 * math.sqrt(2)))
            roc.append(f"{chrom}\t{cov:.2f}\t{min(1.0, frac):.3f}")
    put(os.path.join(i, "indexcov-indexcov.roc"), "\n".join(roc) + "\n")


def sample_qc(d):
    """somalier and VerifyBamID2 files, then step 33's own table writer.

    sex is the sex somalier infers from the reads. original_pedigree_sex
    stays -9 (unknown), as in a real run: step 33 gives somalier no pedigree
    file, so MultiQC's somalier "Sex" column (that pedigree sex) reads -9."""
    q = os.path.join(d, "qc")
    sex = 1 if SEX == "male" else 2
    put(os.path.join(q, "somalier", f"{SAMPLE}.samples.tsv"),
        "#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\toriginal_pedigree_sex\tgt_depth_mean\t"
        "gt_depth_sd\tdepth_mean\tdepth_sd\tab_mean\tab_std\tn_hom_ref\tn_het\tn_hom_alt\tn_unknown\tp_middling_ab\t"
        "X_depth_mean\tX_n\tX_hom_ref\tX_het\tX_hom_alt\tY_depth_mean\tY_n\n"
        f"{SAMPLE}\t{SAMPLE}\t-9\t-9\t{sex}\t-9\t-9\t31.4\t7.2\t31.4\t7.2\t0.41\t0.08\t9817\t6904\t4383\t196\t0.01\t"
        "15.6\t612\t352\t4\t256\t14.9\t41\n")
    put(os.path.join(q, "somalier", f"{SAMPLE}.pairs.tsv"),
        "#sample_a\tsample_b\trelatedness\tibs0\tibs2\thom_concordance\thets_a\thets_b\thets_ab\tshared_hets\t"
        "hom_alts_a\thom_alts_b\tshared_hom_alts\tn\tx_ibs0\tx_ibs2\texpected_relatedness\n")
    put(os.path.join(q, "verifybamid2", f"{SAMPLE}.selfSM"),
        "#SEQ_ID\tRG\tCHIP_ID\t#SNPS\t#READS\tAVG_DP\tFREEMIX\tFREELK1\tFREELK0\tFREE_RH\tFREE_RA\tCHIPMIX\tCHIPLK1\t"
        "CHIPLK0\tCHIP_RH\tCHIP_RA\tDPREF\tRDPHET\tRDPALT\n"
        f"{SAMPLE}\tNA\tNA\t99846\t3087112\t30.92\t0.00412\t-1418823.2\t-1419102.7\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
    if collect_summary.sample_qc_main(["--sample", SAMPLE, "--somalier-samples", os.path.join(q, "somalier", f"{SAMPLE}.samples.tsv"),
                                       "--somalier-pairs", os.path.join(q, "somalier", f"{SAMPLE}.pairs.tsv"),
                                       "--selfsm", os.path.join(q, "verifybamid2", f"{SAMPLE}.selfSM"),
                                       "--declared-sex", SEX, "--marker-check", "passed",
                                       "--out", os.path.join(q, f"{SAMPLE}_sample_qc.tsv")]):
        raise SystemExit("sample-qc failed")


def reads(d, rng):
    """fastp's JSON (step 01b) and samtools flagstat (step 28 writes it from the BAM)."""
    pairs = 412_906_118
    before = pairs * 2
    passed = int(before * 0.9891)
    low_q, too_short, n_reads = int(before * 0.0081), int(before * 0.0026), before - passed - int(before * 0.0081) - int(before * 0.0026)
    length = 151

    def curves(after):
        q = [round(min(37.4, 31.5 + 6 * (1 - math.exp(-i / 4))) - (i / length) ** 3 * (2.2 if after else 3.5)
                   + rng.uniform(-0.15, 0.15), 2) for i in range(length)]
        base = {b: [round(v + rng.uniform(-0.004, 0.004), 4) for v in [p] * length]
                for b, p in (("A", 0.295), ("T", 0.296), ("C", 0.204), ("G", 0.205))}
        content = dict(base, N=[0.0002] * length, GC=[round(base["C"][i] + base["G"][i], 4) for i in range(length)])
        return {"total_reads": (passed if after else before) // 2, "total_bases": (passed if after else before) // 2 * length,
                "q20_bases": int((passed if after else before) // 2 * length * 0.968),
                "q30_bases": int((passed if after else before) // 2 * length * 0.921),
                "total_cycles": length, "quality_curves": dict({b: q for b in "ATCG"}, mean=q),
                "content_curves": content}

    def stats(n, gc):
        return {"total_reads": n, "total_bases": n * length, "q20_bases": int(n * length * 0.968),
                "q30_bases": int(n * length * 0.921), "q20_rate": 0.968, "q30_rate": 0.921,
                "read1_mean_length": length, "read2_mean_length": length, "gc_content": gc}

    hist = [0] * 1000
    for i in range(1000):
        hist[i] = int(4.2e6 * math.exp(-((i - 382) / 118) ** 2 / 2))
    fastp = {
        "summary": {"fastp_version": "1.0.1", "sequencing": "paired end (151 cycles + 151 cycles)",
                    "before_filtering": stats(before, 0.4105), "after_filtering": stats(passed, 0.4098)},
        "filtering_result": {"passed_filter_reads": passed, "low_quality_reads": low_q, "too_many_N_reads": n_reads,
                             "too_short_reads": too_short, "too_long_reads": 0},
        "duplication": {"rate": 0.0712},
        "insert_size": {"peak": 382, "unknown": 5102117, "histogram": hist},
        "adapter_cutting": {"adapter_trimmed_reads": int(before * 0.074), "adapter_trimmed_bases": int(before * 1.9),
                            "read1_adapter_sequence": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                            "read2_adapter_sequence": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"},
        "read1_before_filtering": curves(False), "read2_before_filtering": curves(False),
        "read1_after_filtering": curves(True), "read2_after_filtering": curves(True),
        "command": f"fastp -i /genome/{SAMPLE}/fastq/{SAMPLE}_R1.fastq.gz -I /genome/{SAMPLE}/fastq/{SAMPLE}_R2.fastq.gz "
                   f"-o /genome/{SAMPLE}/fastq_trimmed.part/{SAMPLE}_R1.fastq.gz -O /genome/{SAMPLE}/fastq_trimmed.part/{SAMPLE}_R2.fastq.gz "
                   f"-j /genome/{SAMPLE}/fastq_trimmed.part/{SAMPLE}_fastp.json -h /genome/{SAMPLE}/fastq_trimmed.part/{SAMPLE}_fastp.html",
    }
    put(os.path.join(d, "fastq_trimmed", f"{SAMPLE}_fastp.json"), json.dumps(fastp, indent=1) + "\n")

    primary = passed
    supp = int(primary * 0.0034)
    total = primary + supp
    mapped = int(primary * 0.9962) + supp
    dups = int(primary * 0.0712)
    paired = primary
    proper = int(paired * 0.981)
    single = int(paired * 0.0019)
    both = int(paired * 0.9943)
    diff = int(paired * 0.0061)
    flag = [
        f"{total} + 0 in total (QC-passed reads + QC-failed reads)", f"{primary} + 0 primary",
        "0 + 0 secondary", f"{supp} + 0 supplementary", f"{dups} + 0 duplicates", f"{dups} + 0 primary duplicates",
        f"{mapped} + 0 mapped ({100 * mapped / total:.2f}% : N/A)",
        f"{mapped - supp} + 0 primary mapped ({100 * (mapped - supp) / primary:.2f}% : N/A)",
        f"{paired} + 0 paired in sequencing", f"{paired // 2} + 0 read1", f"{paired // 2} + 0 read2",
        f"{proper} + 0 properly paired ({100 * proper / paired:.2f}% : N/A)",
        f"{both} + 0 with itself and mate mapped", f"{single} + 0 singletons ({100 * single / paired:.2f}% : N/A)",
        f"{diff} + 0 with mate mapped to a different chr", f"{int(diff * 0.62)} + 0 with mate mapped to a different chr (mapQ>=5)",
    ]
    put(os.path.join(d, "aligned", f"{SAMPLE}_flagstat.txt"), "\n".join(flag) + "\n")


# --- structural variants: steps 04, 18, 19, 22 ------------------------------------------

def sv_records(rng, n, caller, pass_rate):
    out = []
    for i in range(n):
        c, length = CHROMS[rng.randrange(len(CHROMS) - 1)]
        t = rng.choices(["DEL", "DUP", "INV", "INS"], weights=[62, 14, 6, 18])[0]
        pos = rng.randint(20000, length - 2_000_000)
        size = int(math.exp(rng.uniform(math.log(60), math.log(400000))))
        end = pos + (1 if t == "INS" else size)
        svlen = -size if t == "DEL" else size
        filt = "PASS" if rng.random() < pass_rate else rng.choice(["MinQUAL", "MinGQ", "LowQual"])
        gt = "1/1" if rng.random() < 0.2 else "0/1"
        out.append((c, pos, f"{caller}{t}:{i}", t, end, svlen, filt, gt, rng.randint(20, 999)))
    out.sort(key=lambda r: ([x[0] for x in CHROMS].index(r[0]), r[1]))
    return out


def structural(d, rng):
    info = ['##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
            '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length of the SV">']

    def write(path, recs, extra=()):
        body = "".join(f"{c}\t{p}\t{i}\tN\t<{t}>\t{q}\t{f}\tEND={e};SVTYPE={t};SVLEN={l}{x}\tGT\t{g}\n"
                       for (c, p, i, t, e, l, f, g, q), x in zip(recs, extra or [""] * len(recs)))
        put(path, vcf_header(info + ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']) + body, gz=True)

    write(os.path.join(d, "manta", "results", "variants", "diploidSV.vcf.gz"),
          sv_records(rng, rng.randint(7600, 8000), "Manta", 0.82))
    write(os.path.join(d, "delly", f"{SAMPLE}_sv.vcf.gz"), sv_records(rng, rng.randint(5900, 6300), "Delly", 0.86))

    # SURVIVOR's merge of Manta, Delly and CNVpytor (step 22): SUPP and SUPP_VEC.
    merged = sv_records(rng, rng.randint(2700, 3100), "SURVIVOR", 1.0)
    vecs = [rng.choices(["110", "101", "011", "111"], weights=[70, 9, 6, 15])[0] for _ in merged]
    write(os.path.join(d, "sv_merged", f"{SAMPLE}_sv_consensus.vcf.gz"), merged,
          [f";SUPP={v.count('1')};SUPP_VEC={v}" for v in vecs])

    lines = []
    for _ in range(rng.randint(1100, 1250)):
        c, length = CHROMS[rng.randrange(len(CHROMS) - 1)]
        t = "deletion" if rng.random() < 0.83 else "duplication"
        size = rng.randrange(1000, 400000, 100)
        start = rng.randrange(10001, length - size, 100)
        cn = rng.uniform(0.02, 0.6) if t == "deletion" else rng.uniform(1.3, 2.1)
        p = 10 ** rng.uniform(-40, -0.5)
        lines.append(f"{t}\t{c}:{start}-{start + size - 1}\t{size}\t{cn:.4f}\t{p:.4e}\t{p * 3:.4e}\t{min(1, p * 5):.4e}\t"
                     f"{min(1, p * 9):.4e}\t{rng.uniform(0, 0.5):.4f}\t{rng.uniform(0, 0.1):.4f}\t{rng.randint(0, 9)}")
    put(os.path.join(d, "cnvpytor", f"{SAMPLE}_cnvs.txt"), "\n".join(lines) + "\n")


# --- repeats, telomeres, ROH, mitochondria, Y: steps 09, 10, 11, 12, 20, 37 --------------

# ExpansionHunter's key loci with the normal range each repeat count is drawn from.
EH_RANGES = [("HTT", "chr4", 3074877, "CAG", 10, 26), ("FMR1", "chrX", 147912051, "CGG", 20, 40),
             ("C9ORF72", "chr9", 27573529, "GGCCCC", 2, 10), ("ATXN1", "chr6", 16327636, "TGC", 25, 35),
             ("DMPK", "chr19", 45770205, "CAG", 5, 30), ("ATXN3", "chr14", 92071011, "GCT", 14, 30),
             ("AR", "chrX", 67545317, "GCA", 17, 26), ("FXN", "chr9", 69037287, "GAA", 6, 30)]


def repeats_and_more(d, rng):
    eh = []
    for rep, c, pos, unit, lo, hi in EH_RANGES:
        hap = c == "chrX" and SEX == "male"
        a = sorted(rng.randint(lo, hi) for _ in range(1 if hap else 2))
        ref = (lo + hi) // 2
        gt = "1" if hap else ("1/2" if a[0] != a[1] else "1/1")
        cn = "/".join(map(str, a))
        alts = ",".join(f"<STR{x}>" for x in dict.fromkeys(a))
        eh.append(f"{c}\t{pos}\t.\t{unit[0]}\t{alts}\t.\tPASS\tEND={pos + ref * len(unit)};REF={ref};RL={ref * len(unit)};"
                  f"RU={unit};VARID={rep};REPID={rep}\tGT:SO:REPCN:LC\t{gt}:SPANNING:{cn}:{rng.uniform(28, 34):.6f}\n")
    put(os.path.join(d, "expansion_hunter", f"{SAMPLE}_eh.vcf"),
        vcf_header(['##ALT=<ID=STR,Description="Short tandem repeat">']) + "".join(eh))

    put(os.path.join(d, "telomere", SAMPLE, f"{SAMPLE}_summary.tsv"),
        "PID\tsample\ttotal_reads\tread_lengths\trepeat_threshold_set\trepeat_threshold_used\tintratelomeric_reads\t"
        "gc_bins_for_correction\ttotal_reads_with_tel_gc\ttel_content\n"
        f"{SAMPLE}\ttumor\t825812236\t151\tn=7 per 100 bp\t10\t61873\t48-52\t185204533\t{rng.uniform(260, 480):.2f}\n")

    segs, total = [], 0
    while total < 44e6:
        c, length = AUTOSOMES[rng.randrange(len(AUTOSOMES))]
        size = int(rng.uniform(0.5e6, 3.9e6))
        start = rng.randint(1_000_000, length - size - 1_000_000)
        segs.append((AUTOSOMES.index((c, length)), c, start, start + size, size))
        total += size
    segs.sort()
    put(os.path.join(d, "vcf", f"{SAMPLE}_roh.txt"),
        "# RG, Regions [2]Sample\t[3]Chromosome\t[4]Start\t[5]End\t[6]Length (bp)\t[7]Number of markers\t[8]Quality\n"
        + "".join(f"RG\t{SAMPLE}\t{c}\t{s}\t{e}\t{n}\t{n // 1400}\t{rng.uniform(40, 90):.1f}\n" for _, c, s, e, n in segs))

    mito = []
    het = {11, 23}
    for i, pos in enumerate(sorted(rng.sample(range(60, 16500), 41))):
        ref = "ACGT"[rng.randrange(4)]
        alt = "ACGT"[("ACGT".index(ref) + 2) % 4]
        if i in het:
            af = rng.uniform(0.11, 0.38)
        else:
            af = rng.uniform(0.985, 1.0) if i < 37 else rng.uniform(0.01, 0.03)
        filt = "PASS" if i < 37 else "weak_evidence"
        mito.append(f"chrM\t{pos}\t.\t{ref}\t{alt}\t.\t{filt}\t.\tGT:AF:DP\t0/1:{af:.3f}:{rng.randint(2400, 3300)}\n")
    put(os.path.join(d, "mito", f"{SAMPLE}_chrM_filtered.vcf.gz"),
        "##fileformat=VCFv4.2\n##contig=<ID=chrM,length=16569>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + SAMPLE + "\n" + "".join(mito), gz=True)
    put(os.path.join(d, "mito", f"{SAMPLE}_haplogroup.txt"),
        '"SampleID"\t"Haplogroup"\t"Rank"\t"Quality"\t"Range"\n'
        f'"{SAMPLE}"\t"H1c3"\t"1"\t"0.9312"\t"1-16569;"\n')
    # Read once wave 11's report has the haplocheck line and the Y card.
    put(os.path.join(d, "mito", f"{SAMPLE}_haplocheck.txt"),
        '"Sample"\t"Contamination Status"\t"Contamination Level"\t"Distance"\t"Sample Coverage"\n'
        f'"{SAMPLE}"\t"NO"\t"ND"\t"0"\t"2861"\n')
    put(os.path.join(d, "y_haplogroup", f"{SAMPLE}_y_haplogroup.txt"),
        "Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\tQC-score\tQC-1\tQC-2\tQC-3\n"
        f"{SAMPLE}_sorted\tI-M253\tM253\t4021876\t1187\t0.982\t1.0\t0.991\t0.982\n")


# --- steps 17, 23, 31 ----------------------------------------------------------------

def clinical(d, rng):
    def gene():
        return f"G{rng.randint(1, 900):04d}"

    impacts = ["HIGH"] * 41 + ["MODERATE"] * 233 + ["LOW"] * 44
    rows = []
    for imp in impacts:
        c, length = AUTOSOMES[rng.randrange(len(AUTOSOMES))]
        rows.append(f"{c}\t{rng.randint(1, length)}\tC\tT\t0/1\t{imp}\t{gene()}\tmissense_variant\t0.0004\t24.1\t0.41\t.\n")
    put(os.path.join(d, "clinical", f"{SAMPLE}_clinical_summary.tsv"),
        "CHROM\tPOS\tREF\tALT\tGT\tIMPACT\tGENE\tConsequence\tMAX_AF\tCADD_PHRED\tREVEL\tAM_CLASS\n" + "".join(rows))

    sl = [f"chr{rng.randint(1, 22)}\t{rng.randint(1, 90_000_000)}\tA\tG\t{rng.choice(['HIGH', 'MODERATE'])}\t{gene()}\n"
          for _ in range(57)]
    put(os.path.join(d, "slivar", f"{SAMPLE}_slivar_summary.tsv"), "CHROM\tPOS\tREF\tALT\tIMPACT\tGENE\n" + "".join(sl))
    ch = []
    for _ in range(3):
        g = gene()
        ch += [f"chr{rng.randint(1, 22)}\t{rng.randint(1, 90_000_000)}\tA\tG\tMODERATE\t{g}\n" for _ in range(2)]
    put(os.path.join(d, "slivar", f"{SAMPLE}_compound_hets.tsv"), "CHROM\tPOS\tREF\tALT\tIMPACT\tGENE\n" + "".join(ch))

    counts = {"Benign": 289, "Likely_Benign": 112, "VUS": 38}
    cls = "".join(f"{SAMPLE}\tchr{rng.randint(1, 22)}:g.{rng.randint(1, 90_000_000)}C>T\t{k}\t{k}\n"
                  for k, n in counts.items() for _ in range(n))
    put(os.path.join(d, "cpsr", f"{SAMPLE}.cpsr.grch38.classification.tsv.gz"),
        "SAMPLE_ID\tGENOMIC_CHANGE\tCPSR_CLASSIFICATION\tCLASSIFICATION\n" + cls, gz=True)
    put(os.path.join(d, "cpsr", f"{SAMPLE}.cpsr.grch38.html"),
        f"<!DOCTYPE html><html><head><title>CPSR {SAMPLE}</title></head><body>Invented sample.</body></html>\n")


def run_status(d, started):
    """logs/run_status.tsv as run-all.sh writes it, every step ok."""
    lines = ["# run-all.sh: when this run started and how each step ended",
             f"meta\tstarted_epoch\t{int(started)}",
             f"meta\tstarted_utc\t{time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime(started))}",
             f"meta\tdeclared_sex\t{SEX}"]
    lines += [f"step\t{s}\tok" for s in sorted(collect_summary.RUN_ALL_STEPS)]
    put(os.path.join(d, "logs", "run_status.tsv"), "\n".join(lines) + "\n")


# --- the two commands ----------------------------------------------------------------

def make(a):
    # A fixed start, so the report's "Latest run" line is the same on every
    # retake; the files written now are newer, so none is marked stale.
    started = RUN_STARTED
    d = os.path.join(a.genome_dir, SAMPLE)
    if os.path.exists(d):
        raise SystemExit(f"{d} exists: give an empty --genome-dir")
    os.makedirs(d)
    run_status(d, started)
    # Each part gets its own generator, so changing one leaves the others' values alone.
    parts = [("variants", lambda r: variants(d, r, a.records)), ("clinvar", lambda r: clinvar(d, a.genome_dir)),
             ("pgx", lambda r: pgx(d, a.genome_dir, r)), ("prs", lambda r: prs(d, a.genome_dir, r)),
             ("coverage", lambda r: coverage(d, r)), ("sample_qc", lambda r: sample_qc(d)),
             ("reads", lambda r: reads(d, r)), ("structural", lambda r: structural(d, r)),
             ("repeats", lambda r: repeats_and_more(d, r)), ("clinical", lambda r: clinical(d, r))]
    for name, fn in parts:
        fn(random.Random(f"{a.seed}:{name}"))
        print(f"  {name}: written")
    print(f"{SAMPLE} written in {d}")
    return 0


def check(a):
    """Fail when a section of step 24's summary is not ok."""
    path = os.path.join(a.genome_dir, SAMPLE, "summary.json")
    with open(path) as f:
        s = json.load(f)
    bad = [(k, v["state"], v.get("note") or "") for k, v in s["sections"].items() if v["state"] != "ok"]
    for k, v in s["sections"].items():
        print(f"  {k:<14} {v['state']}")
    if s["cyp2d6"].get("agree") is not True:
        bad.append(("cyp2d6", "callers do not agree", json.dumps(s["cyp2d6"]["calls"])))
    if bad:
        for k, st, note in bad:
            print(f"FAIL: {k} is {st} {note}".rstrip(), file=sys.stderr)
        print("A card would say 'Not run' or worse: feed that section in tests/demo/make_demo_sample.py.", file=sys.stderr)
        return 1
    print(f"Every section of {path} is ok.")
    return 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ["check"]:
        ap = argparse.ArgumentParser(prog="make_demo_sample.py check")
        ap.add_argument("--genome-dir", required=True)
        return check(ap.parse_args(argv[1:]))
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--genome-dir", required=True, help="an empty folder; DEMO-001 is written inside it")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--records", type=int, default=4_700_000, help="variant records in the VCF (default 4700000)")
    return make(ap.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())
