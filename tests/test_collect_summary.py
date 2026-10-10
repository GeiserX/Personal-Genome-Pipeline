#!/usr/bin/env python3
"""bin/collect_summary.py and bin/render_report.py on synthetic sample folders.

Checks:
  1. a bash output folder: the summary validates against
     tests/schema/summary.schema.json, and the text and HTML reports show the
     same ClinVar count, PRS rows and QC values;
  2. a Nextflow output folder (roh/, coverage/, hla/, pharmcat/): ROH, depth and
     HLA are read (the old reports looked for ROH under vcf/ only);
  3. heteroplasmy counts only chrM calls with an allele fraction of 0.05 to 0.95;
  4. the CPSR breakdown is read from classification.tsv.gz by column name;
  5. a result older than the latest run-all.sh run whose step was skipped is
     marked stale in the summary and in both reports; one whose step was ok, or
     that is newer than the run, is not;
  6. bin/clinvar_hits.awk and collect_summary.py give every review status the
     same number of stars;
  7. a corrupt file (not gzip, or gzip with damaged data) makes its section
     'unreadable' instead of stopping the report;
  8. the secondary-findings tier lists the ClinVar hits and the rare
     HIGH-impact clinical records in ACMG SF v3.3 genes, and nothing else;
  9. haplocheck's contamination status and Yleaf's Y haplogroup (or
     'insufficient markers') reach both reports;
 10. the step 33 verdict with the declared sex as somalier's pedigree:
     somalier keeps the pedigree's sex when the reads cannot tell, which is
     read as unknown (not checked), not as a call that agrees; a sex somalier
     changed from the pedigree's is a mismatch; a call that agrees with
     enough chrX sites is ok; a pedigree female with chrY reads but chrX that
     cannot tell is unknown, not an aneuploidy; without a pedigree the rows
     read as before;
 11. CYP2D6 by step 36's rule: pypgx 'Indeterminate' is not a call and is
     not counted in 'Genes called'; a Cyrius genotype with
     Filter=CYP2D6_depth_unreliable is shown as not usable and is not a
     call; step 36's verdict is shown when its table is there, and ignored
     when the table is older than a caller's result;
 12. repeat expansions: Stranger's STR_STATUS is read, every locus that is
     not normal is listed (also outside the five key loci), a record
     without a status is counted apart, not flagged; a Stranger file older
     than ExpansionHunter's is not read; a stale card says Stale, not
     Complete;
 13. HLA: each allele keeps its own T1K quality; a one-allele row
     ('.', 0, -1) is not low confidence; a quality of 0 or below is, and on
     HLA-A or HLA-B the reports say the gene is not passed to PharmCAT;
 14. a section whose step failed in the latest run says Failed, not Not
     run, in the summary and both reports;
 15. the CPIC card counts the genes called without a function phenotype.

Run: python3 tests/test_collect_summary.py
"""
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.join(REPO, "tests", "schema"))
import collect_summary  # noqa: E402
import render_report  # noqa: E402
import validate  # noqa: E402

FAILS = 0
DAY = 86400


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:400]}")
    if not ok:
        FAILS += 1


def put(path, text, gz=False, age_days=0):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if gz:
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        with open(path, "w") as f:
            f.write(text)
    if age_days:
        t = time.time() - age_days * DAY
        os.utime(path, (t, t))


VCF_HEAD = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
HITS = VCF_HEAD + (
    "chr1\t10\t1\tC\tT\t50\tPASS\tGENEINFO=GENEA:11;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter\tGT\t0/1\n"
    "chr1\t20\t2\tT\tG\t50\tPASS\tGENEINFO=GENEB:22;CLNSIG=Likely_pathogenic;CLNREVSTAT=reviewed_by_expert_panel\tGT\t1/1\n"
    "chr2\t5\t3\tA\tG\t50\tPASS\tGENEINFO=GENEC:33;CLNSIG=Pathogenic;CLNREVSTAT=no_assertion_criteria_provided\tGT\t0|1\n")
MITO = VCF_HEAD.replace("##fileformat=VCFv4.2\n", "") + "".join(
    f"chrM\t{i}\t.\tA\tG\t.\tPASS\t.\tGT:AF\t0/1:{af}\n" for i, af in enumerate(("0.02", "0.30", "0.99", "0.05", "0.949"), 1))


def bash_folder(d, s="S"):
    put(f"{d}/vcf/{s}.vcf.gz", VCF_HEAD + "chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
        "chr1\t2\t.\tAT\tA\t50\tPASS\t.\tGT\t0/1\nchr1\t3\t.\tC\tT\t5\tRefCall\t.\tGT\t0/0\n", gz=True)
    put(f"{d}/clinvar/{s}_clinvar_hits.vcf", HITS)
    put(f"{d}/mito/{s}_chrM_filtered.vcf.gz", MITO, gz=True)
    put(f"{d}/vcf/{s}_roh.txt", "# RG\tsample\tchrom\tstart\tend\tlength\n"
        "RG\tS\tchr1\t100\t7000100\t7000000\t3\t50\nRG\tS\tchrX\t1\t9000000\t9000000\t3\t50\n")
    put(f"{d}/mosdepth/{s}.mosdepth.summary.txt",
        "chrom\tlength\tbases\tmean\tmin\tmax\nchr1\t10\t300\t30.0\t0\t50\ntotal\t10\t300\t29.5\t0\t50\n")
    put(f"{d}/indexcov/indexcov-indexcov.ped", "#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\t"
        "CNchrX\tCNchrY\nS\tS\t-9\t-9\t2\t-9\t2.0\t0.0\n")
    put(f"{d}/prs/{s}_prs_summary.tsv", "Condition\tPGS_ID\tScore_SUM\tVariants_Matched\tVariants_Total\n"
        "Coronary artery disease\tPGS000018\t0.12\t100\t200\nType 2 diabetes\tPGS000014\t-0.4\t50\t90\n")
    put(f"{d}/cpsr/{s}.cpsr.grch38.classification.tsv.gz",
        "SAMPLE_ID\tGENOMIC_CHANGE\tCPSR_CLASSIFICATION\tCLASSIFICATION\n"
        "S\ta\tVUS\tVUS\nS\tb\tVUS\tLikely_Pathogenic\nS\tc\tB\tVUS\n", gz=True)
    put(f"{d}/slivar/{s}_slivar_summary.tsv", "CHROM\tPOS\tREF\tALT\tIMPACT\tSYMBOL\nchr1\t1\tA\tG\tHIGH\tX\n")
    put(f"{d}/hla_t1k/{s}_hla_genotype.tsv", "HLA-A\t2\tA*01:01:01\t30\t60\tA*02:01:01\t20\t55\t\n")
    put(f"{d}/cpic/{s}_phenotypes.tsv", "Gene\tDiplotype\tPhenotype\tStatus\n"
        "CYP2D6\t*1/*4\tIntermediate Metabolizer\tnon-normal\nCYP2C19\t*1/*1\tNormal Metabolizer\tnormal\n"
        "VKORC1\trs9923231 reference (C)/rs9923231 reference (C)\t-1639 GG\tunclassified\n")
    put(f"{d}/cpic/{s}_cpic_recommendations.txt", "Pharmacogenomic Drug Recommendations\n")
    put(f"{d}/pypgx/{s}_pypgx_summary.tsv", "Gene\tDiplotype\tPhenotype\tCNV_call\tSource\nCYP2D6\t*4/*1\tIM\t.\tbam\n")
    put(f"{d}/cyrius/{s}_cyp2d6.tsv", "Sample\tGenotype\tFilter\nS\t*1/*4\tPASS\n")


def render(d, s="S"):
    summ = collect_summary.collect(s, d)
    return summ, render_report.text_report(summ), render_report.html_report(summ)


def main():
    work = tempfile.mkdtemp()
    schema = json.load(open(os.path.join(REPO, "tests", "schema", "summary.schema.json")))
    try:
        # 1. bash folder: schema, and both reports agree
        d = os.path.join(work, "bash", "S")
        bash_folder(d)
        summ, txt, html = render(d)
        errs = validate.validate(schema, json.loads(json.dumps(summ)))
        check("bash folder: the summary validates against the schema", not errs, errs)
        sec = summ["sections"]
        check("CPIC: an unclassified gene (VKORC1 -1639 GG) is not counted as non-normal",
              sec["cpic"]["data"]["non_normal"] == 1, sec["cpic"]["data"])
        check("ClinVar: 3 hits, best-reviewed first", sec["clinvar"]["data"]["count"] == 3
              and [h["gene"] for h in sec["clinvar"]["data"]["hits"]] == ["GENEB", "GENEA", "GENEC"],
              sec["clinvar"]["data"])
        check("ClinVar: count by stars", sec["clinvar"]["data"]["by_stars"] == {"4": 0, "3": 1, "2": 0, "1": 1, "0": 1},
              sec["clinvar"]["data"]["by_stars"])
        check("text report: the ClinVar count line", "  Pathogenic/Likely Pathogenic hits: 3\n" in txt)
        check("HTML report: the same ClinVar count in its badge, yellow (never green) with a hit",
              'class="badge badge-yellow">3</span>' in html and 'badge-green">3<' not in html)
        rows = [l for l in html.splitlines() if l.strip().startswith("<tr><td>chr")]
        check("HTML report: one row per hit, each with a Stars cell",
              len(rows) == 3 and rows[0].strip().endswith("<td>3</td></tr>"), rows)
        for r in sec["prs"]["data"]["scores"]:
            check(f"PRS row {r['pgs_id']} in both reports", r["pgs_id"] in txt and f"<td>{r['pgs_id']}</td>" in html)
        check("QC: mean depth 29.5x in both reports", "Mean depth: 29.5x" in txt and "29.5x" in html,
              sec["coverage"]["data"])
        check("QC: inferred sex female", sec["sex_check"]["data"]["inferred_sex"] == "female")
        check("variants: 3 records, 2 PASS, 2 SNP lines, 1 indel",
              sec["variants"]["data"] == {"total": 3, "pass": 2, "snps": 2, "indels": 1}, sec["variants"]["data"])
        check("CYP2D6: pypgx and Cyrius agree (allele order ignored)", summ["cyp2d6"]["agree"] is True, summ["cyp2d6"])
        # 15. CPIC unclassified
        check("CPIC: the unclassified gene is counted apart (VKORC1 -1639 GG)",
              sec["cpic"]["data"].get("unclassified") == 1, sec["cpic"]["data"])
        check("CPIC: both reports show the unclassified count",
              "Called without a function phenotype (no drug guidance): 1" in txt
              and 'Genes called without a function phenotype (no drug guidance)</span><span class="value">1<' in html,
              [l for l in txt.splitlines() if "Not called" in l])

        # 3. heteroplasmy floor
        check("heteroplasmy: only AF 0.05 to 0.95 counts (0.30, 0.05, 0.949 of 5 PASS)",
              sec["mito"]["data"]["heteroplasmic"] == 3 and sec["mito"]["data"]["pass"] == 5, sec["mito"]["data"])
        # 4. CPSR by column name
        check("CPSR: breakdown read from the CLASSIFICATION column by name",
              sec["cpsr"]["data"].get("classification") == {"VUS": 2, "Likely_Pathogenic": 1}, sec["cpsr"]["data"])

        # 2. Nextflow folder
        n = os.path.join(work, "nf", "S")
        put(f"{n}/roh/S_roh.txt", "RG\tS\tchr2\t1\t6000001\t6000000\t3\t50\n")
        put(f"{n}/roh/S_roh_summary.txt", "summary\n")
        put(f"{n}/coverage/S.mosdepth.summary.txt", "chrom\tlength\tbases\tmean\tmin\tmax\ntotal\t1\t1\t31.25\t0\t9\n")
        put(f"{n}/hla/S_hla_genotype.tsv", "HLA-B\t2\tB*07:02:01\t30\t60\tB*08:01:01\t20\t55\t\n")
        put(f"{n}/pharmcat/S.report.json", json.dumps({"pharmcatVersion": "3.2.0", "genes": {}}))
        summ_n, txt_n, _ = render(n)
        s2 = summ_n["sections"]
        check("Nextflow folder: ROH read from roh/", s2["roh"]["state"] == "ok"
              and s2["roh"]["data"]["autosomal_over_5mb"] == [{"region": "chr2:1-6000001", "mb": 6.0}], s2["roh"])
        check("Nextflow folder: ROH shown in the text report", "Autosomal ROH > 5MB: 1" in txt_n)
        check("Nextflow folder: depth from coverage/", s2["coverage"]["data"]["mean_depth"] == 31.25)
        check("Nextflow folder: HLA from hla/", s2["hla"]["data"]["loci"][0]["alleles"] == ["B*07:02:01", "B*08:01:01"])
        check("Nextflow folder: PharmCAT from pharmcat/", s2["pharmcat"]["data"].get("version") == "3.2.0")
        check("Nextflow folder: the summary validates", not validate.validate(schema, json.loads(json.dumps(summ_n))))

        # 5. stale results
        st = os.path.join(work, "stale", "S")
        bash_folder(st)
        put(f"{st}/slivar/S_slivar_summary.tsv", "CHROM\tPOS\tREF\tALT\tIMPACT\tSYMBOL\nchr1\t1\tA\tG\tHIGH\tX\n", age_days=30)
        put(f"{st}/vcf/S_roh.txt", "RG\tS\tchr1\t1\t6000001\t6000000\t3\t50\n", age_days=30)
        put(f"{st}/logs/run_status.tsv", f"meta\tstarted_epoch\t{time.time() - DAY}\nmeta\tdeclared_sex\tmale\n"
            "step\t31\tskipped (needs VEP, step 13 skipped)\nstep\t11\tok\n")
        summ_s, txt_s, html_s = render(st)
        ss = summ_s["sections"]
        check("stale: step 31 skipped, slivar summary from an earlier run -> stale",
              ss["slivar"]["state"] == "stale" and "skipped" in (ss["slivar"]["note"] or ""), ss["slivar"])
        check("stale: marked in the text report", "[STALE: from an earlier run" in txt_s)
        check("stale: marked in the HTML report", '<div class="stale">STALE: from an earlier run' in html_s)
        check("stale: step 11 ok in that run, an old ROH file is current", ss["roh"]["state"] == "ok", ss["roh"])
        check("stale: a file written after the run started is current", ss["clinvar"]["state"] == "ok")
        check("stale: declared sex from the run status, mismatch flagged",
              summ_s["run"]["declared_sex"] == "male" and "DOES NOT MATCH" in txt_s)
        check("no run status: nothing is stale",
              all(x["state"] != "stale" for x in summ["sections"].values()))

        # 6. stars: awk and python agree
        statuses = list(collect_summary.STARS) + ["no_assertion_criteria_provided", "no_classification_provided",
                                                  "something_new", ""]
        vcf = os.path.join(work, "stars.vcf")
        put(vcf, VCF_HEAD + "".join(f"chr1\t{i}\t.\tA\tG\t50\tPASS\tCLNREVSTAT={st_}\tGT\t0/1\n"
                                    for i, st_ in enumerate(statuses, 1)))
        out = subprocess.run(["awk", "-f", os.path.join(REPO, "bin", "clinvar_hits.awk"), vcf],
                             capture_output=True, text=True).stdout.split("\n")
        awk_stars = [int(l.split("\t")[0]) for l in out if l]
        py_stars = [collect_summary.STARS.get(x, 0) for x in statuses]
        check("bin/clinvar_hits.awk and collect_summary.py give the same stars", awk_stars == py_stars,
              f"awk {awk_stars} python {py_stars}")

        # 7. a corrupt file
        c = os.path.join(work, "corrupt", "S")
        put(f"{c}/vcf/S.vcf.gz", "not gzip\n")
        put(f"{c}/clinvar/S_clinvar_hits.vcf", HITS)
        summ_c, txt_c, _ = render(c)
        check("corrupt VCF: section unreadable, the report still renders",
              summ_c["sections"]["variants"]["state"] == "unreadable" and "hits: 3" in txt_c,
              summ_c["sections"]["variants"])
        # a gzip whose header is fine but whose compressed data is damaged
        # (a truncated copy, a bad disk) raises zlib.error, not OSError
        z = os.path.join(work, "corrupt-deflate", "S")
        os.makedirs(f"{z}/vcf")
        blob = bytearray(gzip.compress((VCF_HEAD + "chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n" * 2000).encode()))
        blob[20:40] = b"\xff" * 20
        with open(f"{z}/vcf/S.vcf.gz", "wb") as f:
            f.write(bytes(blob))
        put(f"{z}/clinvar/S_clinvar_hits.vcf", HITS)
        try:
            summ_z, txt_z, _ = render(z)
            state_z = summ_z["sections"]["variants"]
        except Exception as e:  # the pre-fix behaviour: the whole report dies
            state_z, txt_z = f"raised {type(e).__name__}: {e}", ""
        check("damaged gzip data: section unreadable, the report still renders",
              isinstance(state_z, dict) and state_z["state"] == "unreadable" and "hits: 3" in txt_z, state_z)

        # 8. ACMG SF tier: BRCA2 (on the list) and GENEA (not) in ClinVar; TTN
        # (on the list) and GENEX (not) with HIGH impact; MYH7 on the list
        # but MODERATE only
        a = os.path.join(work, "acmg", "S")
        put(f"{a}/clinvar/S_clinvar_hits.vcf", VCF_HEAD + (
            "chr13\t32340300\t1\tG\tA\t50\tPASS\tGENEINFO=BRCA2:675;CLNSIG=Pathogenic;"
            "CLNREVSTAT=reviewed_by_expert_panel\tGT\t0/1\n"
            "chr1\t10\t2\tC\tT\t50\tPASS\tGENEINFO=GENEA:11;CLNSIG=Pathogenic;"
            "CLNREVSTAT=criteria_provided,_single_submitter\tGT\t0/1\n"))
        csq = ('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
               'Format: Allele|Consequence|IMPACT|SYMBOL">\n')
        put(f"{a}/clinical/S_clinical.vcf.gz", "##fileformat=VCFv4.2\n" + csq + VCF_HEAD.split("\n", 1)[1] + (
            "chr2\t178500000\t.\tC\tT\t50\tPASS\tCSQ=T|stop_gained|HIGH|TTN,T|intron_variant|MODIFIER|TTN-AS1\tGT\t0/1\n"
            "chr3\t100\t.\tA\tG\t50\tPASS\tCSQ=G|frameshift_variant|HIGH|GENEX\tGT\t1/1\n"
            "chr14\t23400000\t.\tG\tA\t50\tPASS\tCSQ=A|missense_variant|MODERATE|MYH7\tGT\t0/1\n"), gz=True)
        summ_a, txt_a, html_a = render(a)
        sf = summ_a["sections"]["clinical"]["data"].get("acmg_sf") or {}
        check("ACMG SF: v3.3, 84 genes", sf.get("version") == "ACMG SF v3.3" and sf.get("genes_on_list") == 84, sf)
        check("ACMG SF: the ClinVar hit in BRCA2 only", [x["gene"] for x in sf.get("clinvar_hits", [])] == ["BRCA2"], sf)
        check("ACMG SF: the HIGH-impact TTN record only (not GENEX, not MODERATE MYH7)",
              [(x["gene"], x["consequence"]) for x in sf.get("high_impact", [])] == [("TTN", "stop_gained")], sf)
        check("ACMG SF: the text report lists BRCA2 and TTN",
              "1 ClinVar P/LP, 1 rare HIGH impact" in txt_a and "    BRCA2" in txt_a and "    TTN" in txt_a, txt_a)
        check("ACMG SF: the HTML report has the card", "Secondary-Findings Genes (ACMG SF v3.3)" in html_a
              and "<td>TTN</td>" in html_a and "<td>GENEX</td>" not in html_a)

        # 9. haplocheck and Yleaf
        h = os.path.join(work, "mito", "S")
        put(f"{h}/mito/S_haplogroup.txt", '"SampleID"\t"Haplogroup"\t"Rank"\n"S"\t"H1a"\t"0.95"\n')
        put(f"{h}/mito/S_haplocheck.txt", '"Sample"\t"Contamination Status"\t"Contamination Level"\t"Distance"\n'
            '"S"\t"YES"\t"0.12"\t"5"\n')
        put(f"{h}/y_haplogroup/S_y_haplogroup.txt", "Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\t"
            "QC-score\tQC-1\tQC-2\tQC-3\nS_sorted\tR-M269\tM269\t100\t42\t0.97\t1\t1\t0.97\n")
        summ_h, txt_h, html_h = render(h)
        hd = summ_h["sections"]["haplogroup"]["data"]
        check("haplocheck: status and level read", hd.get("contamination_status") == "YES"
              and hd.get("contamination_level") == "0.12", hd)
        check("haplocheck: the line in both reports", "Contamination (haplocheck): YES, two mtDNA haplogroups" in txt_h
              and "YES, two mtDNA haplogroups" in html_h, txt_h)
        check("Yleaf: the Y haplogroup in both reports", "R-M269 (42 markers, QC-score 0.97)" in txt_h
              and "Y-Chromosome Haplogroup" in html_h and "R-M269" in html_h, txt_h)
        put(f"{h}/y_haplogroup/S_y_haplogroup.txt", "Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\t"
            "QC-score\tQC-1\tQC-2\tQC-3\nS_sorted\tNA\t\t100\t3\t0\t0\t0\t0\n")
        summ_n, txt_n, _ = render(h)
        check("Yleaf: Hg NA reads as insufficient markers",
              summ_n["sections"]["y_haplogroup"]["data"].get("haplogroup") == "insufficient markers"
              and "Haplogroup: insufficient markers (3 markers" in txt_n, summ_n["sections"]["y_haplogroup"])
        os.remove(f"{h}/mito/S_haplocheck.txt")
        check("no haplocheck file: the report says it was not checked",
              "Contamination (haplocheck): not checked" in render(h)[1])

        # 10. step 33 with the declared sex as somalier's pedigree
        head = ("#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\toriginal_pedigree_sex\t"
                "gt_depth_mean\tn_hom_ref\tn_het\tn_hom_alt\tp_middling_ab\tX_depth_mean\tX_n\tX_hom_ref\t"
                "X_het\tX_hom_alt\tY_depth_mean\tY_n\n")

        def qc(sex, ped, x_n, x_het, x_hom_alt, declared, y_depth="0.0", mid="0.010"):
            path = os.path.join(work, "somalier.samples.tsv")
            put(path, head + f"S\tS\t-9\t-9\t{sex}\t-9\t{ped}\t30.0\t60\t40\t20\t{mid}\t15.0\t{x_n}\t0\t"
                f"{x_het}\t{x_hom_alt}\t{y_depth}\t5\n")
            return dict(collect_summary.sample_qc_table("S", path, declared_sex=declared))

        t = qc("1", "male", 2, 0, 2, "male")
        check("pedigree male kept on 2 chrX sites: unknown, not checked (not ok)",
              t["inferred_sex"] == "unknown" and t["sex_check"] == "not_checked", t)
        t = qc("1", "female", 39, 0, 39, "female")
        check("pedigree female, somalier set male from 39 chrX sites: mismatch",
              t["inferred_sex"] == "male" and t["sex_check"] == "mismatch", t)
        t = qc("1", "male", 39, 0, 39, "male")
        check("pedigree male, 39 homozygous chrX sites agree: ok", t["sex_check"] == "ok", t)
        t = qc("1", "male", 39, 0, 39, "male", mid="0.080")
        check("pedigree male kept on a sample somalier calls low quality (8% middling allele balance): not checked",
              t["sex_check"] == "not_checked", t)
        t = qc("1", "male", 39, 0, 39, "male", mid="nan")
        check("pedigree male kept, middling allele balance nan (no autosomal call): not checked",
              t["sex_check"] == "not_checked", t)
        t = qc("-2", "female", 4, 0, 4, "female", y_depth="14.0")
        check("pedigree female with chrY reads, chrX cannot tell: unknown, not the aneuploidy line",
              t["sex_check"] == "not_checked" and "could not tell" in t["sex_check_reason"], t)
        t = qc("-2", "female", 40, 20, 20, "female", y_depth="14.0")
        check("pedigree female, chrX heterozygous and chrY reads: the aneuploidy line",
              t["sex_check"] == "not_checked" and "chrY has reads" in t["sex_check_reason"], t)
        t = qc("1", "-9", 2, 0, 2, "male")
        check("no pedigree (-9): somalier's sex column is read as before",
              t["inferred_sex"] == "male" and t["sex_check"] == "ok", t)

        # 11. CYP2D6 with a failed depth check: pypgx Indeterminate, Cyrius's genotype filtered
        y = os.path.join(work, "cyp", "S")
        put(f"{y}/pypgx/S_pypgx_summary.tsv", "Gene\tDiplotype\tPhenotype\tCNV_call\tSource\n"
            "CYP2D6\tIndeterminate\tIndeterminate\t.\tbam\nCYP2C19\t*1/*2\tIM\t.\tbam\nCYP2C9\tindeterminate\t.\t.\tbam\n",
            age_days=2)
        put(f"{y}/cyrius/S_cyp2d6.tsv", "Sample\tGenotype\tFilter\nS\t*5/*5\tCYP2D6_depth_unreliable\n", age_days=2)
        summ_y, txt_y, html_y = render(y)
        cy = summ_y["cyp2d6"]
        check("pypgx: Indeterminate (any case) is not counted as called",
              summ_y["sections"]["pypgx"]["data"]["genes_called"] == 1
              and "Genes called: 1/3" in txt_y, summ_y["sections"]["pypgx"]["data"])
        check("CYP2D6: a Cyrius genotype that failed its Filter is shown as not usable",
              cy["calls"]["Cyrius"] == "*5/*5 (not usable: Filter CYP2D6_depth_unreliable)", cy["calls"])
        check("CYP2D6: neither call is usable, so no agreement is claimed (not 'disagree')", cy["agree"] is None, cy)
        check("CYP2D6: the text report says so", "not usable: Filter CYP2D6_depth_unreliable" in txt_y
              and "pypgx and Cyrius agree: fewer than two usable calls" in txt_y, txt_y)
        check("CYP2D6: the HTML card says so", "not usable: Filter CYP2D6_depth_unreliable" in html_y
              and "fewer than two usable calls" in html_y)
        check("CYP2D6: no step 36 table, no verdict", cy.get("consensus") is None and "Step 36" not in txt_y, cy)
        put(f"{y}/pgx_consensus/S_pgx_consensus.tsv", "Gene\tResult\tOutside_call\tReason\tEvidence\n"
            "HLA-A\tnot typed\tno\tHLA typing (step 08) did not run\t-\n"
            "CYP2D6\tindeterminate\tno\tthe CYP2D6 depth check failed\tpypgx: no call (Indeterminate)\n")
        summ_y, txt_y, html_y = render(y)
        v36 = summ_y["cyp2d6"].get("consensus") or {}
        check("CYP2D6: step 36's verdict from its table",
              v36.get("result") == "indeterminate" and v36.get("passed_to_pharmcat") is False, summ_y["cyp2d6"])
        check("CYP2D6: both reports print step 36's verdict and reason",
              "Step 36: indeterminate, not passed to PharmCAT (the CYP2D6 depth check failed)" in txt_y
              and "Not passed to PharmCAT: the CYP2D6 depth check failed" in html_y, txt_y)
        check("CYP2D6 summary validates", not validate.validate(schema, json.loads(json.dumps(summ_y))),
              validate.validate(schema, json.loads(json.dumps(summ_y))))
        put(f"{y}/pgx_consensus/S_pgx_consensus.tsv", open(f"{y}/pgx_consensus/S_pgx_consensus.tsv").read(), age_days=5)
        check("CYP2D6: a step 36 table older than the callers' results is not read",
              render(y)[0]["cyp2d6"].get("consensus") is None)

        # 12. repeat expansions with Stranger
        x = os.path.join(work, "str", "S")
        eh_rec = ("chr4\t3074877\t.\tC\t<STR20>\t.\tPASS\tEND=1;REPID=HTT;VARID=HTT{st}\tGT:REPCN\t0/1:17/20\n"
                  "chr14\t92071011\t.\tG\t<STR75>\t.\tPASS\tEND=1;REPID=ATXN3;VARID=ATXN3{st2}\tGT:REPCN\t0/1:20/75\n"
                  "chr9\t100\t.\tA\t<STR9>\t.\tPASS\tEND=1;REPID=NOTINCAT;VARID=NOTINCAT\tGT:REPCN\t0/1:5/9\n")
        put(f"{x}/expansion_hunter/S_eh.vcf", VCF_HEAD + eh_rec.format(st="", st2=""), age_days=1)
        put(f"{x}/expansion_hunter/S_eh_stranger.vcf",
            VCF_HEAD + eh_rec.format(st=";STR_STATUS=normal", st2=";STR_STATUS=full_mutation"))
        summ_x, txt_x, html_x = render(x)
        ed = summ_x["sections"]["expansions"]["data"]
        check("Stranger: read when it is there", summ_x["sections"]["expansions"]["source"]
              == "expansion_hunter/S_eh_stranger.vcf" and ed.get("stranger") is True, summ_x["sections"]["expansions"])
        check("Stranger: the full_mutation outside the five key loci is listed, the normal one is not",
              ed.get("flagged") == [{"locus": "ATXN3", "repeat_count": "20/75", "status": "full_mutation"}], ed)
        check("Stranger: a record without STR_STATUS is counted apart, not flagged", ed.get("no_status") == 1, ed)
        check("Stranger: the key loci stay", ed["key_loci"][0] == {"locus": "HTT", "repeat_count": "17/20"}, ed)
        check("Stranger: both reports list ATXN3 full_mutation and the short-read caveat",
              "ATXN3    20/75  full_mutation" in txt_x and "<td>ATXN3</td><td>20/75</td><td>full_mutation</td>" in html_x
              and "can be wrong at some loci" in txt_x and "can be wrong at some loci" in html_x, txt_x)
        check("Stranger summary validates", not validate.validate(schema, json.loads(json.dumps(summ_x))))
        put(f"{x}/expansion_hunter/S_eh_stranger.vcf", open(f"{x}/expansion_hunter/S_eh_stranger.vcf").read(), age_days=3)
        summ_x2, txt_x2, _ = render(x)
        check("Stranger: a file older than ExpansionHunter's is not read, and the report says Stranger did not run",
              summ_x2["sections"]["expansions"]["source"] == "expansion_hunter/S_eh.vcf"
              and not summ_x2["sections"]["expansions"]["data"].get("stranger")
              and "Stranger (step 9b) did not run" in txt_x2, summ_x2["sections"]["expansions"])
        put(f"{x}/logs/run_status.tsv", f"meta\tstarted_epoch\t{time.time() - 0.5 * DAY}\nstep\t09\tskipped\n")
        _, _, html_x3 = render(x)
        card_eh = html_x3.split("<h2>Repeat Expansions</h2>", 1)[-1].split("<h2>", 1)[0]
        check("Repeat card of a stale result says Stale, not Complete",
              "Stale" in card_eh and "Complete" not in card_eh, card_eh)

        # 13. HLA quality per allele
        q = os.path.join(work, "hla", "S")
        put(f"{q}/hla_t1k/S_hla_genotype.tsv",
            "HLA-A\t1\tA*01:01:01\t30\t60\t.\t0\t-1\n"
            "HLA-B\t2\tB*57:01:01\t10\t0\tB*08:01:01\t20\t40\n"
            "HLA-C\t2\tC*07:01:01\t5\t-1\tC*07:02:01\t5\t30\n")
        summ_q, txt_q, html_q = render(q)
        loci = {l["gene"]: l for l in summ_q["sections"]["hla"]["data"]["loci"]}
        check("HLA: a one-allele row keeps the allele's quality only and is not low confidence",
              loci["HLA-A"]["alleles"] == ["A*01:01:01"] and loci["HLA-A"]["quality"] == ["60"]
              and loci["HLA-A"].get("low_confidence") is False, loci["HLA-A"])
        check("HLA: HLA-B with a quality 0 allele is low confidence and withheld from PharmCAT",
              loci["HLA-B"]["quality"] == ["0", "40"] and loci["HLA-B"].get("low_confidence") is True
              and loci["HLA-B"].get("withheld_from_pharmcat") is True, loci["HLA-B"])
        check("HLA: HLA-C with a quality -1 allele is low confidence, but PharmCAT never gets HLA-C",
              loci["HLA-C"].get("low_confidence") is True and loci["HLA-C"].get("withheld_from_pharmcat") is False,
              loci["HLA-C"])
        check("HLA: the text report shows each quality and says only HLA-B is not passed to PharmCAT",
              "A*01:01:01 (quality 60)\n" in txt_q and "B*57:01:01 (quality 0) / B*08:01:01 (quality 40); low confidence"
              in txt_q and txt_q.count("low confidence (T1K quality 0 or below), not passed to PharmCAT") == 1
              and "C*07:02:01 (quality 30); low confidence (T1K quality 0 or below)\n" in txt_q, txt_q)
        check("HLA: the HTML card shows the same, and the short-read note",
              "B*57:01:01 (quality 0)" in html_q and html_q.count("not passed to PharmCAT</span>") == 1
              and "HLA typing from short-read WGS is approximate" in html_q)
        check("HLA summary validates", not validate.validate(schema, json.loads(json.dumps(summ_q))))

        # 14. a step that failed in the latest run
        f_ = os.path.join(work, "failed", "S")
        put(f"{f_}/clinvar/S_clinvar_hits.vcf", HITS)
        put(f"{f_}/logs/run_status.tsv", f"meta\tstarted_epoch\t{time.time() - DAY}\nstep\t10\tfailed\nstep\t11\tok\n")
        summ_f, txt_f, html_f = render(f_)
        tel = summ_f["sections"]["telomere"]
        check("failed: no telomere file and step 10 failed -> state failed with a note",
              tel["state"] == "failed" and "step 10 failed" in (tel["note"] or ""), tel)
        check("failed: a step that never ran is still missing", summ_f["sections"]["roh"]["state"] == "missing")
        card_t = html_f.split("<h2>Telomere content (relative)</h2>", 1)[-1].split("<h2>", 1)[0]
        check("failed: the HTML card says Failed, not Not run", 'badge-red">Failed<' in card_t and "Not run" not in card_t,
              card_t)
        check("failed: the text report says FAILED and lists it under Steps Not Run",
              "[FAILED: step 10 failed" in txt_f and "  - Telomere content (TelomereHunter) (step 10 failed)" in txt_f, txt_f)
        check("failed summary validates", not validate.validate(schema, json.loads(json.dumps(summ_f))),
              validate.validate(schema, json.loads(json.dumps(summ_f))))
    finally:
        shutil.rmtree(work)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
