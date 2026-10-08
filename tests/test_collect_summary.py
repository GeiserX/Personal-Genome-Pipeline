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
     'insufficient markers') reach both reports.

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
        "CYP2D6\t*1/*4\tIntermediate Metabolizer\tnon-normal\nCYP2C19\t*1/*1\tNormal Metabolizer\tnormal\n")
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
        check("ClinVar: 3 hits, best-reviewed first", sec["clinvar"]["data"]["count"] == 3
              and [h["gene"] for h in sec["clinvar"]["data"]["hits"]] == ["GENEB", "GENEA", "GENEC"],
              sec["clinvar"]["data"])
        check("ClinVar: count by stars", sec["clinvar"]["data"]["by_stars"] == {"4": 0, "3": 1, "2": 0, "1": 1, "0": 1},
              sec["clinvar"]["data"]["by_stars"])
        check("text report: the ClinVar count line", "  Pathogenic/Likely Pathogenic hits: 3\n" in txt)
        check("HTML report: the same ClinVar count in its badge", 'class="badge badge-green">3</span>' in html)
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
        check("CYP2D6: three callers agree (allele order ignored)", summ["cyp2d6"]["agree"] is True, summ["cyp2d6"])

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
    finally:
        shutil.rmtree(work)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
