#!/usr/bin/env python3
"""The PRS path around pgsc_calc in bin/collect_summary.py, and its card in
bin/render_report.py, on synthetic files shaped like pgsc_calc's own.

Checks:
  1. prs-format writes a harmonised file as pgsc_calc's custom GRCh38 file:
     chr_name/chr_position from hm_chr/hm_pos without "chr", the other allele
     from other_allele or hm_inferOtherAllele (never "A/G"), rows without a
     GRCh38 position dropped, rows off the autosomes (chrX, chrY) dropped
     from both the score and the allele list, the label from
     assets/pgs_scores.tsv; and it
     refuses a GRCh37 file and a non-additive one; with --keep-x (the
     sample's sex is known) the chrX rows stay in both and chrY still goes;
     chrx_rows.tsv records per score how many chrX rows it has and whether
     they were scored or left out;
  2. prs-table without a panel: the sum from aggregated_scores.txt.gz as
     pgsc_calc wrote it but without ".0", matched and total counts from the
     match summary, Percentile and Ancestry_Group NA, ChrX from
     chrx_rows.tsv, Input last; a score pgsc_calc dropped (below its minimum
     overlap) has no sum;
  3. prs-table with a panel: the percentile and group of the target sample
     only (not the reference samples in the same file), and the ancestry table
     with the population and the principal components;
  4. --zero-matches: every score unmatched, no sum; --below-threshold: no
     sum, the matched count unknown and the rate from pgscatalog-match's log;
  5. both reports: a percentile with its group when there is one, the
     "Raw score only" line when there is none, and the not-assessed line about
     percentiles only then; a score whose chrX rows were left out is named,
     one whose chrX rows were scored is not; both captions, with a panel
     and without, carry the ancestry-transfer caveat of docs/25-prs.md;
  6. the ClinVar card: a homozygous 3-star P/LP hit gets no green badge,
     no hit gets a gray one, and both reports ask for every hit's population
     frequency (heterozygous ones too) and say it is unknown without step 13; with step 06's cDNA read check table, both reports print
     the flag beside the flagged hit and "clear" beside the other, the HTML
     hits table gets its column and the card the flagged count; without the
     table (a Nextflow folder) there is no such column;
  5b. the fine-grained ancestry section of step 26b in both reports, rendered
     from the same summary as the PRS card: the groups with their bootstrap
     ranges, a small regional group kept because its range is not small, the
     painted shares, the picture and the "not a nationality" sentence;
  7. the telomere card says relative telomere content with its unit, and
     no report in bin/ says "Telomere Length";
  8. the mosaic card: a man with an estimated loss of Y shows it in both
     reports with its 95% range, the depth ratio and the population-level
     note; a woman shows the allele-fraction line too; a sample without the
     estimate gets a card that says Not run, and no text section.
  9. the paralog card of step 35b is left out when the step did not run, and
     a result older than a run in which step 35b failed is marked Stale.
  10. the blood group card (step 41): every system's row in the HTML table and
     the text report, ABO and RhD on top, the "predicted, not typed" caption;
     the card is absent when step 41 did not run; and neither the reports nor
     bin/bloodtype_report.py say anything about transfusion or donation (the
     table is information).
  11. the traits card (step 42): both reports carry the not-a-diagnosis
      note; a traits table older than the latest run-all.sh run, whose step
      42 was skipped in it, is shown as stale in both reports, and one newer
      than the run is not.
 12. the raw-data export (step 40): both reports show how many chip sites
     got a genotype, the '--' reasons, the two files and the privacy note;
     without the step's table there is no card and no text section.
 13. the couple card (step 39), from a table bin/couple_report.py wrote:
     both reports name the partner, the shared gene with its chance for a
     child and both partners' findings, and point to a genetic counsellor;
     without the table there is no card and no text section.

Run: python3 tests/test_render_report.py
"""
import glob
import gzip
import json
import os
import re
import shutil
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.join(REPO, "tests", "schema"))
import collect_summary  # noqa: E402
import render_report  # noqa: E402
import validate  # noqa: E402

FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS += 1


def put(path, text, gz=False):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with (gzip.open(path, "wt") if gz else open(path, "w")) as f:
        f.write(text)


HM = ("#pgs_id=PGS000001\n#trait_reported=Something else\n#HmPOS_build=GRCh38\n"
      "rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\thm_chr\thm_pos\thm_inferOtherAllele\n"
      "rs1\t1\t100\tA\tG\t0.5\t1\t1100\t\n"
      "rs2\t1\t200\tc\t\t-0.25\tchr2\t2200\tT\n"
      "rs3\t1\t300\tG\t\t1\t3\t3300\tA/C\n"
      "rs4\t1\t400\tT\tC\t2\t\t\t\n"
      "rs5\tX\t500\tA\tG\t1\tX\t51000000\t\n"
      "rs6\tY\t600\tC\tT\t1\tchrY\t2800000\t\n")


def read_rows(path):
    with open(path) as f:
        lines = [l.rstrip("\n").split("\t") for l in f]
    return lines[0], lines[1:]


def pgsc_results(d, adjusted):
    """A pgsc_calc --outdir for sampleset 'sample': two scores, PGS000001 passed, PGS000002 was dropped.
    The reference samples' rows come after the target's, as a reader that
    forgot to skip them would take their values."""
    put(f"{d}/sample/match/sample_summary.csv",
        "dataset,accession,ambiguous,is_multiallelic,match_flipped,duplicate_best_match,duplicate_ID,match_IDs,match_status,count,score_pass,match_rate\n"
        "sample,PGS000001,false,false,false,false,false,NA,matched,2,true,0.67\n"
        "sample,PGS000001,false,false,false,false,false,NA,unmatched,1,true,0.67\n"
        "sample,PGS000002,false,false,false,false,false,NA,unmatched,5,false,0.0\n")
    if not adjusted:
        put(f"{d}/sample/score/aggregated_scores.txt.gz",
            "sampleset\tFID\tIID\tPGS\tSUM\tDENOM\nsample\tHG\tHG\tPGS000001\t62.0\t6\n", gz=True)
        return
    put(f"{d}/sample/score/sample_pgs.txt.gz",
        "sampleset\tFID\tIID\tPGS\tSUM\tZ_MostSimilarPop\tZ_norm1\tZ_norm2\tpercentile_MostSimilarPop\n"
        "sample\tHG\tHG\tPGS000001\t1.25\t1.6\t1.5\t1.4\t94.44\n"
        "reference\tR1\tR1\tPGS000001\t10.0\t0.1\t0.1\t0.1\t12.3\n", gz=True)
    put(f"{d}/sample/score/sample_popsimilarity.txt.gz",
        "sampleset\tFID\tIID\tPC1\tPC2\tPC10\tSuperPop\tUnrelated\tRF_P_AFR\tRF_P_EUR\tMostSimilarPop\tMostSimilarPop_LowConfidence\tREFERENCE\n"
        "reference\tR1\tR1\t-1\t2\t0\tAFR\tTrue\t0.9\t0.1\tAFR\tFalse\tTrue\n"
        "sample\tHG\tHG\t0.012\t-0.034\t0.5\t\t\t0.05\t0.95\tEUR\tFalse\tFalse\n", gz=True)


def main():
    work = tempfile.mkdtemp()
    try:
        # 1. prs-format
        put(f"{work}/in/PGS000001_hmPOS_GRCh38.txt.gz", HM, gz=True)
        put(f"{work}/labels.tsv", "# comment\npgs_id\ttrait_reported\nPGS000001\tCoronary artery disease\n")
        rc = collect_summary.main(["prs-format", "--scores", f"{work}/in", "--labels", f"{work}/labels.tsv",
                                   "--out", f"{work}/pgs", "--alleles", f"{work}/alleles.tsv"])
        check("prs-format exits 0", rc == 0, rc)
        with gzip.open(f"{work}/pgs/PGS000001.txt.gz", "rt") as f:
            text = f.read()
        lines = text.splitlines()
        check("prs-format: the header pgsc_calc reads, labelled from the score list",
              lines[:5] == ["#pgs_id=PGS000001", "#pgs_name=PGS000001", "#trait_reported=Coronary artery disease",
                            "#genome_build=GRCh38", "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight"], lines[:5])
        check("prs-format: GRCh38 positions without chr, other allele from hm_inferOtherAllele, A/C dropped, unplaced and chrX/chrY rows dropped",
              lines[5:] == ["1\t1100\tA\tG\t0.5", "2\t2200\tC\tT\t-0.25", "3\t3300\tG\t\t1"], lines[5:])
        with open(f"{work}/alleles.tsv") as f:
            al = f.read().splitlines()
        check("prs-format: every effect and other allele of the autosomes, chr-prefixed and sorted",
              al == ["chr1\t1100\tA", "chr1\t1100\tG", "chr2\t2200\tC", "chr2\t2200\tT", "chr3\t3300\tG"], al)
        with open(f"{work}/pgs/chrx_rows.tsv") as f:
            cx = f.read().splitlines()
        check("prs-format: chrx_rows.tsv says the score's one chrX row was left out",
              cx == ["pgs_id\tchrx_rows\tchrx", "PGS000001\t1\tleft out"], cx)
        # --keep-x: the sample's sex is known
        rc = collect_summary.main(["prs-format", "--scores", f"{work}/in", "--labels", f"{work}/labels.tsv",
                                   "--out", f"{work}/pgsx", "--alleles", f"{work}/allelesx.tsv", "--keep-x"])
        check("prs-format --keep-x exits 0", rc == 0, rc)
        with gzip.open(f"{work}/pgsx/PGS000001.txt.gz", "rt") as f:
            lx = f.read().splitlines()
        check("prs-format --keep-x: the chrX row stays, the chrY row still goes",
              lx[5:] == ["1\t1100\tA\tG\t0.5", "2\t2200\tC\tT\t-0.25", "3\t3300\tG\t\t1", "X\t51000000\tA\tG\t1"], lx[5:])
        with open(f"{work}/allelesx.tsv") as f:
            alx = f.read().splitlines()
        check("prs-format --keep-x: the chrX alleles are in the allele list, after the autosomes",
              alx[-2:] == ["chrX\t51000000\tA", "chrX\t51000000\tG"] and not any(a.startswith("chrY") for a in alx), alx)
        with open(f"{work}/pgsx/chrx_rows.tsv") as f:
            cxx = f.read().splitlines()
        check("prs-format --keep-x: chrx_rows.tsv says the chrX row was scored",
              cxx == ["pgs_id\tchrx_rows\tchrx", "PGS000001\t1\tscored"], cxx)
        put(f"{work}/bad37/PGS000009.txt.gz", HM.replace("HmPOS_build=GRCh38", "HmPOS_build=GRCh37"), gz=True)
        rc = collect_summary.main(["prs-format", "--scores", f"{work}/bad37", "--out", f"{work}/o37", "--alleles", f"{work}/a37"])
        check("prs-format refuses a file harmonised to GRCh37", rc == 1, rc)
        put(f"{work}/badnonadd/PGS000010.txt.gz",
            "#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\teffect_weight\tis_recessive\n1\t5\tA\t1\tTrue\n", gz=True)
        rc = collect_summary.main(["prs-format", "--scores", f"{work}/badnonadd", "--out", f"{work}/onon", "--alleles", f"{work}/anon"])
        check("prs-format refuses a non-additive score", rc == 1, rc)
        rc = collect_summary.main(["prs-format", "--scores", f"{work}/in", "--ids", "PGS000001,PGS000777",
                                   "--out", f"{work}/omiss", "--alleles", f"{work}/amiss"])
        check("prs-format: an id of the list with no file is an error", rc == 1, rc)

        # A second formatted score, for the table
        put(f"{work}/pgs/PGS000002.txt.gz", "#pgs_id=PGS000002\n#trait_reported=Trait two\n#genome_build=GRCh38\n"
            "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n" + "1\t1\tA\t\t1\n" * 5, gz=True)

        # 2. prs-table, no panel
        pgsc_results(f"{work}/raw", adjusted=False)
        out = f"{work}/S_prs_summary.tsv"
        rc = collect_summary.main(["prs-table", "--sample", "S", "--results", f"{work}/raw", "--sampleset", "sample",
                                   "--scores", f"{work}/pgs", "--input-kind", "gvcf", "--out", out])
        hdr, rows = read_rows(out)
        check("prs-table exits 0", rc == 0, rc)
        check("prs-table: the columns, Input last", hdr == collect_summary.PRS_COLUMNS and hdr[-1] == "Input", hdr)
        check("prs-table: sum 62 (not 62.0), 2 of 3 matched, no percentile, its chrX row left out",
              rows[0] == ["Coronary artery disease", "PGS000001", "62", "2", "3", "66.7", "NA", "NA", "1 left out", "gvcf"], rows[0])
        check("prs-table: a score pgsc_calc dropped has no sum and its counts; no chrX entry for a score prs-format did not write",
              rows[1] == ["Trait two", "PGS000002", "NA", "0", "5", "0.0", "NA", "NA", "", "gvcf"], rows[1])

        # 3. prs-table with a panel
        pgsc_results(f"{work}/adj", adjusted=True)
        anc = f"{work}/S_ancestry.tsv"
        rc = collect_summary.main(["prs-table", "--sample", "S", "--results", f"{work}/adj", "--sampleset", "sample",
                                   "--scores", f"{work}/pgs", "--input-kind", "gvcf", "--panel", "pgsc_1000G_v1",
                                   "--ancestry-out", anc, "--out", out])
        hdr, rows = read_rows(out)
        check("prs-table with a panel exits 0", rc == 0, rc)
        check("prs-table: the target's sum, percentile and group, not the reference sample's",
              rows[0][2:8] == ["1.25", "2", "3", "66.7", "94.4", "EUR"], rows[0])
        kv = dict(l.rstrip("\n").split("\t", 1) for l in open(anc))
        check("ancestry table: population, panel, probabilities and the PCs in order",
              kv.get("population") == "EUR" and kv.get("reference_panel") == "pgsc_1000G_v1"
              and kv.get("probability_EUR") == "0.95" and kv.get("PC1") == "0.012" and kv.get("PC2") == "-0.034"
              and list(k for k in kv if k.startswith("PC")) == ["PC1", "PC2", "PC10"], kv)

        # 4. zero matches
        rc = collect_summary.main(["prs-table", "--sample", "S", "--results", f"{work}/none", "--sampleset", "sample",
                                   "--scores", f"{work}/pgs", "--input-kind", "vcf", "--zero-matches", "--out", out])
        hdr, rows = read_rows(out)
        check("prs-table --zero-matches: every score unmatched, Input vcf",
              rc == 0 and [r[2:6] + r[9:] for r in rows] == [["NA", "0", "3", "0.0", "vcf"], ["NA", "0", "5", "0.0", "vcf"]], rows)
        rc = collect_summary.main(["prs-table", "--sample", "S", "--results", f"{work}/none", "--sampleset", "sample",
                                   "--scores", f"{work}/pgs", "--input-kind", "vcf", "--out", out])
        check("prs-table: no pgsc_calc output and no --zero-matches is an error", rc == 1, rc)

        # 4b. every score under pgsc_calc's minimum overlap: no sum, matched unknown, the rate from the log
        put(f"{work}/below.log", "ERROR Score PGS000001 fails minimum matching threshold (33.33% variants match)\n"
            "ERROR pgscatalog.core.lib.pgsexceptions.ZeroMatchesError: All scores fail to meet match threshold 0.75\n")
        rc = collect_summary.main(["prs-table", "--sample", "S", "--results", f"{work}/none", "--sampleset", "sample",
                                   "--scores", f"{work}/pgs", "--input-kind", "gvcf", "--below-threshold", f"{work}/below.log",
                                   "--out", out])
        hdr, rows = read_rows(out)
        check("prs-table --below-threshold: no sum, matched NA (not 0), the rate from the log, NA without one",
              rc == 0 and [r[2:6] for r in rows] == [["NA", "NA", "3", "33.3"], ["NA", "NA", "5", "NA"]], rows)

        # 5. the reports
        d = f"{work}/sample_raw/S"
        put(f"{d}/prs/S_prs_summary.tsv", "\t".join(collect_summary.PRS_COLUMNS) + "\n"
            "Coronary artery disease\tPGS000001\t62\t2\t3\t66.7\tNA\tNA\t8 left out\tgvcf\n"
            "Trait two\tPGS000002\t1\t5\t5\t100.0\tNA\tNA\t3 scored\tgvcf\n")
        summ = collect_summary.collect("S", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        check("no panel: the raw-score line in both reports", "Raw score only" in txt and "Raw score only" in html)
        check("chrX left out: both reports name the score and its row count, not the scored one",
              "Rows on chrX were left out, because the sample's sex was not given: PGS000001 (8)." in txt
              and "Rows on chrX were left out, because the sample&#x27;s sex was not given: PGS000001 (8)." in html
              and "PGS000002 (3)" not in txt, txt)
        check("no panel: the HTML row says raw score only", "<td>raw score only</td>" in html)
        check("no panel: not-assessed names the missing percentiles",
              collect_summary.PRS_NOT_ADJUSTED in summ["not_assessed"])
        d = f"{work}/sample_adj/S"
        put(f"{d}/prs/S_prs_summary.tsv", "\t".join(collect_summary.PRS_COLUMNS) + "\n"
            "Coronary artery disease\tPGS000001\t1.25\t2\t3\t66.7\t94.4\tEUR\t0\tgvcf\n")
        put(f"{d}/ancestry/S_ancestry.tsv", "key\tvalue\nsample\tS\nreference_panel\tpgsc_1000G_v1\npopulation\tEUR\n")
        summ = collect_summary.collect("S", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        check("panel: the percentile with its group in the text report", "percentile 94.4 (EUR)" in txt, txt)
        check("panel: the percentile with its group in the HTML report", "<td>94.4 (EUR)</td>" in html)
        check("panel: the note names the group and the panel, no raw-score line, no chrX line",
              "among the EUR samples of the pgsc_1000G_v1 reference panel" in txt and "Raw score only" not in txt
              and "Raw score only" not in html and "Rows on chrX" not in txt, txt)
        check("panel: no not-assessed line about percentiles", collect_summary.PRS_NOT_ADJUSTED not in summ["not_assessed"])
        caveat = ("Most PRS scores were developed in European-ancestry populations and predict less well for "
                  "other ancestries, even with an ancestry-adjusted percentile.")
        check("panel: the ancestry-transfer caveat in both reports", caveat in txt and caveat in html, txt)
        raw = collect_summary.collect("S", f"{work}/sample_raw/S")
        check("no panel: the ancestry-transfer caveat in both reports",
              caveat in render_report.text_report(raw) and caveat in render_report.html_report(raw))

        # 5b. the fine-grained ancestry section (step 26b), in the same sample
        # folder as the PRS card: both run in one report, and the PRS block's
        # own local names must not break this one.
        put(f"{d}/ancestry_fine/S_ancestry_fine.tsv",
            "level\tcode\tgroup\tfraction\tse\tlow\thigh\n"
            "continent\tWEA\tWest Eurasia (Europe, Middle East, North Africa)\t0.8200\t0.0300\t0.7600\t0.8800\n"
            "continent\tCSA\tCentral and South Asia\t0.1500\t0.0300\t0.0900\t0.2100\n"
            "continent\tAFR\tAfrica\t0.0300\t0.0100\t0.0100\t0.0500\n"
            "region\tLEV\tLevant\t0.4600\t0.1000\t0.2600\t0.6600\n"
            "region\tSAR\tSardinia\t0.0050\t0.0200\t0.0000\t0.0600\n"
            "region\tOCE\tOceania\t0.0010\t0.0005\t0.0000\t0.0020\n"
            "painting\tWEA\tWest Eurasia (Europe, Middle East, North Africa)\t0.6400\t.\t.\t.\n"
            "painting\tAFR\tAfrica\t0.3600\t.\t.\t.\n")
        with open(f"{d}/ancestry_fine/S_painting.png", "wb") as f_:
            f_.write(b"\x89PNG\r\n\x1a\n" + b"0" * 64)
        summ = collect_summary.collect("S", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        check("26b: both reports still render with the PRS card beside this one",
              "percentile 94.4 (EUR)" in txt and "<td>94.4 (EUR)</td>" in html, txt)
        check("26b: the text section has the group, its fraction and its bootstrap range",
              "## Fine-Grained Ancestry (ADMIXTURE, RFMix)" in txt
              and "West Eurasia (Europe, Middle East, North Africa)     82% (76% to 88%)" in txt, txt)
        check("26b: a regional group of 0.5% whose range reaches 6% is shown, one of 0.1% is not",
              "Sardinia" in txt and "Oceania" not in txt.split("## Fine-Grained")[1].split("##")[0], txt)
        check("26b: the painted shares and the picture are in the HTML card",
              "<h2>Fine-Grained Ancestry</h2>" in html and "West Eurasia (Europe, Middle East, North Africa) 64%" in html
              and 'src="data:image/png;base64,' in html, html[:200])
        check("26b: both reports say a group name is not a nationality",
              "not a nationality" in txt and "not a nationality" in html)

        # 6. ClinVar badge: one homozygous 3-star Pathogenic hit
        d = f"{work}/clinvar_hom/S"
        put(f"{d}/clinvar/S_clinvar_hits.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
            "chr1\t20\t2\tT\tG\t50\tPASS\tGENEINFO=GENEB:22;CLNSIG=Pathogenic;"
            "CLNREVSTAT=reviewed_by_expert_panel\tGT\t1/1\n")
        summ = collect_summary.collect("S", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        badge_line = next((l for l in html.splitlines() if "ClinVar matches" in l), "")
        check("ClinVar: a homozygous 3-star P/LP hit gets no green badge", "badge-green" not in badge_line
              and 'badge-yellow">1<' in badge_line, badge_line)
        freq = "a heterozygous hit in a dominant gene as well as a homozygous one: a hit above 1% in gnomAD"
        check("ClinVar: both reports ask for every hit's population frequency, not only a homozygous one's",
              freq in txt and freq in html and "Look up a homozygous P/LP hit" not in txt + html, txt)
        check("ClinVar: without step 13 both reports say the hit's frequency is unknown",
              "frequency: gnomAD AF unknown" in txt and "<td>gnomAD AF unknown</td>" in html, txt)
        put(f"{d}/clinvar/S_clinvar_hits.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
        html0 = render_report.html_report(collect_summary.collect("S", d))
        check("ClinVar: no hit gets a gray badge, not a green all-clear",
              'ClinVar matches</span><span class="value"><span class="badge badge-gray">0<' in html0)

        # 6b. ClinVar: step 06's cDNA read check beside each hit
        d = f"{work}/clinvar_cdna/S"
        put(f"{d}/clinvar/S_clinvar_hits.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
            "chr1\t20\t2\tT\tG\t50\tPASS\tGENEINFO=GENEB:22;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel\tGT\t0/1\n"
            "chr2\t30\t3\tC\tA\t50\tPASS\tGENEINFO=GENEC:33;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel\tGT\t0/1\n")
        summ = collect_summary.collect("S", d)
        html = render_report.html_report(summ)
        check("cDNA check: no table, no column (a Nextflow folder)", "cDNA read check" not in html
              and "cdna_flagged" not in summ["sections"]["clinvar"]["data"], html)
        flag = "possible cDNA/transgene contamination: check the reads"
        put(f"{d}/clinvar/S_clinvar_cdna_check.tsv",
            "chrom\tpos\tref\talt\talt_reads\talt_clipped\tref_reads\tref_clipped\tcdna_check\tcdna_note\n"
            f"chr1\t20\tT\tG\t9\t8\t11\t0\tflagged\t{flag}\n"
            "chr2\t30\tC\tA\t12\t0\t10\t0\tclear\t\n")
        summ = collect_summary.collect("S", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        shown = f"{flag} (8 of 9 alt reads and 0 of 11 ref reads clipped at an exon edge)"
        check("cDNA check: the text report prints the flag under the flagged hit",
              "chr1:20 T>G" in txt and f"reads: {shown}" in txt and "reads: clear" in txt, txt)
        check("cDNA check: the HTML hits table has the column, the flag and clear",
              "<th>cDNA read check</th>" in html and f"<td>{shown}</td>" in html and "<td>clear</td>" in html, html)
        check("cDNA check: the card counts one flagged hit",
              'Flagged by the cDNA read check</span><span class="value"><span class="badge badge-red">1<' in html, html)
        check("cDNA check: the summary validates against the schema", not validate.validate(
            json.load(open(os.path.join(REPO, "tests", "schema", "summary.schema.json"))),
            json.loads(json.dumps(summ))))

        # 7. telomere card
        d = f"{work}/telomere/S"
        put(f"{d}/telomere/S/S_summary.tsv", "PID\tsample\tintratelomeric_reads\ttotal_reads_with_tel_gc\ttel_content\n"
            "S\ttumor\t61873\t185204533\t334.08\n")
        summ = collect_summary.collect("S", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        check("telomere: the text heading and the value with its unit",
              "## Telomere content (relative, TelomereHunter)" in txt
              and "Telomere content: 334.08 intratelomeric reads per million reads with 48-52% GC" in txt, txt)
        check("telomere: the HTML card title and the unit", "<h2>Telomere content (relative)</h2>" in html
              and "intratelomeric reads per million reads with 48-52% GC" in html and "not a telomere length" in html)
        said = [p_ for p_ in glob.glob(os.path.join(REPO, "bin", "*.py")) if "Telomere Length" in open(p_).read()]
        check("no file in bin/ says 'Telomere Length'", not said, said)

        # 8. mosaic loss of chrY or chrX
        d = f"{work}/mosaic/S"

        def mosaic_table(rows):
            put(f"{d}/qc/S_sample_qc.tsv", "key\tvalue\nsample\tS\ninferred_sex\tmale\n"
                + "".join(f"mosaic_{k}\t{v}\n" for k, v in rows.items()))
            summ_ = collect_summary.collect("S", d)
            return render_report.text_report(summ_), render_report.html_report(summ_)
        male = {"sex": "male", "chrom": "chrY", "depth_ratio": "0.351", "depth_ratio_expected": "0.5",
                "loss_depth": "0.298", "loss_depth_low": "0.293", "loss_depth_high": "0.303", "floor": "0.100",
                "call": "detected", "reason": "the low end of the 95% range is above 10%"}
        txt, html = mosaic_table(male)
        line = ("From depth: estimated share of the cells that lost the Y 29.8% (95% range 29.3% to 30.3%); "
                "chrY to autosome depth ratio 0.351")
        card_m = html.split("<h2>Mosaic Loss of Chromosome Y or X</h2>", 1)[-1].split("<h2>", 1)[0]
        check("mosaic: the text report heading, call and estimate",
              "## Mosaic loss of chromosome Y or X" in txt and "chrY: detected (the low end" in txt and line in txt, txt)
        check("mosaic: the HTML card shows the same", 'badge-yellow">detected<' in card_m and line in card_m, card_m)
        check("mosaic: both say population-level, not a diagnosis, and the 10% floor",
              all("population-level association, not a diagnosis" in t and "Below about 10% of the cells" in t
                  for t in (txt, card_m)), card_m)
        check("mosaic: a man has no allele-fraction line", "allele fractions at" not in txt, txt)
        txt, html = mosaic_table(dict(male, sex="female", chrom="chrX", depth_ratio="0.849",
                                      depth_ratio_expected="1.0", x_het_sites="600", loss_af="0.316",
                                      loss_af_low="0.288", loss_af_high="0.340"))
        check("mosaic: a woman shows one X and the allele-fraction line",
              "of the cells that lost one X" in txt and "From chrX allele fractions at 600 heterozygous sites: 31.6% "
              "(95% range 28.8% to 34.0%)" in txt and "31.6%" in html, txt)
        txt, html = mosaic_table(dict(male, depth_ratio="", loss_depth="0.190", loss_depth_low="-0.002",
                                      loss_depth_high="0.382", call="not_detected"))
        card_m = html.split("<h2>Mosaic Loss of Chromosome Y or X</h2>", 1)[-1].split("<h2>", 1)[0]
        check("mosaic: not detected says 'above 10%', and a missing depth ratio prints '?'",
              "chrY: not detected above 10%" in txt and "below 10%" not in txt and "depth ratio ?," in txt
              and 'badge-gray">not detected above 10%<' in card_m, txt)
        put(f"{d}/qc/S_sample_qc.tsv", "key\tvalue\nsample\tS\ninferred_sex\tmale\n")
        summ = collect_summary.collect("S", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        card_m = html.split("<h2>Mosaic Loss of Chromosome Y or X</h2>", 1)[-1].split("<h2>", 1)[0]
        check("mosaic: no estimate, the card says Not run and the text has no section",
              'badge-gray">Not run<' in card_m and "## Mosaic loss" not in txt, card_m)
        # 9. paralog card (step 35b): hidden when the step did not run; a
        # result older than a run in which step 35b failed is marked Stale
        check("paralogs: no card when step 35b did not run", "Paralog Genes (Paraphase)" not in html)
        d = f"{work}/paralogs/S"
        put(f"{d}/paralogs/S_paraphase_copy_number.tsv",
            "region\tgene\tcopy_number\tstatus\tdetail\ttotal_cn\tregion_depth\tgenome_depth\n"
            "pms2\tPMS2\t2\tcalled\t.\t4\t60\t30\n")
        os.utime(f"{d}/paralogs/S_paraphase_copy_number.tsv", (1, 1))
        put(f"{d}/logs/run_status.tsv", "meta\tstarted_epoch\t1000\nstep\t35b\tfailed\n")
        summ = collect_summary.collect("S", d)
        html = render_report.html_report(summ)
        card_p = html.split("<h2>Paralog Genes (Paraphase)</h2>", 1)[-1].split("<h2>", 1)[0]
        check("paralogs: a result older than a run where step 35b failed is stale, with its PMS2 row",
              summ["sections"]["paralogs"]["state"] == "stale" and 'badge-yellow">Stale<' in card_p
              and "<td>PMS2</td><td>2</td>" in card_p, card_p)
        # 12. raw-data export card
        d = f"{work}/consumer/S"
        summ = collect_summary.collect("S", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        check("raw-data export: no card and no text section without the step's table",
              "Raw-Data Export" not in txt and "Raw-Data Export" not in html)
        put(f"{d}/consumer_export/S_consumer_export.tsv", "key\tvalue\nsites\t200\ncalled\t150\nno_call\t50\n"
            "from_reference_blocks\t120\nfile_23andme\tS_23andme_v5.zip\nfile_ancestrydna\tS_ancestrydna_layout.zip\n"
            "no_call: no gVCF record\t30\nno_call: filtered call\t20\n")
        summ = collect_summary.collect("S", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        want = ["Sites with a genotype: 150 of 200 (75.0%); '--' at 50",
                "Homozygous reference from the gVCF's reference blocks: 120",
                "consumer_export/S_23andme_v5.zip and consumer_export/S_ancestrydna_layout.zip",
                "'--' because: no gVCF record 30; filtered call 20", "Uploading one to a service shares your genotypes"]
        check("raw-data export: the text section has every line", "## Raw-Data Export (23andMe v5 layout)" in txt
              and all(w_ in txt for w_ in want), txt)
        check("raw-data export: the HTML card has every line", "<h2>Raw-Data Export (23andMe v5 layout)</h2>" in html
              and all(render_report.E(w_) in html for w_ in want), html)

        # 10. blood groups
        d = f"{work}/blood/S"
        put(f"{d}/bloodtype/S_bloodtype.tsv", "system\tphenotype\tgenotype\tconfidence\tnote\n"
            "ABO\tO\tO.01/O.01\thigh\t\nRhD\tD+\tRHD/RHD (RHD copies by read depth 1.78)\thigh\t\n"
            "Duffy\tFy(a+b-)\tFY*A/FY*B GATA\thigh\tc.-67T>C\nHPA-1\tHPA-1(a+b+)\t1a/1b\thigh\t\n")
        summ = collect_summary.collect("S", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        check("blood groups: the text section has every system",
              "## Blood Groups (predicted from the genome)" in txt
              and all(x in txt for x in ("ABO       O ", "RhD       D+ ", "Duffy     Fy(a+b-) ", "HPA-1     HPA-1(a+b+) ")), txt)
        check("blood groups: the HTML card has ABO, RhD and the table",
              "<h2>Blood Groups (predicted)</h2>" in html and "<td>Fy(a+b-)</td>" in html
              and "<td>FY*A/FY*B GATA</td>" in html and ">D+</span>" in html, html[-3000:])
        check("blood groups: both reports say it is predicted, not typed from blood",
              txt.count("not typed from blood") == 1 and html.count("not typed from blood") == 1)
        nothing = collect_summary.collect("S", f"{work}/telomere/S")
        check("blood groups: no card without step 41", "Blood Groups" not in render_report.html_report(nothing))
        words = re.compile(r"transfus|donat|donor", re.I)
        said = [n for n, t in (("text report", txt), ("HTML report", html),
                               ("bin/bloodtype_report.py", open(os.path.join(REPO, "bin", "bloodtype_report.py")).read()))
                if words.search(t)]
        check("blood groups: nothing about transfusion or donation", not said, said)
        check("control: the wording check finds the word", bool(words.search("for a transfusion")))

        # 11. traits card
        import time
        import traits_report
        d = f"{work}/traits/S"
        rows = traits_report.read_table(os.path.join(REPO, "assets", "traits.tsv"))
        put(f"{d}/traits/S_traits.tsv", "")
        traits_report.write_tsv(f"{d}/traits/S_traits.tsv", traits_report.report(rows, {}))
        summ = collect_summary.collect("S", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        check("traits: both reports say a reading is not a diagnosis",
              render_report.TRAITS_NOTE in txt and render_report.TRAITS_NOTE in html, txt)
        check("traits: with no reads anywhere every row is not read, and the APOE line too",
              summ["sections"]["traits"]["data"]["called"] == 0
              and summ["sections"]["traits"]["data"]["apoe"][-1]["status"] == "not_covered")
        old = time.time() - 86400
        os.utime(f"{d}/traits/S_traits.tsv", (old, old))
        put(f"{d}/logs/run_status.tsv", f"meta\tstarted_epoch\t{time.time() - 3600}\nstep\t42\tskipped (no gVCF: S/vcf/S.g.vcf.gz)\n")
        summ = collect_summary.collect("S", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        card_tr = html.split("<h2>Traits</h2>", 1)[-1].split("<h2>", 1)[0]
        check("traits: older than a run that skipped step 42 -> stale in both reports",
              summ["sections"]["traits"]["state"] == "stale" and "[STALE: from an earlier run" in txt
              and 'class="stale">STALE: from an earlier run' in card_tr, card_tr)
        os.utime(f"{d}/traits/S_traits.tsv", None)
        check("traits: newer than the run -> not stale",
              collect_summary.collect("S", d)["sections"]["traits"]["state"] == "ok")

        # 13. couple card (step 39)
        import couple_report
        genes = couple_report.read_genes(os.path.join(REPO, "assets", "carrier_genes.tsv"))
        g = genes["NAGA"]
        pos = (g["start"] + g["end"]) // 2
        cv = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
              "chr22\t%d\t1\tA\tG\t.\tPASS\tGENEINFO=NAGA:4668;CLNSIG=Pathogenic;"
              "CLNREVSTAT=reviewed_by_expert_panel\tGT\t0/1\n")
        c = f"{work}/couple"
        put(f"{c}/ped.tsv", "sample\tsex\tfather\tmother\tpartner\nP1\tfemale\t.\t.\tP2\nP2\tmale\t.\t.\tP1\n")
        put(f"{c}/P1.vcf", cv % pos)
        put(f"{c}/P2.vcf", cv % (pos + 5))
        couple_report.main(["screen", "--pedigree", f"{c}/ped.tsv", "--genes", os.path.join(REPO, "assets", "carrier_genes.tsv"),
                            "--clinvar", f"P1={c}/P1.vcf", "--clinvar", f"P2={c}/P2.vcf", "--out-dir", f"{c}/out"])
        d = f"{c}/P2"
        put(f"{d}/couple/P2_couple_carrier.tsv", open(f"{c}/out/P2_couple_carrier.tsv").read())
        summ = collect_summary.collect("P2", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        card_c = html.split("<h2>Couple Carrier Screen</h2>", 1)[-1].split("<h2>", 1)[0]
        want = ("Partner: P1", "Shared: NAGA (Schindler disease): 1 in 4 per pregnancy",
                "This sample: chr22:%d A&gt;G het" % (pos + 5), "P1: chr22:%d A&gt;G het, Pathogenic, 3 stars" % pos,
                "genetic counsellor")
        check("couple: the card names the partner, the shared gene, its chance, both findings and the counsellor",
              all(w in card_c for w in want), card_c)
        check("couple: the text report says the same",
              all(w.replace("&gt;", ">") in txt for w in want) and "## Couple Carrier Screen (step 39)" in txt, txt)
        os.remove(f"{d}/couple/P2_couple_carrier.tsv")
        summ = collect_summary.collect("P2", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        check("couple: without the table there is no card and no text section",
              "Couple Carrier Screen" not in html and "Couple Carrier Screen" not in txt)
    finally:
        shutil.rmtree(work)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
