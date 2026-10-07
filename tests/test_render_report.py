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
     refuses a GRCh37 file and a non-additive one;
  2. prs-table without a panel: the sum from aggregated_scores.txt.gz as
     pgsc_calc wrote it but without ".0", matched and total counts from the
     match summary, Percentile and Ancestry_Group NA, Input last; a score
     pgsc_calc dropped (below its minimum overlap) has no sum;
  3. prs-table with a panel: the percentile and group of the target sample
     only (not the reference samples in the same file), and the ancestry table
     with the population and the principal components;
  4. --zero-matches: every score unmatched, no sum; --below-threshold: no
     sum, the matched count unknown and the rate from pgscatalog-match's log;
  5. both reports: a percentile with its group when there is one, the
     "Raw score only" line when there is none, and the not-assessed line about
     percentiles only then.

Run: python3 tests/test_render_report.py
"""
import gzip
import os
import shutil
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import collect_summary  # noqa: E402
import render_report  # noqa: E402

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
        check("prs-table: sum 62 (not 62.0), 2 of 3 matched, no percentile",
              rows[0] == ["Coronary artery disease", "PGS000001", "62", "2", "3", "66.7", "NA", "NA", "gvcf"], rows[0])
        check("prs-table: a score pgsc_calc dropped has no sum and its counts",
              rows[1] == ["Trait two", "PGS000002", "NA", "0", "5", "0.0", "NA", "NA", "gvcf"], rows[1])

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
              rc == 0 and [r[2:6] + r[8:] for r in rows] == [["NA", "0", "3", "0.0", "vcf"], ["NA", "0", "5", "0.0", "vcf"]], rows)
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
            "Coronary artery disease\tPGS000001\t62\t2\t3\t66.7\tNA\tNA\tgvcf\n")
        summ = collect_summary.collect("S", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        check("no panel: the raw-score line in both reports", "Raw score only" in txt and "Raw score only" in html)
        check("no panel: the HTML row says raw score only", "<td>raw score only</td>" in html)
        check("no panel: not-assessed names the missing percentiles",
              collect_summary.PRS_NOT_ADJUSTED in summ["not_assessed"])
        d = f"{work}/sample_adj/S"
        put(f"{d}/prs/S_prs_summary.tsv", "\t".join(collect_summary.PRS_COLUMNS) + "\n"
            "Coronary artery disease\tPGS000001\t1.25\t2\t3\t66.7\t94.4\tEUR\tgvcf\n")
        put(f"{d}/ancestry/S_ancestry.tsv", "key\tvalue\nsample\tS\nreference_panel\tpgsc_1000G_v1\npopulation\tEUR\n")
        summ = collect_summary.collect("S", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        check("panel: the percentile with its group in the text report", "percentile 94.4 (EUR)" in txt, txt)
        check("panel: the percentile with its group in the HTML report", "<td>94.4 (EUR)</td>" in html)
        check("panel: the note names the group and the panel, no raw-score line",
              "among the EUR samples of the pgsc_1000G_v1 reference panel" in txt and "Raw score only" not in txt
              and "Raw score only" not in html, txt)
        check("panel: no not-assessed line about percentiles", collect_summary.PRS_NOT_ADJUSTED not in summ["not_assessed"])
    finally:
        shutil.rmtree(work)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
