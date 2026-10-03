#!/usr/bin/awk -f
# clinvar_hits.awk: one TSV row per record of a step 6 ClinVar hits VCF.
#
#   awk -f clinvar_hits.awk HITS.vcf
#
# Used by scripts/06-clinvar-screen.sh (host awk) and by the CLINVAR_SCREEN and
# HTML_REPORT Nextflow modules (the bcftools image's awk; bin/ is on the task
# PATH). bin/collect_summary.py reads the same file for the two reports and
# maps review status to stars with the same table; tests/test_collect_summary.py
# checks that the two agree.
#
# Columns: STARS CHROM POS REF ALT GENOTYPE GENE SIGNIFICANCE REVIEW_STATUS
#   STARS          ClinVar's review-status stars: 4 practice guideline,
#                  3 expert panel, 2 multiple submitters with no conflict,
#                  1 single submitter or conflicting classifications,
#                  0 no assertion criteria (and anything unknown).
#   GENOTYPE       het, hom, or the GT as written.
#   GENE           GENEINFO's symbols, comma-joined; "." when absent.
#   SIGNIFICANCE, REVIEW_STATUS  CLNSIG and CLNREVSTAT with "_" as spaces;
#                  "." when absent.
BEGIN { FS = OFS = "\t" }
/^#/ { next }
{
    geneinfo = ""; clnsig = ""; rev = ""
    n = split($8, kv, ";")
    for (i = 1; i <= n; i++) {
        p = index(kv[i], "=")
        if (p == 0) continue
        k = substr(kv[i], 1, p - 1); v = substr(kv[i], p + 1)
        if (k == "GENEINFO") geneinfo = v
        else if (k == "CLNSIG") clnsig = v
        else if (k == "CLNREVSTAT") rev = v
    }
    stars = 0
    if (rev == "practice_guideline") stars = 4
    else if (rev == "reviewed_by_expert_panel") stars = 3
    else if (rev == "criteria_provided,_multiple_submitters,_no_conflicts") stars = 2
    else if (rev ~ /^criteria_provided,_(single_submitter|conflicting_classifications|conflicting_interpretations)$/) stars = 1
    gene = ""
    m = split(geneinfo, g, "|")
    for (i = 1; i <= m; i++) { split(g[i], sym, ":"); gene = gene (i > 1 ? "," : "") sym[1] }
    if (gene == "") gene = "."
    if (clnsig == "") clnsig = "."
    gsub(/_/, " ", clnsig)
    if (rev == "") rev = "."
    gsub(/_/, " ", rev)
    split($10, f, ":"); gt = f[1]
    if (gt == "0/1" || gt == "1/0" || gt == "0|1" || gt == "1|0") zyg = "het"
    else if (gt == "1/1" || gt == "1|1") zyg = "hom"
    else zyg = gt
    print stars, $1, $2, $4, $5, zyg, gene, clnsig, rev
}
