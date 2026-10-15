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
# Columns: STARS CHROM POS REF ALT GENOTYPE GENE SIGNIFICANCE REVIEW_STATUS COMBINED_ALLELE
#   STARS          ClinVar's review-status stars: 4 practice guideline,
#                  3 expert panel, 2 multiple submitters with no conflict,
#                  1 single submitter or conflicting classifications,
#                  0 no assertion criteria (and anything unknown).
#   GENOTYPE       het, hom, or the GT as written.
#   GENE           GENEINFO's symbols, comma-joined; "." when absent.
#   SIGNIFICANCE, REVIEW_STATUS  CLNSIG and CLNREVSTAT with "_" as spaces;
#                  "." when absent.
#   COMBINED_ALLELE  from step 6's tags (bin/clinvar_join_adjacent.awk):
#                  "combined allele of A + B" for a combined record,
#                  "part of a combined allele: C" for each of its parts,
#                  "possibly one allele with D: check ClinVar for a combined
#                  record" for two hits that matched no combined record;
#                  "." otherwise. Variants read CHROM:POS REF>ALT. On an
#                  unphased het call (0/1) the first two notes end with
#                  "; unphased, so the parts may sit on different
#                  haplotypes": the genotype cannot say they share one.
BEGIN { FS = OFS = "\t" }

# "chr1:10:C:T|chr1:15:G:A" -> "chr1:10 C>T<sep>chr1:15 G>A"
function variants(v, sep,    n, x, i, f, out) {
    n = split(v, x, "|"); out = ""
    for (i = 1; i <= n; i++) {
        split(x[i], f, ":")
        out = out (i > 1 ? sep : "") f[1] ":" f[2] " " f[3] ">" f[4]
    }
    return out
}
/^#/ { next }
{
    geneinfo = ""; clnsig = ""; rev = ""; joined = ""; partof = ""; one = ""
    n = split($8, kv, ";")
    for (i = 1; i <= n; i++) {
        p = index(kv[i], "=")
        if (p == 0) continue
        k = substr(kv[i], 1, p - 1); v = substr(kv[i], p + 1)
        if (k == "GENEINFO") geneinfo = v
        else if (k == "CLNSIG") clnsig = v
        else if (k == "CLNREVSTAT") rev = v
        else if (k == "PGP_JOINED_FROM") joined = v
        else if (k == "PGP_PART_OF") partof = v
        else if (k == "PGP_ONE_ALLELE") one = v
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
    note = "."; phase = ""
    if (gt ~ /^[0-9]+\/[0-9]+$/) {
        split(gt, ga, "/")
        if (ga[1] != ga[2]) phase = "; unphased, so the parts may sit on different haplotypes"
    }
    if (joined != "") note = "combined allele of " variants(joined, " + ") phase
    else if (partof != "") note = "part of a combined allele: " variants(partof, ", ") phase
    else if (one != "") note = "possibly one allele with " variants(one, ", ") ": check ClinVar for a combined record"
    print stars, $1, $2, $4, $5, zyg, gene, clnsig, rev, note
}
