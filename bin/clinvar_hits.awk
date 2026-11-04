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
#          DP GQ LOW_CONFIDENCE
#   STARS          ClinVar's review-status stars: 4 practice guideline,
#                  3 expert panel, 2 multiple submitters with no conflict,
#                  1 single submitter or conflicting classifications,
#                  0 no assertion criteria (and anything unknown).
#   GENOTYPE       het, hom, "one copy, other copy not called ref (GT)" for
#                  a GT of one alt allele and one missing one (1/., ./1,
#                  1|.), or the GT as written. Step 6 splits a
#                  multiallelic record with the other alt allele written
#                  as missing, not as ref (bcftools norm --multi-overlaps .),
#                  so a 1/2 call (at a poly-T tract, T5 on one copy and T9
#                  on the other) gives 1/. and ./1 here, never het.
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
#   DP, GQ         the call's FORMAT/DP and FORMAT/GQ; "." when the caller
#                  wrote none.
#   LOW_CONFIDENCE "yes" when GQ or DP is below step 6's floor, the
#                  ##PGP_CALL_FLOOR=GQ=..,DP=.. header line step 6 writes
#                  for a floor other than the default (GQ 30 and DP 10 when
#                  the line is missing); "no" otherwise,
#                  also when neither value is there to compare.
#
# The always-report list (the known deep intronic sites):
#
#   awk -v mode=sites [-v min_gq=N] -f clinvar_hits.awk SITES.tsv VCF_CALLS.txt GVCF_CALLS.txt
#
# SITES.tsv is assets/deep_intronic_sites.tsv (comment lines start with #, then
# a header: chrom pos ref alt gene hgvs clinvar_id). VCF_CALLS.txt holds step
# 6's split PASS records at those positions and GVCF_CALLS.txt the gVCF records
# that cover them (an empty file without a gVCF), both as bcftools query
# '%CHROM\t%POS\t%END\t%REF\t%ALT\t[%GT]\t[%GQ]\t[%DP]\n'; set gvcf=no when
# there is no gVCF. One row per site, every site, carrier or not:
#   chrom pos ref alt gene hgvs clinvar_id call detail
#   call    the genotype words above (het, hom, ...) when the sample carries
#           the site's allele; "ref" for a reference call (a 0/0 record, or a
#           gVCF reference block over the site), "ref (low confidence)" when
#           its GQ is below min_gq (step 6's call floor, default 30);
#           "other allele" when the sample carries another allele there;
#           "no call" otherwise.
#   detail  "GQ n, DP n" of the call ("GQ n (gVCF block)" for a block, with
#           "below the call floor of GQ n" for a low-confidence ref), the
#           other allele, or why there is no call.
BEGIN { FS = OFS = "\t"; if (min_gq == "") min_gq = 30; if (min_dp == "") min_dp = 10 }

# A GT in words: het, hom, one copy of an alt beside a missing allele, or as written.
function zygosity(gt) {
    if (gt == "0/1" || gt == "1/0" || gt == "0|1" || gt == "1|0") return "het"
    if (gt == "1/1" || gt == "1|1") return "hom"
    if (gt ~ /^([1-9][\/|]\.|\.[\/|][1-9])$/) return "one copy, other copy not called ref (" gt ")"
    return gt
}
# A reference call in words: "ref", or "ref (low confidence)" below the GQ floor.
function refcall(gq) { return (gq ~ /^[0-9]+$/ && gq + 0 < min_gq + 0) ? "ref (low confidence)" : "ref" }
function floornote(gq) { return (gq ~ /^[0-9]+$/ && gq + 0 < min_gq + 0) ? ", below the call floor of GQ " min_gq : "" }
function quality(gq, dp,    out) {
    out = (gq ~ /^[0-9]+$/) ? "GQ " gq : ""
    if (dp ~ /^[0-9]+$/) out = out (out != "" ? ", " : "") "DP " dp
    return out
}
mode == "sites" && FILENAME == ARGV[1] {
    if ($0 ~ /^#/ || $1 == "chrom" || NF < 7) next
    ns++; site[ns] = $0
    next
}
mode == "sites" && FILENAME == ARGV[2] { nv++; vcf[nv] = $0; next }
mode == "sites" && FILENAME == ARGV[3] { ng++; gv[ng] = $0; next }
mode == "sites" { next }
END {
    if (mode != "sites") exit
    print "chrom", "pos", "ref", "alt", "gene", "hgvs", "clinvar_id", "call", "detail"
    for (i = 1; i <= ns; i++) {
        split(site[i], st, "\t")
        call = ""; detail = ""; other = ""
        for (j = 1; j <= nv; j++) {
            split(vcf[j], vr, "\t")
            if (vr[1] != st[1] || vr[2] != st[2]) continue
            if (vr[4] == st[3] && vr[5] == st[4]) {
                if (vr[6] ~ /[1-9]/) { call = zygosity(vr[6]); detail = quality(vr[7], vr[8]) }
                else if (vr[6] ~ /^0[\/|]0$/) { call = refcall(vr[7]); detail = quality(vr[7], vr[8]) floornote(vr[7]) }
            } else if (vr[6] ~ /[1-9]/ && other == "") {
                other = vr[4] ">" vr[5] " " zygosity(vr[6])
            }
        }
        if (call == "" && other != "") { call = "other allele"; detail = other }
        if (call == "") {
            for (j = 1; j <= ng; j++) {
                split(gv[j], gr, "\t")
                if (gr[1] != st[1] || gr[2] + 0 > st[2] + 0 || gr[3] + 0 < st[2] + 0) continue
                if (gr[6] ~ /^0[\/|]0$/) {
                    call = refcall(gr[7])
                    detail = quality(gr[7], "") (gr[7] ~ /^[0-9]+$/ ? " " : "") "(gVCF block)" floornote(gr[7])
                } else {
                    call = "no call"; detail = "the gVCF record over the site is not a PASS reference call (GT " gr[6] ")"
                }
                break
            }
        }
        if (call == "") {
            call = "no call"
            detail = gvcf == "no" ? "no record in the VCF, and no gVCF to tell a reference call from no call" \
                                  : "no record in the VCF or the gVCF covers the site"
        }
        print st[1], st[2], st[3], st[4], st[5], st[6], st[7], call, detail
    }
}
/^##PGP_CALL_FLOOR=/ {
    if (match($0, /GQ=[0-9]+/)) min_gq = substr($0, RSTART + 3, RLENGTH - 3) + 0
    if (match($0, /DP=[0-9]+/)) min_dp = substr($0, RSTART + 3, RLENGTH - 3) + 0
    next
}

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
    dp = "."; gq = "."
    nk = split($9, fk, ":")
    for (i = 1; i <= nk; i++) {
        if (fk[i] == "DP" && f[i] ~ /^[0-9]+$/) dp = f[i]
        else if (fk[i] == "GQ" && f[i] ~ /^[0-9]+$/) gq = f[i]
    }
    low = ((gq != "." && gq + 0 < min_gq) || (dp != "." && dp + 0 < min_dp)) ? "yes" : "no"
    zyg = zygosity(gt)
    note = "."; phase = ""
    if (gt ~ /^[0-9]+\/[0-9]+$/) {
        split(gt, ga, "/")
        if (ga[1] != ga[2]) phase = "; unphased, so the parts may sit on different haplotypes"
    }
    if (joined != "") note = "combined allele of " variants(joined, " + ") phase
    else if (partof != "") note = "part of a combined allele: " variants(partof, ", ") phase
    else if (one != "") note = "possibly one allele with " variants(one, ", ") ": check ClinVar for a combined record"
    print stars, $1, $2, $4, $5, zyg, gene, clnsig, rev, note, dp, gq, low
}
