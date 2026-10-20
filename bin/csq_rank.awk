#!/usr/bin/awk -f
# csq_rank.awk: rank each VEP record on one transcript, give it a frequency
# tier and judge its genotype call. A VCF goes in on stdin, the same VCF comes
# out with INFO tags added.
#
# The one ranking rule of the pipeline: steps 23 and 31 run it on the host and
# the CLINICAL_FILTER and SLIVAR_PRIORITIZE Nextflow modules run it in the
# bcftools image (no python there), so the four cannot disagree about a record.
# bin/collect_summary.py applies the same rule to the clinical VCF.
#
#   ... | awk [-v keep=all|noncommon] [-v rare_af=0.01] [-v min_gq=30]
#             [-v min_dp=10] [-v het_lo=0.2] [-v het_hi=0.8] [-v hom_lo=0.8]
#             -f csq_rank.awk
#
# Transcript: the CSQ entry of the MANE Select transcript (MANE_SELECT set),
# else of the canonical one (CANONICAL=YES), else of any transcript, the most
# severe IMPACT within that set. A HIGH consequence on another transcript does
# not make a record HIGH when its MANE Select consequence is intronic.
#   PGP_IMPACT, PGP_CSQ, PGP_GENE  IMPACT, Consequence and SYMBOL of that entry
#   PGP_TX                         MANE, CANONICAL or ANY (NONE: no CSQ)
# Frequency: MAX_AF when the CSQ has that field, else gnomADe_AF and
# gnomADg_AF; the highest value over the record's entries ('&'-joined values
# too).
#   PGP_FREQ  rare (below rare_af), common, unknown (no value: VEP found no
#             known variant for the allele, which is not the same as rare) or
#             none (the CSQ has no frequency field at all)
# Call, from the first sample's FORMAT: GQ below min_gq, DP below min_dp, an
# alt read fraction (AD) outside het_lo..het_hi on a heterozygote or below
# hom_lo on a homozygous alt call makes the call weak.
#   PGP_CALL      ok, weak or unknown (no GQ, DP or AD to judge)
#   PGP_CALL_WHY  what made it weak, e.g. GQ11,AB0.29
# The floors are written to the header as
#   ##PGP_CALL_FLOOR=GQ=..,DP=..,HET_AB=lo-hi,HOM_AB=..
# which bin/collect_summary.py reads. PGP_* tags and that line from an earlier
# pass are replaced, not repeated.
#
# keep=noncommon drops the records whose PGP_FREQ is common.
function add(s, x) { return s == "" ? x : s "," x }
function field(name) { return (name in col) ? v[col[name]] : "" }
BEGIN {
    FS = OFS = "\t"
    if (keep == "") keep = "all"
    if (rare_af == "") rare_af = 0.01
    if (min_gq == "") min_gq = 30
    if (min_dp == "") min_dp = 10
    if (het_lo == "") het_lo = 0.2
    if (het_hi == "") het_hi = 0.8
    if (hom_lo == "") hom_lo = 0.8
    rank["HIGH"] = 4; rank["MODERATE"] = 3; rank["LOW"] = 2; rank["MODIFIER"] = 1
}
/^##INFO=<ID=CSQ,/ {
    s = $0
    sub(/.*Format: /, "", s)
    sub(/".*$/, "", s)
    n = split(s, f, "|")
    for (i = 1; i <= n; i++) col[f[i]] = i
    print
    next
}
/^##INFO=<ID=PGP_(IMPACT|CSQ|GENE|TX|FREQ|CALL|CALL_WHY),/ || /^##PGP_CALL_FLOOR=/ { next }
/^#CHROM/ {
    nfreq = 0
    if ("MAX_AF" in col) freq[++nfreq] = col["MAX_AF"]
    else {
        if ("gnomADe_AF" in col) freq[++nfreq] = col["gnomADe_AF"]
        if ("gnomADg_AF" in col) freq[++nfreq] = col["gnomADg_AF"]
    }
    print "##INFO=<ID=PGP_IMPACT,Number=1,Type=String,Description=\"IMPACT of the ranked CSQ entry: MANE Select, else canonical, else any transcript (bin/csq_rank.awk)\">"
    print "##INFO=<ID=PGP_CSQ,Number=1,Type=String,Description=\"Consequence of the ranked CSQ entry (bin/csq_rank.awk)\">"
    print "##INFO=<ID=PGP_GENE,Number=1,Type=String,Description=\"SYMBOL of the ranked CSQ entry (bin/csq_rank.awk)\">"
    print "##INFO=<ID=PGP_TX,Number=1,Type=String,Description=\"Transcript set of the ranked CSQ entry: MANE, CANONICAL or ANY (bin/csq_rank.awk)\">"
    print "##INFO=<ID=PGP_FREQ,Number=1,Type=String,Description=\"Frequency tier: rare, common, unknown (no frequency for the allele) or none (no frequency field) (bin/csq_rank.awk)\">"
    print "##INFO=<ID=PGP_CALL,Number=1,Type=String,Description=\"Genotype call check: ok, weak or unknown, against ##PGP_CALL_FLOOR (bin/csq_rank.awk)\">"
    print "##INFO=<ID=PGP_CALL_WHY,Number=.,Type=String,Description=\"What made the call weak: GQ, DP or alt read fraction AB (bin/csq_rank.awk)\">"
    printf "##PGP_CALL_FLOOR=GQ=%s,DP=%s,HET_AB=%s-%s,HOM_AB=%s\n", min_gq, min_dp, het_lo, het_hi, hom_lo
    print
    next
}
/^#/ { print; next }
{
    csq = ""; info = ""
    n = split($8, kv, ";")
    for (i = 1; i <= n; i++) {
        if (substr(kv[i], 1, 4) == "CSQ=") csq = substr(kv[i], 5)
        if (kv[i] == "." || kv[i] ~ /^PGP_(IMPACT|CSQ|GENE|TX|FREQ|CALL|CALL_WHY)=/) continue
        info = info (info == "" ? "" : ";") kv[i]
    }
    jm = jc = jx = 0; bm = bc = bx = 0; known = 0; maxaf = -1
    ne = csq == "" ? 0 : split(csq, tr, ",")
    for (j = 1; j <= ne; j++) {
        split(tr[j], v, "|")
        imp = field("IMPACT")
        r = ((imp in rank) ? rank[imp] : 0) + 1
        if (field("MANE_SELECT") != "" && r > bm) { bm = r; jm = j }
        if (field("CANONICAL") == "YES" && r > bc) { bc = r; jc = j }
        if (r > bx) { bx = r; jx = j }
        for (k = 1; k <= nfreq; k++) {
            na = split(v[freq[k]], av, "&")
            for (a = 1; a <= na; a++)
                if (av[a] ~ /^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) {
                    known = 1
                    if (av[a] + 0 > maxaf) maxaf = av[a] + 0
                }
        }
    }
    if (jm) { j = jm; tx = "MANE" } else if (jc) { j = jc; tx = "CANONICAL" } else if (jx) { j = jx; tx = "ANY" } else { j = 0; tx = "NONE" }
    imp = cons = gene = ""
    if (j) {
        split(tr[j], v, "|")
        imp = field("IMPACT"); cons = field("Consequence"); gene = field("SYMBOL")
    }
    if (imp == "") imp = "NONE"
    if (cons == "") cons = "."
    if (gene == "") gene = "."
    if (nfreq == 0) ft = "none"
    else if (!known) ft = "unknown"
    else if (maxaf < rare_af + 0) ft = "rare"
    else ft = "common"
    if (keep == "noncommon" && ft == "common") next

    call = "unknown"; why = ""
    if (NF >= 10) {
        nk = split($9, fk, ":"); split($10, fv, ":")
        gt = gq = dp = ad = ""
        for (i = 1; i <= nk; i++) {
            if (fk[i] == "GT") gt = fv[i]
            else if (fk[i] == "GQ") gq = fv[i]
            else if (fk[i] == "DP") dp = fv[i]
            else if (fk[i] == "AD") ad = fv[i]
        }
        judged = 0
        if (gq ~ /^[0-9]+$/) { judged = 1; if (gq + 0 < min_gq + 0) why = add(why, "GQ" gq) }
        if (dp ~ /^[0-9]+$/) { judged = 1; if (dp + 0 < min_dp + 0) why = add(why, "DP" dp) }
        g = gt; gsub(/\|/, "/", g)
        if (ad ~ /^[0-9]+(,[0-9]+)+$/ && split(g, al, "/") == 2 && al[1] ~ /^[0-9]+$/ && al[2] ~ /^[0-9]+$/) {
            na = split(ad, av, ","); tot = 0; alt = 0
            for (a = 1; a <= na; a++) { tot += av[a]; if (a > 1) alt += av[a] }
            if (tot > 0) {
                judged = 1
                ab = alt / tot
                if ((al[1] + 0 == 0) != (al[2] + 0 == 0)) {
                    if (ab < het_lo + 0 || ab > het_hi + 0) why = add(why, "AB" sprintf("%.2f", ab))
                } else if (al[1] + 0 == al[2] + 0 && al[1] + 0 != 0) {
                    if (ab < hom_lo + 0) why = add(why, "AB" sprintf("%.2f", ab))
                }
            }
        }
        if (why != "") call = "weak"
        else if (judged) call = "ok"
    }
    tags = "PGP_IMPACT=" imp ";PGP_CSQ=" cons ";PGP_GENE=" gene ";PGP_TX=" tx ";PGP_FREQ=" ft ";PGP_CALL=" call
    if (why != "") tags = tags ";PGP_CALL_WHY=" why
    $8 = info == "" ? tags : info ";" tags
    print
}
