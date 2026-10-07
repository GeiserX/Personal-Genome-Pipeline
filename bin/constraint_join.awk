#!/usr/bin/awk -f
# constraint_join.awk: append gnomAD v4.1 gene constraint columns to a TSV.
#
# The one constraint loader of the pipeline: steps 23 and 31 run it on the host
# and the SLIVAR_PRIORITIZE Nextflow module runs it in the bcftools image (no
# python there), so the three cannot disagree about a gene.
#
#   awk -f constraint_join.awk gene_col=NAME [constrained=1] CONSTRAINT.tsv TABLE.tsv
#
# CONSTRAINT.tsv  gnomad.v4.1.constraint_metrics.tsv: columns found by header
#                 name (gene, canonical, transcript, lof.oe_ci.upper, lof.pLI,
#                 mis.z_score). Only canonical rows count; v4.1 lists an Ensembl
#                 and a RefSeq canonical row per gene, and the Ensembl (ENST)
#                 one wins wherever it sits in the file.
# TABLE.tsv       has a header row; gene_col names its gene-symbol column.
#
# Output: TABLE with LOEUF, pLI and mis_z appended (and CONSTRAINED, with
# constrained=1: YES when LOEUF < 0.35 or pLI > 0.9, NO when either value is
# known and neither says so, "." when the gene has neither). A missing value is ".".
# Exit 3: a needed column is missing. Exit 4: rows carry a gene symbol but not
# one matched the table, which means a wrong or broken constraint file.
# A summary line goes to stderr.
BEGIN { FS = OFS = "\t" }
NR == FNR {
    if (FNR == 1) {
        for (i = 1; i <= NF; i++) col[$i] = i
        n = split("gene canonical transcript lof.oe_ci.upper lof.pLI mis.z_score", need, " ")
        for (k = 1; k <= n; k++) if (!(need[k] in col)) {
            print "ERROR: column " need[k] " missing from the constraint table" > "/dev/stderr"
            bad = 3; exit 3
        }
        next
    }
    if ($col["canonical"] != "true") next
    g = $col["gene"]
    if (g == "") next
    ens = ($col["transcript"] ~ /^ENST/)
    if ((g in val) && (src[g] || !ens)) next
    l = $col["lof.oe_ci.upper"]; p = $col["lof.pLI"]; m = $col["mis.z_score"]
    if (l == "NA" || l == "") l = "."
    if (p == "NA" || p == "") p = "."
    if (m == "NA" || m == "") m = "."
    val[g] = l OFS p OFS m
    loeuf[g] = l; pli[g] = p; src[g] = ens
    next
}
FNR == 1 {
    gc = 0
    for (i = 1; i <= NF; i++) if ($i == gene_col) gc = i
    if (!gc) {
        print "ERROR: the table has no column named '" gene_col "'" > "/dev/stderr"
        bad = 3; exit 3
    }
    print $0, "LOEUF", "pLI", "mis_z" (constrained ? OFS "CONSTRAINED" : "")
    next
}
{
    rows++
    g = $gc
    if (g != "." && g != "") with_gene++
    c = "."
    if (g in val) {
        matched[g] = 1
        if (loeuf[g] != "." || pli[g] != ".") c = "NO"
        if ((loeuf[g] != "." && loeuf[g] + 0 < 0.35) || (pli[g] != "." && pli[g] + 0 > 0.9)) c = "YES"
        out = val[g]
    } else {
        out = "." OFS "." OFS "."
    }
    print $0, out (constrained ? OFS c : "")
}
END {
    if (bad) exit bad
    n = 0
    for (g in matched) n++
    printf "gnomAD constraint: %d of %d rows carry a gene symbol; %d distinct genes matched\n", with_gene, rows, n > "/dev/stderr"
    if (with_gene > 0 && n == 0) {
        print "ERROR: no gene in the table matched the constraint file; check that it is the gnomAD v4.1 constraint table" > "/dev/stderr"
        exit 4
    }
}
