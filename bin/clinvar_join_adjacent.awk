#!/usr/bin/awk -f
# clinvar_join_adjacent.awk: the combined allele of two neighbouring records,
# for step 06's ClinVar match.
#
# A caller can write one complex allele (a delins) as two records a few bases
# apart on the same haplotype. Matched one by one, the parts can hit two
# single-submitter ClinVar records and miss the better-reviewed record ClinVar
# holds for the whole allele. bcftools norm -m +any joins only records at the
# same position, so it cannot build that allele. This script does, in three
# modes (-v mode=...), all reading and writing VCF text:
#
#   prepare  the sample's split, left-aligned PASS calls in; out, one record
#            per pair of neighbours that can be one allele: POS and REF span
#            both parts, the gap between them is N, and INFO/PGP_JOIN_SEQ is
#            left for `bcftools +fill-from-fasta -- -c INFO/PGP_JOIN_SEQ` to
#            fill with the reference bases under REF.
#   finish   the filled records in; out, the combined allele itself:
#            REF = the reference from the first part's start to the second
#            part's end, ALT = ALT1 + the reference bases between the parts +
#            ALT2, and INFO/PGP_JOINED_FROM naming the two parts. A record
#            whose parts disagree with the reference is dropped.
#   tag      the step 06 hits VCF in, read twice (awk -f ... HITS HITS); out,
#            the same hits with each part of a matched combined allele marked
#            PGP_PART_OF, and two hits in one gene that could be one allele
#            but matched no combined record marked PGP_ONE_ALLELE.
#
# Two records can be one allele when the second starts at most WINDOW bases
# after the first ends (WINDOW = 10; a record that overlaps the other is never
# joined), and their genotypes allow both ALT alleles on the same haplotype:
# phased (|) with the same GT in the same phase set (FORMAT/PS, or none on
# both), or unphased with the same GT (0/1 and 1/0 are the same het call).
# Two unphased het calls may also sit on different haplotypes; the genotype
# cannot tell, so bin/clinvar_hits.awk and bin/collect_summary.py say so in
# the note of such a hit. A call with a missing allele (./1) is never joined. The gene is not known
# here (the sample VCF is not annotated yet): a combined allele counts only
# when ClinVar holds a record for it, and that record names its gene.
#
# Used by scripts/06-clinvar-screen.sh and the CLINVAR_SCREEN Nextflow module
# (both inside BCFTOOLS_IMAGE); tests/test_clinvar_join_adjacent.py checks it.
BEGIN {
    FS = OFS = "\t"
    if (window == "") window = 10
    window += 0
    if (mode != "prepare" && mode != "finish" && mode != "tag") {
        print "clinvar_join_adjacent.awk: set -v mode=prepare, finish or tag" > "/dev/stderr"
        bad = 1
        exit 2
    }
}

# The phase signature of a record's genotype: two records with the same
# non-empty signature can carry their ALT alleles on the same haplotype
# (they do when phased or homozygous; unphased het calls may not).
function hapsig(fmt, smp,    nf, f, ns, s, i, gt, ps, a, na, phased, sorted) {
    nf = split(fmt, f, ":"); ns = split(smp, s, ":")
    gt = ""; ps = "."
    for (i = 1; i <= nf; i++) {
        if (f[i] == "GT" && i <= ns) gt = s[i]
        if (f[i] == "PS" && i <= ns) ps = s[i]
    }
    if (gt == "" || gt ~ /[^0-9|\/]/ || gt !~ /[1-9]/) return ""
    phased = (index(gt, "|") > 0)
    if (phased && index(gt, "/") > 0) return ""
    if (phased) return "P:" ps ":" gt
    na = split(gt, a, "/")
    if (na == 2 && a[1] + 0 > a[2] + 0) return "U:" a[2] "/" a[1]
    return "U:" gt
}

function genotype_out(sig,    p) {
    if (substr(sig, 1, 2) == "U:") return substr(sig, 3)
    p = substr(sig, 3)
    return substr(p, index(p, ":") + 1)
}

function ps_out(sig,    p) {
    p = substr(sig, 3)
    return substr(p, 1, index(p, ":") - 1)
}

function nstr(n,    s) { s = ""; while (n-- > 0) s = s "N"; return s }

function info_get(info, key,    n, kv, i) {
    n = split(info, kv, ";")
    for (i = 1; i <= n; i++) if (index(kv[i], key "=") == 1) return substr(kv[i], length(key) + 2)
    return ""
}

function info_add(info, tag) { return (info == "." || info == "") ? tag : info ";" tag }

function genes_of(info,    gi, n, g, i, s, out) {
    gi = info_get(info, "GENEINFO"); out = ""
    n = split(gi, g, "|")
    for (i = 1; i <= n; i++) { split(g[i], s, ":"); if (s[1] != "") out = out "," s[1] "," }
    return out
}

function share_gene(a, b,    n, g, i) {
    n = split(a, g, ",")
    for (i = 1; i <= n; i++) if (g[i] != "" && index(b, "," g[i] ",")) return 1
    return 0
}

function plain(allele) { return allele ~ /^[ACGTNacgtn]+$/ }

# Print the pending combined records that start before `limit` ("" = all),
# lowest POS first, in the order they were found when POS is the same.
function flush_pending(limit,    i, best, k) {
    while (1) {
        best = 0
        for (i = 1; i <= np; i++)
            if (!pdone[i] && (limit == "" || ppos[i] < limit + 0) && (best == 0 || ppos[i] < ppos[best])) best = i
        if (best == 0) break
        print pline[best]; pdone[best] = 1
    }
    k = 0
    for (i = 1; i <= np; i++) {
        if (pdone[i]) continue
        k++; ppos[k] = ppos[i]; pline[k] = pline[i]; pdone[k] = 0
    }
    for (i = k + 1; i <= np; i++) delete pdone[i]
    np = k
}

# ---- prepare -------------------------------------------------------------------
mode == "prepare" && /^##INFO=<ID=PGP_(JOIN_SEQ|JOINED_FROM|PART_OF|ONE_ALLELE),/ { next }
mode == "prepare" && /^##/ { print; next }
mode == "prepare" && /^#CHROM/ {
    print "##INFO=<ID=PGP_JOIN_SEQ,Number=1,Type=String,Description=\"Step 06, temporary: the reference bases under REF\">"
    print "##INFO=<ID=PGP_JOINED_FROM,Number=1,Type=String,Description=\"Step 06: the combined allele of these two neighbouring sample records, CHROM:POS:REF:ALT joined by |\">"
    print "##INFO=<ID=PGP_PART_OF,Number=1,Type=String,Description=\"Step 06: this hit is part of a combined allele that also matched ClinVar, CHROM:POS:REF:ALT (several joined by |)\">"
    print "##INFO=<ID=PGP_ONE_ALLELE,Number=1,Type=String,Description=\"Step 06: this hit and the one named, CHROM:POS:REF:ALT, are in one gene and can be one allele, but no combined ClinVar record matched: check ClinVar for one\">"
    print
    next
}
mode == "prepare" {
    if ($1 != chrom) { flush_pending(""); nb = 0; chrom = $1 }
    pos = $2 + 0; ref = $4; alt = $5
    if (!plain(ref) || !plain(alt) || NF < 10) next
    sig = hapsig($9, $10)
    if (sig == "") next
    end = pos + length(ref) - 1
    # Pair with every buffered record that ends at most WINDOW bases before this one.
    for (j = 1; j <= nb; j++) {
        d = pos - bend[j]
        if (d < 1 || d > window || bsig[j] != sig) continue
        gt = genotype_out(sig)
        fmt = "GT"; smp = gt
        if (substr(sig, 1, 2) == "P:" && ps_out(sig) != ".") { fmt = "GT:PS"; smp = gt ":" ps_out(sig) }
        np++
        ppos[np] = bpos[j]
        pline[np] = chrom OFS bpos[j] OFS "." OFS bref[j] nstr(d - 1) ref OFS balt[j] nstr(d - 1) alt OFS "." OFS "PASS" OFS \
              "PGP_JOIN_SEQ=.;PGP_JOINED_FROM=" chrom ":" bpos[j] ":" bref[j] ":" balt[j] "|" chrom ":" pos ":" ref ":" alt OFS \
              fmt OFS smp
    }
    # Keep only records a later one can still pair with (input is sorted by POS).
    k = 0
    for (j = 1; j <= nb; j++) {
        if (bend[j] < pos - window) continue
        k++; bpos[k] = bpos[j]; bend[k] = bend[j]; bref[k] = bref[j]; balt[k] = balt[j]; bsig[k] = bsig[j]
    }
    nb = k + 1
    bpos[nb] = pos; bend[nb] = end; bref[nb] = ref; balt[nb] = alt; bsig[nb] = sig
    # A later pair starts at a buffered record, so every pending record that
    # starts before the first of them can go out, in POS order.
    flush_pending(bpos[1])
    next
}
# ---- finish --------------------------------------------------------------------
mode == "finish" && /^##INFO=<ID=PGP_JOIN_SEQ,/ { next }
mode == "finish" && /^#/ { print; next }
mode == "finish" {
    seq = toupper(info_get($8, "PGP_JOIN_SEQ"))
    from = info_get($8, "PGP_JOINED_FROM")
    if (split(from, part, "|") != 2) next
    split(part[1], a, ":"); split(part[2], b, ":")
    r1 = toupper(a[3]); a1 = toupper(a[4]); r2 = toupper(b[3]); a2 = toupper(b[4])
    p1 = a[2] + 0; p2 = b[2] + 0
    gap = p2 - (p1 + length(r1))
    if (gap < 0 || length(seq) != p2 + length(r2) - p1) next
    if (substr(seq, 1, length(r1)) != r1 || substr(seq, p2 - p1 + 1) != r2) next
    alt = a1 substr(seq, length(r1) + 1, gap) a2
    if (alt == seq) next
    $2 = p1; $4 = seq; $5 = alt; $8 = "PGP_JOINED_FROM=" from
    print
    next
}

# ---- tag -----------------------------------------------------------------------
mode == "tag" && FNR == NR {
    if (/^#/) next
    n++
    key[n] = $1 ":" $2 ":" $4 ":" $5
    hchr[n] = $1; hpos[n] = $2 + 0; hend[n] = $2 + length($4) - 1
    hgenes[n] = genes_of($8)
    hsig[n] = (NF >= 10) ? hapsig($9, $10) : ""
    from = info_get($8, "PGP_JOINED_FROM")
    hcomb[n] = (from != "")
    if (from != "") {
        m = split(from, part, "|")
        for (i = 1; i <= m; i++) {
            # Test before assigning: mawk and busybox awk create the element
            # on the left before they evaluate the right.
            prev = (part[i] in partof) ? partof[part[i]] "|" : ""
            partof[part[i]] = prev key[n]
        }
    }
    next
}
mode == "tag" && /^#/ { print; next }
mode == "tag" {
    r++
    if (key[r] in partof) {
        $8 = info_add($8, "PGP_PART_OF=" partof[key[r]])
    } else if (!hcomb[r] && hsig[r] != "") {
        other = ""
        for (j = 1; j <= n; j++) {
            if (j == r || hcomb[j] || (key[j] in partof) || hchr[j] != hchr[r] || hsig[j] != hsig[r]) continue
            if (hpos[j] > hend[r]) d = hpos[j] - hend[r]
            else if (hpos[r] > hend[j]) d = hpos[r] - hend[j]
            else continue
            if (d > window || !share_gene(hgenes[r], hgenes[j])) continue
            other = other (other == "" ? "" : "|") key[j]
        }
        if (other != "") $8 = info_add($8, "PGP_ONE_ALLELE=" other)
    }
    print
}

END { if (!bad && mode == "prepare") flush_pending("") }
