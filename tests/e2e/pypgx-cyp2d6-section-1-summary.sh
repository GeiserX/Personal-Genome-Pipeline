#!/usr/bin/env bash
# Step 32 and the PYPGX module on the fixture, after cases 38, 51 and
# nextflow-hardening-2-nf, which ran them on the same BAM and VCF:
#   - both summaries have the columns of bin/pgx_outside_calls.py
#     pypgx-summary, and the same gene, source, status and copy-number call
#     per gene;
#   - GSTT1, which GRCh38 has only on chr22_KI270879v1_alt, is 'not callable
#     on this reference' in both (the fixture's reference has no ALT contig),
#     not FAILED, and both reports say it is not assessable on this reference;
#   - the defining variants of every called gene were looked up in the VCF;
#   - both CYP2D6 depth checks carry the MAPQ >= 20 and gene-body shares,
#     measured in the pypgx image, and agree; the verdict stays ok.
. "$(dirname "$0")/lib.sh"

B="${GENOME_DIR}/${SAMPLE}/pypgx"
N="${GENOME_DIR}/nf-hardening/${SAMPLE}/pypgx"
HEAD=$'Gene\tDiplotype\tPhenotype\tCNV_call\tSource\tStatus\tDefining_variants\tMin_GQ\tWeakest_variant\tDefining_filter\tQuality'
for side in "step 32:${B}" "PYPGX module:${N}"; do
  name=${side%%:*} dir=${side#*:}
  tsv="${dir}/${SAMPLE}_pypgx_summary.tsv"
  echo "--- ${name}: ${tsv}"; cat "$tsv" 2>/dev/null
  check_eq "${name}: the summary's header" "$(head -n 1 "$tsv" 2>/dev/null)" "$HEAD"
  check_eq "${name}: GSTT1 is not callable on this reference" \
    "$(awk -F'\t' '$1 == "GSTT1" {print $2 "|" $5 "|" $6}' "$tsv" 2>/dev/null)" "N/A|BAM|not callable on this reference"
  check_eq "${name}: no gene says FAILED for GSTT1" "$(awk -F'\t' '$1 == "GSTT1" && $2 == "FAILED"' "$tsv" 2>/dev/null | grep -c . || true)" 0
  check_eq "${name}: every called gene had its defining variants looked up" \
    "$(awk -F'\t' 'NR > 1 && $6 == "called" && ($7 == "not looked up" || $7 == "N/A")' "$tsv" 2>/dev/null | grep -c . || true)" 0
  check_ge "${name}: called genes" "$(awk -F'\t' 'NR > 1 && $6 == "called"' "$tsv" 2>/dev/null | grep -c . || true)" 1
  chk="${dir}/${SAMPLE}_cyp2d6_depth_check.tsv"
  echo "--- ${name}: ${chk}"; cat "$chk" 2>/dev/null
  val() { awk -F'\t' -v k="$1" '$1 == k {print $2}' "$chk" 2>/dev/null; }
  check_eq "${name}: the depth check verdict stays ok" "$(val status)" ok
  for k in gene_mapq20_share body_mapq1_share body_mapq20_share; do
    check "${name}: ${k} is measured ($(val "$k"))" has '^[01]\.[0-9]{2}$' "$(val "$k")"
  done
done
cols() { awk -F'\t' -v OFS='\t' '{print $1, $4, $5, $6}' "$1" 2>/dev/null; }
check_eq "the two summaries agree on gene, copy-number call, source and status" \
  "$(cols "${N}/${SAMPLE}_pypgx_summary.tsv")" "$(cols "${B}/${SAMPLE}_pypgx_summary.tsv")"
shares() { grep -E '^(gene|body)_' "$1" 2>/dev/null; }
check_eq "the two depth checks agree on the depths and shares" \
  "$(shares "${N}/${SAMPLE}_cyp2d6_depth_check.tsv")" "$(shares "${B}/${SAMPLE}_cyp2d6_depth_check.tsv")"

TXT="${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.txt"
HTML="${GENOME_DIR}/nf-hardening/${SAMPLE}/${SAMPLE}_report.html"
check "the text report (case 51) says GSTT1 is not assessable on this reference" \
  grep -q 'Not assessable on this reference: GSTT1' "$TXT"
check "the Nextflow HTML report says GSTT1 is not assessable on this reference" \
  grep -q 'Not assessable on this reference: GSTT1' "$HTML"
check "the text report has the PGx copy number block" grep -q '^## PGx copy number$' "$TXT"

finish
