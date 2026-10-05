#!/usr/bin/env bash
# Step 33: sample identity and contamination (somalier, VerifyBamID2)
# Input: sorted BAM with .bai (step 02), the somalier sites file and the
#        VerifyBamID2 marker panel (setup.sh installs both)
# Output: ${SAMPLE}/qc/${SAMPLE}_sample_qc.tsv, a key/value table: the sex
#         somalier infers from the reads, the declared sex, VerifyBamID2's
#         FREEMIX (the estimated share of reads from another person) and its
#         verdict; somalier's and VerifyBamID2's own files beside it
#
# Usage: ./scripts/33-sample-qc.sh <sample_name> [declared_sex: male|female]
#
# Two things go wrong with a sample before any result is read: it is not the
# person you think (a swap at the lab or in your own files), or another
# person's DNA is mixed in. With a declared sex, the step exits non-zero when
# somalier's sex, from heterozygosity at chrX sites and from chrY depth,
# differs; SEX_CHECK=warn prints the mismatch and exits 0, as in step 16.
# When somalier cannot tell (too few chrX sites with reads, as on a sliced or
# targeted BAM), it says so and does not stop. FREEMIX above FREEMIX_WARN
# (default 0.03) is a warning only: contamination makes calls less reliable,
# it does not make them someone else's.
#
# Env: SEX_CHECK          fail (default) or warn
#      FREEMIX_WARN       FREEMIX above this is reported as possible contamination (0.03)
#      SOMALIER_SITES     sites VCF (default reference/somalier/sites.hg38.vcf.gz)
#      VERIFYBAMID2_PANEL panel prefix, the path without .UD/.mu/.bed (default
#                         reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat)
#      Both paths must lie inside GENOME_DIR, the only directory containers see.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> [declared_sex: male|female]}
DECLARED_SEX=${2:-}
SEX_CHECK=${SEX_CHECK:-fail}
FREEMIX_WARN=${FREEMIX_WARN:-0.03}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image SOMALIER_IMAGE VERIFYBAMID2_IMAGE PYTHON_IMAGE
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
OUT="${SAMPLE_DIR}/qc"
# setup.sh installs these two under these names.
SITES=${SOMALIER_SITES:-${GENOME_DIR}/reference/somalier/sites.hg38.vcf.gz}
PANEL=${VERIFYBAMID2_PANEL:-${GENOME_DIR}/reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat}
case "$SITES" in /*) ;; *) SITES="${GENOME_DIR}/${SITES}" ;; esac
case "$PANEL" in /*) ;; *) PANEL="${GENOME_DIR}/${PANEL}" ;; esac

echo "=== Step 33: sample identity and contamination: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Somalier sites: ${SITES}"
echo "VerifyBamID2 panel: ${PANEL}"
echo "Output: ${OUT}/"

case "$DECLARED_SEX" in
  ""|male|female) ;;
  *) echo "ERROR: declared sex must be 'male' or 'female', got '${DECLARED_SEX}'" >&2; exit 1 ;;
esac
case "$SEX_CHECK" in
  fail|warn) ;;
  *) echo "ERROR: SEX_CHECK must be 'fail' or 'warn', got '${SEX_CHECK}'" >&2; exit 1 ;;
esac
if ! awk -v v="$FREEMIX_WARN" 'BEGIN { exit !(v ~ /^[0-9]*\.?[0-9]+$/ && v + 0 < 1) }'; then
  echo "ERROR: FREEMIX_WARN must be a fraction between 0 and 1, got '${FREEMIX_WARN}'" >&2
  exit 1
fi

for f in "$BAM" "${BAM}.bai" "$REF_FASTA" "${REF_FASTA}.fai" "$SITES" "${PANEL}.UD" "${PANEL}.mu" "${PANEL}.bed"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    case "$f" in
      "$SITES"|"$PANEL".*) echo "  Install it with: ./scripts/setup.sh ${GENOME_DIR}" >&2 ;;
    esac
    exit 1
  fi
done
SITES_C=$(cpath "$SITES") || exit 1
PANEL_C=$(cpath "$PANEL") || exit 1

mkdir -p "${OUT}/somalier" "${OUT}/verifybamid2"
# Never read a file an earlier run left: each tool writes its own anew.
rm -f "${OUT}/somalier/"*.somalier "${OUT}/somalier/${SAMPLE}".* "${OUT}/verifybamid2/${SAMPLE}".* \
  "${OUT}/${SAMPLE}_sample_qc.tsv"

# --- 1. somalier: genotypes and depth at known sites -------------------------------
# extract names its file after the @RG SM tag of the BAM, which is the name
# relate reports the sample under.
echo ""
echo "--- somalier extract"
run_in --cpus 2 --memory 4g "$SOMALIER_IMAGE" \
  somalier extract \
    -d "/genome/${SAMPLE}/qc/somalier" \
    --sites "$SITES_C" \
    -f "$REF_FASTA_C" \
    "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
EXTRACTED=("${OUT}/somalier/"*.somalier)
if [ "${#EXTRACTED[@]}" -ne 1 ] || [ ! -s "${EXTRACTED[0]}" ]; then
  echo "ERROR: somalier extract wrote ${#EXTRACTED[@]} files in ${OUT}/somalier/, expected one" >&2
  exit 1
fi
SOMALIER_ID=$(basename "${EXTRACTED[0]}" .somalier)
echo "somalier reads the sample as '${SOMALIER_ID}' (the BAM's @RG SM)"

echo "--- somalier relate"
run_in --cpus 1 --memory 2g "$SOMALIER_IMAGE" \
  somalier relate --infer \
    --sites "$SITES_C" \
    -o "/genome/${SAMPLE}/qc/somalier/${SAMPLE}" \
    "/genome/${SAMPLE}/qc/somalier/${SOMALIER_ID}.somalier"
for f in samples pairs; do
  [ -s "${OUT}/somalier/${SAMPLE}.${f}.tsv" ] || { echo "ERROR: somalier relate wrote no ${SAMPLE}.${f}.tsv" >&2; exit 1; }
done

# --- 2. VerifyBamID2: contamination ----------------------------------------------
# VerifyBamID2 refuses to estimate when fewer than 1,000 panel markers have
# reads (a targeted or sliced BAM). Then it runs again without that check,
# and the table says on how many markers FREEMIX rests.
echo ""
echo "--- VerifyBamID2"
VB_OUT="/genome/${SAMPLE}/qc/verifybamid2/${SAMPLE}"
VB_LOG="${OUT}/verifybamid2/${SAMPLE}.log"
vb2() {
  run_in --cpus "$THREADS" --memory 4g "$VERIFYBAMID2_IMAGE" \
    verifybamid2 \
      --SVDPrefix "$PANEL_C" \
      --Reference "$REF_FASTA_C" \
      --BamFile "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" \
      --Output "$VB_OUT" \
      --NumThread "$THREADS" "$@" > "$VB_LOG" 2>&1
}
SANITY_CHECK=passed
if ! vb2; then
  if grep -q 'Insufficient Available markers' "$VB_LOG"; then
    SANITY_CHECK=skipped
    echo "Fewer than 1,000 panel markers have reads; running again without VerifyBamID2's marker check."
    vb2 --DisableSanityCheck || { tail -n 20 "$VB_LOG" >&2; echo "ERROR: VerifyBamID2 failed" >&2; exit 1; }
  else
    tail -n 20 "$VB_LOG" >&2
    echo "ERROR: VerifyBamID2 failed (log: ${VB_LOG})" >&2
    exit 1
  fi
fi
[ -s "${OUT}/verifybamid2/${SAMPLE}.selfSM" ] || { echo "ERROR: VerifyBamID2 wrote no ${SAMPLE}.selfSM" >&2; exit 1; }
tail -n 3 "$VB_LOG"

# --- 3. The verdict ---------------------------------------------------------------
# bin/collect_summary.py holds the rule, shared with the Nextflow SAMPLE_QC
# process and read again by the reports.
echo ""
echo "--- verdict"
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYTHON_IMAGE" \
  python3 /pgp-bin/collect_summary.py sample-qc \
    --sample "$SAMPLE" \
    --somalier-id "$SOMALIER_ID" \
    --somalier-samples "/genome/${SAMPLE}/qc/somalier/${SAMPLE}.samples.tsv" \
    --somalier-pairs "/genome/${SAMPLE}/qc/somalier/${SAMPLE}.pairs.tsv" \
    --selfsm "/genome/${SAMPLE}/qc/verifybamid2/${SAMPLE}.selfSM" \
    --declared-sex "$DECLARED_SEX" \
    --freemix-warn "$FREEMIX_WARN" \
    --marker-check "$SANITY_CHECK" \
    --out "/genome/${SAMPLE}/qc/${SAMPLE}_sample_qc.tsv"
TABLE="${OUT}/${SAMPLE}_sample_qc.tsv"
[ -s "$TABLE" ] || { echo "ERROR: no ${TABLE}" >&2; exit 1; }
val() { awk -F'\t' -v k="$1" '$1 == k { print $2; exit }' "$TABLE"; }

echo ""
echo "=== Step 33 complete ==="
echo "Results: ${TABLE}"
echo "  somalier:     ${OUT}/somalier/${SAMPLE}.samples.tsv, ${SAMPLE}.html"
echo "  VerifyBamID2: ${OUT}/verifybamid2/${SAMPLE}.selfSM"
echo ""
echo "Contamination: FREEMIX $(val freemix) (VerifyBamID2's marker check ${SANITY_CHECK})"
if [ "$(val contamination)" = warn ]; then
  echo "!!! FREEMIX $(val freemix) is above ${FREEMIX_WARN}: about $(awk -v f="$(val freemix)" 'BEGIN { printf "%d", f * 100 + 0.5 }')% of the reads may come" >&2
  echo "!!! from another person. Calls, above all heterozygous ones, are less reliable; see docs/33-sample-qc.md." >&2
fi
echo "Sex from the reads (somalier): $(val inferred_sex)"
if [ -n "$(val same_person_as)" ]; then
  echo "!!! The same person as: $(val same_person_as)" >&2
fi
case "$(val sex_check)" in
  ok) echo "  Declared sex:  ${DECLARED_SEX}"; echo "  Sex check: OK" ;;
  not_checked) echo "  Sex check: not done ($(val sex_check_reason))" ;;
  mismatch)
    echo "" >&2
    echo "!!! SEX CHECK MISMATCH: declared ${DECLARED_SEX}, somalier infers $(val inferred_sex) from the reads" >&2
    echo "!!! (chrX sites $(val x_sites): $(val x_het) heterozygous, $(val x_hom_alt) homozygous ALT; chrY depth ratio $(val y_depth_ratio))." >&2
    echo "!!! Either the sample is not who you think it is, the declared sex is wrong, or the" >&2
    echo "!!! sample has a sex-chromosome aneuploidy. Steps that use the sex (DeepVariant, ExpansionHunter) would get the wrong value." >&2
    if [ "$SEX_CHECK" = warn ]; then
      echo "!!! SEX_CHECK=warn: continuing." >&2
    else
      echo "!!! Set SEX_CHECK=warn to continue anyway." >&2
      exit 1
    fi ;;
  *) echo "ERROR: ${TABLE} has no sex_check value" >&2; exit 1 ;;
esac
