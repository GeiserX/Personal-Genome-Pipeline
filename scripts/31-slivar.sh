#!/usr/bin/env bash
# 31-slivar.sh — Variant prioritization and compound heterozygote detection
# Usage: ./scripts/31-slivar.sh <sample_name>
#
# Prioritizes clinically interesting variants using tiered filters and detects
# compound heterozygote candidates using slivar. Optionally annotates results
# with gnomAD gene constraint metrics (LOEUF, pLI).
#
# Rarity is the rule of step 23 (bin/csq_rank.awk, PGP_FREQ): a known MAX_AF
# (else gnomADe_AF and gnomADg_AF) below 1%. A HIGH or MODERATE variant with
# no frequency at all is not rare: it goes into its own freq_unknown tier,
# which the prioritized VCF and the compound-het search still include, and
# the TSVs say FREQ_TIER=unknown. The tiers take the worst consequence; the
# TSVs add the MANE Select consequence (MANE_Consequence) and the gnomAD AF.
#
# Each compound-het pair is phased from the sample's BAM (bin/comphet_phase.py,
# pysam from PYPGX_IMAGE): cis, trans or unknown with the reads that say so,
# in <sample>_compound_het_pairs.tsv and the PHASE column of the per-variant
# table. Without the BAM every pair is unknown.
#
# Input: vcfanno-enriched VCF (step 30) or VEP-annotated VCF (step 13), and
#   aligned/<sample>_sorted.bam (ALIGN_DIR) with its .bai for the phase
# Output: prioritized VCF, compound het TSV and pairs TSV, summary TSV
# Requires: VEP-annotated VCF. Step 30 (vcfanno) recommended for full filtering.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"


SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
OUTDIR="${SAMPLE_DIR}/slivar"

# Input VCF: prefer vcfanno-enriched (step 30), fall back to VEP (step 13)
ANNOTATED_VCF="${SAMPLE_DIR}/vep/${SAMPLE}_annotated.vcf.gz"
VEP_VCF_GZ="${SAMPLE_DIR}/vep/${SAMPLE}_vep.vcf.gz"
VEP_VCF="${SAMPLE_DIR}/vep/${SAMPLE}_vep.vcf"

# Optional gene constraint data
CONSTRAINT_TSV="${GENOME_DIR}/annotations/gnomad_v4.1_constraint.tsv"

# A derived file is used only when it is newer than what it was built from,
# so a re-run of step 13 is never hidden behind an older _vep.vcf.gz or an
# older step 30 output. An empty file (an interrupted write) never counts.
VEP_SRC=""
if [ -s "$VEP_VCF" ] && { [ ! -s "$VEP_VCF_GZ" ] || [ "$VEP_VCF" -nt "$VEP_VCF_GZ" ]; }; then
  VEP_SRC="$VEP_VCF"
elif [ -s "$VEP_VCF_GZ" ]; then
  VEP_SRC="$VEP_VCF_GZ"
fi
INPUT=""
HAS_VCFANNO=0
if [ -s "$ANNOTATED_VCF" ] && { [ -z "$VEP_SRC" ] || [ "$ANNOTATED_VCF" -nt "$VEP_SRC" ]; }; then
  INPUT="$ANNOTATED_VCF"
  HAS_VCFANNO=1
elif [ -n "$VEP_SRC" ]; then
  INPUT="$VEP_SRC"
  if [ -f "$ANNOTATED_VCF" ]; then
    echo "NOTICE: ${ANNOTATED_VCF} is older than ${VEP_SRC}; ignoring it. Run step 30 again for the score filters."
  fi
else
  echo "ERROR: No annotated VCF found. Run step 13 (VEP) first."
  echo "  Expected: ${ANNOTATED_VCF} or ${VEP_VCF_GZ} or ${VEP_VCF}"
  exit 1
fi

CONTAINER_INPUT="/genome/${SAMPLE}/vep/$(basename "$INPUT")"

echo "============================================"
echo "  Step 31: Variant Prioritization (slivar)"
echo "  Sample: ${SAMPLE}"
echo "  Input:  ${INPUT}"
echo "  vcfanno annotations: $([ "$HAS_VCFANNO" -eq 1 ] && echo 'yes (full filtering)' || echo 'no (VEP-only mode)')"
echo "  Output: ${OUTDIR}/"
echo "============================================"
echo ""

mkdir -p "$OUTDIR"
# The report reads the summary; a failed run must not leave an older one behind.
rm -f "${OUTDIR}/${SAMPLE}_slivar_summary.tsv"

# ── Step 1: Ensure VCF is bgzipped and indexed ────────────────────────
echo "[1/5] Preparing input VCF..."
if [ "$INPUT" = "$VEP_VCF" ]; then
  echo "  Compressing VEP VCF..."
  rm -f "${VEP_VCF_GZ:?}" "${VEP_VCF_GZ:?}.tbi"
  run_in --cpus 2 --memory 2g \
    "${BCFTOOLS_IMAGE}" \
    bash -c "bcftools view /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf -Oz \
      -o /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz.tmp && \
      mv /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz.tmp /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz && \
      bcftools index -f -t /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz"
  INPUT="$VEP_VCF_GZ"
  CONTAINER_INPUT="/genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz"
  echo "  Done."
elif [ ! -f "${INPUT}.tbi" ] || [ "$INPUT" -nt "${INPUT}.tbi" ]; then
  echo "  Indexing VCF..."
  run_in --cpus 2 --memory 2g \
    "${BCFTOOLS_IMAGE}" \
    bcftools index -f -t "$CONTAINER_INPUT"
  echo "  Done."
else
  echo "  VCF already compressed and indexed."
fi

# ── Step 2: Detect available annotation fields ────────────────────────
echo ""
echo "[2/5] Detecting annotation fields..."

VEP_FIELDS=$(run_in \
  --cpus 2 --memory 2g \
  "${BCFTOOLS_IMAGE}" \
  bcftools +split-vep -l "$CONTAINER_INPUT" 2>/dev/null || echo "")

if [ -z "$VEP_FIELDS" ]; then
  echo "ERROR: No CSQ annotation found in VCF. Was VEP step 13 run correctly?"
  exit 1
fi

HAS_CLINVAR=0
printf '%s\n' "$VEP_FIELDS" | awk -F'\t' '$2 == "CLIN_SIG" {f = 1} END {exit !f}' && HAS_CLINVAR=1

# Rarity: VEP's MAX_AF (highest frequency in any 1000 Genomes or gnomAD
# exome/genome population), or gnomADe_AF and gnomADg_AF when MAX_AF is
# absent. Exome frequency alone would call a variant common in genomes but
# absent from exomes rare. bin/csq_rank.awk applies the rule, the same as
# step 23; this only names the fields.
has_csq() { printf '%s\n' "$VEP_FIELDS" | awk -F'\t' -v f="$1" '$2 == f {x = 1} END {exit !x}'; }
FREQ_NAME=""
if has_csq MAX_AF; then
  FREQ_NAME="MAX_AF"
else
  for f in gnomADe_AF gnomADg_AF; do
    if has_csq "$f"; then FREQ_NAME="${FREQ_NAME:+${FREQ_NAME} and }${f}"; fi
  done
fi

# Check for vcfanno INFO fields
HAS_CADD=0
HAS_CADD_INDEL=0
HAS_REVEL=0
HAS_AM=0
HAS_SPLICEAI=0
HAS_SPLICEAI_INDEL=0
if [ "$HAS_VCFANNO" -eq 1 ]; then
  INFO_HEADER=$(run_in \
    --cpus 2 --memory 2g \
    "${BCFTOOLS_IMAGE}" \
    bcftools view -h "$CONTAINER_INPUT" 2>/dev/null | grep '^##INFO' || echo "")
  grep -q 'ID=CADD_PHRED,' <<< "$INFO_HEADER" && HAS_CADD=1
  grep -q 'ID=CADD_PHRED_indel,' <<< "$INFO_HEADER" && HAS_CADD_INDEL=1
  grep -q 'ID=REVEL' <<< "$INFO_HEADER" && HAS_REVEL=1
  grep -q 'ID=AM_class' <<< "$INFO_HEADER" && HAS_AM=1
  grep -q 'ID=SpliceAI,' <<< "$INFO_HEADER" && HAS_SPLICEAI=1
  grep -q 'ID=SpliceAI_indel,' <<< "$INFO_HEADER" && HAS_SPLICEAI_INDEL=1
fi

echo "  VEP CSQ fields: frequency=${FREQ_NAME:-none}, ClinVar=$([ "$HAS_CLINVAR" -eq 1 ] && echo 'yes' || echo 'no')"
if [ -z "$FREQ_NAME" ]; then
  echo "  NOTICE: no MAX_AF, gnomADe_AF or gnomADg_AF in the VEP output: the rare tiers are not filtered by frequency."
fi
echo "  vcfanno INFO fields: CADD=$([ "$HAS_CADD" -eq 1 ] && echo 'yes' || echo 'no'), REVEL=$([ "$HAS_REVEL" -eq 1 ] && echo 'yes' || echo 'no'), AlphaMissense=$([ "$HAS_AM" -eq 1 ] && echo 'yes' || echo 'no'), SpliceAI=$([ "$HAS_SPLICEAI" -eq 1 ] && echo 'yes' || echo 'no')"
echo ""

# ── Step 3: Variant prioritization with bcftools ──────────────────────
# Uses bcftools +split-vep for CSQ fields and bcftools view -i for INFO fields.
# Three filter tiers: rare_high, rare_moderate_deleterious, clinvar_pathogenic.
echo "[3/5] Prioritizing variants..."

# rank KEEP: bin/csq_rank.awk on stdin, on the host (PGP_FREQ, the MANE
# Select consequence PGP_CSQ and PGP_TX); KEEP=noncommon drops the records
# with a known frequency of 1% or more.
rank() { LC_ALL=C awk -v keep="$1" -f "${PGP_ROOT}/bin/csq_rank.awk"; }
S="/genome/${SAMPLE}/slivar"
RANKED="${S}/${SAMPLE}_ranked.vcf.gz"
run_in --cpus 2 --memory 2g "${BCFTOOLS_IMAGE}" bcftools view -f PASS "${CONTAINER_INPUT}" \
  | rank noncommon \
  | run_in -i --cpus 2 --memory 2g "${BCFTOOLS_IMAGE}" bcftools view -Oz -o "$RANKED"
run_in --cpus 2 --memory 2g "${BCFTOOLS_IMAGE}" bcftools index -f -t "$RANKED"
# Rare: a known frequency below 1%, or any record when the VEP output has no
# frequency field at all.
RARE_FREQ='(INFO/PGP_FREQ="rare" || INFO/PGP_FREQ="none")'

# --- Filter 1: rare_high ---
# PASS + HIGH VEP impact + rare
echo "  [a] rare_high: PASS + HIGH impact + ${FREQ_NAME:-no frequency filter}..."
run_in --cpus 2 --memory 4g \
  "${BCFTOOLS_IMAGE}" \
  bash -o pipefail -c "bcftools +split-vep ${RANKED} -c IMPACT -s worst \
      -i 'IMPACT=\"HIGH\"' | \
    bcftools view -i '${RARE_FREQ}' \
      -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_rare_high.vcf.gz && \
    bcftools index -f -t /genome/${SAMPLE}/slivar/${SAMPLE}_rare_high.vcf.gz"

RARE_HIGH_COUNT=$(run_in \
  --cpus 2 --memory 2g \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -H "/genome/${SAMPLE}/slivar/${SAMPLE}_rare_high.vcf.gz" 2>/dev/null | wc -l || echo 0)
echo "      Found: ${RARE_HIGH_COUNT} variants"

# --- Filter 2: rare_moderate_deleterious ---
# PASS + MODERATE impact + gnomAD AF < 0.01 + at least one deleterious predictor
echo "  [b] rare_moderate_deleterious: PASS + MODERATE + rare + deleterious predictors..."

# The rare MODERATE records; the predictors below narrow them.
MODERATE_SEL="bcftools +split-vep ${RANKED} -c IMPACT -s worst -i 'IMPACT=\"MODERATE\"' | bcftools view -i '${RARE_FREQ}'"

# If vcfanno annotations are available, add predictor thresholds as a second pass
if [ "$HAS_VCFANNO" -eq 1 ]; then
  # First pass: extract rare MODERATE via split-vep, then second pass: filter on INFO fields

  # Build INFO-level predictor filter — only reference tags that exist in the header
  PREDICTOR_PARTS=()
  if [ "$HAS_CADD" -eq 1 ] && [ "$HAS_CADD_INDEL" -eq 1 ]; then
    PREDICTOR_PARTS+=('INFO/CADD_PHRED>=20 || INFO/CADD_PHRED_indel>=20')
  elif [ "$HAS_CADD" -eq 1 ]; then
    PREDICTOR_PARTS+=('INFO/CADD_PHRED>=20')
  elif [ "$HAS_CADD_INDEL" -eq 1 ]; then
    PREDICTOR_PARTS+=('INFO/CADD_PHRED_indel>=20')
  fi
  [ "$HAS_REVEL" -eq 1 ] && PREDICTOR_PARTS+=('INFO/REVEL>=0.5')
  [ "$HAS_AM" -eq 1 ] && PREDICTOR_PARTS+=('INFO/AM_class="likely_pathogenic"')
  # Note: SpliceAI is a pipe-delimited string, not a numeric field.
  # bcftools cannot numerically compare sub-fields, so we use presence check only.
  # Variants with SpliceAI annotations are included; threshold filtering happens in step 23.
  if [ "$HAS_SPLICEAI" -eq 1 ] && [ "$HAS_SPLICEAI_INDEL" -eq 1 ]; then
    PREDICTOR_PARTS+=('INFO/SpliceAI!="." || INFO/SpliceAI_indel!="."')
  elif [ "$HAS_SPLICEAI" -eq 1 ]; then
    PREDICTOR_PARTS+=('INFO/SpliceAI!="."')
  elif [ "$HAS_SPLICEAI_INDEL" -eq 1 ]; then
    PREDICTOR_PARTS+=('INFO/SpliceAI_indel!="."')
  fi

  if [ ${#PREDICTOR_PARTS[@]} -gt 0 ]; then
    PREDICTOR_EXPR=$(printf ' || %s' "${PREDICTOR_PARTS[@]}")
    PREDICTOR_EXPR="${PREDICTOR_EXPR:4}"  # strip leading ' || '

    run_in --cpus 2 --memory 4g \
      "${BCFTOOLS_IMAGE}" \
      bash -o pipefail -c "${MODERATE_SEL} | \
        bcftools view -i '${PREDICTOR_EXPR}' \
          -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz && \
        bcftools index -t /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz"
  else
    # No predictors available despite vcfanno — fall back to all rare MODERATE
    run_in --cpus 2 --memory 4g \
      "${BCFTOOLS_IMAGE}" \
      bash -o pipefail -c "${MODERATE_SEL} \
          -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz && \
        bcftools index -t /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz"
  fi
else
  # No vcfanno — include all rare MODERATE (same as step 23 behavior)
  echo "    WARNING: No vcfanno annotations — including all rare MODERATE variants."
  echo "    Run step 30 (vcfanno) for CADD/REVEL/AlphaMissense/SpliceAI filtering."

  run_in --cpus 2 --memory 4g \
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "${MODERATE_SEL} \
        -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz && \
      bcftools index -t /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz"
fi

MODERATE_DEL_COUNT=$(run_in \
  --cpus 2 --memory 2g \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -H "/genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz" 2>/dev/null | wc -l || echo 0)
echo "      Found: ${MODERATE_DEL_COUNT} variants"

# --- Filter 2b: freq_unknown ---
# PASS + HIGH or MODERATE + no known frequency: not rare (an indel VEP could
# not match to gnomAD may be common), kept for the compound-het search.
FREQ_UNKNOWN_COUNT=0
FREQ_UNKNOWN_FILE=""
if [ -n "$FREQ_NAME" ]; then
  echo "  [b2] freq_unknown: PASS + HIGH or MODERATE + no known frequency (not rare)..."
  run_in --cpus 2 --memory 4g \
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools +split-vep ${RANKED} -c IMPACT -s worst \
        -i 'IMPACT=\"HIGH\" || IMPACT=\"MODERATE\"' | \
      bcftools view -i 'INFO/PGP_FREQ=\"unknown\"' \
        -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_freq_unknown.vcf.gz && \
      bcftools index -f -t /genome/${SAMPLE}/slivar/${SAMPLE}_freq_unknown.vcf.gz"
  FREQ_UNKNOWN_COUNT=$(run_in \
    --cpus 2 --memory 2g \
    "${BCFTOOLS_IMAGE}" \
    bcftools view -H "/genome/${SAMPLE}/slivar/${SAMPLE}_freq_unknown.vcf.gz" 2>/dev/null | wc -l || echo 0)
  echo "      Found: ${FREQ_UNKNOWN_COUNT} variants"
  FREQ_UNKNOWN_FILE="/genome/${SAMPLE}/slivar/${SAMPLE}_freq_unknown.vcf.gz"
fi

# --- Filter 3: clinvar_pathogenic ---
CLINVAR_COUNT=0
CLINVAR_FILE=""
if [ "$HAS_CLINVAR" -eq 1 ]; then
  echo "  [c] clinvar_pathogenic: ClinVar pathogenic/likely_pathogenic..."
  # Ranked too, so the TSVs show its frequency tier and MANE consequence.
  run_in --cpus 2 --memory 4g \
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
      bcftools +split-vep - -c CLIN_SIG \
        -i 'CLIN_SIG~\"pathogenic\" && CLIN_SIG!~\"conflicting\"'" \
    | rank all \
    | run_in -i --cpus 2 --memory 2g "${BCFTOOLS_IMAGE}" \
      bcftools view -Oz -o "/genome/${SAMPLE}/slivar/${SAMPLE}_clinvar_path.vcf.gz"
  run_in --cpus 2 --memory 2g "${BCFTOOLS_IMAGE}" \
    bcftools index -f -t "/genome/${SAMPLE}/slivar/${SAMPLE}_clinvar_path.vcf.gz"

  CLINVAR_COUNT=$(run_in \
    --cpus 2 --memory 2g \
    "${BCFTOOLS_IMAGE}" \
    bcftools view -H "/genome/${SAMPLE}/slivar/${SAMPLE}_clinvar_path.vcf.gz" 2>/dev/null | wc -l || echo 0)
  echo "      Found: ${CLINVAR_COUNT} variants"
  CLINVAR_FILE="/genome/${SAMPLE}/slivar/${SAMPLE}_clinvar_path.vcf.gz"
else
  echo "  [c] clinvar_pathogenic: skipped (CLIN_SIG not in VEP annotations)"
fi

# Merge all tiers into a single prioritized VCF
echo ""
echo "  Merging filter tiers into prioritized VCF..."
MERGE_FILES="/genome/${SAMPLE}/slivar/${SAMPLE}_rare_high.vcf.gz /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz"
[ -n "$FREQ_UNKNOWN_FILE" ] && MERGE_FILES="${MERGE_FILES} ${FREQ_UNKNOWN_FILE}"
[ -n "$CLINVAR_FILE" ] && MERGE_FILES="${MERGE_FILES} ${CLINVAR_FILE}"

run_in --cpus 2 --memory 4g \
  "${BCFTOOLS_IMAGE}" \
  bash -o pipefail -c "bcftools concat -a -D \
    ${MERGE_FILES} | \
    bcftools sort -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_prioritized.vcf.gz && \
    bcftools index -t /genome/${SAMPLE}/slivar/${SAMPLE}_prioritized.vcf.gz"

PRIORITIZED_COUNT=$(run_in \
  --cpus 2 --memory 2g \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -H "/genome/${SAMPLE}/slivar/${SAMPLE}_prioritized.vcf.gz" 2>/dev/null | wc -l || echo 0)
echo "  Total prioritized variants: ${PRIORITIZED_COUNT}"

# Columns both TSVs end with: the frequency fields VEP wrote ('.' for a field
# it did not), the frequency tier, the transcript set step 23 ranks on and
# the MANE Select consequence ('.' when the record has no MANE transcript).
EXTRA_HEAD=$'MAX_AF\tgnomADe_AF\tgnomADg_AF\tFREQ_TIER\tTRANSCRIPT\tMANE_Consequence'
EXTRA_FMT=""
for f in MAX_AF gnomADe_AF gnomADg_AF; do
  if has_csq "$f"; then EXTRA_FMT+="\\t%${f}"; else EXTRA_FMT+="\\t."; fi
done
EXTRA_FMT+='\t%INFO/PGP_FREQ\t%INFO/PGP_TX\t%INFO/PGP_CSQ'
# mane_col: the last column (PGP_CSQ) is the MANE consequence only when the
# column before it says MANE.
mane_col() { awk -F'\t' 'BEGIN {OFS = "\t"} $(NF - 1) != "MANE" {$NF = "."} {print}'; }

# ── Step 4: Compound heterozygote detection with slivar ───────────────
echo ""
echo "[4/5] Detecting compound heterozygote candidates..."

# slivar compound-hets requires:
#   --vcf: input VCF with CSQ annotations
#   --ped: PED file describing sample relationships
#   --allow-non-trios: required for singleton/duo samples (no trio structure)
# It groups heterozygous variants by gene (from CSQ) and outputs VCF to stdout
# with pairs of variants per gene. From one sample's genotypes these are
# candidates: below, each pair is phased from the reads.
COMPHET_VCF="${OUTDIR}/${SAMPLE}_compound_hets.vcf.gz"
COMPHET_TSV="${OUTDIR}/${SAMPLE}_compound_hets.tsv"
COMPHET_PAIRS_TSV="${OUTDIR}/${SAMPLE}_compound_het_pairs.tsv"
COMPHET_PED="${OUTDIR}/${SAMPLE}.ped"

# Generate minimal PED file for single sample (no parents, unaffected)
# Format: family_id sample_id father mother sex phenotype
echo -e "${SAMPLE}\t${SAMPLE}\t0\t0\t0\t-9" > "$COMPHET_PED"

COMPHET_LOG="${OUTDIR}/${SAMPLE}_compound_hets.log"
# Remove stale outputs from a previous run so a failure cannot leave old results behind
rm -f "$COMPHET_VCF" "$COMPHET_TSV" "$COMPHET_PAIRS_TSV" "${COMPHET_TSV}.unphased"

# A failure here (wrong image, slivar crash, broken input) stops the step with
# slivar's own error. It must never be reported as "no candidates found".
if ! run_in --cpus 2 --memory 4g \
  "${SLIVAR_IMAGE}" \
  slivar compound-hets \
    --allow-non-trios \
    --vcf "/genome/${SAMPLE}/slivar/${SAMPLE}_prioritized.vcf.gz" \
    --ped "/genome/${SAMPLE}/slivar/${SAMPLE}.ped" \
  2>"$COMPHET_LOG" | \
  run_in -i    --cpus 2 --memory 2g \
    "${BCFTOOLS_IMAGE}" \
    bcftools view -Oz -o "/genome/${SAMPLE}/slivar/${SAMPLE}_compound_hets.vcf.gz"; then
  echo "ERROR: slivar compound-hets failed (image: ${SLIVAR_IMAGE})." >&2
  echo "  Output of the failed command (also in ${COMPHET_LOG}):" >&2
  sed 's/^/    /' "$COMPHET_LOG" >&2 || true
  rm -f "$COMPHET_VCF"
  exit 1
fi

# Count pairs from slivar_comphet INFO field (not VCF record count).
# slivar compound-hets outputs one record per unique VARIANT, with the
# slivar_comphet INFO field listing all partner variants. Format:
#   sample/GENE/PAIR_ID/chrom/pos/ref/alt  (comma-separated when multiple)
# A gene with N variants produces C(N,2) pairs but only N VCF records.
COMPHET_PAIRS=0
COMPHET_GENES=0
COMPHET_RECORDS=$(run_in \
  --cpus 2 --memory 2g \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -H "/genome/${SAMPLE}/slivar/${SAMPLE}_compound_hets.vcf.gz" | wc -l | tr -d ' ')
if [ "$COMPHET_RECORDS" -gt 0 ]; then
  COMPHET_STATS=$(run_in \
    --cpus 2 --memory 2g \
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools query -f '%INFO/slivar_comphet\n' \
      /genome/${SAMPLE}/slivar/${SAMPLE}_compound_hets.vcf.gz \
    | tr ',' '\n' \
    | awk -F'/' '{pairs[\$3]=1; genes[\$2]=1} END{print length(pairs), length(genes)}'")
  COMPHET_PAIRS=$(echo "$COMPHET_STATS" | awk '{print $1}')
  COMPHET_GENES=$(echo "$COMPHET_STATS" | awk '{print $2}')

  # Export to TSV sorted by gene for human review, with the population
  # frequency and the MANE Select consequence (EXTRA_COLS below).
  run_in --cpus 2 --memory 2g \
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools +split-vep \
        /genome/${SAMPLE}/slivar/${SAMPLE}_compound_hets.vcf.gz \
        -f '%SYMBOL\t%CHROM\t%POS\t%REF\t%ALT\t%IMPACT\t%Consequence[\t%GT]${EXTRA_FMT}\n' \
        -s worst -d \
      | sort -t\$'\t' -k1,1 -k2,2V -k3,3n \
      > /genome/${SAMPLE}/slivar/${SAMPLE}_compound_hets_raw.tsv"
  {
    printf 'GENE\tCHROM\tPOS\tREF\tALT\tIMPACT\tConsequence\tGT\t%s\n' "$EXTRA_HEAD"
    mane_col < "${OUTDIR}/${SAMPLE}_compound_hets_raw.tsv"
  } > "${COMPHET_TSV}.unphased"
  rm -f "${OUTDIR}/${SAMPLE}_compound_hets_raw.tsv"
  echo "  Found: ${COMPHET_PAIRS} compound het candidate pairs across ${COMPHET_GENES} genes"
else
  echo "  slivar ran and returned no records: no compound heterozygote candidates found."
fi

if [ ! -s "${COMPHET_TSV}.unphased" ]; then
  printf 'GENE\tCHROM\tPOS\tREF\tALT\tIMPACT\tConsequence\tGT\t%s\n' "$EXTRA_HEAD" > "${COMPHET_TSV}.unphased"
fi

# Phase of each pair from the reads (bin/comphet_phase.py): a read, or a read
# and its mate, that covers both positions shows the two alt alleles together
# (cis, one copy of the gene hit) or apart (trans, both copies hit). Pairs
# further apart than a read pair spans stay unknown. The pairs table holds the
# verdict and the read counts; the per-variant table gets a PHASE column after
# GT. pysam is in PYPGX_IMAGE, as for step 06's read check.
BAM="${SAMPLE_DIR}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam"
if [ -f "$BAM" ] && [ -f "${BAM}.bai" ]; then
  PHASE_ARGS=(--bam "$(cpath "$BAM")")
else
  PHASE_ARGS=(--skip "no BAM (${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam with its .bai)")
fi
if ! run_in --cpus 1 --memory 2g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYPGX_IMAGE" \
  python3 /pgp-bin/comphet_phase.py --vcf "$(cpath "$COMPHET_VCF")" --out "$(cpath "${COMPHET_PAIRS_TSV}.part")" \
    --variants "$(cpath "${COMPHET_TSV}.unphased")" --variants-out "$(cpath "${COMPHET_TSV}.part")" "${PHASE_ARGS[@]}"; then
  echo "ERROR: the phase of the compound-het pairs failed (bin/comphet_phase.py, image: ${PYPGX_IMAGE})." >&2
  echo "  Input: ${COMPHET_VCF} and ${PHASE_ARGS[*]}; would write ${COMPHET_PAIRS_TSV} and ${COMPHET_TSV}." >&2
  rm -f "${COMPHET_PAIRS_TSV}.part" "${COMPHET_TSV}.part"
  exit 1
fi
mv "${COMPHET_PAIRS_TSV}.part" "$COMPHET_PAIRS_TSV"
mv "${COMPHET_TSV}.part" "$COMPHET_TSV"
rm -f "${COMPHET_TSV}.unphased"
PHASE_TRANS=$(awk -F'\t' 'NR > 1 && $10 == "trans"' "$COMPHET_PAIRS_TSV" | wc -l | tr -d ' ')
PHASE_CIS=$(awk -F'\t' 'NR > 1 && $10 == "cis"' "$COMPHET_PAIRS_TSV" | wc -l | tr -d ' ')
PHASE_UNKNOWN=$(awk -F'\t' 'NR > 1 && $10 == "unknown"' "$COMPHET_PAIRS_TSV" | wc -l | tr -d ' ')

# ── Step 5: Generate summary TSV with gene constraint ─────────────────
echo ""
echo "[5/5] Generating summary with gene constraint annotations..."

SUMMARY_TSV="${OUTDIR}/${SAMPLE}_slivar_summary.tsv"

# Extract variant info from prioritized VCF into a TSV
if ! run_in --cpus 2 --memory 4g \
  "${BCFTOOLS_IMAGE}" \
  bash -o pipefail -c "bcftools +split-vep \
    /genome/${SAMPLE}/slivar/${SAMPLE}_prioritized.vcf.gz \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%IMPACT\t%SYMBOL\t%Consequence\t%Existing_variation[\t%GT]${EXTRA_FMT}\n' \
    -s worst -d \
  > /genome/${SAMPLE}/slivar/${SAMPLE}_variants_raw.tsv"; then
  echo "ERROR: Summary extraction from ${OUTDIR}/${SAMPLE}_prioritized.vcf.gz failed." >&2
  rm -f "${OUTDIR}/${SAMPLE}_variants_raw.tsv"
  exit 1
fi

# Add header and optional gene constraint columns (bin/constraint_join.awk,
# the loader step 23 and the Nextflow slivar module use: canonical rows,
# mis.z_score, the Ensembl row over the RefSeq one; it exits non-zero when rows
# carry genes and not one matches the file).
{
  printf 'CHROM\tPOS\tREF\tALT\tIMPACT\tSYMBOL\tConsequence\tExisting_variation\tGT\t%s\n' "$EXTRA_HEAD"
  mane_col < "${OUTDIR}/${SAMPLE}_variants_raw.tsv"
} > "${SUMMARY_TSV}.raw"
if [ -f "$CONSTRAINT_TSV" ]; then
  echo "  Joining with gnomAD gene constraint metrics..."
  awk -f "${PGP_ROOT}/bin/constraint_join.awk" gene_col=SYMBOL constrained=1 \
    "$CONSTRAINT_TSV" "${SUMMARY_TSV}.raw" > "${SUMMARY_TSV}.tmp"
  rm -f "${SUMMARY_TSV:?}.raw"
else
  echo "  Gene constraint file not found (optional): ${CONSTRAINT_TSV}"
  echo "  Generating summary without constraint annotations."
  mv "${SUMMARY_TSV}.raw" "${SUMMARY_TSV}.tmp"
fi
mv "${SUMMARY_TSV}.tmp" "$SUMMARY_TSV"

# Clean up intermediate files
rm -f "${OUTDIR}/${SAMPLE}_variants_raw.tsv" "${OUTDIR}/${SAMPLE}_ranked.vcf.gz" "${OUTDIR}/${SAMPLE}_ranked.vcf.gz.tbi"

# ── Summary ───────────────────────────────────────────────────────────
echo ""
echo "============================================"
echo "  Step 31 complete: ${SAMPLE}"
echo ""
echo "  Filter results:"
echo "    rare_high (HIGH, rare):            ${RARE_HIGH_COUNT}"
echo "    rare_moderate_deleterious:          ${MODERATE_DEL_COUNT}"
if [ -n "$FREQ_NAME" ]; then
echo "    freq_unknown (HIGH/MODERATE, not rare): ${FREQ_UNKNOWN_COUNT}"
fi
if [ "$HAS_CLINVAR" -eq 1 ]; then
echo "    clinvar_pathogenic:                 ${CLINVAR_COUNT}"
fi
echo "    ────────────────────────────────────"
echo "    Total prioritized (deduplicated):   ${PRIORITIZED_COUNT}"
echo ""
echo "  Compound heterozygote candidates:     ${COMPHET_PAIRS} pairs (${COMPHET_GENES} genes)"
echo "    phase from the reads: ${PHASE_TRANS} trans, ${PHASE_CIS} cis, ${PHASE_UNKNOWN} unknown"
echo ""

# Highlight constrained genes if constraint data was used
if [ -f "$CONSTRAINT_TSV" ] && [ -f "$SUMMARY_TSV" ]; then
  CONSTRAINED_COUNT=$(tail -n +2 "$SUMMARY_TSV" 2>/dev/null | awk -F'\t' '$NF=="YES"' | wc -l | tr -d ' ')
  if [ "$CONSTRAINED_COUNT" -gt 0 ]; then
    echo "  Variants in constrained genes:        ${CONSTRAINED_COUNT}"
    echo "  (LOEUF < 0.35 or pLI > 0.9 — loss-of-function intolerant)"
    echo ""
    echo "  Top constrained gene hits:"
    tail -n +2 "$SUMMARY_TSV" | awk -F'\t' '$NF=="YES" {print "    "$6" ("$5") "$7}' | sort -u | head -10
    echo ""
  fi
fi

echo "  Output files:"
echo "    ${OUTDIR}/${SAMPLE}_prioritized.vcf.gz         (all prioritized variants)"
echo "    ${OUTDIR}/${SAMPLE}_compound_hets.vcf.gz       (compound het VCF)"
echo "    ${OUTDIR}/${SAMPLE}_compound_hets.tsv           (compound het candidates, PHASE column)"
echo "    ${OUTDIR}/${SAMPLE}_compound_het_pairs.tsv      (each pair: cis, trans or unknown, with the reads)"
echo "    ${OUTDIR}/${SAMPLE}_slivar_summary.tsv          (summary + gene constraint)"
echo "    ${OUTDIR}/${SAMPLE}_rare_high.vcf.gz            (HIGH impact tier)"
echo "    ${OUTDIR}/${SAMPLE}_rare_moderate_del.vcf.gz    (MODERATE + deleterious)"
if [ -n "$FREQ_UNKNOWN_FILE" ]; then
echo "    ${OUTDIR}/${SAMPLE}_freq_unknown.vcf.gz         (HIGH/MODERATE, frequency unknown)"
fi
if [ "$HAS_CLINVAR" -eq 1 ]; then
echo "    ${OUTDIR}/${SAMPLE}_clinvar_path.vcf.gz         (ClinVar P/LP)"
fi
echo "============================================"
echo ""
echo "NOTE: A compound-het pair is two-hit (both copies of the gene) only when"
echo "  it is in trans. cis means one copy carries both variants. unknown means"
echo "  no read or read pair covers both positions: a parent's genome (step 38)"
echo "  or long reads can phase it."
echo ""
echo "Next: Review ${SAMPLE}_slivar_summary.tsv or load prioritized VCF in IGV/gene.iobio"
