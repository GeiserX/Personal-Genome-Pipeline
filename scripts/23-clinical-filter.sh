#!/usr/bin/env bash
# 23-clinical-filter.sh — Extract clinically interesting variants from annotated VCF
# Usage: ./scripts/23-clinical-filter.sh <sample_name>
#
# Produces a small VCF of PASS variants that are:
#   - rare (MAX_AF < 1%, or no frequency) AND HIGH or MODERATE VEP impact
#   - OR in the ClinVar screen of step 6 (P/LP alleles; any frequency)
#   - OR rare with a high CADD score (>= 20) outside HIGH/MODERATE
#   - OR rare with a high SpliceAI delta score (>= 0.2), any gene of the value
#   - OR rare with REVEL >= 0.644 or AlphaMissense >= 0.564
# The rarity filter uses VEP's MAX_AF (highest frequency in any population of
# 1000 Genomes, gnomAD exomes and gnomAD genomes), or gnomADe_AF and gnomADg_AF
# when MAX_AF is absent. A variant common in genomes but absent from exomes is
# not rare. With no frequency field at all every tier keeps all variants, and
# the step says so.
#
# The summary TSV takes the gene, impact and consequence of the worst
# consequence from `bcftools +split-vep`, and gnomAD constraint columns from
# bin/constraint_join.awk (the loader step 31 and the Nextflow slivar module use).
# Requires: VEP-annotated VCF from step 13 (step 30 vcfanno enrichment recommended)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"

ANNOTATED_VCF="${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_annotated.vcf.gz"
VEP_VCF="${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_vep.vcf"
VEP_VCF_GZ="${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz"
OUTDIR="${GENOME_DIR}/${SAMPLE}/clinical"
CONSTRAINT_TSV="${GENOME_DIR}/annotations/gnomad_v4.1_constraint.tsv"
CLINVAR_HITS="${GENOME_DIR}/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.vcf"
C="/genome/${SAMPLE}/clinical"
mkdir -p "$OUTDIR"

# Input: a derived file is used only when it is newer than what it was built
# from, so a re-run of step 13 is never hidden behind an older _vep.vcf.gz or
# an older step 30 output.
VEP_SRC=""
if [ -f "$VEP_VCF" ] && { [ ! -f "$VEP_VCF_GZ" ] || [ "$VEP_VCF" -nt "$VEP_VCF_GZ" ]; }; then
  VEP_SRC="$VEP_VCF"
elif [ -f "$VEP_VCF_GZ" ]; then
  VEP_SRC="$VEP_VCF_GZ"
fi
INPUT=""
if [ -f "$ANNOTATED_VCF" ] && { [ -z "$VEP_SRC" ] || [ "$ANNOTATED_VCF" -nt "$VEP_SRC" ]; }; then
  INPUT="$ANNOTATED_VCF"
elif [ -n "$VEP_SRC" ]; then
  INPUT="$VEP_SRC"
  if [ -f "$ANNOTATED_VCF" ]; then
    echo "NOTICE: ${ANNOTATED_VCF} is older than ${VEP_SRC}; ignoring it. Run step 30 again for the score tiers."
  fi
else
  echo "ERROR: No annotated VCF found. Run step 13 (VEP) first."
  echo "  Expected: ${ANNOTATED_VCF} or ${VEP_VCF_GZ} or ${VEP_VCF}"
  exit 1
fi

echo "============================================"
echo "  Step 23: Clinical Variant Filter"
echo "  Sample: ${SAMPLE}"
echo "  Input:  ${INPUT}"
echo "  Output: ${OUTDIR}/"
echo "============================================"
echo ""

# Step 1: compress and index if needed (bcftools filters need an indexed .vcf.gz)
if [ "$INPUT" = "$VEP_VCF" ]; then
  echo "[1] Compressing VEP VCF (required for bcftools filtering)..."
  rm -f "${VEP_VCF_GZ:?}" "${VEP_VCF_GZ:?}.tbi"
  run_in --cpus 4 --memory 4g \
    "${BCFTOOLS_IMAGE}" \
    bash -c "bcftools view /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf -Oz \
      -o /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz.tmp && \
      mv /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz.tmp /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz && \
      bcftools index -f -t /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz"
  INPUT="$VEP_VCF_GZ"
elif [ ! -f "${INPUT}.tbi" ] || [ "$INPUT" -nt "${INPUT}.tbi" ]; then
  echo "[1] Indexing $(basename "$INPUT")..."
  run_in --cpus 2 --memory 2g "${BCFTOOLS_IMAGE}" bcftools index -f -t "$(cpath "$INPUT")"
else
  echo "[1] Input already compressed and indexed."
fi
CONTAINER_INPUT=$(cpath "$INPUT")

# Detect available CSQ subfields and INFO tags (header only)
VEP_FIELDS=$(run_in "${BCFTOOLS_IMAGE}" bcftools +split-vep -l "$CONTAINER_INPUT" 2>/dev/null | cut -f2 || true)
if [ -z "$VEP_FIELDS" ]; then
  echo "ERROR: No CSQ/BCSQ annotation found in VEP VCF."
  echo "  Was VEP step 13 run correctly? The VCF must contain a CSQ INFO field."
  exit 1
fi
has_field() { printf '%s\n' "$VEP_FIELDS" | grep -qx "$1"; }
for f in IMPACT SYMBOL Consequence; do
  has_field "$f" || { echo "ERROR: the CSQ annotation has no ${f} field." >&2; exit 1; }
done
VCF_HEADER=$(run_in "${BCFTOOLS_IMAGE}" bcftools view -h "$CONTAINER_INPUT" 2>/dev/null || true)
has_info() { printf '%s\n' "$VCF_HEADER" | grep -q "^##INFO=<ID=$1,"; }

# Frequency: MAX_AF, or the gnomAD exome and genome fields VEP was asked for.
FREQ_COLS=""
FREQ_EXPR=""
FREQ_NAME=""
if has_field MAX_AF; then
  FREQ_COLS="MAX_AF:Float"
  FREQ_EXPR='(MAX_AF<0.01 || MAX_AF=".")'
  FREQ_NAME="MAX_AF"
else
  for f in gnomADe_AF gnomADg_AF; do
    has_field "$f" || continue
    FREQ_COLS="${FREQ_COLS:+${FREQ_COLS},}${f}:Float"
    FREQ_EXPR="${FREQ_EXPR:+${FREQ_EXPR} && }(${f}<0.01 || ${f}=\".\")"
    FREQ_NAME="${FREQ_NAME:+${FREQ_NAME} and }${f}"
  done
fi
HAS_CLINSIG=0; has_field CLIN_SIG && HAS_CLINSIG=1
HAS_CADD=0; has_info CADD_PHRED && HAS_CADD=1
HAS_CADD_INDEL=0; has_info CADD_PHRED_indel && HAS_CADD_INDEL=1
HAS_SPLICEAI=0; has_info SpliceAI && HAS_SPLICEAI=1
HAS_SPLICEAI_INDEL=0; has_info SpliceAI_indel && HAS_SPLICEAI_INDEL=1
HAS_REVEL=0; has_info REVEL && HAS_REVEL=1
HAS_ALPHAMISSENSE=0; has_info AM_pathogenicity && HAS_ALPHAMISSENSE=1
HAS_CONSTRAINT=0; [ -f "$CONSTRAINT_TSV" ] && HAS_CONSTRAINT=1
HAS_HITS=0; [ -f "$CLINVAR_HITS" ] && HAS_HITS=1

if [ "$HAS_HITS" -eq 1 ]; then
  CLINVAR_PLAN="step 6 hits (${CLINVAR_HITS})"
elif [ "$HAS_CLINSIG" -eq 1 ]; then
  CLINVAR_PLAN="VEP's cached CLIN_SIG (run step 6 to use the current ClinVar file)"
else
  CLINVAR_PLAN="none (run step 6)"
fi
echo "  Population frequency: ${FREQ_NAME:-none in the VEP output (no tier is filtered by frequency)}"
echo "  ClinVar tier: ${CLINVAR_PLAN}"
echo "  CADD scores: $([ "$HAS_CADD" -eq 1 ] && echo 'available' || echo 'not annotated (run step 30)')$([ "$HAS_CADD_INDEL" -eq 1 ] && echo ' (+indels)')"
echo "  SpliceAI scores: $([ "$HAS_SPLICEAI" -eq 1 ] && echo 'available' || echo 'not annotated (run step 30)')$([ "$HAS_SPLICEAI_INDEL" -eq 1 ] && echo ' (+indels)')"
echo "  REVEL scores: $([ "$HAS_REVEL" -eq 1 ] && echo 'available' || echo 'not annotated (run step 30)')"
echo "  AlphaMissense: $([ "$HAS_ALPHAMISSENSE" -eq 1 ] && echo 'available' || echo 'not annotated (run step 30)')"
echo "  gnomAD constraint: $([ "$HAS_CONSTRAINT" -eq 1 ] && echo 'available' || echo 'not downloaded')"
echo ""
if [ -z "$FREQ_EXPR" ]; then
  echo "  NOTICE: the VEP output has no MAX_AF, gnomADe_AF or gnomADg_AF field, so no tier"
  echo "  is filtered by frequency and every MODERATE variant is kept. Run step 13 with"
  echo "  --everything (or --max_af) for population frequencies."
  echo ""
fi

# count FILE: records in a tier file.
count() {
  run_in "${BCFTOOLS_IMAGE}" bcftools view -H "${C}/$1" 2>/dev/null | wc -l | tr -d ' ' || echo 0
}

# Step 2: the rare PASS records every non-ClinVar tier starts from.
echo "[2] Selecting PASS records$([ -n "$FREQ_EXPR" ] && echo " with ${FREQ_NAME} < 1% or missing")..."
if [ -n "$FREQ_EXPR" ]; then
  run_in --cpus 4 --memory 4g "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
      bcftools +split-vep - -c '${FREQ_COLS}' -s worst -i '${FREQ_EXPR}' \
        -Oz -o ${C}/${SAMPLE}_rare_pass.vcf.gz && \
      bcftools index -f -t ${C}/${SAMPLE}_rare_pass.vcf.gz"
else
  run_in --cpus 4 --memory 4g "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} -Oz -o ${C}/${SAMPLE}_rare_pass.vcf.gz && \
      bcftools index -f -t ${C}/${SAMPLE}_rare_pass.vcf.gz"
fi
RARE="${C}/${SAMPLE}_rare_pass.vcf.gz"

# Step 3: HIGH impact (stop-gain, frameshift, splice site), rare
echo "[3] Extracting rare HIGH impact variants (stop-gain, frameshift, splice)..."
run_in --cpus 4 --memory 4g "${BCFTOOLS_IMAGE}" \
  bash -o pipefail -c "bcftools +split-vep ${RARE} -c IMPACT -s worst -i 'IMPACT=\"HIGH\"' \
      -Oz -o ${C}/${SAMPLE}_high_impact.vcf.gz && \
    bcftools index -f -t ${C}/${SAMPLE}_high_impact.vcf.gz"
HIGH_COUNT=$(count "${SAMPLE}_high_impact.vcf.gz")
echo "  Found: ${HIGH_COUNT} HIGH impact variants"

# Step 4: MODERATE impact (missense, in-frame indel), rare
echo "[4] Extracting rare MODERATE impact variants (missense, in-frame indel)..."
run_in --cpus 4 --memory 4g "${BCFTOOLS_IMAGE}" \
  bash -o pipefail -c "bcftools +split-vep ${RARE} -c IMPACT -s worst -i 'IMPACT=\"MODERATE\"' \
      -Oz -o ${C}/${SAMPLE}_rare_moderate.vcf.gz && \
    bcftools index -f -t ${C}/${SAMPLE}_rare_moderate.vcf.gz"
MODERATE_COUNT=$(count "${SAMPLE}_rare_moderate.vcf.gz")
echo "  Found: ${MODERATE_COUNT} MODERATE impact variants"
MERGE_FILES="${C}/${SAMPLE}_high_impact.vcf.gz ${C}/${SAMPLE}_rare_moderate.vcf.gz"

# Step 5: ClinVar pathogenic/likely pathogenic, at any frequency (a common
# pathogenic allele such as HFE p.C282Y must stay). Step 6 screens the sample
# against the ClinVar file in clinvar/, which setup.sh refreshes; VEP's
# CLIN_SIG comes from its cache and is only the fallback.
CLINVAR_COUNT=0
CLINVAR_SOURCE="none"
CLINVAR_TIER="${OUTDIR}/${SAMPLE}_clinvar_pathogenic.vcf.gz"
rm -f "${CLINVAR_TIER:?}" "${CLINVAR_TIER:?}.tbi"
if [ "$HAS_HITS" -eq 1 ]; then
  echo "[5] Extracting the step 6 ClinVar hits..."
  CLINVAR_SOURCE="step 6 ClinVar screen"
  TARGETS="${OUTDIR}/${SAMPLE}_clinvar_targets.tsv"
  awk -F'\t' 'BEGIN {OFS = "\t"} !/^#/ {print $1, $2}' "$CLINVAR_HITS" | sort -u > "$TARGETS"
  if [ -s "$TARGETS" ]; then
    run_in --cpus 4 --memory 4g "${BCFTOOLS_IMAGE}" \
      bash -o pipefail -c "bcftools view -f PASS -T ${C}/${SAMPLE}_clinvar_targets.tsv ${CONTAINER_INPUT} \
          -Oz -o ${C}/${SAMPLE}_clinvar_pathogenic.vcf.gz && \
        bcftools index -f -t ${C}/${SAMPLE}_clinvar_pathogenic.vcf.gz"
  else
    echo "  Step 6 found no ClinVar hit."
  fi
  rm -f "${TARGETS:?}"
elif [ "$HAS_CLINSIG" -eq 1 ]; then
  echo "[5] Extracting VEP CLIN_SIG pathogenic/likely pathogenic (VEP's cached ClinVar; run step 6 for the current file)..."
  CLINVAR_SOURCE="VEP CLIN_SIG (cache release)"
  run_in --cpus 4 --memory 4g "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
      bcftools +split-vep - -c CLIN_SIG \
        -i 'CLIN_SIG~\"pathogenic\" && CLIN_SIG!~\"conflicting\"' \
        -Oz -o ${C}/${SAMPLE}_clinvar_pathogenic.vcf.gz && \
      bcftools index -f -t ${C}/${SAMPLE}_clinvar_pathogenic.vcf.gz"
else
  echo "[5] Skipping the ClinVar tier (no step 6 hits and no CLIN_SIG in the VEP output)."
fi
if [ -f "$CLINVAR_TIER" ]; then
  CLINVAR_COUNT=$(count "${SAMPLE}_clinvar_pathogenic.vcf.gz")
  echo "  Found: ${CLINVAR_COUNT} ClinVar pathogenic/likely pathogenic variants (${CLINVAR_SOURCE})"
  MERGE_FILES="${MERGE_FILES} ${C}/${SAMPLE}_clinvar_pathogenic.vcf.gz"
fi

# Step 6: high CADD outside HIGH/MODERATE, rare (only tags present in the header)
CADD_COUNT=0
if [ "$HAS_CADD" -eq 1 ] || [ "$HAS_CADD_INDEL" -eq 1 ]; then
  CADD_EXPR=''
  [ "$HAS_CADD" -eq 1 ] && CADD_EXPR='INFO/CADD_PHRED>=20'
  [ "$HAS_CADD_INDEL" -eq 1 ] && CADD_EXPR="${CADD_EXPR:+${CADD_EXPR} || }INFO/CADD_PHRED_indel>=20"
  echo "[6] Extracting rare high-CADD variants (PHRED >= 20, non-HIGH/MODERATE)..."
  run_in --cpus 4 --memory 4g "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools +split-vep ${RARE} -c IMPACT -s worst \
        -i 'IMPACT!=\"HIGH\" && IMPACT!=\"MODERATE\" && (${CADD_EXPR})' \
        -Oz -o ${C}/${SAMPLE}_cadd_high.vcf.gz && \
      bcftools index -f -t ${C}/${SAMPLE}_cadd_high.vcf.gz"
  CADD_COUNT=$(count "${SAMPLE}_cadd_high.vcf.gz")
  echo "  Found: ${CADD_COUNT} high-CADD non-coding variants (PHRED >= 20)"
  MERGE_FILES="${MERGE_FILES} ${C}/${SAMPLE}_cadd_high.vcf.gz"
fi

# Step 7: SpliceAI cryptic splice variants, rare
SPLICEAI_COUNT=0
if [ "$HAS_SPLICEAI" -eq 1 ] || [ "$HAS_SPLICEAI_INDEL" -eq 1 ]; then
  SPLICEAI_PREFILTER=''
  [ "$HAS_SPLICEAI" -eq 1 ] && SPLICEAI_PREFILTER='INFO/SpliceAI!="."'
  [ "$HAS_SPLICEAI_INDEL" -eq 1 ] && SPLICEAI_PREFILTER="${SPLICEAI_PREFILTER:+${SPLICEAI_PREFILTER} || }INFO/SpliceAI_indel!=\".\""
  echo "[7] Extracting rare cryptic splice variants (SpliceAI delta >= 0.2)..."
  # A SpliceAI value is ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL,
  # one per gene joined by ','. Every gene's four delta scores are tested.
  run_in --cpus 4 --memory 4g "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -i '${SPLICEAI_PREFILTER}' ${RARE} | \
      awk -F'\t' 'BEGIN{OFS=\"\t\"} /^#/{print;next} {
        hit=0
        n=split(\$8, kv, \";\")
        for(i=1;i<=n;i++){
          if(kv[i] !~ /^SpliceAI(_indel)?=/) continue
          v=kv[i]; sub(/^[^=]*=/,\"\",v)
          na=split(v, genes, \",\")
          for(a=1;a<=na;a++){
            split(genes[a], sp, \"|\")
            for(j=3;j<=6;j++) if(sp[j]!=\"\" && sp[j]!=\".\" && sp[j]+0>=0.2) hit=1
          }
        }
        if(hit) print
      }' | \
      bcftools view -Oz -o ${C}/${SAMPLE}_spliceai_high.vcf.gz.tmp - && \
    mv ${C}/${SAMPLE}_spliceai_high.vcf.gz.tmp ${C}/${SAMPLE}_spliceai_high.vcf.gz && \
    bcftools index -f -t ${C}/${SAMPLE}_spliceai_high.vcf.gz"
  SPLICEAI_COUNT=$(count "${SAMPLE}_spliceai_high.vcf.gz")
  echo "  Found: ${SPLICEAI_COUNT} cryptic splice variants (SpliceAI >= 0.2)"
  MERGE_FILES="${MERGE_FILES} ${C}/${SAMPLE}_spliceai_high.vcf.gz"
fi

# Step 8: high-confidence deleterious missense (REVEL/AlphaMissense), rare
MISSENSE_COUNT=0
if [ "$HAS_REVEL" -eq 1 ] || [ "$HAS_ALPHAMISSENSE" -eq 1 ]; then
  MISSENSE_FILTER=""
  MISSENSE_LABEL=""
  if [ "$HAS_REVEL" -eq 1 ]; then
    MISSENSE_FILTER="INFO/REVEL>=0.644"
    MISSENSE_LABEL="REVEL >= 0.644 (ClinGen PP3 Supporting)"
  fi
  if [ "$HAS_ALPHAMISSENSE" -eq 1 ]; then
    MISSENSE_FILTER="${MISSENSE_FILTER:+${MISSENSE_FILTER} || }INFO/AM_pathogenicity>=0.564"
    MISSENSE_LABEL="${MISSENSE_LABEL:+${MISSENSE_LABEL} or }AlphaMissense >= 0.564 (its likely_pathogenic class boundary, not an ACMG evidence level)"
  fi
  echo "[8] Extracting rare deleterious missense variants: ${MISSENSE_LABEL}..."
  run_in --cpus 4 --memory 4g "${BCFTOOLS_IMAGE}" \
    bash -c "bcftools view -i '${MISSENSE_FILTER}' ${RARE} \
      -Oz -o ${C}/${SAMPLE}_missense_deleterious.vcf.gz && \
    bcftools index -f -t ${C}/${SAMPLE}_missense_deleterious.vcf.gz"
  MISSENSE_COUNT=$(count "${SAMPLE}_missense_deleterious.vcf.gz")
  echo "  Found: ${MISSENSE_COUNT} deleterious missense variants"
  MERGE_FILES="${MERGE_FILES} ${C}/${SAMPLE}_missense_deleterious.vcf.gz"
fi

# Step 9: merge into the combined clinical VCF
echo "[9] Merging into combined clinical VCF..."
run_in --cpus 2 --memory 2g "${BCFTOOLS_IMAGE}" \
  bash -o pipefail -c "bcftools concat -a -D ${MERGE_FILES} | \
    bcftools sort -Oz -o ${C}/${SAMPLE}_clinical.vcf.gz && \
    bcftools index -f -t ${C}/${SAMPLE}_clinical.vcf.gz"
RARE_HOST="${OUTDIR}/${SAMPLE}_rare_pass.vcf.gz"
rm -f "${RARE_HOST:?}" "${RARE_HOST:?}.tbi"
TOTAL_COUNT=$(count "${SAMPLE}_clinical.vcf.gz")

# Summary TSV: one row per variant, from its worst consequence. A score or
# frequency the input does not carry is written as '.'.
echo ""
echo "Generating human-readable summary..."
COL_CADD='.'; [ "$HAS_CADD" -eq 1 ] && COL_CADD='%INFO/CADD_PHRED'
COL_REVEL='.'; [ "$HAS_REVEL" -eq 1 ] && COL_REVEL='%INFO/REVEL'
COL_AM='.'; has_info AM_class && COL_AM='%INFO/AM_class'
COL_FREQ='.'; has_field MAX_AF && COL_FREQ='%MAX_AF'
SUMMARY="${OUTDIR}/${SAMPLE}_clinical_summary.tsv"
printf 'CHROM\tPOS\tREF\tALT\tGT\tIMPACT\tGENE\tConsequence\tMAX_AF\tCADD_PHRED\tREVEL\tAM_CLASS\n' > "${SUMMARY}.tmp"
run_in --cpus 2 --memory 2g "${BCFTOOLS_IMAGE}" \
  bcftools +split-vep "${C}/${SAMPLE}_clinical.vcf.gz" -s worst \
    -f "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\t%IMPACT\t%SYMBOL\t%Consequence\t${COL_FREQ}\t${COL_CADD}\t${COL_REVEL}\t${COL_AM}\n" \
  >> "${SUMMARY}.tmp"
# An empty SYMBOL (an intergenic worst consequence) is written as '.'.
awk -F'\t' 'BEGIN {OFS = "\t"} NR > 1 && $7 == "" {$7 = "."} {print}' "${SUMMARY}.tmp" > "${SUMMARY}.tmp2"
mv "${SUMMARY}.tmp2" "${SUMMARY}.tmp"

# gnomAD gene constraint columns (bin/constraint_join.awk exits non-zero when
# rows carry genes and not one matches: a wrong file, never a silent '.').
if [ "$HAS_CONSTRAINT" -eq 1 ]; then
  echo "Adding gnomAD gene constraint metrics..."
  awk -f "${PGP_ROOT}/bin/constraint_join.awk" gene_col=GENE "$CONSTRAINT_TSV" "${SUMMARY}.tmp" > "${SUMMARY}.tmp2"
  mv "${SUMMARY}.tmp2" "${SUMMARY}.tmp"
fi
mv "${SUMMARY}.tmp" "$SUMMARY"

echo ""
echo "============================================"
echo "  Clinical filter complete: ${SAMPLE}"
echo "  Total clinically interesting variants: ${TOTAL_COUNT}"
echo "    HIGH impact (LoF):           ${HIGH_COUNT}"
echo "    Rare MODERATE impact:        ${MODERATE_COUNT}"
echo "    ClinVar pathogenic/LP:       ${CLINVAR_COUNT} (${CLINVAR_SOURCE})"
{ [ "$HAS_CADD" -eq 1 ] || [ "$HAS_CADD_INDEL" -eq 1 ]; } && \
echo "    High CADD non-coding:        ${CADD_COUNT}"
{ [ "$HAS_SPLICEAI" -eq 1 ] || [ "$HAS_SPLICEAI_INDEL" -eq 1 ]; } && \
echo "    Cryptic splice (SpliceAI):   ${SPLICEAI_COUNT}"
{ [ "$HAS_REVEL" -eq 1 ] || [ "$HAS_ALPHAMISSENSE" -eq 1 ]; } && \
echo "    Deleterious missense:        ${MISSENSE_COUNT}"
echo "  Frequency filter: ${FREQ_NAME:-none (no frequency field in the VEP output)}"
echo ""
echo "  Output files:"
echo "    ${OUTDIR}/${SAMPLE}_clinical.vcf.gz              (combined)"
echo "    ${OUTDIR}/${SAMPLE}_clinical_summary.tsv         (human-readable table)"
echo "    ${OUTDIR}/${SAMPLE}_high_impact.vcf.gz           (HIGH only)"
echo "    ${OUTDIR}/${SAMPLE}_rare_moderate.vcf.gz         (rare MODERATE only)"
[ -f "$CLINVAR_TIER" ] && \
echo "    ${CLINVAR_TIER}    (ClinVar P/LP only)"
{ [ "$HAS_CADD" -eq 1 ] || [ "$HAS_CADD_INDEL" -eq 1 ]; } && \
echo "    ${OUTDIR}/${SAMPLE}_cadd_high.vcf.gz             (CADD >= 20 non-coding)"
{ [ "$HAS_SPLICEAI" -eq 1 ] || [ "$HAS_SPLICEAI_INDEL" -eq 1 ]; } && \
echo "    ${OUTDIR}/${SAMPLE}_spliceai_high.vcf.gz         (SpliceAI >= 0.2)"
{ [ "$HAS_REVEL" -eq 1 ] || [ "$HAS_ALPHAMISSENSE" -eq 1 ]; } && \
echo "    ${OUTDIR}/${SAMPLE}_missense_deleterious.vcf.gz  (REVEL/AlphaMissense)"
echo "============================================"
echo ""
echo "Next: Review ${SAMPLE}_clinical_summary.tsv or load the VCF in IGV/gene.iobio"
