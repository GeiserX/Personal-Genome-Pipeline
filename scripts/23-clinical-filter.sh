#!/usr/bin/env bash
# 23-clinical-filter.sh — Extract clinically interesting variants from annotated VCF
# Usage: ./scripts/23-clinical-filter.sh <sample_name>
#
# Produces a small VCF of PASS variants that are:
#   - rare (MAX_AF < 1%) AND HIGH or MODERATE VEP impact
#   - OR HIGH or MODERATE with no known frequency: their own tier
#     (_freq_unknown.vcf.gz), never counted as rare
#   - OR in the ClinVar screen of step 6 (P/LP alleles; any frequency), else
#     P/LP in ClinVar_CLNSIG (step 13's --custom annotation with the same
#     ClinVar file), else in the cache's CLIN_SIG
#   - OR rare or of unknown frequency with a high CADD score (>= 20) outside
#     HIGH/MODERATE
#   - OR rare or of unknown frequency with a high SpliceAI delta score
#     (>= 0.2), any gene of the value
#   - OR rare or of unknown frequency with REVEL >= 0.644 or AlphaMissense >= 0.564
# The rarity rule uses VEP's MAX_AF (highest frequency in any population of
# 1000 Genomes, gnomAD exomes and gnomAD genomes), or gnomADe_AF and gnomADg_AF
# when MAX_AF is absent. A variant common in genomes but absent from exomes is
# not rare. A variant with no value at all is "frequency unknown", not rare:
# VEP finds no known variant for an indel written differently from gnomAD
# (step 13 normalises its input for that reason). With no frequency field at
# all every tier keeps all variants, and the step says so.
#
# Impact: bin/csq_rank.awk ranks each record on its MANE Select transcript,
# else the canonical one, else any transcript (the most severe IMPACT in that
# set). A HIGH consequence on another transcript does not make a record HIGH
# when its MANE Select consequence is intronic; the summary shows both.
# Genotype: the same pass marks a call weak (CALL_CHECK) when GQ is below
# CLINICAL_MIN_GQ (default 30, step 6's floor), DP below CLINICAL_MIN_DP
# (default 10), a heterozygote's alt read fraction outside
# CLINICAL_HET_AB_MIN..CLINICAL_HET_AB_MAX (0.2..0.8) or a homozygous alt
# call's below CLINICAL_HOM_AB_MIN (0.8). Weak calls stay in the tiers and are
# counted apart.
#
# The summary TSV takes the gene, impact and consequence of the ranked
# transcript, the worst consequence of `bcftools +split-vep` beside it, and
# gnomAD constraint columns from bin/constraint_join.awk (the loader step 31
# and the Nextflow slivar module use).
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
# an older step 30 output. An empty file (an interrupted write) never counts.
VEP_SRC=""
if [ -s "$VEP_VCF" ] && { [ ! -s "$VEP_VCF_GZ" ] || [ "$VEP_VCF" -nt "$VEP_VCF_GZ" ]; }; then
  VEP_SRC="$VEP_VCF"
elif [ -s "$VEP_VCF_GZ" ]; then
  VEP_SRC="$VEP_VCF_GZ"
fi
INPUT=""
if [ -s "$ANNOTATED_VCF" ] && { [ -z "$VEP_SRC" ] || [ "$ANNOTATED_VCF" -nt "$VEP_SRC" ]; }; then
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
# Here-strings, not a pipe: grep -q exits at the first match, and under
# pipefail the writer's SIGPIPE on a large header would read as "no match".
has_field() { grep -qx -- "$1" <<< "$VEP_FIELDS"; }
for f in IMPACT SYMBOL Consequence; do
  has_field "$f" || { echo "ERROR: the CSQ annotation has no ${f} field." >&2; exit 1; }
done
VCF_HEADER=$(run_in "${BCFTOOLS_IMAGE}" bcftools view -h "$CONTAINER_INPUT" 2>/dev/null || true)
has_info() { grep -q "^##INFO=<ID=$1," <<< "$VCF_HEADER"; }

# Frequency: MAX_AF, or the gnomAD exome and genome fields VEP was asked for.
# bin/csq_rank.awk applies the rule (PGP_FREQ); this only names the fields.
FREQ_NAME=""
if has_field MAX_AF; then
  FREQ_NAME="MAX_AF"
else
  for f in gnomADe_AF gnomADg_AF; do
    if has_field "$f"; then FREQ_NAME="${FREQ_NAME:+${FREQ_NAME} and }${f}"; fi
  done
fi
HAS_CLINSIG=0; has_field CLIN_SIG && HAS_CLINSIG=1
HAS_CUSTOM=0; has_field ClinVar_CLNSIG && HAS_CUSTOM=1
HAS_CADD=0; has_info CADD_PHRED && HAS_CADD=1
HAS_CADD_INDEL=0; has_info CADD_PHRED_indel && HAS_CADD_INDEL=1
HAS_SPLICEAI=0; has_info SpliceAI && HAS_SPLICEAI=1
HAS_SPLICEAI_INDEL=0; has_info SpliceAI_indel && HAS_SPLICEAI_INDEL=1
HAS_REVEL=0; has_info REVEL && HAS_REVEL=1
HAS_ALPHAMISSENSE=0; has_info AM_pathogenicity && HAS_ALPHAMISSENSE=1
HAS_CONSTRAINT=0; [ -f "$CONSTRAINT_TSV" ] && HAS_CONSTRAINT=1
HAS_HITS=0; [ -f "$CLINVAR_HITS" ] && HAS_HITS=1
has_format() { grep -q "^##FORMAT=<ID=$1," <<< "$VCF_HEADER"; }

# Call floors of the genotype check (bin/csq_rank.awk writes them into the
# clinical VCF as ##PGP_CALL_FLOOR, which the reports read).
MIN_GQ=${CLINICAL_MIN_GQ:-30}
MIN_DP=${CLINICAL_MIN_DP:-10}
HET_AB_MIN=${CLINICAL_HET_AB_MIN:-0.2}
HET_AB_MAX=${CLINICAL_HET_AB_MAX:-0.8}
HOM_AB_MIN=${CLINICAL_HOM_AB_MIN:-0.8}
for v in "CLINICAL_MIN_GQ=${MIN_GQ}" "CLINICAL_MIN_DP=${MIN_DP}"; do
  [[ "${v#*=}" =~ ^[0-9]+$ ]] || { echo "ERROR: ${v%%=*} must be a whole number, got '${v#*=}'" >&2; exit 1; }
done
for v in "CLINICAL_HET_AB_MIN=${HET_AB_MIN}" "CLINICAL_HET_AB_MAX=${HET_AB_MAX}" "CLINICAL_HOM_AB_MIN=${HOM_AB_MIN}"; do
  [[ "${v#*=}" =~ ^(0(\.[0-9]+)?|1(\.0+)?)$ ]] || { echo "ERROR: ${v%%=*} must be a fraction from 0 to 1, got '${v#*=}'" >&2; exit 1; }
done

if [ "$HAS_HITS" -eq 1 ]; then
  CLINVAR_PLAN="step 6 hits (${CLINVAR_HITS})"
elif [ "$HAS_CUSTOM" -eq 1 ]; then
  CLINVAR_PLAN="ClinVar_CLNSIG from step 13's --custom ClinVar file"
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
echo "  Genotype check: GQ < ${MIN_GQ}, DP < ${MIN_DP}, het alt fraction outside ${HET_AB_MIN}-${HET_AB_MAX}, hom alt fraction < ${HOM_AB_MIN} is weak"
echo ""
if [ -z "$FREQ_NAME" ]; then
  echo "  NOTICE: the VEP output has no MAX_AF, gnomADe_AF or gnomADg_AF field, so no tier"
  echo "  is filtered by frequency and every MODERATE variant is kept. Run step 13 with"
  echo "  --everything (or --max_af) for population frequencies."
  echo ""
fi

# count FILE: records in a tier file.
count() {
  run_in "${BCFTOOLS_IMAGE}" bcftools view -H "${C}/$1" 2>/dev/null | wc -l | tr -d ' ' || echo 0
}

# rank_vcf OUT KEEP VIEW_ARGS...: `bcftools view VIEW_ARGS` through
# bin/csq_rank.awk on the host (ranked transcript, frequency tier, call check
# as PGP_* INFO tags), bgzipped and indexed as clinical/OUT. KEEP=noncommon
# drops the records with a known frequency of 1% or more.
rank_vcf() {
  local out=$1 keep=$2
  shift 2
  run_in --cpus 2 --memory 2g "${BCFTOOLS_IMAGE}" bcftools view "$@" \
    | LC_ALL=C awk -v keep="$keep" -v min_gq="$MIN_GQ" -v min_dp="$MIN_DP" \
        -v het_lo="$HET_AB_MIN" -v het_hi="$HET_AB_MAX" -v hom_lo="$HOM_AB_MIN" \
        -f "${PGP_ROOT}/bin/csq_rank.awk" \
    | run_in -i --cpus 2 --memory 2g "${BCFTOOLS_IMAGE}" bcftools view -Oz -o "${C}/${out}"
  run_in --cpus 2 --memory 2g "${BCFTOOLS_IMAGE}" bcftools index -f -t "${C}/${out}"
}
# tier OUT EXPR: the records of RARE that match EXPR, as clinical/OUT.
tier() {
  run_in --cpus 2 --memory 2g "${BCFTOOLS_IMAGE}" bcftools view -i "$2" -Oz -o "${C}/$1" "$RARE"
  run_in --cpus 2 --memory 2g "${BCFTOOLS_IMAGE}" bcftools index -f -t "${C}/$1"
}

# Step 2: the PASS records every non-ClinVar tier starts from, ranked, with
# the common ones left out. RARE_FREQ is "rare" for the HIGH and MODERATE
# tiers: a known frequency below 1%, or any record when the VEP output has no
# frequency field at all.
echo "[2] Ranking PASS records on their MANE Select transcript$([ -n "$FREQ_NAME" ] && echo ", leaving out those with ${FREQ_NAME} of 1% or more")..."
rank_vcf "${SAMPLE}_rare_pass.vcf.gz" noncommon -f PASS "${CONTAINER_INPUT}"
RARE="${C}/${SAMPLE}_rare_pass.vcf.gz"
RARE_FREQ='(INFO/PGP_FREQ="rare" || INFO/PGP_FREQ="none")'

# Step 3: HIGH impact (stop-gain, frameshift, splice site), rare
echo "[3] Extracting rare HIGH impact variants (stop-gain, frameshift, splice)..."
tier "${SAMPLE}_high_impact.vcf.gz" "INFO/PGP_IMPACT=\"HIGH\" && ${RARE_FREQ}"
HIGH_COUNT=$(count "${SAMPLE}_high_impact.vcf.gz")
echo "  Found: ${HIGH_COUNT} HIGH impact variants"

# Step 4: MODERATE impact (missense, in-frame indel), rare
echo "[4] Extracting rare MODERATE impact variants (missense, in-frame indel)..."
tier "${SAMPLE}_rare_moderate.vcf.gz" "INFO/PGP_IMPACT=\"MODERATE\" && ${RARE_FREQ}"
MODERATE_COUNT=$(count "${SAMPLE}_rare_moderate.vcf.gz")
echo "  Found: ${MODERATE_COUNT} MODERATE impact variants"
MERGE_FILES="${C}/${SAMPLE}_high_impact.vcf.gz ${C}/${SAMPLE}_rare_moderate.vcf.gz"

# Step 4b: HIGH and MODERATE with no known frequency. Not rare: an indel VEP
# could not match to gnomAD may be common. Its own tier and its own count.
FREQ_UNKNOWN_COUNT=0
FREQ_UNKNOWN="${OUTDIR}/${SAMPLE}_freq_unknown.vcf.gz"
rm -f "${FREQ_UNKNOWN:?}" "${FREQ_UNKNOWN:?}.tbi"
if [ -n "$FREQ_NAME" ]; then
  echo "[4b] Extracting HIGH and MODERATE variants with no known frequency (their own tier, not rare)..."
  tier "${SAMPLE}_freq_unknown.vcf.gz" '(INFO/PGP_IMPACT="HIGH" || INFO/PGP_IMPACT="MODERATE") && INFO/PGP_FREQ="unknown"'
  FREQ_UNKNOWN_COUNT=$(count "${SAMPLE}_freq_unknown.vcf.gz")
  echo "  Found: ${FREQ_UNKNOWN_COUNT} HIGH or MODERATE variants of unknown frequency"
  MERGE_FILES="${MERGE_FILES} ${C}/${SAMPLE}_freq_unknown.vcf.gz"
fi

# Step 5: ClinVar pathogenic/likely pathogenic, at any frequency (a common
# pathogenic allele such as HFE p.C282Y must stay). Step 6 screens the sample
# against the ClinVar file in clinvar/, which setup.sh refreshes; VEP's
# CLIN_SIG comes from its cache and is only the fallback.
CLINVAR_COUNT=0
CLINVAR_SOURCE="none"
CLINVAR_TIER="${OUTDIR}/${SAMPLE}_clinvar_pathogenic.vcf.gz"
CLINVAR_UNRANKED="${OUTDIR}/${SAMPLE}_clinvar_unranked.vcf.gz"
rm -f "${CLINVAR_TIER:?}" "${CLINVAR_TIER:?}.tbi" "${CLINVAR_UNRANKED:?}"
if [ "$HAS_HITS" -eq 1 ]; then
  echo "[5] Extracting the step 6 ClinVar hits..."
  CLINVAR_SOURCE="step 6 ClinVar screen"
  TARGETS="${OUTDIR}/${SAMPLE}_clinvar_targets.tsv"
  awk -F'\t' 'BEGIN {OFS = "\t"} !/^#/ {print $1, $2}' "$CLINVAR_HITS" | sort -u > "$TARGETS"
  if [ -s "$TARGETS" ]; then
    run_in --cpus 4 --memory 4g "${BCFTOOLS_IMAGE}" \
      bcftools view -f PASS -T "${C}/${SAMPLE}_clinvar_targets.tsv" -Oz -o "${C}/${SAMPLE}_clinvar_unranked.vcf.gz" "${CONTAINER_INPUT}"
  else
    echo "  Step 6 found no ClinVar hit."
  fi
  rm -f "${TARGETS:?}"
elif [ "$HAS_CUSTOM" -eq 1 ]; then
  echo "[5] Extracting ClinVar_CLNSIG pathogenic/likely pathogenic (step 13's --custom ClinVar file)..."
  CLINVAR_SOURCE="ClinVar_CLNSIG (step 13 --custom)"
  run_in --cpus 4 --memory 4g "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
      bcftools +split-vep - -c ClinVar_CLNSIG \
        -i 'ClinVar_CLNSIG~\"athogenic\" && ClinVar_CLNSIG!~\"onflicting\"' \
        -Oz -o ${C}/${SAMPLE}_clinvar_unranked.vcf.gz"
elif [ "$HAS_CLINSIG" -eq 1 ]; then
  echo "[5] Extracting VEP CLIN_SIG pathogenic/likely pathogenic (VEP's cached ClinVar; run step 6 for the current file)..."
  CLINVAR_SOURCE="VEP CLIN_SIG (cache release)"
  run_in --cpus 4 --memory 4g "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
      bcftools +split-vep - -c CLIN_SIG \
        -i 'CLIN_SIG~\"pathogenic\" && CLIN_SIG!~\"conflicting\"' \
        -Oz -o ${C}/${SAMPLE}_clinvar_unranked.vcf.gz"
else
  echo "[5] Skipping the ClinVar tier (no step 6 hits, no ClinVar_CLNSIG and no CLIN_SIG in the VEP output)."
fi
# The ClinVar records are ranked like the others, so the summary shows their
# MANE consequence and call check too.
if [ -f "$CLINVAR_UNRANKED" ]; then
  rank_vcf "${SAMPLE}_clinvar_pathogenic.vcf.gz" all "${C}/${SAMPLE}_clinvar_unranked.vcf.gz"
  rm -f "${CLINVAR_UNRANKED:?}"
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
  tier "${SAMPLE}_cadd_high.vcf.gz" "INFO/PGP_IMPACT!=\"HIGH\" && INFO/PGP_IMPACT!=\"MODERATE\" && (${CADD_EXPR})"
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

# Summary TSV: one row per variant. IMPACT, GENE and Consequence are the
# ranked transcript's (TRANSCRIPT: MANE, CANONICAL or ANY), WORST_IMPACT and
# WORST_Consequence the worst over every transcript. FREQ_TIER, CALL_CHECK and
# CALL_WHY come from bin/csq_rank.awk. A score, frequency or FORMAT field the
# input does not carry is written as '.'.
echo ""
echo "Generating human-readable summary..."
COL_CADD='.'; [ "$HAS_CADD" -eq 1 ] && COL_CADD='%INFO/CADD_PHRED'
COL_REVEL='.'; [ "$HAS_REVEL" -eq 1 ] && COL_REVEL='%INFO/REVEL'
COL_AM='.'; has_info AM_class && COL_AM='%INFO/AM_class'
COL_FREQ='.'; has_field MAX_AF && COL_FREQ='%MAX_AF'
COL_GQ='.'; has_format GQ && COL_GQ='%GQ'
COL_DP='.'; has_format DP && COL_DP='%DP'
COL_AD='.'; has_format AD && COL_AD='%AD'
SUMMARY="${OUTDIR}/${SAMPLE}_clinical_summary.tsv"
printf 'CHROM\tPOS\tREF\tALT\tGT\tIMPACT\tGENE\tConsequence\tMAX_AF\tCADD_PHRED\tREVEL\tAM_CLASS\tTRANSCRIPT\tWORST_IMPACT\tWORST_Consequence\tFREQ_TIER\tGQ\tDP\tAD\tCALL_CHECK\tCALL_WHY\n' > "${SUMMARY}.tmp"
run_in --cpus 2 --memory 2g "${BCFTOOLS_IMAGE}" \
  bcftools +split-vep "${C}/${SAMPLE}_clinical.vcf.gz" -s worst \
    -f "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\t%INFO/PGP_IMPACT\t%INFO/PGP_GENE\t%INFO/PGP_CSQ\t${COL_FREQ}\t${COL_CADD}\t${COL_REVEL}\t${COL_AM}\t%INFO/PGP_TX\t%IMPACT\t%Consequence\t%INFO/PGP_FREQ[\t${COL_GQ}\t${COL_DP}\t${COL_AD}]\t%INFO/PGP_CALL\t%INFO/PGP_CALL_WHY\n" \
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
WEAK_COUNT=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $c["CALL_CHECK"] == "weak" {n++} END {print n + 0}' "$SUMMARY")

echo ""
echo "============================================"
echo "  Clinical filter complete: ${SAMPLE}"
echo "  Total clinically interesting variants: ${TOTAL_COUNT}"
echo "    HIGH impact (LoF):           ${HIGH_COUNT}"
echo "    Rare MODERATE impact:        ${MODERATE_COUNT}"
[ -n "$FREQ_NAME" ] && \
echo "    Frequency unknown (not rare): ${FREQ_UNKNOWN_COUNT}"
echo "    ClinVar pathogenic/LP:       ${CLINVAR_COUNT} (${CLINVAR_SOURCE})"
{ [ "$HAS_CADD" -eq 1 ] || [ "$HAS_CADD_INDEL" -eq 1 ]; } && \
echo "    High CADD non-coding:        ${CADD_COUNT}"
{ [ "$HAS_SPLICEAI" -eq 1 ] || [ "$HAS_SPLICEAI_INDEL" -eq 1 ]; } && \
echo "    Cryptic splice (SpliceAI):   ${SPLICEAI_COUNT}"
{ [ "$HAS_REVEL" -eq 1 ] || [ "$HAS_ALPHAMISSENSE" -eq 1 ]; } && \
echo "    Deleterious missense:        ${MISSENSE_COUNT}"
echo "  Frequency filter: ${FREQ_NAME:-none (no frequency field in the VEP output)}"
echo "  Weak genotype calls (CALL_CHECK=weak, check the reads): ${WEAK_COUNT}"
echo ""
echo "  Output files:"
echo "    ${OUTDIR}/${SAMPLE}_clinical.vcf.gz              (combined)"
echo "    ${OUTDIR}/${SAMPLE}_clinical_summary.tsv         (human-readable table)"
echo "    ${OUTDIR}/${SAMPLE}_high_impact.vcf.gz           (HIGH only)"
echo "    ${OUTDIR}/${SAMPLE}_rare_moderate.vcf.gz         (rare MODERATE only)"
[ -f "$FREQ_UNKNOWN" ] && \
echo "    ${FREQ_UNKNOWN}  (HIGH/MODERATE, frequency unknown)"
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
