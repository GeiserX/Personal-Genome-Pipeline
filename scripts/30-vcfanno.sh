#!/usr/bin/env bash
# vcfanno — Annotate VEP output with CADD, SpliceAI, REVEL, and AlphaMissense scores
# Input: VEP-annotated VCF from step 13
# Output: Fully annotated VCF with additional pathogenicity scores in INFO field
# Runtime: ~5-15 minutes depending on variant count
#
# Handles chromosome naming mismatch: CADD files use bare names (1, 2, 3)
# while the VCF and other databases use chr-prefixed names (chr1, chr2, chr3).
# Solved via two-pass annotation with chromosome renaming between passes.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=../versions.env
. "${SCRIPT_DIR}/../versions.env"

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}

# Validate sample name to prevent shell injection in bash -c strings
if [[ "$SAMPLE" =~ [^a-zA-Z0-9._-] ]]; then
  echo "ERROR: Sample name contains invalid characters. Use only a-z, A-Z, 0-9, ., _, -" >&2
  exit 1
fi

SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
VEP_DIR="${SAMPLE_DIR}/vep"
ANNOT_DIR="${GENOME_DIR}/annotations"
OUTPUT_FILE="${VEP_DIR}/${SAMPLE}_annotated.vcf.gz"

echo "=== vcfanno Annotation: ${SAMPLE} ==="

# Skip if output already exists
if [ -f "$OUTPUT_FILE" ]; then
  echo "Output already exists: ${OUTPUT_FILE}"
  echo "Skipping. Delete the file to re-run."
  exit 0
fi

# --- Find VEP input (compressed or uncompressed) ---
VEP_VCF=""
if [ -f "${VEP_DIR}/${SAMPLE}_vep.vcf.gz" ]; then
  VEP_VCF="${VEP_DIR}/${SAMPLE}_vep.vcf.gz"
elif [ -f "${VEP_DIR}/${SAMPLE}_vep.vcf" ]; then
  VEP_VCF="${VEP_DIR}/${SAMPLE}_vep.vcf"
else
  echo "ERROR: VEP-annotated VCF not found." >&2
  echo "  Expected: ${VEP_DIR}/${SAMPLE}_vep.vcf[.gz]" >&2
  echo "  Run step 13 (VEP annotation) first." >&2
  exit 1
fi

echo "Input VCF: ${VEP_VCF}"

# --- bgzip and index if uncompressed ---
if [[ "$VEP_VCF" == *.vcf ]] && [ ! -f "${VEP_VCF}.gz" ]; then
  echo "Compressing VEP output with bgzip..."
  docker run --rm --user root \
    --cpus 2 --memory 2g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bgzip -c "/genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf" \
    > "${VEP_DIR}/${SAMPLE}_vep.vcf.gz"
  VEP_VCF="${VEP_DIR}/${SAMPLE}_vep.vcf.gz"
fi

if [ ! -f "${VEP_VCF}.tbi" ]; then
  echo "Indexing VEP VCF..."
  docker run --rm --user root \
    --cpus 2 --memory 2g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    tabix -p vcf "/genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz"
fi

echo "Input VCF (indexed): ${VEP_VCF}"

# --- Discover available annotation tracks ---
# CADD (no chr prefix)
CADD_SNV="${ANNOT_DIR}/whole_genome_SNVs.tsv.gz"
CADD_INDEL="${ANNOT_DIR}/gnomad.genomes.r4.0.indel.tsv.gz"
# SpliceAI (chr prefix)
SPLICEAI_SNV="${ANNOT_DIR}/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL="${ANNOT_DIR}/spliceai_scores.raw.indel.hg38.vcf.gz"
# REVEL (chr prefix)
REVEL="${ANNOT_DIR}/revel_grch38.tsv.gz"
# AlphaMissense (chr prefix)
ALPHAMISSENSE="${ANNOT_DIR}/AlphaMissense_hg38.tsv.gz"

NOCHR_TRACKS=()   # Tracks that use bare chromosome names (1, 2, 3)
CHR_TRACKS=()     # Tracks that use chr-prefixed names (chr1, chr2, chr3)
APPLIED_NAMES=()  # Human-readable names for summary

# Check CADD SNVs
if [ -f "$CADD_SNV" ] && [ -f "${CADD_SNV}.tbi" ]; then
  NOCHR_TRACKS+=("cadd_snv")
  APPLIED_NAMES+=("CADD SNVs")
  echo "  Found: CADD SNVs (whole_genome_SNVs.tsv.gz)"
else
  echo "  Skipping: CADD SNVs (not found at ${CADD_SNV})"
fi

# Check CADD indels
if [ -f "$CADD_INDEL" ] && [ -f "${CADD_INDEL}.tbi" ]; then
  NOCHR_TRACKS+=("cadd_indel")
  APPLIED_NAMES+=("CADD indels")
  echo "  Found: CADD indels (gnomad.genomes.r4.0.indel.tsv.gz)"
else
  echo "  Skipping: CADD indels (not found at ${CADD_INDEL})"
fi

# Check SpliceAI SNVs
if [ -f "$SPLICEAI_SNV" ] && [ -f "${SPLICEAI_SNV}.tbi" ]; then
  CHR_TRACKS+=("spliceai_snv")
  APPLIED_NAMES+=("SpliceAI SNVs")
  echo "  Found: SpliceAI SNVs"
else
  echo "  Skipping: SpliceAI SNVs (not found at ${SPLICEAI_SNV})"
fi

# Check SpliceAI indels
if [ -f "$SPLICEAI_INDEL" ] && [ -f "${SPLICEAI_INDEL}.tbi" ]; then
  CHR_TRACKS+=("spliceai_indel")
  APPLIED_NAMES+=("SpliceAI indels")
  echo "  Found: SpliceAI indels"
else
  echo "  Skipping: SpliceAI indels (not found at ${SPLICEAI_INDEL})"
fi

# Check REVEL
if [ -f "$REVEL" ] && [ -f "${REVEL}.tbi" ]; then
  CHR_TRACKS+=("revel")
  APPLIED_NAMES+=("REVEL")
  echo "  Found: REVEL"
else
  echo "  Skipping: REVEL (not found at ${REVEL})"
fi

# Check AlphaMissense
if [ -f "$ALPHAMISSENSE" ] && [ -f "${ALPHAMISSENSE}.tbi" ]; then
  CHR_TRACKS+=("alphamissense")
  APPLIED_NAMES+=("AlphaMissense")
  echo "  Found: AlphaMissense"
else
  echo "  Skipping: AlphaMissense (not found at ${ALPHAMISSENSE})"
fi

# Exit gracefully if nothing to annotate
if [ ${#NOCHR_TRACKS[@]} -eq 0 ] && [ ${#CHR_TRACKS[@]} -eq 0 ]; then
  echo ""
  echo "No annotation databases found in ${ANNOT_DIR}/."
  echo "Download them as described in docs/00-reference-setup.md and docs/30-vcfanno.md."
  echo "Copying VEP output as-is."
  cp "${VEP_VCF}" "${OUTPUT_FILE}"
  cp "${VEP_VCF}.tbi" "${OUTPUT_FILE}.tbi"
  exit 0
fi

echo ""

# --- Helper: Generate TOML config ---
# vcfanno uses TOML config to define annotation sources.
# We generate it dynamically based on which files are present.

generate_nochr_toml() {
  local toml=""
  for track in "${NOCHR_TRACKS[@]}"; do
    case "$track" in
      cadd_snv)
        toml+='[[annotation]]
file="/genome/annotations/whole_genome_SNVs.tsv.gz"
columns=[6]
names=["CADD_PHRED"]
ops=["self"]

'
        ;;
      cadd_indel)
        toml+='[[annotation]]
file="/genome/annotations/gnomad.genomes.r4.0.indel.tsv.gz"
columns=[6]
names=["CADD_PHRED_indel"]
ops=["self"]

'
        ;;
    esac
  done
  echo "$toml"
}

generate_chr_toml() {
  local toml=""
  for track in "${CHR_TRACKS[@]}"; do
    case "$track" in
      spliceai_snv)
        toml+='[[annotation]]
file="/genome/annotations/spliceai_scores.raw.snv.hg38.vcf.gz"
fields=["SpliceAI"]
names=["SpliceAI"]
ops=["self"]

'
        ;;
      spliceai_indel)
        toml+='[[annotation]]
file="/genome/annotations/spliceai_scores.raw.indel.hg38.vcf.gz"
fields=["SpliceAI"]
names=["SpliceAI_indel"]
ops=["self"]

'
        ;;
      revel)
        toml+='[[annotation]]
file="/genome/annotations/revel_grch38.tsv.gz"
columns=[5]
names=["REVEL"]
ops=["self"]

'
        ;;
      alphamissense)
        toml+='[[annotation]]
file="/genome/annotations/AlphaMissense_hg38.tsv.gz"
columns=[9,10]
names=["AM_pathogenicity","AM_class"]
ops=["self","self"]

'
        ;;
    esac
  done
  echo "$toml"
}

# --- Create working directory for TOML configs and temp files ---
WORK_DIR="${VEP_DIR}/vcfanno_tmp"
mkdir -p "$WORK_DIR"
trap 'rm -rf "$WORK_DIR"' EXIT

# --- Pass 1: CADD annotation (no-chr tracks) ---
# CADD files use bare chromosome names (1, 2, 3).
# Strategy: strip chr prefix from VCF -> annotate -> re-add chr prefix.
CURRENT_VCF="/genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz"

if [ ${#NOCHR_TRACKS[@]} -gt 0 ]; then
  echo "=== Pass 1/2: Annotating with CADD (bare chromosome names) ==="

  # Generate chromosome rename files
  # strip_chr.txt: chr1 1, chr2 2, ... (for bcftools annotate --rename-chrs)
  # add_chr.txt:   1 chr1, 2 chr2, ... (reverse)
  for i in $(seq 1 22) X Y M; do
    echo "chr${i} ${i}" >> "${WORK_DIR}/strip_chr.txt"
    echo "${i} chr${i}" >> "${WORK_DIR}/add_chr.txt"
  done

  # Write TOML for CADD tracks
  generate_nochr_toml > "${WORK_DIR}/nochr.toml"

  # Run: strip chr -> vcfanno -> add chr back
  # All within a single Docker invocation using bcftools + vcfanno
  # Step 1a: Strip chr prefix
  echo "  Stripping chr prefix from VCF..."
  docker run --rm --user root \
    --cpus 2 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bash -c "
      bcftools annotate --rename-chrs /genome/${SAMPLE}/vep/vcfanno_tmp/strip_chr.txt \
        ${CURRENT_VCF} \
        -Oz -o /genome/${SAMPLE}/vep/vcfanno_tmp/nochr_input.vcf.gz && \
      tabix -p vcf /genome/${SAMPLE}/vep/vcfanno_tmp/nochr_input.vcf.gz
    "

  # Step 1b: Run vcfanno with CADD
  echo "  Running vcfanno with CADD tracks..."
  docker run --rm --user root \
    --cpus 4 --memory 8g \
    -v "${GENOME_DIR}:/genome" \
    "${VCFANNO_IMAGE}" \
    vcfanno -p 4 \
      "/genome/${SAMPLE}/vep/vcfanno_tmp/nochr.toml" \
      "/genome/${SAMPLE}/vep/vcfanno_tmp/nochr_input.vcf.gz" \
    > "${WORK_DIR}/nochr_annotated.vcf"

  # Step 1c: Re-add chr prefix and compress
  echo "  Re-adding chr prefix..."
  docker run --rm --user root \
    --cpus 2 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bash -c "
      bgzip -c /genome/${SAMPLE}/vep/vcfanno_tmp/nochr_annotated.vcf \
        > /genome/${SAMPLE}/vep/vcfanno_tmp/nochr_annotated.vcf.gz && \
      bcftools annotate --rename-chrs /genome/${SAMPLE}/vep/vcfanno_tmp/add_chr.txt \
        /genome/${SAMPLE}/vep/vcfanno_tmp/nochr_annotated.vcf.gz \
        -Oz -o /genome/${SAMPLE}/vep/vcfanno_tmp/pass1_output.vcf.gz && \
      tabix -p vcf /genome/${SAMPLE}/vep/vcfanno_tmp/pass1_output.vcf.gz
    "

  CURRENT_VCF="/genome/${SAMPLE}/vep/vcfanno_tmp/pass1_output.vcf.gz"
  echo "  CADD annotation complete."
  echo ""
fi

# --- Pass 2: chr-prefixed tracks (SpliceAI, REVEL, AlphaMissense) ---
if [ ${#CHR_TRACKS[@]} -gt 0 ]; then
  PASS_LABEL="2/2"
  if [ ${#NOCHR_TRACKS[@]} -eq 0 ]; then
    PASS_LABEL="1/1"
  fi
  echo "=== Pass ${PASS_LABEL}: Annotating with SpliceAI, REVEL, AlphaMissense ==="

  # Write TOML for chr-prefixed tracks
  generate_chr_toml > "${WORK_DIR}/chr.toml"

  # Run vcfanno
  echo "  Running vcfanno with chr-prefixed tracks..."
  docker run --rm --user root \
    --cpus 4 --memory 8g \
    -v "${GENOME_DIR}:/genome" \
    "${VCFANNO_IMAGE}" \
    vcfanno -p 4 \
      "/genome/${SAMPLE}/vep/vcfanno_tmp/chr.toml" \
      "${CURRENT_VCF}" \
    > "${WORK_DIR}/pass2_output.vcf"

  CURRENT_VCF="/genome/${SAMPLE}/vep/vcfanno_tmp/pass2_output.vcf"
  echo "  Annotation complete."
  echo ""
elif [ ${#NOCHR_TRACKS[@]} -gt 0 ]; then
  # Only CADD was annotated, pass1 output is the final
  # Need to decompress for the final bgzip step below
  docker run --rm --user root \
    --cpus 2 --memory 2g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bcftools view "${CURRENT_VCF}" \
    > "${WORK_DIR}/pass2_output.vcf"
  CURRENT_VCF="/genome/${SAMPLE}/vep/vcfanno_tmp/pass2_output.vcf"
fi

# --- Compress and index final output ---
echo "=== Compressing and indexing output ==="
docker run --rm --user root \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}:/genome" \
  "${BCFTOOLS_IMAGE}" \
  bash -c "
    bgzip -c /genome/${SAMPLE}/vep/vcfanno_tmp/pass2_output.vcf \
      > /genome/${SAMPLE}/vep/${SAMPLE}_annotated.vcf.gz && \
    tabix -p vcf /genome/${SAMPLE}/vep/${SAMPLE}_annotated.vcf.gz
  "

# --- Summary ---
echo ""
echo "=== vcfanno annotation complete ==="
echo "Output: ${OUTPUT_FILE}"
echo ""
echo "Annotation tracks applied:"
for name in "${APPLIED_NAMES[@]}"; do
  echo "  - ${name}"
done

# Count variants in output
VARIANT_COUNT=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -H "/genome/${SAMPLE}/vep/${SAMPLE}_annotated.vcf.gz" \
  2>/dev/null | wc -l || echo "0")
echo ""
echo "Total variants in output: ${VARIANT_COUNT}"
echo ""
echo "Query examples:"
echo "  # Variants with CADD PHRED >= 20 (top 1% most deleterious):"
echo "  bcftools view -i 'INFO/CADD_PHRED>=20' ${OUTPUT_FILE} | head"
echo ""
echo "  # Variants with high SpliceAI score (>= 0.5):"
echo "  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SpliceAI\n' -i 'INFO/SpliceAI!=\".\"' ${OUTPUT_FILE} | head"
echo ""
echo "  # AlphaMissense likely pathogenic variants:"
echo "  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AM_pathogenicity\t%INFO/AM_class\n' -i 'INFO/AM_class=\"likely_pathogenic\"' ${OUTPUT_FILE} | head"
