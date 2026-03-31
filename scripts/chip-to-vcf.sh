#!/usr/bin/env bash
# Chip-to-VCF Converter — converts consumer genotyping array data to GRCh38 VCF
# Supports 23andMe, AncestryDNA, and MyHeritage raw data formats
#
# This script uses bcftools convert --tsv2vcf (NOT plink) because plink's binary
# format cannot represent both alleles for monomorphic single-sample sites,
# silently corrupting all homozygous ALT genotypes.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> [format]}
FORMAT=${2:-auto}  # auto, 23andme, myheritage, ancestrydna
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}

RAW_DIR="${GENOME_DIR}/${SAMPLE}/raw"
VCF_DIR="${GENOME_DIR}/${SAMPLE}/vcf"
REF_HG19="${GENOME_DIR}/reference_hg19/human_g1k_v37.fasta"
REF_HG38="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
CHAIN="${GENOME_DIR}/liftover/hg19ToHg38.over.chain.gz"

echo "=== Chip-to-VCF Converter: ${SAMPLE} ==="
echo "Format: ${FORMAT}"

# --- Validate prerequisites ---
for f in "$REF_HG19" "${REF_HG19}.fai" "$REF_HG38" "$CHAIN"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: Required file not found: ${f}" >&2
    echo "  Run the chip data prerequisite downloads first." >&2
    echo "  See docs/chip-data-guide.md for instructions." >&2
    exit 1
  fi
done

mkdir -p "$RAW_DIR" "$VCF_DIR"

# --- Step 0: Detect and normalize input format ---
# All formats need to become a TSV with columns: rsid chr pos genotype
RAW_TSV="${RAW_DIR}/${SAMPLE}_raw.txt"

if [ "$FORMAT" = "auto" ]; then
  # Auto-detect based on available files
  if [ -f "${RAW_DIR}/MyHeritage_raw_dna_data.csv" ]; then
    FORMAT="myheritage"
  elif [ -f "${RAW_TSV}" ]; then
    # Check if it's 23andMe or AncestryDNA by peeking at the header
    if head -1 "$RAW_TSV" | grep -qi "ancestrydna"; then
      FORMAT="ancestrydna"
    else
      FORMAT="23andme"
    fi
  else
    echo "ERROR: No raw data file found in ${RAW_DIR}/" >&2
    echo "  Expected one of:" >&2
    echo "    ${RAW_DIR}/MyHeritage_raw_dna_data.csv" >&2
    echo "    ${RAW_DIR}/${SAMPLE}_raw.txt" >&2
    exit 1
  fi
  echo "Auto-detected format: ${FORMAT}"
fi

case "$FORMAT" in
  myheritage)
    INPUT="${RAW_DIR}/MyHeritage_raw_dna_data.csv"
    if [ ! -f "$INPUT" ]; then
      echo "ERROR: MyHeritage file not found: ${INPUT}" >&2
      exit 1
    fi
    echo "Converting MyHeritage CSV to TSV..."
    grep -v "^#" "$INPUT" | \
      grep -v "^RSID" | \
      sed 's/"//g' | \
      awk -F',' '{print $1"\t"$2"\t"$3"\t"$4}' \
      > "$RAW_TSV"
    ;;
  23andme|ancestrydna)
    if [ ! -f "$RAW_TSV" ]; then
      echo "ERROR: Raw data file not found: ${RAW_TSV}" >&2
      echo "  Place your 23andMe/AncestryDNA file at this path." >&2
      exit 1
    fi
    ;;
  *)
    echo "ERROR: Unknown format '${FORMAT}'. Use: auto, 23andme, myheritage, ancestrydna" >&2
    exit 1
    ;;
esac

VARIANT_COUNT=$(grep -c -v "^#" "$RAW_TSV" || true)
echo "Input: ${VARIANT_COUNT} genotyped positions"

# --- Step 1: Convert to hg19 VCF with proper REF/ALT ---
echo ""
echo "--- Stage 1: Converting to hg19 VCF (bcftools convert --tsv2vcf) ---"
echo "  This looks up the reference allele at each position from the FASTA."
echo "  Homozygous ALT genotypes will be correctly encoded as GT 1/1."

docker run --rm --user root --cpus 2 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools convert --tsv2vcf "/genome/${SAMPLE}/raw/${SAMPLE}_raw.txt" \
    -f /genome/reference_hg19/human_g1k_v37.fasta \
    -s "${SAMPLE}" \
    -c ID,CHROM,POS,AA \
    -Oz -o "/genome/${SAMPLE}/raw/${SAMPLE}_hg19.vcf.gz"

# --- Step 2: Add chr prefix for liftover chain compatibility ---
echo ""
echo "--- Adding chr prefix to chromosome names ---"

CHR_RENAME="${GENOME_DIR}/reference_hg19/chr_rename.txt"
if [ ! -f "$CHR_RENAME" ]; then
  printf '%s\n' $(seq 1 22) X Y MT | \
    awk '{print $1" chr"$1}' > "$CHR_RENAME"
fi

docker run --rm --user root --cpus 2 --memory 2g \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools annotate \
    --rename-chrs /genome/reference_hg19/chr_rename.txt \
    "/genome/${SAMPLE}/raw/${SAMPLE}_hg19.vcf.gz" \
    -Oz -o "/genome/${SAMPLE}/raw/${SAMPLE}_hg19_chr.vcf.gz"

docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools index -t "/genome/${SAMPLE}/raw/${SAMPLE}_hg19_chr.vcf.gz"

# --- Step 3: Liftover to GRCh38 ---
echo ""
echo "--- Stage 2: Liftover to GRCh38 (Picard LiftoverVcf) ---"

docker run --rm --user root --cpus 2 --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  broadinstitute/picard:latest \
  java -jar /usr/picard/picard.jar LiftoverVcf \
    I="/genome/${SAMPLE}/raw/${SAMPLE}_hg19_chr.vcf.gz" \
    O="/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
    CHAIN=/genome/liftover/hg19ToHg38.over.chain.gz \
    R=/genome/reference/Homo_sapiens_assembly38.fasta \
    REJECT="/genome/${SAMPLE}/raw/${SAMPLE}_liftover_rejected.vcf.gz" \
    WARN_ON_MISSING_CONTIG=true

# --- Step 4: Index the final VCF ---
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools index -t -f "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"

# --- Summary ---
echo ""
echo "=== Conversion complete ==="
echo "  Output VCF: ${VCF_DIR}/${SAMPLE}.vcf.gz"

docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools stats "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" 2>/dev/null | \
  grep "^SN" | sed 's/^SN\t0\t/  /'

REJECTED_COUNT=$(docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools view -H "/genome/${SAMPLE}/raw/${SAMPLE}_liftover_rejected.vcf.gz" 2>/dev/null | wc -l || echo "0")
echo "  Liftover rejected: ${REJECTED_COUNT} variants"
echo ""
echo "You can now run pipeline steps 6, 7, 11, 25, and 27 on this VCF."
echo "See docs/chip-data-guide.md for which steps work and their limitations."
