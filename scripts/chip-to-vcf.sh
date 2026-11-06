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
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"

RAW_DIR="${GENOME_DIR}/${SAMPLE}/raw"
VCF_DIR="${GENOME_DIR}/${SAMPLE}/vcf"
REF_HG19="${GENOME_DIR}/reference_hg19/human_g1k_v37.fasta"
REF_HG38="$REF_FASTA"
CHAIN="${GENOME_DIR}/liftover/hg19ToHg38.over.chain.gz"

echo "=== Chip-to-VCF Converter: ${SAMPLE} ==="
echo "Format: ${FORMAT}"

# --- Validate prerequisites ---
# Picard LiftoverVcf reads the GRCh38 sequence dictionary (setup.sh creates it);
# the chrM records take their contig lines from the GRCh38 .fai.
for f in "$REF_HG19" "${REF_HG19}.fai" "$REF_HG38" "${REF_HG38}.fai" "$REF_DICT" "$CHAIN"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: Required file not found: ${f}" >&2
    echo "  Run the chip data prerequisite downloads first." >&2
    echo "  See docs/chip-data-guide.md for instructions." >&2
    exit 1
  fi
done

mkdir -p "$RAW_DIR" "$VCF_DIR"

# --- Step 0: Detect and normalize input format ---
# All formats become one TSV with the columns bcftools reads (-c ID,CHROM,POS,AA):
# rsid, chromosome, position, genotype as two letters.
RAW_TSV="${RAW_DIR}/${SAMPLE}_raw.txt"
CHIP_TSV="${RAW_DIR}/${SAMPLE}_tsv2vcf.tsv"

if [ "$FORMAT" = "auto" ]; then
  # Auto-detect based on available files
  if [ -f "${RAW_DIR}/MyHeritage_raw_dna_data.csv" ]; then
    FORMAT="myheritage"
  elif [ -f "${RAW_TSV}" ]; then
    # AncestryDNA names itself in its comment header and has a column header
    # with allele1 and allele2; 23andMe has neither.
    if head -n 40 "$RAW_TSV" | grep -qiE 'ancestrydna|^rsid[[:space:]]+chromosome[[:space:]]+position[[:space:]]+allele1'; then
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
    # MyHeritage: a comment header, a column header row (RSID, CHROMOSOME,
    # POSITION, RESULT), then quoted comma-separated rows, often with
    # Windows (CRLF) line endings. The carriage return is dropped before
    # anything else, or it would stay on the genotype column.
    echo "Converting MyHeritage CSV to TSV..."
    tr -d '\r' < "$INPUT" | sed 's/"//g' | awk -F',' -v OFS='\t' '
      /^#/ || NF == 0 || tolower($1) == "rsid" { next }
      { print $1, $2, $3, $4 }' > "$CHIP_TSV"
    ;;
  23andme)
    if [ ! -f "$RAW_TSV" ]; then
      echo "ERROR: Raw data file not found: ${RAW_TSV}" >&2
      echo "  Place your 23andMe file at this path." >&2
      exit 1
    fi
    # rsid, chromosome, position, genotype already; drop the comment header.
    grep -v '^#' "$RAW_TSV" | tr -d '\r' > "$CHIP_TSV"
    # 23andMe numbers its MT sites on hg19's chrM (the Yoruba sequence,
    # 16,571 bases), not on the rCRS that stage 1 reads REF from: its v5 site
    # list puts m.3243 at 3244. Each MT position moves to the rCRS through the
    # chrM chain of the hg19ToHg38 chain file (hg38's chrM is the rCRS). A
    # site inside one of that chain's gaps has no rCRS position and is
    # dropped, with its count printed.
    MT_BLOCKS="${RAW_DIR}/${SAMPLE}_chrM_blocks.tsv"
    gzip -dc "$CHAIN" | awk -v OFS='\t' '
      $1 == "chain" { inm = ($3 == "chrM" && $8 == "chrM" && $5 == "+" && $10 == "+"); t = $6; q = $11; next }
      inm && $1 ~ /^[0-9]+$/ { print t, q, $1; t += $1 + $2; q += $1 + $3 }' > "$MT_BLOCKS"
    awk -F'\t' -v OFS='\t' -v blocks="$MT_BLOCKS" '
      BEGIN { while ((getline l < blocks) > 0) { split(l, b, "\t"); n++; T[n] = b[1]; Q[n] = b[2]; S[n] = b[3] } }
      $2 == "MT" {
        mt++
        p = $3 + 0; to = 0
        for (i = 1; i <= n; i++) if (p > T[i] && p <= T[i] + S[i]) { to = Q[i] + p - T[i]; break }
        if (!to) { gap++; next }
        $3 = to
      }
      { print }
      END {
        if (mt && !n) printf "  23andMe MT: no chrM chain in the chain file; %d MT site(s) dropped\n", mt > "/dev/stderr"
        else if (mt) printf "  23andMe MT: %d site(s) moved from hg19 chrM to rCRS positions, %d in a chain gap dropped\n", mt - gap, gap > "/dev/stderr"
      }' "$CHIP_TSV" > "${CHIP_TSV}.rcrs"
    mv -f "${CHIP_TSV}.rcrs" "$CHIP_TSV"
    ;;
  ancestrydna)
    if [ ! -f "$RAW_TSV" ]; then
      echo "ERROR: Raw data file not found: ${RAW_TSV}" >&2
      echo "  Place your AncestryDNA file at this path." >&2
      exit 1
    fi
    # AncestryDNA: a comment header, a column header row, then five columns
    # (rsid, chromosome, position, allele1, allele2). Chromosomes are numbers:
    # 23 is X, 24 is Y, 25 the X pseudoautosomal regions (X positions) and 26
    # the mitochondrion. A no-call is allele 0. The two alleles are joined into
    # the two-letter genotype bcftools reads; a no-call becomes "--", which
    # bcftools writes as a missing genotype (./.), as it does for 23andMe.
    # A data row without exactly five columns or with an empty allele (a cut
    # or damaged file) stops the conversion with its line number: skipping it
    # would drop a genotype, and one allele would become a haploid call.
    echo "Converting AncestryDNA (five columns, numeric chromosomes) to TSV..."
    tr -d '\r' < "$RAW_TSV" | awk -F'\t' -v OFS='\t' -v file="$RAW_TSV" '
      /^#/ || NF == 0 || tolower($1) == "rsid" { next }
      NF != 5 || $4 == "" || $5 == "" {
        printf "ERROR: %s line %d: want five tab-separated columns (rsid, chromosome, position, allele1, allele2) with both alleles, found %d column(s): %s\n", file, NR, NF, $0 > "/dev/stderr"
        exit 1
      }
      {
        c = $2
        if (c == "23" || c == "25") c = "X"
        else if (c == "24") c = "Y"
        else if (c == "26") c = "MT"
        a = $4; b = $5
        if (a == "0" || b == "0") { a = "-"; b = "-" }
        print $1, c, $3, a b
      }' > "$CHIP_TSV"
    ;;
  *)
    echo "ERROR: Unknown format '${FORMAT}'. Use: auto, 23andme, myheritage, ancestrydna" >&2
    exit 1
    ;;
esac

# Sorted by chromosome, then position: bcftools writes rows in input order and
# the index needs each chromosome in one block. AncestryDNA lists its X
# pseudoautosomal rows (25) after Y (24).
LC_ALL=C sort -t "$(printf '\t')" -k2,2 -k3,3n "$CHIP_TSV" > "${CHIP_TSV}.sorted"
mv -f "${CHIP_TSV}.sorted" "$CHIP_TSV"

VARIANT_COUNT=$(grep -c . "$CHIP_TSV" || true)
echo "Input: ${VARIANT_COUNT} genotyped positions"

# --- Step 1: Convert to hg19 VCF with proper REF/ALT ---
echo ""
echo "--- Stage 1: Converting to hg19 VCF (bcftools convert --tsv2vcf) ---"
echo "  This looks up the reference allele at each position from the FASTA."
echo "  Homozygous ALT genotypes will be correctly encoded as GT 1/1."

# bcftools writes the hg19 FASTA's .fai next to it when it is missing, so that
# directory is writable here.
run_in --rw "$(dirname "$REF_HG19")" --cpus 2 --memory 4g \
  "${BCFTOOLS_IMAGE}" \
  bcftools convert --tsv2vcf "/genome/${SAMPLE}/raw/${SAMPLE}_tsv2vcf.tsv" \
    -f /genome/reference_hg19/human_g1k_v37.fasta \
    -s "${SAMPLE}" \
    -c ID,CHROM,POS,AA \
    -Oz -o "/genome/${SAMPLE}/raw/${SAMPLE}_hg19.vcf.gz"

# --- Step 2: Add chr prefix for liftover chain compatibility ---
echo ""
echo "--- Adding chr prefix to chromosome names ---"

# Build chr-prefix rename map. GRCh38 uses chrM (not chrMT) for mitochondria.
CHR_RENAME="${GENOME_DIR}/reference_hg19/chr_rename.txt"
{
  seq 1 22 | awk '{print $1" chr"$1}'
  echo "X chrX"
  echo "Y chrY"
  echo "MT chrM"
} > "$CHR_RENAME"

run_in --cpus 2 --memory 2g \
  "${BCFTOOLS_IMAGE}" \
  bcftools annotate \
    --rename-chrs /genome/reference_hg19/chr_rename.txt \
    "/genome/${SAMPLE}/raw/${SAMPLE}_hg19.vcf.gz" \
    -Oz -o "/genome/${SAMPLE}/raw/${SAMPLE}_hg19_chr.vcf.gz"

run_in "${BCFTOOLS_IMAGE}" \
  bcftools index -f -t "/genome/${SAMPLE}/raw/${SAMPLE}_hg19_chr.vcf.gz"

# The mitochondrion is not lifted. GRCh37's MT (the human_g1k_v37 FASTA that
# stage 1 read REF from) and GRCh38's chrM are the same sequence, the rCRS,
# so its positions stay as they are. The UCSC chain maps hg19's chrM instead
# (the Yoruba sequence, 16,571 bases): it would move rCRS positions one or two
# bases (m.1555 to 1553) and keep the record where the base happens to match.
# So chrM goes around Picard and is merged back below, with GRCh38's contigs.
run_in "${BCFTOOLS_IMAGE}" \
  bcftools view -t ^chrM "/genome/${SAMPLE}/raw/${SAMPLE}_hg19_chr.vcf.gz" \
    -Oz -o "/genome/${SAMPLE}/raw/${SAMPLE}_hg19_nuclear.vcf.gz"
run_in "${BCFTOOLS_IMAGE}" \
  bcftools view -t chrM "/genome/${SAMPLE}/raw/${SAMPLE}_hg19_chr.vcf.gz" \
    -Oz -o "/genome/${SAMPLE}/raw/${SAMPLE}_chrM_unlifted.vcf.gz"
run_in "${BCFTOOLS_IMAGE}" \
  bcftools reheader --fai "${REF_FASTA_C}.fai" \
    -o "/genome/${SAMPLE}/raw/${SAMPLE}_chrM.vcf.gz" \
    "/genome/${SAMPLE}/raw/${SAMPLE}_chrM_unlifted.vcf.gz"

# --- Step 3: Liftover to GRCh38 ---
echo ""
echo "--- Stage 2: Liftover to GRCh38 (Picard LiftoverVcf, chrM kept in place) ---"

run_in --cpus 2 --memory 8g \
  "${PICARD_IMAGE}" \
  java -jar /usr/picard/picard.jar LiftoverVcf \
    I="/genome/${SAMPLE}/raw/${SAMPLE}_hg19_nuclear.vcf.gz" \
    O="/genome/${SAMPLE}/raw/${SAMPLE}_lifted.vcf.gz" \
    CHAIN=/genome/liftover/hg19ToHg38.over.chain.gz \
    R="${REF_FASTA_C}" \
    REJECT="/genome/${SAMPLE}/raw/${SAMPLE}_liftover_rejected.vcf.gz" \
    WARN_ON_MISSING_CONTIG=true

# The lifted records and chrM in one file, sorted in the GRCh38 contig order.
run_in --cpus 2 --memory 2g "${BCFTOOLS_IMAGE}" \
  bcftools concat "/genome/${SAMPLE}/raw/${SAMPLE}_lifted.vcf.gz" \
    "/genome/${SAMPLE}/raw/${SAMPLE}_chrM.vcf.gz" \
    -Oz -o "/genome/${SAMPLE}/raw/${SAMPLE}_merged.vcf.gz"
run_in --cpus 2 --memory 2g "${BCFTOOLS_IMAGE}" \
  bcftools sort "/genome/${SAMPLE}/raw/${SAMPLE}_merged.vcf.gz" \
    -Oz -o "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
MT_COUNT=$(run_in "${BCFTOOLS_IMAGE}" \
  bcftools view -H "/genome/${SAMPLE}/raw/${SAMPLE}_chrM.vcf.gz" | wc -l | tr -d ' ')

# --- Step 4: Index the final VCF ---
run_in "${BCFTOOLS_IMAGE}" \
  bcftools index -t -f "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"

# --- Summary ---
echo ""
echo "=== Conversion complete ==="
echo "  Output VCF: ${VCF_DIR}/${SAMPLE}.vcf.gz"

run_in "${BCFTOOLS_IMAGE}" \
  bcftools stats "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" 2>/dev/null | \
  grep "^SN" | sed 's/^SN\t0\t/  /'

REJECTED_COUNT=$(run_in "${BCFTOOLS_IMAGE}" \
  bcftools view -H "/genome/${SAMPLE}/raw/${SAMPLE}_liftover_rejected.vcf.gz" 2>/dev/null | wc -l || echo "0")
echo "  Liftover rejected: ${REJECTED_COUNT} variants"
echo "  chrM: ${MT_COUNT} variants kept at their rCRS positions (not lifted)"
echo ""
echo "You can now run pipeline steps 6, 7, 11, 25, and 27 on this VCF."
echo "See docs/chip-data-guide.md for which steps work and their limitations."
