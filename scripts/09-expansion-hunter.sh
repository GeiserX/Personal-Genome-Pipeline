#!/usr/bin/env bash
# ExpansionHunter v5 — Screen for pathogenic short tandem repeat (STR) expansions
# Detects: Huntington's, Fragile X, Friedreich's ataxia, ALS/FTD, SCAs, myotonic dystrophy, etc.
# Input: sorted BAM + GRCh38 reference FASTA
# Output: JSON + VCF with repeat genotypes in $GENOME_DIR/<sample>/expansion_hunter/
#
# Upgraded from v2.5.5 to v5.0.0:
#   - --bam → --reads, --ref-fasta → --reference, --repeat-specs → --variant-catalog
#   - --vcf/--json/--log → --output-prefix (auto-generates .vcf, .json)
#   - Multithreading support (--threads)
#   - Bundled GRCh38 variant catalog (31 pathogenic loci) inside the container
# Catalog: STRchive's 82 disease loci (reference/strchive-<version>/catalog.json,
# which setup.sh builds with bin/strchive_convert.py), or the bundled 31 when
# that file is missing. EH_CATALOG=bundled picks the bundled one;
# EH_CATALOG=FILE another catalog (a JSON file under GENOME_DIR).
# ExpansionHunter stops when a catalog locus is on a contig the reference
# lacks, so a reduced reference needs a reduced catalog. The catalog that ran
# is copied beside the VCF (<sample>_eh_catalog.json): step 9b grades with it
# when it carries Stranger's thresholds, so the two steps read one catalog.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> <male|female>}
SEX=${2:?Usage: $0 <sample_name> <male|female>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}   # common.sh defaults to 8
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
ALIGN_DIR=${ALIGN_DIR:-aligned}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
OUTPUT_DIR="${SAMPLE_DIR}/expansion_hunter"

echo "=== ExpansionHunter v5: ${SAMPLE} (${SEX}) ==="

for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

# ExpansionHunter v5.0.0 via biocontainer
# The variant catalog (31 pathogenic GRCh38 loci) is bundled inside the container
# at /usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json
CATALOG_C=/usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json
STRCHIVE_CATALOG="${GENOME_DIR}/reference/strchive-${STRCHIVE_VERSION}/catalog.json"
CATALOG=""
case "${EH_CATALOG:-}" in
  bundled) echo "Variant catalog: the 31 loci inside the image (EH_CATALOG=bundled)" ;;
  "")
    if [ -s "$STRCHIVE_CATALOG" ]; then
      CATALOG=$STRCHIVE_CATALOG
      echo "Variant catalog: STRchive ${STRCHIVE_VERSION}, ${CATALOG}"
    else
      echo "Variant catalog: the 31 loci inside the image (${STRCHIVE_CATALOG#"${GENOME_DIR}"/} not found; setup.sh --strchive-data installs STRchive's 82)"
    fi ;;
  *)
    [ -f "$EH_CATALOG" ] || { echo "ERROR: EH_CATALOG not found: ${EH_CATALOG}" >&2; exit 1; }
    CATALOG=$EH_CATALOG
    echo "Variant catalog: ${EH_CATALOG} (EH_CATALOG)" ;;
esac
if [ -n "$CATALOG" ]; then
  CATALOG_C=$(cpath "$CATALOG") || exit 2
  # bin/strchive_convert.py leaves Stranger's thresholds off a locus whose
  # repeat short reads cannot size, or that STRchive gives no pathogenic size
  # (NoThresholds, one field per line as it writes them): its count is not a call.
  awk -F'"' '$2 == "LocusId" {id = $4} $2 == "NoThresholds" {print "  Not graded (no call from this count): " id ": " $4}' "$CATALOG"
fi
run_in \
  --cpus "${THREADS}" --memory 4g \
  "${EXPANSIONHUNTER_IMAGE}" \
  ExpansionHunter \
    --reads "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
    --reference "${REF_FASTA_C}" \
    --variant-catalog "$CATALOG_C" \
    --output-prefix "/genome/${SAMPLE}/expansion_hunter/${SAMPLE}_eh" \
    --threads "${THREADS}" \
    --sex "$SEX"

# The catalog that ran, for step 9b; none for the bundled one, so 9b falls
# back to Stranger's bundled GRCh38 catalog as before.
if [ -n "$CATALOG" ]; then
  cp "$CATALOG" "${OUTPUT_DIR}/${SAMPLE}_eh_catalog.json.part"
  mv -f "${OUTPUT_DIR}/${SAMPLE}_eh_catalog.json.part" "${OUTPUT_DIR}/${SAMPLE}_eh_catalog.json"
else
  rm -f "${OUTPUT_DIR}/${SAMPLE}_eh_catalog.json"
fi

echo "=== ExpansionHunter complete ==="
echo "VCF:  ${OUTPUT_DIR}/${SAMPLE}_eh.vcf"
echo "JSON: ${OUTPUT_DIR}/${SAMPLE}_eh.json"
echo ""
echo "Key disease thresholds (alleles >= threshold = pathogenic):"
echo "  HTT(Huntington): 36  FMR1(FragileX): 55  FXN(Friedreich): 66"
echo "  C9ORF72(ALS/FTD): 30  DMPK(MyotonicDystrophy): 50  ATXN1(SCA1): 39"
