#!/usr/bin/env bash
# 42-traits.sh — Traits panel: lactase persistence, caffeine, alcohol flush, bitter taste, earwax, APOE...
# Usage: ./scripts/42-traits.sh <sample_name>
#
# Reads the sample's genotype at each position of assets/traits.tsv, a
# hand-curated table of well-replicated variants with the reading of each
# genotype and the primary papers behind it (no SNPedia text), and writes the
# reading per row. After the two APOE rows comes an APOE line with the
# e2/e3/e4 pair.
#
# Reads the gVCF of step 3 (vcf/<sample>.g.vcf.gz), never the variant-only
# VCF: a position where the sample matches the reference is then a real
# ref/ref from a reference block, and a position with no reads is "not
# covered" instead of being taken for the reference. The same commands run in
# the Nextflow TRAITS_SITES and TRAITS processes.
#
# Output: traits/<sample>_traits.tsv (rsid, gene, chrom, pos, trait, genotype,
# status, reading, ancestry, note, pmids). status: called, not_covered, or
# other_allele (an allele the table does not describe).
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image BCFTOOLS_IMAGE PYTHON_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
GVCF="${S}/vcf/${SAMPLE}.g.vcf.gz"
TABLE="${PGP_ROOT}/assets/traits.tsv"
OUTDIR="${S}/traits"
OUT="${OUTDIR}/${SAMPLE}_traits.tsv"
SITES="${OUTDIR}/.${SAMPLE}_traits_sites.tsv"
GENOTYPES="${OUTDIR}/${SAMPLE}_traits_genotypes.tsv"

echo "============================================"
echo "  Step 42: Traits panel"
echo "  Sample: ${SAMPLE}"
echo "  Input:  ${GVCF}"
echo "  Table:  assets/traits.tsv"
echo "  Output: ${OUT}"
echo "============================================"
for f in "$GVCF" "${GVCF}.tbi" "$REF_FASTA" "${REF_FASTA}.fai"; do
  [ -f "$f" ] || { echo "ERROR: File not found: ${f}. Step 42 reads the gVCF of step 3 (a variant-only VCF cannot tell a reference match from a position without reads)." >&2; exit 1; }
done
[ -f "$TABLE" ] || { echo "ERROR: File not found: ${TABLE}. The traits table ships with the repository; restore it with git checkout -- assets/traits.tsv." >&2; exit 1; }

mkdir -p "$OUTDIR"
rm -f "$OUT"
# The table's positions (comment lines and the header left out), sorted.
grep -v '^#' "$TABLE" | awk -F'\t' 'NR > 1 { print $1 "\t" $2 }' | LC_ALL=C sort -k1,1 -k2,2n -u > "$SITES"

# gvcf2vcf expands each reference block over a position into a record with
# the reference base, -T keeps the table's positions, --trim-alt-alleles drops
# <*>, and a no-call is dropped, so a position without reads has no line.
echo "[1/2] Genotyping the table's positions from the gVCF..."
# shellcheck disable=SC2016  # $1 to $3 belong to the inner bash
run_in --cpus 1 --memory 2g "$BCFTOOLS_IMAGE" \
  bash -euo pipefail -c '
    bcftools convert --gvcf2vcf -f "$1" -R "$2" -Ou "$3" \
      | bcftools view -T "$2" --trim-alt-alleles -i "GT!=\"mis\"" -Ou \
      | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t[%TGT]\n"' \
  _ "${REF_FASTA_C}" "$(cpath "$SITES")" "$(cpath "$GVCF")" > "${GENOTYPES}.tmp" \
  || { rm -f "${GENOTYPES}.tmp" "$SITES"; echo "ERROR: bcftools could not read ${GVCF} at the table's positions." >&2; exit 1; }
mv -f "${GENOTYPES}.tmp" "$GENOTYPES"
rm -f "$SITES"

echo "[2/2] Reading each genotype against the table..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" -v "${PGP_ROOT}/assets:/pgp-assets:ro" --cpus 1 --memory 1g "$PYTHON_IMAGE" \
  python3 /pgp-bin/traits_report.py --table /pgp-assets/traits.tsv \
    --genotypes "$(cpath "$GENOTYPES")" --sample "$SAMPLE" --out "$(cpath "$OUT")"

echo ""
echo "Traits (${OUT}):"
awk -F'\t' 'NR > 1 { printf "  %-45s %-12s %s\n", substr($5, 1, 45), $6, $7 }' "$OUT"
echo ""
echo "A reading is what the cited studies found on average, not a diagnosis. APOE is a risk factor"
echo "for Alzheimer's disease: docs/42-traits.md says what it means before you read it."
