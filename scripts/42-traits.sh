#!/usr/bin/env bash
# 42-traits.sh — Traits panel: lactase persistence, caffeine, alcohol flush, bitter taste, earwax, APOE...
# Usage: ./scripts/42-traits.sh <sample_name>
#
# Reads the sample's genotype at each position of assets/traits.tsv, a
# hand-curated table of well-replicated variants with the reading of each
# genotype and the primary papers behind it (no SNPedia text), and writes the
# reading per row. assets/genotype_rules.tsv reads several sites together:
# APOE's e2/e3/e4 pair and TAS2R38's taster haplotypes (PAV, AVI). The APOE
# line comes after the two APOE rows of the traits table; every rule has a
# line in the rules table.
#
# Reads the gVCF of step 3 (vcf/<sample>.g.vcf.gz), never the variant-only
# VCF: a position where the sample matches the reference is then a real
# ref/ref from a reference block, and a position with no reads is "not
# covered" instead of being taken for the reference. A rule that needs phase
# (TAS2R38) with two or more heterozygous sites is phased from the reads with
# WhatsHap, when the BAM of step 2 is there (aligned/<sample>_sorted.bam);
# sites no read or read pair joins stay unphased and the rule is not read.
# The same commands run in the Nextflow TRAITS_SITES, TRAITS_PHASE and TRAITS
# processes.
#
# Output: traits/<sample>_traits.tsv (rsid, gene, chrom, pos, trait, genotype,
# status, reading, ancestry, note, pmids). status: called, not_covered, or
# other_allele (an allele the table does not describe).
# traits/<sample>_genotype_rules.tsv: one line per rule (docs/42-traits.md).
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image BCFTOOLS_IMAGE PYTHON_IMAGE WHATSHAP_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
GVCF="${S}/vcf/${SAMPLE}.g.vcf.gz"
BAM="${S}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam"
TABLE="${PGP_ROOT}/assets/traits.tsv"
RULES="${PGP_ROOT}/assets/genotype_rules.tsv"
OUTDIR="${S}/traits"
OUT="${OUTDIR}/${SAMPLE}_traits.tsv"
RULES_OUT="${OUTDIR}/${SAMPLE}_genotype_rules.tsv"
PHASE_IN="${OUTDIR}/.${SAMPLE}_traits_phase_sites.vcf"
PHASED="${OUTDIR}/.${SAMPLE}_traits_phased.vcf"
SITES="${OUTDIR}/.${SAMPLE}_traits_sites.tsv"
GENOTYPES="${OUTDIR}/${SAMPLE}_traits_genotypes.tsv"

echo "============================================"
echo "  Step 42: Traits panel"
echo "  Sample: ${SAMPLE}"
echo "  Input:  ${GVCF}"
echo "  Tables: assets/traits.tsv, assets/genotype_rules.tsv"
echo "  Output: ${OUT}, ${RULES_OUT}"
echo "============================================"
for f in "$GVCF" "${GVCF}.tbi" "$REF_FASTA" "${REF_FASTA}.fai"; do
  [ -f "$f" ] || { echo "ERROR: File not found: ${f}. Step 42 reads the gVCF of step 3 (a variant-only VCF cannot tell a reference match from a position without reads)." >&2; exit 1; }
done
for f in "$TABLE" "$RULES"; do
  [ -f "$f" ] || { echo "ERROR: File not found: ${f}. The table ships with the repository; restore it with git checkout -- assets/$(basename "$f")." >&2; exit 1; }
done

mkdir -p "$OUTDIR"
rm -f "$OUT" "$RULES_OUT" "$PHASE_IN" "$PHASED"
# The positions of the traits table (comment lines and the header left out)
# and of the rule sites, sorted.
{ grep -v '^#' "$TABLE" | awk -F'\t' 'NR > 1 { print $1 "\t" $2 }'
  awk -F'\t' '$1 == "site" { print $4 "\t" $5 }' "$RULES"; } | LC_ALL=C sort -k1,1 -k2,2n -u > "$SITES"

# gvcf2vcf expands each reference block over a position into a record with
# the reference base (with the block's GQ), -T keeps the table's positions,
# --trim-alt-alleles drops <*>, and a no-call is dropped, so a position
# without reads has no line.
echo "[1/3] Genotyping the table's positions from the gVCF..."
# shellcheck disable=SC2016  # $1 to $3 belong to the inner bash
run_in --cpus 1 --memory 2g "$BCFTOOLS_IMAGE" \
  bash -euo pipefail -c '
    n=$(bcftools query -l "$3" | wc -l | tr -d " ")
    [ "$n" -eq 1 ] || { echo "ERROR: $3 holds $n samples; step 42 reads a one-sample gVCF (keep one with bcftools view -s <sample>)." >&2; exit 1; }
    bcftools convert --gvcf2vcf -f "$1" -R "$2" -Ou "$3" \
      | bcftools view -T "$2" --trim-alt-alleles -i "GT!=\"mis\"" -Ou \
      | bcftools query -u -f "%CHROM\t%POS\t%REF\t%ALT\t[%TGT]\t[%GQ]\n"' \
  _ "${REF_FASTA_C}" "$(cpath "$SITES")" "$(cpath "$GVCF")" > "${GENOTYPES}.tmp" \
  || { rm -f "${GENOTYPES}.tmp" "$SITES"; echo "ERROR: step 42 could not genotype ${GVCF} at the table's positions (the line above says why)." >&2; exit 1; }
mv -f "${GENOTYPES}.tmp" "$GENOTYPES"
rm -f "$SITES"

# WhatsHap reads the BAM at the heterozygous sites of the rules that need
# phase (none when no such rule has two). A read, or a read pair, that shows
# alleles at two of them puts them in one phase block; the rest stay unphased.
echo "[2/3] Phasing the rule sites that need it from the reads (WhatsHap)..."
if [ -f "$BAM" ] && [ -f "${BAM}.bai" ]; then
  # shellcheck disable=SC2016  # $1 to $6 belong to the inner bash
  run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" -v "${PGP_ROOT}/assets:/pgp-assets:ro" --cpus 1 --memory 2g "$WHATSHAP_IMAGE" \
    bash -euo pipefail -c '
      python3 /pgp-bin/traits_report.py --rules /pgp-assets/genotype_rules.tsv --genotypes "$1" --sample "$2" --phase-sites "$3"
      if grep -qv "^#" "$3"; then
        whatshap phase --reference "$4" --ignore-read-groups -o "$5" "$3" "$6"
      else
        : > "$5"
      fi' \
    _ "$(cpath "$GENOTYPES")" "$SAMPLE" "$(cpath "$PHASE_IN")" "${REF_FASTA_C}" "$(cpath "$PHASED")" "$(cpath "$BAM")" \
    && PHASE_OK=1 || PHASE_OK=0
  rm -f "$PHASE_IN"
  if [ "$PHASE_OK" = 1 ]; then
    PHASE_ARGS=(--phased "$(cpath "$PHASED")")
  else
    # Same as Nextflow, where TRAITS_PHASE may fail and TRAITS still runs:
    # the traits table is still written, and the rules read unphased.
    rm -f "$PHASED"
    echo "WARNING: step 42 could not phase the rule sites from ${BAM} (the lines above say why); the rules read unphased." >&2
    PHASE_ARGS=()
  fi
else
  echo "  No BAM (${BAM} with its .bai): a rule that needs phase reads only when at most one of its sites is heterozygous."
  PHASE_ARGS=()
fi

echo "[3/3] Reading each genotype and each rule against the tables..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" -v "${PGP_ROOT}/assets:/pgp-assets:ro" --cpus 1 --memory 1g "$PYTHON_IMAGE" \
  python3 /pgp-bin/traits_report.py --table /pgp-assets/traits.tsv --rules /pgp-assets/genotype_rules.tsv \
    --genotypes "$(cpath "$GENOTYPES")" ${PHASE_ARGS[@]+"${PHASE_ARGS[@]}"} --sample "$SAMPLE" \
    --out "$(cpath "$OUT")" --rules-out "$(cpath "$RULES_OUT")"
rm -f "$PHASED"

echo ""
echo "A reading is what the cited studies found on average, not a diagnosis. APOE is a risk factor"
echo "for Alzheimer's disease: docs/42-traits.md says what it means before you read it."
echo ""
echo "Traits (${OUT}):"
awk -F'\t' 'NR > 1 { printf "  %-45s %-12s %s\n", substr($5, 1, 45), $6, $8 }' "$OUT"
echo ""
echo "Genotype rules (${RULES_OUT}):"
awk -F'\t' 'NR > 1 && $1 != "APOE" { printf "  %-45s %-12s %s\n", substr($3, 1, 45), ($11 == "called" ? $10 : $11), $12 }' "$RULES_OUT"
