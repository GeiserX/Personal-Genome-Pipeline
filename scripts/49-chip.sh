#!/usr/bin/env bash
# 49-chip.sh — Clonal haematopoiesis (CHIP): acquired blood-cell mutations in DNMT3A, TET2, ASXL1, JAK2 and others
# Usage: ./scripts/49-chip.sh <sample_name>
#
# Calls low-VAF changes with GATK Mutect2 in tumour-only mode on the coding
# exons (10 bases wider on each side) of the genes in assets/chip_genes.tsv,
# each on its MANE Select transcript from Ensembl's gene annotation (setup.sh
# --chip-data). gnomAD's AF-only germline resource and the panel of normals of
# step 29 are used when they are in GENOME_DIR/somatic/. bcftools csq names
# each change on that transcript, and bin/chip_report.py keeps the changes in
# the table's CHIP classes, leaves out the ones whose VAF looks inherited (35%
# and up) and states the detection limit at the sample's depth over the same
# exons. The same commands run in the Nextflow CHIP_* processes.
#
# Not a diagnosis: CHIP is common with age and most people with it never get
# a blood cancer. A negative result does not rule out a clone smaller than the
# detection limit. docs/49-chip.md says how to read it.
#
# Input:  aligned/<sample>_sorted.bam (ALIGN_DIR), the reference with .fai and .dict,
#         reference/ensembl/Homo_sapiens.GRCh38.<ENSEMBL_GFF_RELEASE>.chr.gff3.gz
# Output: chip/<sample>_chip.tsv (one row per protein change, with its status),
#         chip/<sample>_chip_depth.tsv (depth and detection limit per gene),
#         chip/<sample>_chip_filtered.vcf.gz (Mutect2's calls)
# CHIP_GENES=FILE uses another gene table with the columns of assets/chip_genes.tsv.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}   # common.sh defaults to 8
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image GATK_IMAGE BCFTOOLS_IMAGE PYTHON_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
ALIGN_DIR=${ALIGN_DIR:-aligned}
BAM="${S}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
GFF="${GENOME_DIR}/reference/ensembl/Homo_sapiens.GRCh38.${ENSEMBL_GFF_RELEASE}.chr.gff3.gz"
GENES=${CHIP_GENES:-${PGP_ROOT}/assets/chip_genes.tsv}
OUTDIR="${S}/chip"
P="${OUTDIR}/${SAMPLE}_chip"
GNOMAD_VCF="${GENOME_DIR}/somatic/af-only-gnomad.hg38.vcf.gz"
PON_VCF="${GENOME_DIR}/somatic/1000g_pon.hg38.vcf.gz"

echo "============================================"
echo "  Step 49: Clonal haematopoiesis (CHIP)"
echo "  Sample: ${SAMPLE}"
echo "  Input:  ${BAM}"
echo "  Genes:  ${GENES}"
echo "  Output: ${P}.tsv"
echo "============================================"
for f in "$BAM" "${BAM}.bai" "$REF_FASTA" "${REF_FASTA}.fai" "$GENES"; do
  [ -f "$f" ] || { echo "ERROR: File not found: ${f}" >&2; exit 1; }
done
[ -f "$REF_DICT" ] || { echo "ERROR: Sequence dictionary not found: ${REF_DICT}. Mutect2 needs it; setup.sh creates it." >&2; exit 1; }
[ -f "$GFF" ] || { echo "ERROR: File not found: ${GFF}. Install Ensembl's gene annotation with: scripts/setup.sh --chip-data \"\${GENOME_DIR}\"" >&2; exit 1; }

mkdir -p "$OUTDIR"
rm -f "${P}.tsv" "${P}_depth.tsv"
# Containers see GENOME_DIR only: the table is copied beside the outputs.
cp "$GENES" "${OUTDIR}/chip_genes.tsv"
PY=(run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 2g "$PYTHON_IMAGE" python3 /pgp-bin/chip_report.py)

echo "[1/4] The CHIP genes' MANE Select transcripts and coding exons (Ensembl ${ENSEMBL_GFF_RELEASE})..."
"${PY[@]}" prepare --gff "$(cpath "$GFF")" --genes "$(cpath "${OUTDIR}/chip_genes.tsv")" --fai "${REF_FASTA_C}.fai" \
  --out-gff "$(cpath "${P}.gff3")" --out-bed "$(cpath "${P}.bed")"
[ -s "${P}.bed" ] || { echo "ERROR: none of the CHIP genes is on a contig of ${REF_FASTA}." >&2; exit 1; }

echo "[2/4] Mutect2, tumour-only, on $(wc -l < "${P}.bed" | tr -d ' ') coding exons..."
M2=(gatk Mutect2 -R "${REF_FASTA_C}" -I "$(cpath "$BAM")" -L "$(cpath "${P}.bed")"
  -O "$(cpath "${P}_unfiltered.vcf.gz")" --f1r2-tar-gz "$(cpath "${P}_f1r2.tar.gz")"
  --native-pair-hmm-threads "$THREADS" --max-mnp-distance 0)
if [ -f "$GNOMAD_VCF" ] && [ -f "${GNOMAD_VCF}.tbi" ]; then
  M2+=(--germline-resource "$(cpath "$GNOMAD_VCF")"); echo "  germline resource: ${GNOMAD_VCF}"
else
  echo "  no gnomAD germline resource (${GNOMAD_VCF}): inherited variants are left out by their VAF only (docs/29-mutect2-somatic.md)."
fi
if [ -f "$PON_VCF" ] && [ -f "${PON_VCF}.tbi" ]; then
  M2+=(-pon "$(cpath "$PON_VCF")"); echo "  panel of normals: ${PON_VCF}"
fi
run_in --cpus "$THREADS" --memory 8g "$GATK_IMAGE" "${M2[@]}"
run_in --cpus 2 --memory 4g "$GATK_IMAGE" gatk LearnReadOrientationModel -I "$(cpath "${P}_f1r2.tar.gz")" \
  -O "$(cpath "${P}_orientation.tar.gz")"
run_in --cpus 2 --memory 4g "$GATK_IMAGE" gatk FilterMutectCalls -R "${REF_FASTA_C}" -V "$(cpath "${P}_unfiltered.vcf.gz")" \
  --ob-priors "$(cpath "${P}_orientation.tar.gz")" -O "$(cpath "${P}_filtered.vcf.gz")"

echo "[3/4] Naming each change on its transcript (bcftools csq) and counting the depth (bcftools mpileup)..."
# norm splits multi-allelic records, so AD and AF hold one ALT each; csq
# --local-csq reads each record alone (a clone's changes need not share a
# haplotype). mpileup counts what Mutect2 reads by default: mapping quality
# 20 and up, base quality 10 and up, duplicates and secondary reads left out.
# shellcheck disable=SC2016  # $1 to $7 belong to the inner bash
run_in --cpus 1 --memory 2g "$BCFTOOLS_IMAGE" bash -euo pipefail -c '
    bcftools norm -m -any -f "$1" -Ou "$2" \
      | bcftools csq -v 0 -f "$1" -g "$3" -s - --local-csq -Ou \
      | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/BCSQ\t[%AD]\t[%AF]\t[%DP]\n" > "$5"
    bcftools mpileup -f "$1" -R "$4" -q 20 -Q 10 -B -d 100000 -a FORMAT/DP -Ou "$6" \
      | bcftools query -f "%CHROM\t%POS\t[%DP]\n" > "$7"' \
  _ "${REF_FASTA_C}" "$(cpath "${P}_filtered.vcf.gz")" "$(cpath "${P}.gff3")" "$(cpath "${P}.bed")" \
  "$(cpath "${P}_calls.tsv")" "$(cpath "$BAM")" "$(cpath "${P}_pileup.tsv")"

echo "[4/4] CHIP classes, germline exclusion and detection limit..."
"${PY[@]}" report --genes "$(cpath "${OUTDIR}/chip_genes.tsv")" --gff "$(cpath "${P}.gff3")" --bed "$(cpath "${P}.bed")" \
  --calls "$(cpath "${P}_calls.tsv")" --depth "$(cpath "${P}_pileup.tsv")" --sample "$SAMPLE" \
  --out "$(cpath "${P}.tsv")" --depth-out "$(cpath "${P}_depth.tsv")"
rm -f "${P}_pileup.tsv" "${P}_f1r2.tar.gz" "${P}_orientation.tar.gz"

echo ""
echo "CHIP calls (${P}.tsv):"
awk -F'\t' 'NR > 1 && $13 == "chip" { printf "  %-8s %-10s VAF %-6s %s/%s reads  %s\n", $1, $7, $9, $10, $11, $15; n++ }
  END { if (!n) print "  none" }' "${P}.tsv"
echo "Not a diagnosis. docs/49-chip.md says what a call means and what to ask a haematologist."
