#!/usr/bin/env bash
# 48-exomiser.sh — Rank the variants that could explain a person's symptoms (Exomiser and HPO terms)
# Usage: HPO=HP:0001250,HP:0001263 ./scripts/48-exomiser.sh <sample_name>
#
# Opt-in: runs only with HPO terms, the person's symptoms written as Human
# Phenotype Ontology ids (find them at https://hpo.jax.org; docs/48-exomiser.md
# says how to pick them). Exomiser (EXOMISER_IMAGE) scores each rare coding or
# splice variant of the VCF by how damaging it looks and each gene by how well
# its known diseases and animal models match the terms, and ranks the genes
# under each mode of inheritance. bin/exomiser_report.py checks the terms
# against the HPO release of Exomiser's data, writes the analysis (Exomiser's
# exome preset) and reads the result: the top genes, their variants, scores,
# mode of inheritance and Exomiser's ACMG classification, in words too. The
# EXOMISER_* Nextflow processes run the same code. A ranking to bring to a
# clinical geneticist, not a diagnosis.
#
# Input:  vcf/<sample>.vcf.gz (one sample), and Exomiser's data
#         (setup.sh --exomiser-data: reference/exomiser-<EXOMISER_DATA_VERSION>/
#         with <version>_hg38 and <version>_phenotype: 37 GB to download, 58 GB unpacked).
# Env:    HPO              the terms, comma or space separated (required)
#         EXOMISER_FAMILY  the folder of a step 38 family run (its family.vcf.gz
#                          and pedigree.tsv): Exomiser then reads the members'
#                          genotypes, with this sample as the affected proband
#         EXOMISER_AFFECTED, EXOMISER_UNAFFECTED
#                          members also affected, and members known not to be;
#                          every other member is unknown, which Exomiser does not
#                          filter on
#         EXOMISER_MEM     memory of Exomiser's container in GB (default 16; Java
#                          gets three quarters of it)
#         EXOMISER_DATA, EXOMISER_ASSEMBLY
#                          another Exomiser data folder (inside GENOME_DIR) and its
#                          assembly (default hg38). Only CI sets them: it runs
#                          Exomiser's cut-down hg19 test release.
# Output: exomiser/<sample>_exomiser_top.tsv   the top 10 ranks, one row per variant
#         exomiser/<sample>_exomiser.txt       the same in words
#         exomiser/<sample>_exomiser_terms.tsv the terms searched, with their names
#         exomiser/<sample>_exomiser.{variants.tsv,genes.tsv,html}  Exomiser's own files (top 50 genes)
set -euo pipefail

SAMPLE=${1:?Usage: HPO=HP:0001250,... $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image EXOMISER_IMAGE PYTHON_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
OUTDIR="${S}/exomiser"
DATA=${EXOMISER_DATA:-${GENOME_DIR}/reference/exomiser-${EXOMISER_DATA_VERSION}}
case "$DATA" in /*) ;; *) DATA="${GENOME_DIR}/${DATA}" ;; esac
ASSEMBLY=${EXOMISER_ASSEMBLY:-hg38}
MEM=${EXOMISER_MEM:-16}
P="${OUTDIR}/${SAMPLE}_exomiser"

echo "=== Step 48: symptom-driven variant ranking (Exomiser): ${SAMPLE} ==="
if [ -z "${HPO:-}" ]; then
  echo "ERROR: no HPO terms: set HPO to the symptoms as HPO ids, e.g. HPO=HP:0001250,HP:0001263 (https://hpo.jax.org)." >&2
  exit 2
fi
[[ "$ASSEMBLY" =~ ^(hg38|hg19)$ ]] || { echo "ERROR: EXOMISER_ASSEMBLY must be hg38 or hg19, got '${ASSEMBLY}'" >&2; exit 2; }
[[ "$MEM" =~ ^[0-9]+$ ]] && [ "$MEM" -ge 4 ] || { echo "ERROR: EXOMISER_MEM must be a whole number of GB, 4 or more, got '${MEM}'" >&2; exit 2; }
# Exomiser's data folder holds <version>_<assembly> and <version>_phenotype; the
# version is read from the one phenotype folder there.
shopt -s nullglob
PHENO=("${DATA}"/*_phenotype)
shopt -u nullglob
if [ "${#PHENO[@]}" -ne 1 ] || [ ! -d "${PHENO[0]%_phenotype}_${ASSEMBLY}" ]; then
  echo "ERROR: Exomiser's data is not installed in ${DATA} (one <version>_phenotype and its <version>_${ASSEMBLY} folder):" >&2
  echo "  run scripts/setup.sh --exomiser-data ${GENOME_DIR}" >&2
  exit 1
fi
VERSION=$(basename "${PHENO[0]}" _phenotype)
OBO=""
for o in "${PHENO[0]}/hp.obo" "${PHENO[0]}/phenix/hp.obo"; do [ -s "$o" ] && { OBO=$o; break; }; done

ARGS=()
if [ -n "${EXOMISER_FAMILY:-}" ]; then
  validate_sample "$EXOMISER_FAMILY"
  F="${GENOME_DIR}/${EXOMISER_FAMILY}"
  VCF="${F}/family.vcf.gz"
  for f in "$VCF" "${VCF}.tbi" "${F}/pedigree.tsv"; do
    [ -f "$f" ] || { echo "ERROR: File not found: ${f} (run scripts/38-family.sh ${EXOMISER_FAMILY} first)" >&2; exit 1; }
  done
  ARGS+=(--pedigree "$(cpath "${F}/pedigree.tsv")" --ped "$(cpath "${P}.ped")"
         --affected "${EXOMISER_AFFECTED:-}" --unaffected "${EXOMISER_UNAFFECTED:-}")
  echo "Family: ${F} (proband ${SAMPLE})"
else
  [ -z "${EXOMISER_AFFECTED:-}${EXOMISER_UNAFFECTED:-}" ] \
    || { echo "ERROR: EXOMISER_AFFECTED and EXOMISER_UNAFFECTED name family members: set EXOMISER_FAMILY too." >&2; exit 2; }
  VCF="${S}/vcf/${SAMPLE}.vcf.gz"
  [ -f "$VCF" ] || { echo "ERROR: File not found: ${VCF}" >&2; exit 1; }
fi
echo "VCF: ${VCF}"
echo "Exomiser data: ${DATA} (release ${VERSION}, ${ASSEMBLY})"
echo "HPO ontology: ${OBO:-not found in the phenotype data; the terms are checked for their form only}"

mkdir -p "$OUTDIR"
rm -f "${P}"* "${OUTDIR}/${SAMPLE}_analysis.yml"
O=$(cpath "$OUTDIR")
echo "[1/3] HPO terms and the analysis (bin/exomiser_report.py)..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 1g "$PYTHON_IMAGE" \
  python3 /pgp-bin/exomiser_report.py prepare --hpo "$HPO" ${OBO:+--obo "$(cpath "$OBO")"} --sample "$SAMPLE" \
    --vcf "$(cpath "$VCF")" --assembly "$ASSEMBLY" --out-dir "$O" --analysis "${O}/${SAMPLE}_analysis.yml" \
    --terms "$(cpath "${P}_terms.tsv")" ${ARGS[@]+"${ARGS[@]}"}
echo "[2/3] Exomiser..."
# The image's entrypoint is bash; its exomiser wrapper starts Java with JAVA_OPTS.
# Exomiser reads its data folder and versions from these variables (its
# installation docs) and needs no network.
run_in --cpus "$THREADS" --memory "${MEM}g" --entrypoint /usr/local/bin/exomiser \
  -e JAVA_OPTS="-Xmx$((MEM * 3 / 4))g" -e EXOMISER_DATA_DIRECTORY="$(cpath "$DATA")" \
  -e "EXOMISER_$(echo "$ASSEMBLY" | tr "[:lower:]" "[:upper:]")_DATA_VERSION=${VERSION}" -e EXOMISER_PHENOTYPE_DATA_VERSION="$VERSION" \
  "$EXOMISER_IMAGE" analyse --analysis "${O}/${SAMPLE}_analysis.yml"
[ -f "${P}.variants.tsv" ] || { echo "ERROR: Exomiser wrote no ${P}.variants.tsv" >&2; exit 1; }
echo "[3/3] The top genes in words..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 1g "$PYTHON_IMAGE" \
  python3 /pgp-bin/exomiser_report.py report --sample "$SAMPLE" --variants "${O}/${SAMPLE}_exomiser.variants.tsv" \
    --terms "$(cpath "${P}_terms.tsv")" --out "$(cpath "${P}_top.tsv")" --text "$(cpath "${P}.txt")"

echo ""
echo "Top genes: ${P}_top.tsv"
echo "In words: ${P}.txt"
echo "Exomiser's report: ${P}.html"
