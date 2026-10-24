#!/usr/bin/env bash
# Step 38: family analysis of two or more processed samples (opt-in)
# Usage: ./scripts/38-family.sh <family_name> [pedigree.tsv]
# Input: a pedigree (assets/pedigree.example.tsv: sample, sex, father, mother,
#        partner, tab separated; default $GENOME_DIR/<family_name>/pedigree.tsv)
#        naming two or more samples that have a DeepVariant gVCF in
#        $GENOME_DIR/<sample>/vcf/<sample>.g.vcf.gz (step 03), and somalier's
#        sites (setup.sh --sample-qc-data)
# Output: $GENOME_DIR/<family_name>/, a folder of its own (the name may not be
#   a sample's):
#   pedigree.tsv           the pedigree the step read
#   family.vcf.gz          the members joint-called by GLnexus
#   somalier/              somalier relate over the members' gVCFs
#   relationships.tsv      each pair's declared and measured relatedness, verdict
#   family.phased.vcf.gz   WhatsHap's phasing by inheritance (complete trios only)
#   <child>_inheritance.tsv.gz  the allele from each parent, site by site
#   de_novo.tsv            de novo candidates (children with both parents sequenced)
#   compound_het.tsv       pairs of a child's clinical-filter candidates (step 23)
#                          in one gene, in trans or in cis
#   family_report.md       all of it in words
#
# The step stops (exit 1) when a pair's measured relatedness disagrees with
# the pedigree, before any phasing: a sample swap or a wrong pedigree line.
# relationships.tsv says which pair and why. The rules are bin/family_report.py's,
# the code the Nextflow FAMILY workflow (--pedigree) runs too.
#
# Env: SOMALIER_SITES  sites VCF (default reference/somalier/sites.hg38.vcf.gz),
#                      inside GENOME_DIR
#      THREADS         CPUs of GLnexus (default 8)
#      GLNEXUS_MEM     memory of GLnexus's container in GB (default 16; GLnexus's
#                      own --mem-gbytes budget is half of it)
set -euo pipefail

FAMILY=${1:?Usage: $0 <family_name> [pedigree.tsv]}
PEDIGREE=${2:-}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# No one sample's folder is writable here: run_in mounts only the family folder.
unset SAMPLE
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
# The family name names its folder, so it follows the sample-name rule.
validate_sample "$FAMILY"
require_image GLNEXUS_IMAGE BCFTOOLS_IMAGE SOMALIER_IMAGE WHATSHAP_IMAGE PYTHON_IMAGE
GLNEXUS_MEM=${GLNEXUS_MEM:-16}
OUT="${GENOME_DIR}/${FAMILY}"
OUT_C="/genome/${FAMILY}"
PEDIGREE=${PEDIGREE:-${OUT}/pedigree.tsv}
SITES=${SOMALIER_SITES:-${GENOME_DIR}/reference/somalier/sites.hg38.vcf.gz}
case "$SITES" in /*) ;; *) SITES="${GENOME_DIR}/${SITES}" ;; esac

echo "=== Step 38: family analysis ==="
echo "Pedigree: ${PEDIGREE}"
echo "Somalier sites: ${SITES}"
echo "Output: ${OUT}/"

for f in "$PEDIGREE" "$REF_FASTA" "${REF_FASTA}.fai" "$SITES"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    [ "$f" = "$SITES" ] && echo "  Install it with: ./scripts/setup.sh --sample-qc-data ${GENOME_DIR}" >&2
    exit 1
  fi
done
SITES_C=$(cpath "$SITES") || exit 1

# The samples the pedigree names (its sample column; the first non-blank
# line is the header), and which of them have a gVCF here. Read as
# family_report.py reads it: CRLF line ends, whitespace-only lines and spaces
# around a column name are allowed. family_report.py prepare checks
# everything else.
mapfile -t NAMED < <(awk -F'\t' '
  { sub(/\r+$/, "") }
  /^[ \t]*$/ { next }
  !h { sub(/^#+/, ""); for (i = 1; i <= NF; i++) { f = tolower($i); gsub(/^[ \t]+|[ \t]+$/, "", f); if (f == "sample") c = i }; h = 1; next }
  /^#/ || !c { next }
  { s = $c; gsub(/^[ \t]+|[ \t]+$/, "", s); if (s != "") print s }' "$PEDIGREE")
MEMBERS=()
for s in ${NAMED[@]+"${NAMED[@]}"}; do
  validate_sample "$s"
  g="${GENOME_DIR}/${s}/vcf/${s}.g.vcf.gz"
  if [ -f "$g" ] && [ -f "${g}.tbi" ]; then
    MEMBERS+=("$s")
  else
    echo "  ${s}: no gVCF with its index at ${g#"${GENOME_DIR}"/}; left out (run step 03 for it)"
  fi
done
for s in ${NAMED[@]+"${NAMED[@]}"}; do
  if [ "$s" = "$FAMILY" ]; then
    echo "ERROR: the family name '${FAMILY}' is also a sample of the pedigree; its folder would be the sample's." >&2
    echo "  Give the family another name." >&2
    exit 2
  fi
done
if [ -d "${OUT}/vcf" ] || [ -d "${OUT}/aligned" ]; then
  echo "ERROR: ${OUT} holds a sample's vcf/ or aligned/ folder; give the family a name of its own." >&2
  exit 2
fi
if [ "${#MEMBERS[@]}" -lt 2 ]; then
  echo "ERROR: the pedigree names ${#MEMBERS[@]} sample(s) with a gVCF (${MEMBERS[*]:-none}); the family analysis" >&2
  echo "  needs two or more. Run steps 02 and 03 for the others first." >&2
  exit 1
fi
echo "Members: ${MEMBERS[*]}"

mkdir -p "$OUT"
[ "$PEDIGREE" -ef "${OUT}/pedigree.tsv" ] || cp -f "$PEDIGREE" "${OUT}/pedigree.tsv"
PY=(run_in --rw "$OUT" --cpus 2 --memory 4g -v "${PGP_ROOT}/bin:/pgp-bin:ro" -v "${PGP_ROOT}/assets/par_grch38.bed:/pgp/par_grch38.bed:ro"
    "$PYTHON_IMAGE" python3 /pgp-bin/family_report.py)

# --- 1. The pedigree, and each gVCF's sample name --------------------------------------
GVCF_ARGS=()
GVCFS_C=()
for s in "${MEMBERS[@]}"; do
  GVCF_ARGS+=(--gvcf "${s}=/genome/${s}/vcf/${s}.g.vcf.gz")
  GVCFS_C+=("/genome/${s}/vcf/${s}.g.vcf.gz")
done
"${PY[@]}" prepare --pedigree "${OUT_C}/pedigree.tsv" "${GVCF_ARGS[@]}" --out-dir "$OUT_C"

# --- 2. Joint call (GLnexus) --------------------------------------------------------------
echo ""
echo "--- joint call (GLnexus, DeepVariantWGS)"
rm -rf "${OUT}/glnexus_db"
# --mem-gbytes sizes only GLnexus's RocksDB (block cache and bulk-load write
# buffers); the alleles it discovers and the unifier's sites are held outside
# it, so a budget near the limit leaves no room for them. GLnexus gets half
# of GLNEXUS_MEM (the Nextflow GLNEXUS rule too).
run_in --rw "$OUT" --cpus "$THREADS" --memory "${GLNEXUS_MEM}g" "$GLNEXUS_IMAGE" \
  glnexus_cli --config DeepVariantWGS --dir "${OUT_C}/glnexus_db" --threads "$THREADS" \
    --mem-gbytes "$((GLNEXUS_MEM / 2 > 1 ? GLNEXUS_MEM / 2 : 1))" "${GVCFS_C[@]}" > "${OUT}/family.bcf.tmp"
rm -rf "${OUT}/glnexus_db"
mv -f "${OUT}/family.bcf.tmp" "${OUT}/family.bcf"
run_in --rw "$OUT" --cpus 2 --memory 4g "$BCFTOOLS_IMAGE" bcftools view -Oz -o "${OUT_C}/family.vcf.gz" "${OUT_C}/family.bcf"
run_in --rw "$OUT" --cpus 2 --memory 4g "$BCFTOOLS_IMAGE" bcftools index -f -t "${OUT_C}/family.vcf.gz"
rm -f "${OUT}/family.bcf"
have_output "${OUT}/family.vcf.gz" || { echo "ERROR: no ${OUT}/family.vcf.gz" >&2; exit 1; }

# --- 3. Relationships (somalier) ----------------------------------------------------------
echo ""
echo "--- relationships (somalier)"
rm -rf "${OUT}/somalier"
mkdir -p "${OUT}/somalier/extract"
EXTRACTED=()
for s in "${MEMBERS[@]}"; do
  run_in --rw "$OUT" --cpus 2 --memory 4g "$SOMALIER_IMAGE" \
    somalier extract -d "${OUT_C}/somalier/extract" --sites "$SITES_C" -f "$REF_FASTA_C" "/genome/${s}/vcf/${s}.g.vcf.gz"
  EXTRACTED+=("${OUT_C}/somalier/extract/${s}.somalier")
done
run_in --rw "$OUT" --cpus 1 --memory 2g "$SOMALIER_IMAGE" \
  somalier relate --ped "${OUT_C}/family.ped" -o "${OUT_C}/somalier/family" "${EXTRACTED[@]}"
rm -rf "${OUT}/somalier/extract"
"${PY[@]}" check --members "${OUT_C}/members.tsv" --pairs "${OUT_C}/somalier/family.pairs.tsv" \
  --out "${OUT_C}/relationships.tsv"

# --- 4. Phasing by inheritance (WhatsHap), complete trios only -------------------------------
echo ""
echo "--- phasing by inheritance (WhatsHap)"
rm -f "${OUT}/family.phased.vcf.gz" "${OUT}/family.phased.vcf.gz.tbi"
PHASED=()
if [ -s "${OUT}/trios.ped" ]; then
  run_in --rw "$OUT" --cpus 4 --memory 8g "$WHATSHAP_IMAGE" \
    whatshap phase --ped "${OUT_C}/trios.ped" --use-ped-samples -o "${OUT_C}/family.phased.vcf.gz" "${OUT_C}/family.vcf.gz"
  run_in --rw "$OUT" --cpus 1 --memory 2g "$WHATSHAP_IMAGE" \
    python -c "import pysam; pysam.tabix_index('${OUT_C}/family.phased.vcf.gz', preset='vcf', force=True)"
  PHASED=(--phased "${OUT_C}/family.phased.vcf.gz")
else
  echo "No child has both parents sequenced: nothing to phase by inheritance."
fi

# --- 5. The report ------------------------------------------------------------------------
echo ""
echo "--- report"
CANDS=()
for s in "${MEMBERS[@]}"; do
  c="${GENOME_DIR}/${s}/clinical/${s}_clinical.vcf.gz"
  [ -f "$c" ] && CANDS+=(--candidates "${s}=/genome/${s}/clinical/${s}_clinical.vcf.gz")
done
rm -f "${OUT}"/*_inheritance.tsv.gz
"${PY[@]}" report --members "${OUT_C}/members.tsv" --relationships "${OUT_C}/relationships.tsv" \
  --vcf "${OUT_C}/family.vcf.gz" ${PHASED[@]+"${PHASED[@]}"} ${CANDS[@]+"${CANDS[@]}"} \
  --par-bed /pgp/par_grch38.bed --out-dir "$OUT_C"

echo ""
echo "=== Step 38 complete ==="
echo "Report: ${OUT}/family_report.md"
echo "Relationships: ${OUT}/relationships.tsv"
