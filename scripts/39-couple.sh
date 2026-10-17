#!/usr/bin/env bash
# Step 39: couple carrier screen for two partners (opt-in)
# Usage: ./scripts/39-couple.sh <family_name> [pedigree.tsv]
# Input: the pedigree of step 38 (assets/pedigree.example.tsv: sample, sex,
#        father, mother, partner, tab separated; default
#        $GENOME_DIR/<family_name>/pedigree.tsv) with two samples marked as each
#        other's partner, and for each of them step 06's ClinVar hits
#        ($GENOME_DIR/<sample>/clinvar/<sample>_clinvar_hits.vcf). Read when
#        they are there: Delly's calls (step 19, <sample>/delly/<sample>_sv.vcf.gz),
#        ExpansionHunter's (step 09, <sample>/expansion_hunter/<sample>_eh.vcf),
#        Parascopy's SMN1 table (step 35, <sample>/paralogs/<sample>_smn_copy_number.tsv)
#        and GENCODE's exons (setup.sh: reference/gencode.v<release>.basic.exons.bed).
#        Step 38 does not have to run first.
# Output: $GENOME_DIR/<family_name>/, the family's folder, for each couple A and B:
#   couple_A_B.tsv   one row per gene of assets/carrier_genes.tsv where either
#                    partner has a finding: shared (both carry a change in one
#                    recessive gene), x_linked, or one_partner, with the chance
#                    for a child and a plain explanation
#   couple_A_B.md    the same in words, what was not assessed, and the residual risk
# and $GENOME_DIR/<sample>/couple/<sample>_couple_carrier.tsv, the same table,
# for each partner's own reports (steps 24 and generate-report.sh).
#
# A screen, not a diagnosis: docs/39-couple.md says what it does not see. The
# rules are bin/couple_report.py's, the code the Nextflow COUPLE_CARRIER
# process (--pedigree with partners) runs too.
#
# Env: CARRIER_GENES  the carrier list (default assets/carrier_genes.tsv)
#      CARRIER_EXONS  exon BED (default the GENCODE exons setup.sh installs),
#                     inside GENOME_DIR; without it deletions are not used
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
require_image PYTHON_IMAGE
OUT="${GENOME_DIR}/${FAMILY}"
OUT_C="/genome/${FAMILY}"
PEDIGREE=${PEDIGREE:-${OUT}/pedigree.tsv}
GENES=${CARRIER_GENES:-${PGP_ROOT}/assets/carrier_genes.tsv}
EXONS=${CARRIER_EXONS:-${GENOME_DIR}/reference/gencode.v${GENCODE_RELEASE}.basic.exons.bed}
case "$EXONS" in /*) ;; *) EXONS="${GENOME_DIR}/${EXONS}" ;; esac

echo "=== Step 39: couple carrier screen ==="
echo "Pedigree: ${PEDIGREE}"
echo "Carrier list: ${GENES}"
echo "Output: ${OUT}/"

for f in "$PEDIGREE" "$GENES"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done
if [ -d "${OUT}/vcf" ] || [ -d "${OUT}/aligned" ]; then
  echo "ERROR: ${OUT} holds a sample's vcf/ or aligned/ folder; give the family a name of its own." >&2
  exit 2
fi

# The samples the pedigree names, read as step 38 reads them (CRLF line ends,
# whitespace-only lines and spaces around a column name are allowed);
# couple_report.py reads the partners and checks everything else.
mapfile -t NAMED < <(awk -F'\t' '
  { sub(/\r+$/, "") }
  /^[ \t]*$/ { next }
  !h { sub(/^#+/, ""); for (i = 1; i <= NF; i++) { f = tolower($i); gsub(/^[ \t]+|[ \t]+$/, "", f); if (f == "sample") c = i }; h = 1; next }
  /^#/ || !c { next }
  { s = $c; gsub(/^[ \t]+|[ \t]+$/, "", s); if (s != "") print s }' "$PEDIGREE")
ARGS=()
for s in ${NAMED[@]+"${NAMED[@]}"}; do
  validate_sample "$s"
  if [ "$s" = "$FAMILY" ]; then
    echo "ERROR: the family name '${FAMILY}' is also a sample of the pedigree; its folder would be the sample's." >&2
    exit 2
  fi
  [ -d "${GENOME_DIR}/${s}" ] || continue
  for kind in "clinvar:clinvar/${s}_clinvar_hits.vcf:06" "delly:delly/${s}_sv.vcf.gz:19" \
              "eh:expansion_hunter/${s}_eh.vcf:09" "smn:paralogs/${s}_smn_copy_number.tsv:35"; do
    IFS=: read -r k rel step <<< "$kind"
    if [ -f "${GENOME_DIR}/${s}/${rel}" ]; then
      ARGS+=("--${k}" "${s}=/genome/${s}/${rel}")
      echo "  ${s}: ${rel} (step ${step})"
    fi
  done
done
if [ -f "$EXONS" ]; then
  EXONS_C=$(cpath "$EXONS") || exit 1
  ARGS+=(--exons "$EXONS_C")
  echo "Exons: ${EXONS}"
else
  echo "No exon list at ${EXONS}: deletions are not used (./scripts/setup.sh installs it)."
fi

mkdir -p "$OUT"
rm -f "${OUT}"/couple_*.tsv "${OUT}"/couple_*.md "${OUT}"/*_couple_carrier.tsv
# and each named sample's copy, so a failed screen or a sample no longer
# paired leaves no stale table for its own reports to read.
for s in ${NAMED[@]+"${NAMED[@]}"}; do rm -f "${GENOME_DIR}/${s}/couple/${s}_couple_carrier.tsv"; done
[ "$PEDIGREE" -ef "${OUT}/pedigree.tsv" ] || cp -f "$PEDIGREE" "${OUT}/pedigree.tsv"
run_in --rw "$OUT" --cpus 1 --memory 2g -v "${PGP_ROOT}/bin:/pgp-bin:ro" \
  -v "${PGP_ROOT}/assets/par_grch38.bed:/pgp/par_grch38.bed:ro" -v "${GENES}:/pgp/carrier_genes.tsv:ro" \
  "$PYTHON_IMAGE" python3 /pgp-bin/couple_report.py screen --pedigree "${OUT_C}/pedigree.tsv" \
    --genes /pgp/carrier_genes.tsv --par-bed /pgp/par_grch38.bed ${ARGS[@]+"${ARGS[@]}"} --out-dir "$OUT_C"

# Each partner's copy, where its own reports read it.
for t in "${OUT}"/*_couple_carrier.tsv; do
  [ -f "$t" ] || continue
  s=$(basename "$t" _couple_carrier.tsv)
  mkdir -p "${GENOME_DIR}/${s}/couple"
  mv -f "$t" "${GENOME_DIR}/${s}/couple/${s}_couple_carrier.tsv"
done

echo ""
echo "=== Step 39 complete ==="
for m in "${OUT}"/couple_*.md; do echo "Report: ${m}"; done
