#!/usr/bin/env bash
# 43-archaic.sh — Neanderthal and Denisovan ancestry: the share and the segments you carry
# Usage: ./scripts/43-archaic.sh <sample_name>
#
# Opt-in. IBDmix (Chen et al. 2020, Cell) compares your genome with three
# high-coverage archaic genomes, the Altai and Vindija 33.19 Neanderthals and
# the Altai Denisovan, and calls the stretches you share with one of them.
# It needs the archaic data setup.sh --archaic builds (reference/archaic/):
#   1. bin/archaic_report.py genotypes writes IBDmix's genotype file of each
#      chromosome: the archaic genotypes, yours from the VCF, and 100
#      stand-ins carrying the 1000 Genomes allele count of ARCHAIC_POP at each
#      site (IBDmix takes the modern frequency from the modern columns, and
#      one genome alone gives it none; docs/43-archaic.md says why this is
#      the same score);
#   2. IBDmix scores your column against each archaic genome (LOD 4);
#   3. bin/archaic_report.py report keeps the segments of 50 kb or more,
#      merges the two Neanderthals, and writes the share, the segments and
#      the ten longest with their protein-coding genes.
#
# Without the archaic data the step says so in one line and exits 0.
# Env: ARCHAIC_DIR (another data folder inside GENOME_DIR, or none to skip),
# ARCHAIC_POP (the 1000 Genomes group whose frequencies are used: ALL, the
# default, or AFR, AMR, EAS, EUR, SAS), THREADS (default 4).
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image PYTHON_IMAGE IBDMIX_IMAGE

DATA=${ARCHAIC_DIR:-${GENOME_DIR}/reference/archaic/${ARCHAIC_DATA_VERSION}}
POP=${ARCHAIC_POP:-ALL}
S="${GENOME_DIR}/${SAMPLE}"
VCF="${S}/vcf/${SAMPLE}.vcf.gz"
OUT="${S}/archaic"
WORK="${OUT}/work"

if [ "${ARCHAIC_DIR:-}" = none ] || [ ! -s "${DATA}/DONE" ]; then
  echo "Step 43 skipped: no archaic data at ${DATA}; install it with scripts/setup.sh --archaic ${GENOME_DIR}"
  exit 0
fi
case "$POP" in ALL|AFR|AMR|EAS|EUR|SAS) ;; *) echo "ERROR: ARCHAIC_POP must be ALL, AFR, AMR, EAS, EUR or SAS, got '${POP}'" >&2; exit 2 ;; esac
[ -f "$VCF" ] || { echo "ERROR: File not found: ${VCF}" >&2; exit 1; }
D=$(cpath "$DATA")

echo "============================================"
echo "  Step 43: Archaic ancestry (IBDmix)"
echo "  Archaic genomes: Altai and Vindija 33.19 Neanderthals, Altai Denisovan"
echo "  Modern frequencies: 1000 Genomes ${POP}"
echo "  Sample: ${SAMPLE}"
echo "  Output: ${OUT}"
echo "============================================"
echo ""

rm -rf "$WORK" "${OUT}/${SAMPLE}_archaic".* "${OUT}/${SAMPLE}_archaic_segments.tsv"
mkdir -p "${WORK}/ibdmix"
W=$(cpath "$WORK")

echo "[1/3] Writing IBDmix's genotype files (your VCF, the archaic genomes, the stand-ins)..."
run_in --cpus 1 --memory 4g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYTHON_IMAGE" \
  python3 /pgp-bin/archaic_report.py genotypes --table-dir "$D" --vcf "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
    --pop "$POP" --out-dir "${W}/gt"

# One IBDmix run per chromosome and archaic genome, THREADS at a time, with the
# options of Chen et al. 2020 (their j.cell.2020.01.012-workflow/config.yaml).
echo "[2/3] IBDmix against each archaic genome..."
# shellcheck disable=SC2016  # $1 and $2 belong to the inner bash
run_in --cpus "$THREADS" --memory 8g "$IBDMIX_IMAGE" bash -euo pipefail -c '
    one() {
      c=$(basename "$1" .gt.gz)
      gzip -dc "$1" > "$1.txt"
      for a in Altai Vindija Denisova; do
        ibdmix -g "$1.txt" -s "$2/gt/samples.txt" -n "$a" -o "$2/ibdmix/$c.$a.txt.part" \
          -d 4 -m 1 -a 0.01 -e 0.002 -c 2 > "$2/ibdmix/$c.$a.log" 2>&1 \
          || { tail -n 20 "$2/ibdmix/$c.$a.log"; exit 1; }
        mv "$2/ibdmix/$c.$a.txt.part" "$2/ibdmix/$c.$a.txt"
      done
      rm -f "$1.txt"
    }
    n=0
    for f in "$2"/gt/chr*.gt.gz; do
      one "$f" "$2" &
      n=$((n + 1))
      if [ "$n" -ge "$1" ]; then wait -n; n=$((n - 1)); fi
    done
    while [ "$n" -gt 0 ]; do wait -n; n=$((n - 1)); done' _ "$THREADS" "$W" \
  || { echo "ERROR: IBDmix failed (log above)" >&2; exit 1; }

echo "[3/3] Writing the segments and the summary..."
echo ""
GENES=$(data_file gencode_genes || true)
GARG=()
if [ -s "$GENES" ]; then
  GARG=(--genes "$(cpath "$GENES")")
else
  echo "  No GENCODE gene file (${GENES}): the segments are listed without their genes (setup.sh installs it)."
fi
run_in --cpus 1 --memory 2g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYTHON_IMAGE" \
  python3 /pgp-bin/archaic_report.py report --sample "$SAMPLE" --ibdmix-dir "${W}/ibdmix" \
    --spans "${W}/gt/spans.tsv" --pop "$POP" ${GARG[@]+"${GARG[@]}"} --out-dir "$(cpath "$OUT")"
for f in "${OUT}/${SAMPLE}_archaic.tsv" "${OUT}/${SAMPLE}_archaic_segments.tsv"; do
  [ -s "$f" ] || { echo "ERROR: step 43 wrote no ${f}" >&2; exit 1; }
done
# The genotype files are large (the stand-ins); IBDmix's raw output stays.
cp "${WORK}/gt/spans.tsv" "${WORK}/ibdmix/spans.tsv"
rm -rf "${WORK}/gt"
echo ""
echo "  Summary:  ${OUT}/${SAMPLE}_archaic.tsv and ${SAMPLE}_archaic.txt"
echo "  Segments: ${OUT}/${SAMPLE}_archaic_segments.tsv"
echo "  IBDmix's own output: ${WORK}/ibdmix"
echo "  See docs/43-archaic.md for what the numbers mean and where they end."
