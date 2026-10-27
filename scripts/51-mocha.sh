#!/usr/bin/env bash
# 51-mocha.sh — Mosaic chromosomal alterations (large gains, losses and CN-LOH) with MoChA
# Usage: ./scripts/51-mocha.sh <sample_name>
#
# Opt-in. A mosaic chromosomal alteration is a large gain, loss or copy-neutral
# loss of heterozygosity (CN-LOH) in a share of the blood cells. MoChA
# (freeseek/mocha, a bcftools plugin) finds one as an allele balance that leans
# the same way along a phased haplotype. The VCF's heterozygous SNVs are phased
# with Beagle against the 1000 Genomes high-coverage panel (setup.sh
# --mocha-data), MoChA calls on chr1 to chr22, and bin/mocha_report.py keeps the
# calls MoChA's UK Biobank filter keeps whose cell fraction is at or above the
# floor it prints for the sample's depth. Mosaic loss of Y and X is step 33's.
# The same commands as the MOCHA_* Nextflow processes. Research result, not a
# diagnosis (docs/51-mocha.md).
#
# Input:  vcf/<sample>.vcf.gz with .tbi, with GT, AD, DP and GQ (DeepVariant
#         writes all four); reference/mocha-<MOCHA_PANEL_VERSION>/ (setup.sh
#         --mocha-data; MOCHA_DATA overrides, inside GENOME_DIR), whose
#         cnps.bed (the panel's common germline CNVs) is MoChA's --cnp list
#         MOCHA_PHASED_VCF=<file inside GENOME_DIR>: a VCF already phased (GT
#         with '|'), read instead of the sample's; no Beagle run, no panel
#         needed (without the folder's cnps.bed MoChA runs with no --cnp list)
# Output: mosaic/<sample>_mocha.tsv          one row per MoChA call, reported or not and why
#         mosaic/<sample>_mocha_summary.tsv  status, depth, heterozygous sites, the floor
#         mosaic/<sample>_mocha_calls_raw.tsv, _mocha_stats.tsv  MoChA's own tables
#         mosaic/<sample>_mocha_sites.tsv     per chromosome: span, sites, hets, phased hets
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image MOCHA_IMAGE BEAGLE_IMAGE PYTHON_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
OUTDIR="${S}/mosaic"
W="${OUTDIR}/work"
P="${OUTDIR}/${SAMPLE}"
VCF="${S}/vcf/${SAMPLE}.vcf.gz"
PHASING=beagle
if [ -n "${MOCHA_PHASED_VCF:-}" ]; then
  VCF=$MOCHA_PHASED_VCF; PHASING=input
  case "$VCF" in /*) ;; *) VCF="${GENOME_DIR}/${VCF}" ;; esac
fi
DATA=${MOCHA_DATA:-${GENOME_DIR}/reference/mocha-${MOCHA_PANEL_VERSION}}
case "$DATA" in /*) ;; *) DATA="${GENOME_DIR}/${DATA}" ;; esac
# MoChA's GRCh38 MHC and KIR regions (README of freeseek/mocha): left out of the model, kept in the output.
MHC=chr6:27518932-33480487 KIR=chr19:54071493-54992731

echo "=== Mosaic chromosomal alterations (MoChA): ${SAMPLE} ==="
for f in "$VCF" "${VCF}.tbi"; do [ -f "$f" ] || { echo "ERROR: File not found: ${f}" >&2; exit 1; }; done
if [ "$PHASING" = beagle ] && { [ ! -s "${DATA}/DONE" ] || [ ! -f "${DATA}/cnps.bed" ]; }; then
  echo "ERROR: MoChA's phasing panel not found or older than its CNV list: ${DATA}/DONE and cnps.bed (run: ./scripts/setup.sh --mocha-data ${GENOME_DIR})" >&2
  exit 1
fi
CNP=""
if [ -f "${DATA}/cnps.bed" ]; then
  CNP=$(cpath "${DATA}/cnps.bed")
  echo "Common germline CNVs (--cnp): ${DATA}/cnps.bed ($(grep -c . "${DATA}/cnps.bed" || true) regions)"
else
  echo "[WARN] No ${DATA}/cnps.bed: MoChA runs with no --cnp list, so a short common germline CNV can be reported as mosaic"
fi
echo "VCF: ${VCF}$([ "$PHASING" = input ] && echo ' (already phased: no Beagle run)')"
[ "$PHASING" = input ] || echo "Phasing panel: ${DATA}"
rm -rf "$W"
mkdir -p "${W}/phase"
rm -f "${P}_mocha.tsv" "${P}_mocha_summary.tsv" "${P}_mocha_calls_raw.tsv" "${P}_mocha_stats.tsv" "${P}_mocha_sites.tsv"
V=$(cpath "$VCF") WC=$(cpath "$W") D=$(cpath "$DATA") O=$(cpath "$P")

# 1. The biallelic PASS SNVs of chr1 to chr22 with GT and AD; a genotype under
# 10 reads or GQ 20 set missing (MoChA's README); INFO/GC from the reference
# (mochatools). phase_in.vcf.gz holds the genotypes Beagle phases.
echo "[1/4] Biallelic SNVs of chr1 to chr22, GC content (bcftools +mochatools)..."
AUTO=$(printf 'chr%s,' $(seq 1 22)); AUTO=${AUTO%,}
# shellcheck disable=SC2016  # $1 to $4 belong to the inner bash
run_in --cpus 2 --memory 4g "$MOCHA_IMAGE" bash -euo pipefail -c '
    bcftools view --no-version -f PASS,. -m2 -M2 -v snps -t "$4" -Ou "$1" \
      | bcftools filter --no-version -e "FMT/DP<10 | FMT/GQ<20" --set-GT . -Ou \
      | bcftools annotate --no-version -x ID,QUAL,FILTER,INFO,^FMT/GT,^FMT/AD -Ou \
      | bcftools +mochatools --no-version -Ob -o "$2/prep.bcf" - -- -t GC -f "$3"
    bcftools index -f "$2/prep.bcf"
    bcftools view --no-version -e "GT=\"mis\"" -Oz -o "$2/phase_in.vcf.gz" "$2/prep.bcf"
    bcftools index -f -t "$2/phase_in.vcf.gz"
    bcftools index -s "$2/prep.bcf" | cut -f1 > "$2/chroms.txt"' _ "$V" "$WC" "$REF_FASTA_C" "$AUTO"
CHROMS=$(tr '\n' ' ' < "${W}/chroms.txt")
[ -n "${CHROMS// /}" ] || { echo "ERROR: no biallelic PASS SNV on chr1 to chr22 in ${VCF}" >&2; exit 1; }

# 2. Beagle against the panel, one chromosome at a time; a chromosome the panel
# lacks (a panel cut to regions) is left out.
if [ "$PHASING" = beagle ]; then
  PHASE=""
  for c in $CHROMS; do [ -s "${DATA}/panel/${c}.vcf.gz" ] && PHASE+="${c} "; done
  [ -n "$PHASE" ] || { echo "ERROR: the panel has none of the VCF's chromosomes (${CHROMS})" >&2; exit 1; }
  echo "[2/4] Phasing with Beagle against the panel (${PHASE% })..."
  # shellcheck disable=SC2016  # $1 to $4 belong to the inner bash
  run_in --cpus "$THREADS" --memory 16g "$BEAGLE_IMAGE" bash -euo pipefail -c '
      for c in $4; do
        beagle -Xmx14g gt="$2/phase_in.vcf.gz" chrom="$c" ref="$1/panel/$c.vcf.gz" map="$1/maps/$c.map" \
          out="$2/phase/$c" impute=false seed=1 nthreads="$3" > "$2/phase/$c.out" 2>&1 \
          || { tail -n 20 "$2/phase/$c.out"; exit 1; }
      done' _ "$D" "$WC" "$THREADS" "$PHASE" \
    || { echo "ERROR: Beagle failed (log above)" >&2; exit 1; }
  tr ' ' '\n' <<< "${PHASE% }" > "${W}/phased_chroms.txt"
else
  echo "[2/4] The VCF is already phased: no Beagle run."
  cp "${W}/chroms.txt" "${W}/phased_chroms.txt"
fi

# 3. The phased genotypes onto the prepared sites, then MoChA (WGS mode: AD).
# SITES_AWK sums each chromosome: first and last position, sites, hets, phased hets.
# shellcheck disable=SC2016  # awk's fields, not shell variables
SITES_AWK='$1 != c { if (c != "") print c, a, b, n, h, ph; c = $1; a = $2; n = h = ph = 0 }
  { b = $2; n++; if ($3 ~ /^(0[|\/]1|1[|\/]0)$/) { h++; if ($3 ~ /[|]/) ph++ } }
  END { if (c != "") print c, a, b, n, h, ph }'
echo "[3/4] MoChA (bcftools +mocha, ${MOCHA_VERSION})..."
# shellcheck disable=SC2016  # $1 to $7 belong to the inner bash
run_in --cpus 2 --memory 8g "$MOCHA_IMAGE" bash -euo pipefail -c '
    w=$1; chroms=$(paste -sd, "$w/phased_chroms.txt")
    if [ "$2" = beagle ]; then
      # Beagle writes no ##contig lines and concat refuses records on an
      # undeclared contig: give each Beagle VCF the contig lines of prep.bcf.
      bcftools view -h "$w/prep.bcf" | grep "^##contig=" > "$w/contigs.hdr"
      f=(); for c in $(cat "$w/phased_chroms.txt"); do
        { bcftools view -h "$w/phase/$c.vcf.gz" | grep "^##" | grep -v "^##contig="
          cat "$w/contigs.hdr"; bcftools view -h "$w/phase/$c.vcf.gz" | grep "^#CHROM"; } > "$w/phase/$c.hdr"
        bcftools reheader -h "$w/phase/$c.hdr" -o "$w/phase/$c.ctg.vcf.gz" "$w/phase/$c.vcf.gz"
        f+=("$w/phase/$c.ctg.vcf.gz")
      done
      bcftools concat --no-version -Ob -o "$w/phased_gt.bcf" "${f[@]}"
      bcftools index -f "$w/phased_gt.bcf"
      bcftools annotate --no-version -a "$w/phased_gt.bcf" -c FMT/GT -r "$chroms" -Ob -o "$w/phased.bcf" "$w/prep.bcf"
    else
      bcftools view --no-version -r "$chroms" -Ob -o "$w/phased.bcf" "$w/prep.bcf"
    fi
    bcftools index -f "$w/phased.bcf"
    cnp=(); [ -z "$7" ] || cnp=(--cnp "$7")
    bcftools +mocha --genome GRCh38 --no-version --mhc "$4" --kir "$5" "${cnp[@]}" \
      --calls "$3_mocha_calls_raw.tsv" --stats "$3_mocha_stats.tsv" --log "$w/mocha.log" "$w/phased.bcf"
    bcftools query -f "%CHROM\t%POS[\t%GT]\n" "$w/phased.bcf" | awk -F "\t" -v OFS="\t" "$6" > "$3_mocha_sites.tsv"
    grep -i -m 1 "^mocha [0-9]" "$w/mocha.log" > "$w/mocha_version.txt" || true' \
  _ "$WC" "$PHASING" "$O" "$MHC" "$KIR" "$SITES_AWK" "$CNP"
echo "  $(cat "${W}/mocha_version.txt")"

echo "[4/4] Calls above the floor (bin/mocha_report.py)..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 2g "$PYTHON_IMAGE" \
  python3 /pgp-bin/mocha_report.py --sample "$SAMPLE" --calls "${O}_mocha_calls_raw.tsv" --stats "${O}_mocha_stats.tsv" \
    --sites "${O}_mocha_sites.tsv" --phasing "$PHASING" --out-calls "${O}_mocha.tsv" --out-summary "${O}_mocha_summary.tsv"
cp "${W}/mocha.log" "${P}_mocha.log"
rm -rf "$W"

[ -s "${P}_mocha_summary.tsv" ] || { echo "ERROR: no summary written (${P}_mocha_summary.tsv)" >&2; exit 1; }
echo ""
echo "Calls: ${P}_mocha.tsv"
echo "Summary: ${P}_mocha_summary.tsv"
