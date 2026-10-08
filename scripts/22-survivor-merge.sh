#!/usr/bin/env bash
# 22-survivor-merge.sh — SV consensus: the calls two or more SV callers agree on
# Usage: ./scripts/22-survivor-merge.sh <sample_name>
#
# SURVIVOR merge (SURVIVOR_IMAGE) pairs the calls of the callers below when
# both breakpoints lie within 1,000 bp of each other, the SV type and strands
# agree and the event is at least 50 bp long (`SURVIVOR merge LIST 1000 2 1 1
# 0 50`), and keeps the events at least two callers support. Each record
# carries SUPP (how many callers) and SUPP_VEC (which: one digit per caller,
# in the order of sv_files.txt).
#
# Inputs, each used when it exists (two or more are needed): Manta (step 04),
# Delly (19), CNVpytor (18), TIDDIT (04a) and Sniffles2 (04c). GRIDSS (04b) is
# left out: it reports every event as a pair of breakends (SVTYPE=BND), which
# never match the DEL, DUP and INV records of the other callers.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image SURVIVOR_IMAGE BCFTOOLS_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
OUTDIR="${S}/sv_merged"
IN_DIR="${OUTDIR}/inputs"
OUT="${OUTDIR}/${SAMPLE}_sv_consensus.vcf.gz"
mkdir -p "$IN_DIR"

echo "============================================"
echo "  Step 22: SV consensus (SURVIVOR merge)"
echo "  Sample: ${SAMPLE}"
echo "  Output: ${OUT}"
echo "============================================"
echo ""

# CNVpytor's own table, when step 18 left no VCF: one DEL or DUP record per
# CNV, with a sample column, so SURVIVOR reads it as it reads the others.
CNV_TXT="${S}/cnvpytor/${SAMPLE}_cnvs.txt"
CNV_VCF="${S}/cnvpytor/${SAMPLE}_cnvs.vcf.gz"
if [ ! -f "$CNV_VCF" ] && [ -f "$CNV_TXT" ]; then
  echo "CNVpytor: converting ${CNV_TXT} to a VCF..."
  {
    echo "##fileformat=VCFv4.2"
    echo "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
    echo "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">"
    echo "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">"
    echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    awk '{printf "##contig=<ID=%s,length=%s>\n", $1, $2}' "${REF_FASTA}.fai"
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$SAMPLE"
    awk '$1 == "deletion" || $1 == "duplication" {
      split($2, a, ":"); split(a[2], b, "-")
      t = ($1 == "duplication") ? "DUP" : "DEL"
      len = (t == "DEL") ? -(b[2] - b[1]) : b[2] - b[1]
      printf "%s\t%s\t.\tN\t<%s>\t.\tPASS\tSVTYPE=%s;END=%s;SVLEN=%d\tGT\t./.\n", a[1], b[1], t, t, b[2], len
    }' "$CNV_TXT"
  } > "${IN_DIR}/cnvpytor_txt.vcf"
  run_in --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" bash -euo pipefail -c \
    'bcftools sort "$1" -Oz -o "$2.tmp" && mv "$2.tmp" "$2" && bcftools index -f -t "$2"' \
    _ "$(cpath "${IN_DIR}/cnvpytor_txt.vcf")" "$(cpath "$CNV_VCF")"
  rm -f "${IN_DIR}/cnvpytor_txt.vcf"
fi

# name<TAB>VCF, in the order SUPP_VEC reports them.
CALLERS=(
  "manta	${S}/manta/results/variants/diploidSV.vcf.gz"
  "delly	${S}/delly/${SAMPLE}_sv.vcf.gz"
  "cnvpytor	${CNV_VCF}"
  "tiddit	${S}/sv_tiddit/${SAMPLE}_sv.vcf.gz"
  "sniffles2	${S}/sv_sniffles/${SAMPLE}_sv.vcf.gz"
)
USED=()
rm -f "${IN_DIR}"/*.vcf "${IN_DIR}/sv_files.txt"
for c in "${CALLERS[@]}"; do
  name=${c%%	*} vcf=${c#*	}
  if [ ! -f "$vcf" ]; then
    echo "  [--] ${name}: no ${vcf#"$S"/}"
    continue
  fi
  # PASS (or unfiltered) records as plain text, which SURVIVOR reads, with
  # one sample column named after the caller: SURVIVOR names its output
  # columns after the input samples, and three columns all named after the
  # sample would be one name three times.
  run_in --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" bash -euo pipefail -c '
    bcftools view -f PASS,. "$1" | awk -v s="$3" '\''BEGIN { FS = OFS = "\t" }
      /^##/ { print; next }
      /^#CHROM/ {
        if (NF < 10) { print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"; sites = 1; print $1, $2, $3, $4, $5, $6, $7, $8, "FORMAT", s }
        else print $1, $2, $3, $4, $5, $6, $7, $8, $9, s
        next
      }
      { if (sites) print $1, $2, $3, $4, $5, $6, $7, $8, "GT", "./."; else print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10 }'\'' > "$2"' \
    _ "$(cpath "$vcf")" "$(cpath "${IN_DIR}/${name}.vcf")" "$name"
  n=$(grep -vc '^#' "${IN_DIR}/${name}.vcf" || true)
  echo "  [OK] ${name}: ${n} PASS records (${vcf#"$S"/})"
  echo "$(cpath "${IN_DIR}/${name}.vcf")" >> "${IN_DIR}/sv_files.txt"
  USED+=("$name")
done
if [ -f "${S}/sv_gridss/${SAMPLE}_gridss.vcf.gz" ]; then
  echo "  [--] gridss: left out, its breakend (BND) records do not match the other callers' DEL, DUP and INV"
fi
echo ""

if [ "${#USED[@]}" -lt 2 ]; then
  echo "ERROR: the consensus needs at least 2 SV callers; found ${#USED[@]} (${USED[*]:-none})." >&2
  echo "  Run at least two of: steps 04 (Manta), 19 (Delly), 18 (CNVpytor), 04a (TIDDIT), 04c (Sniffles2)." >&2
  exit 1
fi

echo "[1/2] SURVIVOR merge of ${USED[*]}: breakpoints within 1,000 bp, same type and strands, >= 50 bp, 2+ callers..."
MERGED="${OUTDIR}/${SAMPLE}_sv_merged.vcf"
rm -f "$MERGED"
run_in --cpus 1 --memory 4g "$SURVIVOR_IMAGE" \
  SURVIVOR merge "$(cpath "${IN_DIR}/sv_files.txt")" 1000 2 1 1 0 50 "$(cpath "$MERGED")"
# SURVIVOR exits 0 when it cannot open an input, so its output is the check.
if ! wrote_vcf "$MERGED"; then
  echo "ERROR: SURVIVOR wrote no VCF (${MERGED})." >&2
  exit 1
fi

echo "[2/2] Sorting, compressing and indexing..."
rm -f "${OUT}.tbi"
run_in --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" bash -euo pipefail -c \
  'bcftools sort "$1" -Oz -o "$2.tmp" && mv "$2.tmp" "$2" && bcftools index -f -t "$2"' \
  _ "$(cpath "$MERGED")" "$(cpath "$OUT")"
rm -f "$MERGED"

COUNT=$(run_in "$BCFTOOLS_IMAGE" bcftools view -H "$(cpath "$OUT")" | wc -l | tr -d ' ')
echo ""
echo "============================================"
echo "  SV consensus complete: ${SAMPLE}"
echo "  Callers (SUPP_VEC order): ${USED[*]}"
echo "  Consensus SVs (2+ callers): ${COUNT}"
echo "  Output: ${OUT}"
echo "============================================"
echo "By support (SUPP_VEC: one digit per caller above):"
run_in "$BCFTOOLS_IMAGE" bcftools query -f '%INFO/SVTYPE\t%INFO/SUPP_VEC\n' "$(cpath "$OUT")" \
  | sort | uniq -c | sort -rn | head -20
