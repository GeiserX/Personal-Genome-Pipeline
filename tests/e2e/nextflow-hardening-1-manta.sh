#!/usr/bin/env bash
# Step 04 (Manta) turns inversion breakend pairs into SVTYPE=INV records:
#   - on a planted Manta-style VCF with two inversion BND pairs (one open 3',
#     one open 5') and a deletion, the step converts Manta's results it finds
#     without calling again: two INV records, the deletion unchanged;
#   - on the fixture BAM, with THREADS=4 and MANTA_CALL_REGIONS set to the
#     fixture's regions, it calls only inside them, keeps Manta's own file as
#     diploidSV.raw.vcf.gz and logs that the conversion ran; a second run
#     skips the finished step.
. "$(dirname "$0")/lib.sh"

# bgzip and tabix from Manta's libexec, as the calling user, GENOME_DIR at /genome.
manta_lib() {
  docker run --rm -i -u "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -w /genome "$MANTA_IMAGE" \
    bash -c 'L="$(dirname "$(readlink -f "$(command -v configManta.py)")")/../libexec"; t=$1; shift; exec "$L/$t" "$@"' _ "$@"
}
records() { gzip -cd "$1" 2>/dev/null | grep -v '^#'; }

# The fixture's regions as a call-regions BED (Nextflow case 2 reads it too).
REGIONS="${GENOME_DIR}/reference/fixture_regions.bed.gz"
if [ ! -s "${REGIONS}.tbi" ]; then
  sort -k1,1 -k2,2n "${FIXTURE_DIR}/regions.bed" | manta_lib bgzip -c > "$REGIONS"
  manta_lib tabix -f -p bed /genome/reference/fixture_regions.bed.gz
fi
check "the call-regions BED is indexed" test -s "${REGIONS}.tbi"

# --- 1. Planted inversions ---------------------------------------------------------
INV="${SAMPLE}inv"
V="${GENOME_DIR}/${INV}/manta/results/variants"
mkdir -p "${GENOME_DIR}/${INV}/aligned" "$V"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${GENOME_DIR}/${INV}/aligned/${INV}_sorted.bam"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" "${GENOME_DIR}/${INV}/aligned/${INV}_sorted.bam.bai"
CHR20_LEN=$(awk '$1 == "chr20" {print $2}' "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta.fai")
{
  printf '##fileformat=VCFv4.1\n##source=GenerateSVCandidates 1.6.0\n##contig=<ID=chr20,length=%s>\n' "$CHR20_LEN"
  cat <<'HEADER'
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##ALT=<ID=DEL,Description="Deletion">
HEADER
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$INV"
  b=MantaBND:900:0:1:0:0:0
  c=MantaBND:901:0:1:0:0:0
  printf 'chr20\t10100000\t%s:0\tN\tN]chr20:10200000]\t500\tPASS\tSVTYPE=BND;MATEID=%s:1;IMPRECISE;CIPOS=-50,50;EVENT=%s:0\tGT\t0/1\n' "$b" "$b" "$b"
  printf 'chr20\t10150000\tMantaDEL:902:0:0:0:0:0\tN\t<DEL>\t500\tPASS\tEND=10151000;SVTYPE=DEL;SVLEN=-1000;IMPRECISE;CIPOS=-50,50;CIEND=-50,50\tGT\t0/1\n'
  printf 'chr20\t10200000\t%s:1\tN\tN]chr20:10100000]\t500\tPASS\tSVTYPE=BND;MATEID=%s:0;IMPRECISE;CIPOS=-40,40;EVENT=%s:0\tGT\t0/1\n' "$b" "$b" "$b"
  printf 'chr20\t10250000\t%s:0\tN\t[chr20:10300100[N\t500\tPASS\tSVTYPE=BND;MATEID=%s:1;IMPRECISE;CIPOS=-50,50;EVENT=%s:0\tGT\t0/1\n' "$c" "$c" "$c"
  printf 'chr20\t10300100\t%s:1\tN\t[chr20:10250000[N\t500\tPASS\tSVTYPE=BND;MATEID=%s:0;IMPRECISE;CIPOS=-50,50;EVENT=%s:0\tGT\t0/1\n' "$c" "$c" "$c"
} | manta_lib bgzip -c > "${V}/diploidSV.vcf.gz"
manta_lib tabix -f -p vcf "/genome/${INV}/manta/results/variants/diploidSV.vcf.gz"

run_step 04-manta.sh "$INV"
check_step_exit 04-manta.sh
LOG=$(cat "$STEP_LOG")
check "the step converts the results it finds without calling again" has 'converting the inversions only' "$LOG"
check "the log says the conversion ran (4 breakends in, 2 INV out)" \
  has 'Inversion conversion: 4 inversion breakend records in, 2 SVTYPE=INV records out' "$LOG"
check_eq "records in diploidSV.raw.vcf.gz (Manta's own file)" "$(records "${V}/diploidSV.raw.vcf.gz" | grep -c . || true)" 5
OUT=$(records "${V}/diploidSV.vcf.gz")
printf '%s\n' "$OUT"
check_eq "SVTYPE=INV records" "$(grep -c 'SVTYPE=INV' <<< "$OUT" || true)" 2
check_eq "breakend records left" "$(grep -c 'SVTYPE=BND' <<< "$OUT" || true)" 0
check "the 3' pair is one INV from 10100000 to its mate" \
  has $'^chr20\t10100000\tMantaINV[^\t]*\t[^\t]*\t<INV>\t.*END=10200000;SVTYPE=INV;SVLEN=100000' "$OUT"
check "the 5' pair is one INV, moved one base left as Manta's script does" \
  has $'^chr20\t10249999\tMantaINV[^\t]*\t[ACGTN]\t<INV>\t.*END=10300099;SVTYPE=INV' "$OUT"
check "the deletion passes through unchanged" has $'^chr20\t10150000\tMantaDEL:902:0:0:0:0:0\tN\t<DEL>' "$OUT"
check "the converted file is indexed" test -s "${V}/diploidSV.vcf.gz.tbi"
rm -rf "${GENOME_DIR:?}/${INV}"

# --- 2. The fixture BAM ------------------------------------------------------------
export MANTA_CALL_REGIONS="$REGIONS"
run_step 04-manta.sh "$SAMPLE"
check_step_exit 04-manta.sh
LOG=$(cat "$STEP_LOG")
V="${GENOME_DIR}/${SAMPLE}/manta/results/variants"
check "the log names the call regions" has "Call regions: ${REGIONS}" "$LOG"
check "the log says the conversion step ran" has 'Inversion conversion: [0-9]+ inversion breakend records in' "$LOG"
RAW_N=$(records "${V}/diploidSV.raw.vcf.gz" | grep -c . || true)
OUT_N=$(records "${V}/diploidSV.vcf.gz" | grep -c . || true)
echo "Manta records: ${RAW_N} as Manta wrote them, ${OUT_N} after the conversion"
check_ge "records in diploidSV.raw.vcf.gz" "$RAW_N" 1
check "the converted file is indexed" test -s "${V}/diploidSV.vcf.gz.tbi"
# Each inversion pair (ALT [p[t or t]p] with the mate on the same contig) becomes one record.
read -r BND_IN INV_OUT LEFT < <(
  { records "${V}/diploidSV.raw.vcf.gz" | sed 's/^/raw\t/'; records "${V}/diploidSV.vcf.gz" | sed 's/^/out\t/'; } \
    | awk -F'\t' '{ alt = $6; m = alt; sub(/^[^][]*[][]/, "", m); sub(/:.*/, "", m)
                    inv = (alt ~ /^\[/ || alt ~ /\]$/) && m == $2 }
                  $1 == "raw" && inv { b++ } $1 == "out" && $9 ~ /(^|;)SVTYPE=INV(;|$)/ { i++ } $1 == "out" && inv { l++ }
                  END { print b + 0, i + 0, l + 0 }')
echo "inversion breakends in: ${BND_IN}, INV records out: ${INV_OUT}, inversion breakends left: ${LEFT}"
check_eq "inversion breakends left after the conversion" "$LEFT" 0
check_eq "INV records out (two breakends each)" "$((INV_OUT * 2))" "$BND_IN"
check_eq "records after the conversion" "$OUT_N" "$((RAW_N - BND_IN / 2))"
OUTSIDE=$(records "${V}/diploidSV.raw.vcf.gz" | awk -F'\t' 'NR == FNR { s[NR] = $1; b[NR] = $2; e[NR] = $3; n = NR; next }
  { ok = 0; for (k = 1; k <= n; k++) if ($1 == s[k] && $2 > b[k] && $2 <= e[k]) { ok = 1; break }; if (!ok) o++ }
  END { print o + 0 }' "${FIXTURE_DIR}/regions.bed" -)
check_eq "calls outside the call regions" "$OUTSIDE" 0

run_step 04-manta.sh "$SAMPLE"
check_step_exit 04-manta.sh
check "a second run skips the finished step" has 'Manta already done' "$(cat "$STEP_LOG")"
unset MANTA_CALL_REGIONS

finish
