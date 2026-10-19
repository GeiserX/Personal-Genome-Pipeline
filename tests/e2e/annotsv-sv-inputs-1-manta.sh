#!/usr/bin/env bash
# shellcheck disable=SC2016  # the $ expressions in bash -c bodies expand inside the Manta container
# Step 04 (Manta) keeps a whole-span inversion only where Manta saw both
# junctions. On a planted Manta-style VCF on the fixture's chr20, converted
# with the real Manta image (convertInversion.py, then bin/manta_inv_filter.py
# on the image's Python 2.7):
#   - a reciprocal inversion chr20:10.1-10.6 Mb (an INV3 and an INV5
#     junction, 101 bp apart) stays two INV records, INV_JUNCTIONS=2;
#   - a lone INV3 junction whose breakends lie 50 Mb apart is no INV record:
#     its two breakends stay, INV_JUNCTIONS=1, and no INV record is longer
#     than the paired one;
#   - a lone INV5 junction 400 bp long stays an INV, INV_JUNCTIONS=1 and
#     INVDUP_CANDIDATE;
#   - a deletion with REF 137 bp and ALT 2 bp passes unchanged (step 05
#     rewrites its length for AnnotSV);
#   - the log counts each kind;
#   - the control: convertInversion.py alone, on the same file and image,
#     writes the 50 Mb INV record.
. "$(dirname "$0")/lib.sh"

manta_run() {
  docker run --rm -i -u "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -w /genome "$MANTA_IMAGE" "$@"
}
records() { gzip -cd "$1" 2>/dev/null | grep -v '^#'; }

S="${SAMPLE}invj"
V="${GENOME_DIR}/${S}/manta/results/variants"
mkdir -p "${GENOME_DIR}/${S}/aligned" "$V"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${GENOME_DIR}/${S}/aligned/${S}_sorted.bam"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" "${GENOME_DIR}/${S}/aligned/${S}_sorted.bam.bai"
CHR20_LEN=$(awk '$1 == "chr20" {print $2}' "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta.fai")
REF137=$(printf 'A%.0s' $(seq 1 137))
b() {  # b POS ID ALT MATE: one breakend record
  printf 'chr20\t%s\t%s\tN\t%s\t500\tPASS\tSVTYPE=BND;MATEID=%s;IMPRECISE;CIPOS=-50,50\tGT\t0/1\n' "$1" "$2" "$3" "$4"
}
{
  printf '##fileformat=VCFv4.1\n##source=GenerateSVCandidates 1.6.0\n##contig=<ID=chr20,length=%s>\n' "$CHR20_LEN"
  cat <<'HEADER'
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##ALT=<ID=DEL,Description="Deletion">
HEADER
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$S"
  b 10100000 MantaBND:800:0:1:0:0:0:0 'N]chr20:10600000]' MantaBND:800:0:1:0:0:0:1
  b 10100101 MantaBND:800:1:2:0:0:0:0 '[chr20:10600101[N' MantaBND:800:1:2:0:0:0:1
  printf 'chr20\t10300000\tMantaDEL:801:0:1:0:0:0\t%s\tAC\t500\tPASS\tEND=10300136;SVTYPE=DEL;SVLEN=-136\tGT\t0/1\n' "$REF137"
  b 10600000 MantaBND:800:0:1:0:0:0:1 'N]chr20:10100000]' MantaBND:800:0:1:0:0:0:0
  b 10600101 MantaBND:800:1:2:0:0:0:1 '[chr20:10100101[N' MantaBND:800:1:2:0:0:0:0
  b 12000000 MantaBND:802:0:1:0:0:0:0 'N]chr20:62000000]' MantaBND:802:0:1:0:0:0:1
  b 30000000 MantaBND:803:0:1:0:0:0:0 '[chr20:30000400[N' MantaBND:803:0:1:0:0:0:1
  b 30000400 MantaBND:803:0:1:0:0:0:1 '[chr20:30000000[N' MantaBND:803:0:1:0:0:0:0
  b 62000000 MantaBND:802:0:1:0:0:0:1 'N]chr20:12000000]' MantaBND:802:0:1:0:0:0:0
} | manta_run bash -c 'L="$(dirname "$(readlink -f "$(command -v configManta.py)")")/../libexec"; "$L/bgzip" -c' \
  > "${V}/diploidSV.vcf.gz"
manta_run bash -c 'L="$(dirname "$(readlink -f "$(command -v configManta.py)")")/../libexec"; "$L/tabix" -f -p vcf "$1"' \
  _ "/genome/${S}/manta/results/variants/diploidSV.vcf.gz"

# The control: Manta's script alone writes the 50 Mb inversion.
CONTROL=$(manta_run bash -c 'L="$(dirname "$(readlink -f "$(command -v configManta.py)")")/../libexec"
  "$L/convertInversion.py" "$L/samtools" "$1" "$2"' _ "/genome/reference/GRCh38_no_alt_analysis_set.fasta" \
  "/genome/${S}/manta/results/variants/diploidSV.vcf.gz" | grep -v '^#')
check "control: convertInversion.py alone writes a 50 Mb INV" \
  has $'^chr20\t12000000\tMantaINV[^\t]*\t[^\t]*\t<INV>\t.*END=62000000;SVTYPE=INV;SVLEN=50000000' "$CONTROL"

run_step 04-manta.sh "$S"
check_step_exit 04-manta.sh
LOG=$(cat "$STEP_LOG")
check "the step converts the results it finds" has 'converting the inversions only' "$LOG"
check "the log counts each kind" \
  has 'Inversion check: 2 INV records with both junctions, 1 single-junction INV records kept \(1 shorter than 1000 bp marked INVDUP_CANDIDATE\), 1 single-junction inversions of 1000000 bp or more left as their breakend pairs' "$LOG"
OUT=$(records "${V}/diploidSV.vcf.gz")
printf '%s\n' "$OUT"
check "the reciprocal inversion stays two INV records with INV_JUNCTIONS=2" \
  test "$(grep -cE $'^chr20\t(10100000|10100100)\tMantaINV:800:[01]:[12]:0:0:0\t[^\t]*\t<INV>\t.*INV_JUNCTIONS=2' <<< "$OUT" || true)" -eq 2
LONGEST=$(awk -F'\t' '$8 ~ /(^|;)SVTYPE=INV(;|$)/ { e = $8; sub(/.*(^|;)END=/, "", e); sub(/;.*/, "", e); if (e - $2 > m) m = e - $2 } END { print m + 0 }' <<< "$OUT")
check "no INV record longer than the paired inversion (500 kb)" test "$LONGEST" -le 500001
check "the 50 Mb junction keeps its two breakends, INV_JUNCTIONS=1" \
  test "$(grep -cE $'^chr20\t(12000000|62000000)\tMantaBND:802:0:1:0:0:0:[01]\t.*SVTYPE=BND;.*INV_JUNCTIONS=1' <<< "$OUT" || true)" -eq 2
check "the 400 bp junction stays an INV, INV_JUNCTIONS=1 and INVDUP_CANDIDATE" \
  has $'^chr20\t29999999\tMantaINV:803:0:1:0:0:0\t.*<INV>\t.*INV_JUNCTIONS=1;INVDUP_CANDIDATE' "$OUT"
check "the REF 137 / ALT 2 deletion passes unchanged" \
  has $'^chr20\t10300000\tMantaDEL:801:0:1:0:0:0\tA{137}\tAC\t500\tPASS\tEND=10300136;SVTYPE=DEL;SVLEN=-136\tGT\t0/1$' "$OUT"
check "the header declares INV_JUNCTIONS" has '^##INFO=<ID=INV_JUNCTIONS,' "$(gzip -cd "${V}/diploidSV.vcf.gz" | grep '^##')"
check "the converted file is indexed" test -s "${V}/diploidSV.vcf.gz.tbi"
rm -rf "${GENOME_DIR:?}/${S}"

finish
