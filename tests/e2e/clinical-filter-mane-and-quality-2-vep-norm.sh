#!/usr/bin/env bash
# Step 13 normalises the sample before VEP. VEP itself cannot run here (its
# cache is 26 GB), so a docker wrapper stands in for the VEP image only: it
# keeps the file VEP was given and copies it to VEP's output. bcftools runs
# for real, on a 60-base reference written below:
#   chrT 17 GCAG>G  the last CAG of a CAG repeat deleted, written right-shifted
#                   -> VEP gets chrT 8 ACAG>A, the left-aligned form gnomAD uses
#   chrT 30 G>A,C   a 1/2 call, AD 0,15,13 -> two records, one per ALT, each
#                   AD with the other ALT's reads on REF (13,15 and 15,13)
#   chrT 40 C>A,G   a 0/1 call -> the A record only; the G half is 0/0 and
#                   never reaches VEP
#   chrZ 5 A>G      a contig the reference lacks -> left out, and the step says so
# The Nextflow VEP_NORM process runs the same normalisation (docs/13-vep-annotation.md).
. "$(dirname "$0")/lib.sh"

S=NORM1
GD="${CASE_TMP}/genome"
rm -rf "$GD"
mkdir -p "${GD}/reference" "${GD}/${S}/vcf" "${GD}/${S}/vep" \
  "${GD}/vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38"
echo "species homo_sapiens" > "${GD}/vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt"
SEQ=ACGTTGCACAGCAGCAGCAGTTGACCTAGGATCCGATTACGGCATTAGCCATGGTACCAA
printf '>chrT\n%s\n' "$SEQ" > "${GD}/reference/mini.fa"
printf 'chrT\t%s\t6\t%s\t%s\n' "${#SEQ}" "${#SEQ}" "$(( ${#SEQ} + 1 ))" > "${GD}/reference/mini.fa.fai"
{
  printf '##fileformat=VCFv4.2\n##contig=<ID=chrT,length=%s>\n##contig=<ID=chrZ,length=100>\n' "${#SEQ}"
  printf '##FILTER=<ID=PASS,Description="All filters passed">\n'
  printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
  printf '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$S"
  printf 'chrT\t17\t.\tGCAG\tG\t50\tPASS\t.\tGT:AD\t0/1:12,10\n'
  printf 'chrT\t30\t.\tG\tA,C\t50\tPASS\t.\tGT:AD\t1/2:0,15,13\n'
  printf 'chrT\t40\t.\tC\tA,G\t50\tPASS\t.\tGT:AD\t0/1:10,9,0\n'
  printf 'chrZ\t5\t.\tA\tG\t50\tPASS\t.\tGT:AD\t0/1:8,7\n'
} > "${GD}/${S}/vcf/${S}.vcf"
docker run --rm --user "$(id -u):$(id -g)" -v "${GD}:/genome" -w /genome "$BCFTOOLS_IMAGE" \
  sh -c "bcftools view -Oz -o ${S}/vcf/${S}.vcf.gz ${S}/vcf/${S}.vcf && bcftools index -t ${S}/vcf/${S}.vcf.gz"
check "the sample VCF is indexed" test -s "${GD}/${S}/vcf/${S}.vcf.gz.tbi"

# The wrapper: the VEP image only; every other call goes to the docker the
# PATH without the wrapper finds (the e2e shim, which looks the PATH up again:
# with the wrapper still on it, the two would call each other forever).
WRAP="${CASE_TMP}/wrap"
mkdir -p "$WRAP"
cat > "${WRAP}/docker" <<'WRAPPER'
#!/usr/bin/env bash
set -euo pipefail
self=$(cd "$(dirname "$0")" && pwd)
rest=""
IFS=: read -r -a dirs <<< "$PATH"
for d in "${dirs[@]}"; do
  [ "$d" = "$self" ] || rest+="${rest:+:}${d}"
done
export PATH="$rest"
args=("$@") in="" out="" vep=0
for i in "${!args[@]}"; do
  case "${args[$i]}" in
    "$WRAP_VEP_IMAGE") vep=1 ;;
    --input_file) in=${args[$((i + 1))]} ;;
    -o) out=${args[$((i + 1))]} ;;
  esac
done
if [ "$vep" = 1 ] && [ -n "$in" ] && [ -n "$out" ]; then
  cp "${WRAP_GENOME}${in#/genome}" "$WRAP_SEEN"
  cp "${WRAP_GENOME}${in#/genome}" "${WRAP_GENOME}${out#/genome}"
  exit 0
fi
exec docker "$@"
WRAPPER
chmod +x "${WRAP}/docker"
export WRAP_VEP_IMAGE="$VEP_IMAGE" WRAP_GENOME="$GD" WRAP_SEEN="${CASE_TMP}/vep_input.vcf.gz"
rm -f "$WRAP_SEEN"

GENOME_DIR="$GD" REF_FASTA="${GD}/reference/mini.fa" PATH="${WRAP}:${PATH}" \
  run_step 13-vep-annotation.sh "$S"
check_step_exit 13-vep-annotation.sh
check "step 13 says a contig is not in the reference" has 'NOTICE: 1 contig\(s\) of .* are not in ' "$(cat "$STEP_LOG")"
SEEN=$(gzip -dc "$WRAP_SEEN" 2>/dev/null | grep -v '^#' | cut -f1,2,4,5 | tr '\t' ' ' | paste -sd';' -)
echo "VEP was given: ${SEEN:-nothing}"
check_eq "VEP gets the deletion left-aligned, the multiallelic records split, no 0/0 half, chrZ left out" \
  "${SEEN:-nothing}" "chrT 8 ACAG A;chrT 30 G A;chrT 30 G C;chrT 40 C A"
SEEN_AD=$(gzip -dc "$WRAP_SEEN" 2>/dev/null | awk -F'\t' '!/^#/ && $2 == 30 {split($9, f, ":"); split($10, v, ":"); for (i in f) if (f[i] == "AD") a = v[i]; print $5 "=" a}' | paste -sd';' -)
echo "chrT 30 AD after the split: ${SEEN_AD:-nothing}"
check_eq "the split 1/2 keeps the other ALT's reads on REF" "${SEEN_AD:-nothing}" "A=13,15;C=15,13"
check "step 13 wrote and indexed its output" test -s "${GD}/${S}/vep/${S}_vep.vcf.gz.tbi"
check "step 13 removed its temporary files" \
  test -z "$(find "${GD}/${S}/vep" -name '*.tmp*' -print -quit)"

finish
