#!/usr/bin/env bash
# scripts/chip-to-vcf.sh on the two ten-row vendor fixtures in tests/fixtures/chip/,
# with the real bcftools and Picard images, on a made-up 3 kb "hg19" and
# "GRCh38" (the same sequence under both naming styles) and a one-to-one chain.
#
# Checks, for each fixture:
#   1. the script exits 0 and writes vcf/<sample>.vcf.gz;
#   2. ten records, the no-call row as a missing genotype (./. or .);
#   3. the records sit on chr1, chrX, chrY and chrM (AncestryDNA codes 23-26
#      mapped to X, Y, X and MT, then renamed for GRCh38);
#   4. AncestryDNA only: every genotype is diploid (its two allele columns are
#      joined). The script before this fix read the file as 23andMe, so it
#      took the column header as a row, the first allele as the genotype and
#      the numbers 23-26 as chromosome names.
# Then a cut AncestryDNA file: the script stops with the line number of the
# damaged row and writes no VCF.
# Needs docker.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
# shellcheck source=../versions.env
. "${REPO}/versions.env"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

command -v docker >/dev/null 2>&1 || { echo "FAIL: docker is needed to run this test"; exit 1; }

# write_fasta FILE NAME...: one 3,000-base contig per NAME, the same sequence
# for each NAME at the same index, plus the .fai.
write_fasta() {
  local out=$1
  shift
  awk -v names="$*" -v fa="$out" -v fai="${out}.fai" 'BEGIN {
    n = split(names, name, " ")
    split("A C G T", base, " ")
    s = 7
    off = 0
    for (c = 1; c <= n; c++) {
      hdr = ">" name[c]
      print hdr > fa
      off += length(hdr) + 1
      seq = ""
      for (i = 1; i <= 3000; i++) { s = (s * 75 + 74) % 65537; seq = seq base[s % 4 + 1] }
      lines = 0
      for (i = 1; i <= 3000; i += 60) { print substr(seq, i, 60) > fa; lines++ }
      printf "%s\t3000\t%d\t60\t61\n", name[c], off > fai
      off += 3000 + lines
    }
  }'
}

# setup_genome DIR: the made-up references, dictionary and chain.
setup_genome() {
  local g=$1 c
  mkdir -p "${g}/reference_hg19" "${g}/reference" "${g}/liftover"
  write_fasta "${g}/reference_hg19/human_g1k_v37.fasta" 1 X Y MT
  write_fasta "${g}/reference/test38.fasta" chr1 chrX chrY chrM
  {
    printf '@HD\tVN:1.6\tSO:unsorted\n'
    for c in chr1 chrX chrY chrM; do printf '@SQ\tSN:%s\tLN:3000\n' "$c"; done
  } > "${g}/reference/test38.dict"
  # After the rename the hg19 side is chr1, chrX, chrY, chrM too.
  for c in chr1 chrX chrY chrM; do
    printf 'chain 3000 %s 3000 + 0 3000 %s 3000 + 0 3000 1\n3000\n\n' "$c" "$c"
  done | gzip -c > "${g}/liftover/hg19ToHg38.over.chain.gz"
}

bcf() { docker run --rm -v "${1}:/g:ro" "$BCFTOOLS_IMAGE" bcftools "${@:2}"; }

# run_fixture NAME FIXTURE FORMAT
run_fixture() {
  local name=$1 fixture=$2 format=$3 g="${WORK}/${1}" rc=0 vcf
  setup_genome "$g"
  mkdir -p "${g}/S1/raw"
  cp "${REPO}/tests/fixtures/chip/${fixture}" "${g}/S1/raw/S1_raw.txt"
  echo "--- ${name}: chip-to-vcf.sh S1 ${format}"
  GENOME_DIR="$g" REF_FASTA="${g}/reference/test38.fasta" \
    "${REPO}/scripts/chip-to-vcf.sh" S1 "$format" > "${WORK}/${name}.log" 2>&1 || rc=$?
  sed 's/^/    /' "${WORK}/${name}.log"
  if [ "$rc" -ne 0 ]; then fail "${name}: chip-to-vcf.sh exited ${rc}"; return; fi
  pass "${name}: chip-to-vcf.sh exited 0"
  vcf="${g}/S1/vcf/S1.vcf.gz"
  if [ ! -s "$vcf" ]; then fail "${name}: no ${vcf}"; return; fi
  bcf "$g" query -f '%CHROM\t%POS\t[%GT]\n' /g/S1/vcf/S1.vcf.gz > "${WORK}/${name}.tsv"
  sed 's/^/    /' "${WORK}/${name}.tsv"
  local n chroms
  n=$(grep -c . "${WORK}/${name}.tsv" || true)
  if [ "$n" -eq 10 ]; then pass "${name}: 10 records"; else fail "${name}: ${n} records, want 10"; fi
  local nocall
  nocall=$(awk -F'\t' '$1 == "chr1" && $2 == 404 {print $3}' "${WORK}/${name}.tsv")
  case "$nocall" in
    ./.|.) pass "${name}: the no-call row is a missing genotype (${nocall})" ;;
    *) fail "${name}: the no-call row at chr1:404 is '${nocall}', want a missing genotype" ;;
  esac
  chroms=$(cut -f1 "${WORK}/${name}.tsv" | sort -u | paste -sd' ' -)
  if [ "$chroms" = "chr1 chrM chrX chrY" ]; then pass "${name}: chromosomes ${chroms}"
  else fail "${name}: chromosomes '${chroms}', want 'chr1 chrM chrX chrY'"; fi
  if [ "$name" = ancestrydna ]; then
    local haploid
    haploid=$(awk -F'\t' '$3 !~ /[\/|]/' "${WORK}/${name}.tsv" | grep -c . || true)
    if [ "$haploid" -eq 0 ]; then pass "${name}: every genotype is diploid"
    else fail "${name}: ${haploid} genotype(s) not diploid"; fi
  fi
}

run_fixture ancestrydna ancestrydna.txt auto
run_fixture 23andme 23andme.txt auto

# A cut AncestryDNA file (its last row lost the second allele) stops with the
# line number and writes no VCF. The script before this check skipped the row
# and converted the rest.
cut_g="${WORK}/ancestrydna-cut"
setup_genome "$cut_g"
mkdir -p "${cut_g}/S1/raw"
{ cat "${REPO}/tests/fixtures/chip/ancestrydna.txt"; printf 'rs999\t1\t500\tA\n'; } > "${cut_g}/S1/raw/S1_raw.txt"
cut_line=$(wc -l < "${cut_g}/S1/raw/S1_raw.txt" | tr -d ' ')
echo "--- ancestrydna-cut: chip-to-vcf.sh S1 ancestrydna (line ${cut_line} has four columns)"
cut_rc=0
GENOME_DIR="$cut_g" REF_FASTA="${cut_g}/reference/test38.fasta" \
  "${REPO}/scripts/chip-to-vcf.sh" S1 ancestrydna > "${WORK}/ancestrydna-cut.log" 2>&1 || cut_rc=$?
sed 's/^/    /' "${WORK}/ancestrydna-cut.log"
if [ "$cut_rc" -ne 0 ]; then pass "ancestrydna-cut: chip-to-vcf.sh exited ${cut_rc}"
else fail "ancestrydna-cut: chip-to-vcf.sh exited 0 on a row with four columns"; fi
if grep -q "line ${cut_line}: want five tab-separated columns" "${WORK}/ancestrydna-cut.log"; then
  pass "ancestrydna-cut: the error names line ${cut_line}"
else fail "ancestrydna-cut: no error naming line ${cut_line}"; fi
if [ ! -e "${cut_g}/S1/vcf/S1.vcf.gz" ]; then pass "ancestrydna-cut: no VCF"
else fail "ancestrydna-cut: wrote a VCF from a cut file"; fi

if [ "$FAILS" -gt 0 ]; then
  echo "test_chip_to_vcf: ${FAILS} check(s) failed"
  exit 1
fi
echo "test_chip_to_vcf: all checks passed"
