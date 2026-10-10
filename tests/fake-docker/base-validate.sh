#!/usr/bin/env bash
# validate-setup.sh with Docker up, every required reference file present and
# a sample with a BAM and a VCF: all critical checks pass, exit 0.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1

run_expect 0 validate "${SCRIPTS}/validate-setup.sh" sample1
output_lacks validate 'unbound variable'
output_has validate 'Docker images are pulled'
output_has validate 'Whole pipeline:  \./scripts/run-all\.sh sample1 <male\|female>'
output_has validate 'or step by step, starting with:  \./scripts/06-clinvar-screen\.sh sample1'

# --- free disk space, through a df in front of the fake one --------------------------
mkdir -p "${CASE_WORK}/df-bin"
cat > "${CASE_WORK}/df-bin/df" <<'SHIM'
#!/usr/bin/env bash
echo "Filesystem 1024-blocks Used Available Capacity Mounted on"
echo "fakefs 4294967296 0 $(( ${FAKE_FREE_GB:?} * 1048576 )) 0% /"
SHIM
chmod +x "${CASE_WORK}/df-bin/df"
with_free() { local gb=$1; shift; PATH="${CASE_WORK}/df-bin:${PATH}" FAKE_FREE_GB=$gb "$@"; }

# Under 200 GB on a fresh run: a FAIL, its text matches its threshold.
with_free 150 run_expect 1 disk-low "${SCRIPTS}/validate-setup.sh" sample1
output_has disk-low '\[FAIL\].*Free disk space: 150 GB, under the 200 GB minimum'
output_lacks disk-low 'Need 500 GB'

# The same space with the sample's Nextflow work directory: a resume, so a WARN.
mkdir -p "${GENOME_DIR}/sample1/nextflow/work/ab"
with_free 150 run_expect 0 disk-resume "${SCRIPTS}/validate-setup.sh" sample1
output_has disk-resume "\[WARN\].*Free disk space: 150 GB, under the 200 GB minimum, but .*/sample1/nextflow/work exists"
output_has disk-resume 'a redone alignment'
output_lacks disk-resume '\[FAIL\].*Free disk space'
rm -rf "${GENOME_DIR}/sample1/nextflow"

# --- a FASTQ-only sample: a FASTQ start needs more space than a BAM start ------------
mkdir -p "${GENOME_DIR}/fq1/fastq"
for r in R1 R2; do printf 'placeholder\n' | gzip -c > "${GENOME_DIR}/fq1/fastq/fq1_${r}.fastq.gz"; done
with_free 300 run_expect 0 fastq-low "${SCRIPTS}/validate-setup.sh" fq1
output_has fastq-low '\[WARN\].*FASTQ start with 300 GB free: alignment alone can write about 450 GB'
output_has fastq-low 'Whole pipeline:  \./scripts/run-all\.sh fq1 <male\|female>'
output_has fastq-low 'starting with:  \./scripts/02-alignment\.sh fq1'
with_free 600 run_expect 0 fastq-ok "${SCRIPTS}/validate-setup.sh" fq1
output_lacks fastq-ok 'FASTQ start with'
# A BAM start with the same 300 GB gets no FASTQ warning.
with_free 300 run_expect 0 bam-300 "${SCRIPTS}/validate-setup.sh" sample1
output_lacks bam-300 'FASTQ start with'

# --- the VCF's genome build when its header does not show it -------------------------
# The hook answers `bcftools view -h` with header.txt and the REF spot-check
# (an `sh -c` that runs bcftools norm) with spot.txt: "<records> <kept>".
cat > "${CASE_WORK}/vcf-hook" <<'HOOK'
#!/usr/bin/env bash
shift   # the image
if [ "${1:-}" = bcftools ] && [ "${2:-}" = view ] && [ "${3:-}" = -h ]; then
  cat "${CASE_WORK}/header.txt"
elif [ "${1:-}" = sh ] && [[ "${3:-}" == *"bcftools norm -c x -f "* ]]; then
  echo "$3" > "${CASE_WORK}/spot.cmd"
  [ -f "${CASE_WORK}/spot.txt" ] || exit 1
  cat "${CASE_WORK}/spot.txt"
fi
exit 0
HOOK
chmod +x "${CASE_WORK}/vcf-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/vcf-hook"
COLS=$'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1'
vcf_header() { { echo '##fileformat=VCFv4.2'; [ -z "${1:-}" ] || echo "$1"; echo "$COLS"; } > "${CASE_WORK}/header.txt"; }

# GRCh38 chr1 length: OK, no spot-check.
vcf_header '##contig=<ID=chr1,length=248956422>'
rm -f "${CASE_WORK}/spot.cmd" "${CASE_WORK}/spot.txt"
run_expect 0 build-38 "${SCRIPTS}/validate-setup.sh" sample1
output_has build-38 '\[OK\].*VCF genome build: GRCh38'
[ ! -f "${CASE_WORK}/spot.cmd" ] || fail "build-38: the REF spot-check ran although the header shows GRCh38"

# GRCh37 chr1 length: FAIL.
vcf_header '##contig=<ID=chr1,length=249250621>'
run_expect 1 build-37 "${SCRIPTS}/validate-setup.sh" sample1
output_has build-37 '\[FAIL\].*VCF genome build: GRCh37/hg19'

# An unknown chr1 length: WARN, then the spot-check decides. All REF bases match: OK.
vcf_header '##contig=<ID=chr1,length=12345>'
echo "1000 1000" > "${CASE_WORK}/spot.txt"
run_expect 0 build-unknown "${SCRIPTS}/validate-setup.sh" sample1
output_has build-unknown "\[WARN\].*VCF chr1 length \(12345\) is neither GRCh38's"
output_has build-unknown '\[OK\].*VCF REF bases match the reference in all of the first 1000 records'
grep -q "bcftools norm -c x -f '/genome/reference/GRCh38_no_alt_analysis_set\.fasta'" "${CASE_WORK}/spot.cmd" \
  || fail "build-unknown: the spot-check does not compare with the reference: $(cat "${CASE_WORK}/spot.cmd")"
grep -q "head -n 1000" "${CASE_WORK}/spot.cmd" || fail "build-unknown: the spot-check does not read the first 1000 records"

# No ##contig lines and most REF bases differ: FAIL, the VCF is not on GRCh38.
vcf_header ''
echo "1000 262" > "${CASE_WORK}/spot.txt"
run_expect 1 build-none-wrong "${SCRIPTS}/validate-setup.sh" sample1
output_has build-none-wrong '\[WARN\].*VCF header has no ##contig lines'
output_has build-none-wrong '\[FAIL\].*VCF REF bases differ from the reference in 738 of the first 1000 records: the VCF is not on GRCh38'

# A few differ (5% or less): WARN only.
echo "1000 960" > "${CASE_WORK}/spot.txt"
run_expect 0 build-none-few "${SCRIPTS}/validate-setup.sh" sample1
output_has build-none-few '\[WARN\].*VCF REF bases differ from the reference in 40 of the first 1000 records$'

# Contig lines without chr1: WARN; the spot-check cannot run: WARN, never a silent pass.
vcf_header '##contig=<ID=chr20,length=64444167>'
rm -f "${CASE_WORK}/spot.txt"
run_expect 0 build-nochr1 "${SCRIPTS}/validate-setup.sh" sample1
output_has build-nochr1 '\[WARN\].*VCF header has no chr1 length'
output_has build-nochr1 "\[WARN\].*Could not compare the VCF's REF bases with the reference"
