#!/usr/bin/env bash
# VCF_PRECHECK on a VCF whose header gives no contig lengths: it compares the
# REF bases of the first 1000 records with the reference, with real bcftools
# norm on the fixture reference.
#   - tests/fixtures/vcf/grch37_no_contigs.vcf: 60 invented SNVs on chr1 and
#     chr20 at GRCh37 positions, REF from GRCh37, no ##contig lines. 45 of
#     their REF bases are not GRCh38's at the same positions, so the run stops
#     before any analysis and names the count and the rate. Without the check
#     it ran ROH on the wrong build and exited 0.
#   - the fixture's DeepVariant call with its ##contig lines removed: its REF
#     bases are the reference's, so the run goes on and writes the ROH file.
. "$(dirname "$0")/lib.sh"
. "$(dirname "$0")/vendor-vcf-intake.inc"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

SPOT="${G}/spotcheck"
mkdir -p "$SPOT"

# only_precheck: true when VCF_PRECHECK is the only task the run started.
only_precheck() {
  local others
  others=$(trace_names | grep -v 'VCF_PRECHECK' | paste -sd' ' -)
  [ -z "$others" ] || { echo "    tasks besides VCF_PRECHECK: ${others}"; return 1; }
  trace_names | grep -q VCF_PRECHECK
}

# --- GRCh37 records, no ##contig lines ---------------------------------------------
cp "${REPO}/tests/fixtures/vcf/grch37_no_contigs.vcf" "${SPOT}/grch37.vcf"
in_genome "$BCFTOOLS_IMAGE" sh -c "set -e
  bcftools view -Oz -o spotcheck/${SAMPLE}.vcf.gz spotcheck/grch37.vcf
  bcftools index -f -t spotcheck/${SAMPLE}.vcf.gz"
check "the GRCh37 fixture has no ##contig line" lacks '^##contig' "$(gzip -dc "${SPOT}/${SAMPLE}.vcf.gz" | grep '^#')"
check_eq "the GRCh37 fixture has 60 records" "$(vcf_count "spotcheck/${SAMPLE}.vcf.gz")" 60

sheet "${CASE_TMP}/grch37.csv" "${SPOT}/${SAMPLE}.vcf.gz"
nf_run spot-grch37 "${CASE_TMP}/grch37.csv" roh
LOG=$(cat "$NF_LOG")
check "GRCh37 without contig lengths: the run fails" test "$NF_RC" -ne 0
check "GRCh37 without contig lengths: the message names the sample and the build" \
  has "Sample '${SAMPLE}': ${SAMPLE}\.vcf\.gz is not on GRCh38\. Its header gives no contig lengths" "$LOG"
check "GRCh37 without contig lengths: the message gives the count and the rate" \
  has 'they differ in 45 of the first 60 records \(75%\)' "$LOG"
check "GRCh37 without contig lengths: the message says what to do" has 'start from the FASTQ' "$LOG"
check "GRCh37 without contig lengths: no analysis task started" only_precheck
check "GRCh37 without contig lengths: no ROH file published" test ! -e "${NF_OUT}/${SAMPLE}/roh"

# --- GRCh38 records, no ##contig lines ---------------------------------------------
mkdir -p "${SPOT}/grch38"
in_genome "$BCFTOOLS_IMAGE" sh -c "set -e
  bcftools view -h ${SAMPLE}/vcf/${SAMPLE}.vcf.gz | grep -v '^##contig=' > spotcheck/nocontig_header.txt
  bcftools reheader -h spotcheck/nocontig_header.txt -o spotcheck/grch38/${SAMPLE}.vcf.gz ${SAMPLE}/vcf/${SAMPLE}.vcf.gz
  bcftools index -f -t spotcheck/grch38/${SAMPLE}.vcf.gz"
# Read the file itself: with a .tbi beside it, bcftools adds ID-only ##contig
# lines (no length) from the index.
check "the GRCh38 copy has no ##contig line" lacks '^##contig' "$(gzip -dc "${SPOT}/grch38/${SAMPLE}.vcf.gz" | grep '^#')"
N38=$(vcf_count "spotcheck/grch38/${SAMPLE}.vcf.gz")
check_ge "the GRCh38 copy has records" "$N38" 100
WANT=$(( N38 < 1000 ? N38 : 1000 ))

sheet "${CASE_TMP}/grch38.csv" "${SPOT}/grch38/${SAMPLE}.vcf.gz"
nf_run spot-grch38 "${CASE_TMP}/grch38.csv" roh
check_eq "GRCh38 without contig lengths: nextflow run exits 0" "$NF_RC" 0
check "GRCh38 without contig lengths: no build stop" lacks 'is not on GRCh38' "$(cat "$NF_LOG")"
check_ge "GRCh38 without contig lengths: ROH file has segments" \
  "$(grep -c '^RG' "${NF_OUT}/${SAMPLE}/roh/${SAMPLE}_roh.txt" 2>/dev/null || true)" 1
# The task log of this run's VCF_PRECHECK: the spot-check ran and found no difference.
TASK_OUT=$(grep -rl --include=.command.out "^${SAMPLE}: the header gives no contig lengths; REF bases differ" \
  "${E2E_WORK}/nf-work-intake" 2>/dev/null | xargs -r grep -h "REF bases differ" | sort -u)
echo "    VCF_PRECHECK said: ${TASK_OUT}"
check "GRCh38 without contig lengths: VCF_PRECHECK compared ${WANT} records and none differ" \
  has "REF bases differ from the reference in 0 of the first ${WANT} records \(0%\)" "$TASK_OUT"

finish
