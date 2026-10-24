#!/usr/bin/env bash
# Which samplesheet run-all.sh writes when one from an earlier run is there.
#
# run-all.sh keeps <sample>/nextflow/samplesheet.csv from the last run so a
# rerun gives Nextflow the same inputs and -resume finds its tasks. But a run
# from FASTQ publishes aligned/<sample>_sorted.bam and vcf/<sample>.vcf.gz,
# and a run from a BAM publishes the VCF. Keeping that run's row then starts
# the next run from the FASTQ (minimap2 and DeepVariant again) or from the BAM
# (DeepVariant again), although the first-match rule in the script's header
# picks the BAM and the VCF. The old row is replaced only when the run that
# wrote it ended (main.nf's process_status.tsv with no failed, ignored or
# unfinished task); a stopped or failed run keeps its row so -resume goes on.
# Cases:
#   FASTQ row left over, BAM, VCF and gVCF now,
#     its run ended                             -> BAM+VCF row with the gVCF
#   BAM row left over, VCF now, its run ended   -> BAM+VCF row
#   BAM+VCF row, rerun                          -> the same file, unchanged
#   FASTQ row naming a pair outside fastq/, no
#     other input                               -> kept, unchanged
#   FASTQ row, BAM and VCF now, its run stopped
#     with a task unfinished                    -> kept, unchanged
#   the same, its run ended                     -> BAM+VCF row
#
#   bash tests/test_run_all_samplesheet.sh [RUN_ALL]
#       No Docker, no Nextflow, no Java (guard.yml runs it on every pull
#       request): run-all.sh runs under env -i with the fakes of
#       scripts/ci/fake-docker first on PATH, and the test reads the
#       samplesheet it wrote. RUN_ALL defaults to scripts/run-all.sh; the
#       script before this fix fails cases 1, 2 and 6, and the script
#       that ignores process_status.tsv fails case 5.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
RUN_ALL=${1:-${REPO}/scripts/run-all.sh}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
G="${WORK}/genome" FAILS=0
HEADER=sample,fastq_1,fastq_2,bam,bam_index,vcf,vcf_index,sex
mkdir -p "${G}/reference"
: > "${G}/reference/GRCh38_no_alt_analysis_set.fasta"

touchf() { local f; for f in "$@"; do mkdir -p "$(dirname "$f")"; printf 'placeholder\n' > "$f"; done; }
# run SAMPLE: run-all.sh SAMPLE male; its output in WORK/SAMPLE.out. The exit
# code is not checked: the samplesheet is written before nextflow starts.
run() {
  env -i PATH="${REPO}/scripts/ci/fake-docker:${PATH}" HOME="$WORK" \
    GENOME_DIR="$G" SKIP_VALIDATION=true FAKE_DOCKER_LOG="${WORK}/docker.log" \
    bash "$RUN_ALL" "$1" male > "${WORK}/$1.out" 2>&1 || true
}
sheet() { cat "${G}/$1/nextflow/samplesheet.csv" 2>/dev/null || echo "(no samplesheet)"; }
# leftover SAMPLE ROW: the samplesheet an earlier run wrote, with ROW's six file columns
leftover() { mkdir -p "${G}/$1/nextflow"; printf '%s\n%s,%s,male\n' "$HEADER" "$1" "$2" > "${G}/$1/nextflow/samplesheet.csv"; }
# ended SAMPLE [UNFINISHED]: the process_status.tsv main.nf writes as a run ends, one task left unfinished when given
ended() { mkdir -p "${G}/$1/nextflow"; printf 'process\tsucceeded\tcached\tfailed\tignored\tunfinished\nALIGN_MINIMAP2\t1\t0\t0\t0\t0\nDEEPVARIANT\t%s\t0\t0\t0\t%s\n' \
  "$([ -n "${2:-}" ] && echo 0 || echo 1)" "${2:-0}" > "${G}/$1/nextflow/process_status.tsv"; }
expect() {  # expect SAMPLE WHAT WANT
  if [ "$(sheet "$1")" = "$3" ]; then echo "PASS ${2}"
  else FAILS=$((FAILS + 1)); printf 'FAIL %s\n  want: %s\n  got:  %s\n' "$2" "$3" "$(sheet "$1")"; tail -n 20 "${WORK}/$1.out" | sed 's/^/  | /'; fi
}

# 1. A FASTQ run's row, after that run published the BAM, the VCF and the gVCF
S="${G}/s1" B="${G}/s1/aligned/s1_sorted.bam" V="${G}/s1/vcf/s1.vcf.gz" GV="${G}/s1/vcf/s1.g.vcf.gz"
touchf "${S}/fastq/s1_R1.fastq.gz" "${S}/fastq/s1_R2.fastq.gz" "$B" "${B}.bai" "$V" "${V}.tbi" "$GV" "${GV}.tbi"
leftover s1 "${S}/fastq/s1_R1.fastq.gz,${S}/fastq/s1_R2.fastq.gz,,,,"; ended s1
run s1
expect s1 "a FASTQ row gives way to the BAM, VCF and gVCF its run published" \
  "$(printf '%s,gvcf,gvcf_index\ns1,,,%s,%s,%s,%s,male,%s,%s' "$HEADER" "$B" "${B}.bai" "$V" "${V}.tbi" "$GV" "${GV}.tbi")"

# 2. A BAM run's row, after that run published the VCF (no gVCF here)
S="${G}/s2" B="${G}/s2/aligned/s2_sorted.bam" V="${G}/s2/vcf/s2.vcf.gz"
touchf "$B" "${B}.bai" "$V" "${V}.tbi"
leftover s2 ",,${B},${B}.bai,,"; ended s2
run s2
expect s2 "a BAM row gives way to the BAM and the VCF" \
  "$(printf '%s\ns2,,,%s,%s,%s,%s,male' "$HEADER" "$B" "${B}.bai" "$V" "${V}.tbi")"

# 3. What the kept samplesheet is for: a rerun from the same starting point
#    writes the same file, so -resume finds the tasks again
S="${G}/s3" B="${G}/s3/aligned/s3_sorted.bam" V="${G}/s3/vcf/s3.vcf.gz"
touchf "$B" "${B}.bai" "$V" "${V}.tbi"
run s3
before=$(sheet s3)
run s3
expect s3 "a BAM+VCF rerun writes the same samplesheet" "$before"
[ "$before" = "$(printf '%s\ns3,,,%s,%s,%s,%s,male' "$HEADER" "$B" "${B}.bai" "$V" "${V}.tbi")" ] \
  || { FAILS=$((FAILS + 1)); echo "FAIL the first BAM+VCF run wrote: ${before}"; }

# 4. A row naming files outside the usual folders, with nothing else there, is kept
S="${G}/s4"
touchf "${G}/elsewhere/a_R1.fastq.gz" "${G}/elsewhere/a_R2.fastq.gz"
leftover s4 "${G}/elsewhere/a_R1.fastq.gz,${G}/elsewhere/a_R2.fastq.gz,,,,"
run s4
expect s4 "a row whose files exist, with no other input, is kept" \
  "$(printf '%s\ns4,%s,%s,,,,,male' "$HEADER" "${G}/elsewhere/a_R1.fastq.gz" "${G}/elsewhere/a_R2.fastq.gz")"

# 5. A FASTQ run stopped while DeepVariant ran (the BAM and an older VCF are there): the row stays, so -resume goes on
S="${G}/s5" B="${G}/s5/aligned/s5_sorted.bam" V="${G}/s5/vcf/s5.vcf.gz"
touchf "${S}/fastq/s5_R1.fastq.gz" "${S}/fastq/s5_R2.fastq.gz" "$B" "${B}.bai" "$V" "${V}.tbi"
leftover s5 "${S}/fastq/s5_R1.fastq.gz,${S}/fastq/s5_R2.fastq.gz,,,,"; ended s5 1
run s5
expect s5 "the row of a run that did not end is kept" \
  "$(printf '%s\ns5,%s,%s,,,,,male' "$HEADER" "${S}/fastq/s5_R1.fastq.gz" "${S}/fastq/s5_R2.fastq.gz")"

# 6. The same, after a run that ended
S="${G}/s6" B="${G}/s6/aligned/s6_sorted.bam" V="${G}/s6/vcf/s6.vcf.gz"
touchf "${S}/fastq/s6_R1.fastq.gz" "${S}/fastq/s6_R2.fastq.gz" "$B" "${B}.bai" "$V" "${V}.tbi"
leftover s6 "${S}/fastq/s6_R1.fastq.gz,${S}/fastq/s6_R2.fastq.gz,,,,"; ended s6
run s6
expect s6 "the row of a run that ended gives way to its BAM and VCF" \
  "$(printf '%s\ns6,,,%s,%s,%s,%s,male' "$HEADER" "$B" "${B}.bai" "$V" "${V}.tbi")"

[ "$FAILS" -eq 0 ] || { echo "${FAILS} case(s) failed"; exit 1; }
echo "run-all.sh starts from the most processed input and keeps a samplesheet that still names it."
