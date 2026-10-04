#!/usr/bin/env bash
# INDEXCOV stops the run, before any BAM step, when the samplesheet's sex
# disagrees with the sex indexcov infers, and names both. On the fixture's
# slices indexcov reads HG002 as female (CNchrX and CNchrY near 2, case 34),
# so the row that disagrees is the true one, male; the female row is the
# control that passes. Both use case 21's VCF and case 20's BAM through the
# VCF+BAM entry with mosdepth, a BAM step that must wait for the check.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

G="$GENOME_DIR"

# nf_run SEX: run the pipeline on a one-row samplesheet declaring SEX. Sets
# RC, LOG (its output) and TRACE (its trace file).
nf_run() {
  local sex=$1 dir="${CASE_TMP}/${1}"
  rm -rf "$dir"
  mkdir -p "$dir"
  printf 'sample,vcf,vcf_index,bam,bam_index,sex\n%s,%s,%s,%s,%s,%s\n' "$SAMPLE" \
    "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" \
    "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" \
    "$sex" > "${dir}/samplesheet.csv"
  echo "+ nextflow run main.nf -profile docker (declared ${sex})"
  ( cd "$dir" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
      -work-dir "${dir}/work" \
      --input "${dir}/samplesheet.csv" \
      --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" \
      --tools mosdepth \
      --outdir "${dir}/out" \
      --max_cpus 4 --max_memory 14.GB ) > "${dir}/run.log" 2>&1
  RC=$?
  LOG="${dir}/run.log"
  cat "$LOG"
  TRACE=$(find "${dir}/out/pipeline_info" -name 'trace_*.txt' 2>/dev/null | LC_ALL=C sort | awk 'END {print}')
}
# ran NAME: how many tasks of process NAME the trace lists, in any state
ran() {
  awk -F'\t' -v p="$1" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    { n = $c["name"]; sub(/ \(.*/, "", n); sub(/.*:/, "", n); if (n == p) k++ }
    END {print k + 0}' "${TRACE:-/dev/null}" 2>/dev/null
}

nf_run male
check "declared male: the run stops (exit ${RC})" test "$RC" -ne 0
check "the message names the sample, the declared and the inferred sex" \
  has "Sample '${SAMPLE}': the samplesheet says sex male, but indexcov infers female from the BAM index \\(CNchrX=[0-9.]+, CNchrY=[0-9.]+\\)" "$(cat "$LOG")"
check "the message says how to go on" has 'rerun with --sex_check warn' "$(cat "$LOG")"
check_eq "declared male: INDEXCOV tasks" "$(ran INDEXCOV)" 1
check_eq "declared male: MOSDEPTH tasks (none may start before the check)" "$(ran MOSDEPTH)" 0

nf_run female
check_eq "declared female (what indexcov infers here): the run exits 0" "$RC" 0
check "the log says indexcov agrees" has "Sample '${SAMPLE}': indexcov infers female .*declared female" "$(cat "$LOG")"
check_eq "declared female: MOSDEPTH tasks" "$(ran MOSDEPTH)" 1

finish
