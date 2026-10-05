#!/usr/bin/env bash
# Step 34 (CRAM archive) and CRAM input, on a copy of case 20's BAM under
# another sample name (HG002cram), so the cases that follow keep their BAM:
#   - a CRAM cut short after it was written, and a valid CRAM with half the
#     reads (both planted by a docker wrapper around the real samtools call):
#     the step refuses each, keeps the BAM and leaves no CRAM;
#   - the real CRAM: quickcheck passes, the read count equals the BAM's, it is
#     smaller; --delete-bam then removes the BAM;
#   - --restore writes the BAM back, flagstat equal to the original, and step
#     16b on it reports the depth it reports on the original;
#   - Nextflow: a cram,crai row runs mosdepth through CRAM_TO_BAM with the
#     depth of the BAM row beside it, whose CRAM_ARCHIVE writes a checked CRAM.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
C="${SAMPLE}cram"
A="${G}/${C}/aligned"
mkdir -p "$A"
ln -f "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${A}/${C}_sorted.bam"
ln -f "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" "${A}/${C}_sorted.bam.bai"
BAM_SUM=$(md5sum < "${A}/${C}_sorted.bam")
reads() { sam view -c "$@" 2>/dev/null; }
N_BAM=$(reads "${C}/aligned/${C}_sorted.bam")
echo "reads in the BAM: ${N_BAM}"

# --- 1. the two planted faults ---------------------------------------------------------
# A docker in front of the e2e shim: it runs every call, and after the call
# that writes <sample>_sorted.part.cram it damages that file as FAULT says.
# It takes its own folder off PATH first: the shim looks for the next docker
# on PATH too, and would otherwise call this wrapper back, without end.
mkdir -p "${CASE_TMP}/fault"
cat > "${CASE_TMP}/fault/docker" <<'WRAP'
#!/usr/bin/env bash
self_dir="$(cd "$(dirname "$0")" && pwd)"
real=""
IFS=: read -r -a dirs <<< "$PATH"
keep=()
for d in "${dirs[@]}"; do
  [ "$d" = "$self_dir" ] && continue
  keep+=("$d")
  if [ -z "$real" ] && [ -x "${d}/docker" ]; then real="${d}/docker"; fi
done
PATH=$(IFS=:; echo "${keep[*]}")
export PATH
[ -n "$real" ] || { echo "fault wrapper: no docker on PATH" >&2; exit 127; }
"$real" "$@" || exit $?
args=" $* "
case "$args" in
  *" samtools view "*" -C "*".part.cram "*) ;;
  *) exit 0 ;;
esac
out=$(grep -oE '/genome/[^ ]+\.part\.cram' <<< "$args" | head -n 1)
host="${GENOME_DIR}${out#/genome}"
case "${FAULT:-}" in
  truncate)
    size=$(stat -c %s "$host")
    truncate -s $((size / 2)) "$host"
    echo "fault: ${host} cut to $((size / 2)) bytes" >&2 ;;
  fewer)
    "$real" run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${GENOME_DIR}:/genome" "$SAMTOOLS_IMAGE" \
      samtools view -C -s 3.5 --reference /genome/reference/GRCh38_no_alt_analysis_set.fasta \
        -o "$out" "${out%.part.cram}.bam"
    echo "fault: ${host} rewritten with half the reads" >&2 ;;
esac
WRAP
chmod +x "${CASE_TMP}/fault/docker"
export SAMTOOLS_IMAGE

for fault in truncate fewer; do
  FAULT=$fault PATH="${CASE_TMP}/fault:${PATH}" run_step 34-cram-archive.sh "$C" --delete-bam
  check "${fault}: the step refuses (exit ${STEP_RC})" test "$STEP_RC" -ne 0
  check "${fault}: the planted fault happened" has "^fault: " "$(cat "$STEP_LOG")"
  check "${fault}: the BAM is kept, unchanged" test "$(md5sum < "${A}/${C}_sorted.bam" 2>/dev/null)" = "$BAM_SUM"
  check "${fault}: no CRAM is left" test ! -e "${A}/${C}_sorted.cram" -a ! -e "${A}/${C}_sorted.part.cram"
  cp "$STEP_LOG" "${CASE_TMP}/fault-${fault}.log"
done
check "truncate: quickcheck is what refused it" has 'fails samtools quickcheck' "$(cat "${CASE_TMP}/fault-truncate.log")"
check "fewer: the flagstat comparison refused it" has 'does not hold the same reads as the BAM' "$(cat "${CASE_TMP}/fault-fewer.log")"

# --- 2. the real CRAM ---------------------------------------------------------------------
run_step 34-cram-archive.sh "$C"
check_step_exit 34-cram-archive.sh
CRAM="${A}/${C}_sorted.cram"
check "the CRAM passes quickcheck" sam quickcheck "${C}/aligned/${C}_sorted.cram"
check_eq "the CRAM holds every read of the BAM" \
  "$(reads -T reference/GRCh38_no_alt_analysis_set.fasta "${C}/aligned/${C}_sorted.cram")" "$N_BAM"
BAM_B=$(stat -c %s "${A}/${C}_sorted.bam") CRAM_B=$(stat -c %s "$CRAM" 2>/dev/null || echo 0)
echo "BAM ${BAM_B} bytes, CRAM ${CRAM_B} bytes" | tee -a "$E2E_NOTES"
check "the CRAM is smaller than the BAM" test "$CRAM_B" -gt 0 -a "$CRAM_B" -lt "$BAM_B"
check "the BAM is kept without --delete-bam" test -s "${A}/${C}_sorted.bam"

run_step 34-cram-archive.sh "$C" --delete-bam
check_step_exit 34-cram-archive.sh
check "--delete-bam removed the BAM and its index" test ! -e "${A}/${C}_sorted.bam" -a ! -e "${A}/${C}_sorted.bam.bai"
check "the CRAM is there" test -s "$CRAM" -a -s "${CRAM}.crai"

# --- 3. back to BAM, and one step on it --------------------------------------------------------
run_step 34-cram-archive.sh "$C" --restore
check_step_exit 34-cram-archive.sh
check "the restored BAM passes quickcheck" sam quickcheck "${C}/aligned/${C}_sorted.bam"
check_eq "flagstat of the restored BAM equals the original's" \
  "$(sam flagstat "${C}/aligned/${C}_sorted.bam" 2>/dev/null | md5sum)" \
  "$(sam flagstat "${SAMPLE}/aligned/${SAMPLE}_sorted.bam" 2>/dev/null | md5sum)"
run_step 34-cram-archive.sh "$C" --restore
check "--restore refuses to write over a BAM (exit ${STEP_RC})" test "$STEP_RC" -ne 0
run_step 16b-mosdepth.sh "$C"
check_step_exit 16b-mosdepth.sh
check_eq "step 16b on the restored BAM gives the original's depth summary" \
  "$(md5sum < "${G}/${C}/mosdepth/${C}.mosdepth.summary.txt" 2>/dev/null)" \
  "$(md5sum < "${G}/${SAMPLE}/mosdepth/${SAMPLE}.mosdepth.summary.txt" 2>/dev/null)"

# --- 4. Nextflow: a CRAM row, and CRAM_ARCHIVE on a BAM row -----------------------------------
if command -v nextflow >/dev/null; then
  D="${CASE_TMP}/nf"
  mkdir -p "$D"
  V="${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
  {
    echo 'sample,vcf,vcf_index,bam,bam_index,cram,crai,sex'
    echo "${C},${V},${V}.tbi,,,${CRAM},${CRAM}.crai,female"
    echo "${SAMPLE},${V},${V}.tbi,${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam,${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai,,,female"
  } > "${D}/samplesheet.csv"
  cat "${D}/samplesheet.csv"
  ( cd "$D" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -work-dir "${D}/work" \
      --input "${D}/samplesheet.csv" --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" \
      --tools mosdepth,cram_archive --outdir "${D}/out" --max_cpus 4 --max_memory 14.GB ) > "${D}/run.log" 2>&1
  RC=$?
  grep -vE 'Pulling|Waiting|Verifying|Download complete|Pull complete|Already exists' "${D}/run.log"
  check_eq "the run with a CRAM row exits 0" "$RC" 0
  TRACE=$(find "${D}/out/pipeline_info" -name 'trace_*.txt' 2>/dev/null | LC_ALL=C sort | awk 'END {print}')
  ran() {
    awk -F'\t' -v p="$1" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
      { n = $c["name"]; sub(/ \(.*/, "", n); sub(/.*:/, "", n); if (n == p) k++ }
      END {print k + 0}' "${TRACE:-/dev/null}" 2>/dev/null
  }
  check_eq "CRAM_TO_BAM ran for the CRAM row" "$(ran CRAM_TO_BAM)" 1
  check_eq "CRAM_ARCHIVE ran for the BAM row only" "$(ran CRAM_ARCHIVE)" 1
  check_eq "MOSDEPTH ran for both rows" "$(ran MOSDEPTH)" 2
  check_eq "the CRAM row's depth summary equals the BAM row's" \
    "$(md5sum < "${D}/out/${C}/coverage/${C}.mosdepth.summary.txt" 2>/dev/null)" \
    "$(md5sum < "${D}/out/${SAMPLE}/coverage/${SAMPLE}.mosdepth.summary.txt" 2>/dev/null)"
  NC="${D}/out/${SAMPLE}/aligned/${SAMPLE}_sorted.cram"
  check "CRAM_ARCHIVE published the CRAM and its index" test -s "$NC" -a -s "${NC}.crai"
  check_eq "that CRAM holds every read of the BAM" \
    "$(docker run --rm -v "${G}:/genome:ro" -v "${D}/out:/out:ro" "$SAMTOOLS_IMAGE" \
         samtools view -c -T /genome/reference/GRCh38_no_alt_analysis_set.fasta "/out/${SAMPLE}/aligned/${SAMPLE}_sorted.cram" 2>/dev/null)" \
    "$N_BAM"
  check "no CRAM was written for the CRAM row" test ! -e "${D}/out/${C}/aligned/${C}_sorted.cram"
else
  fail "nextflow is not on PATH"
fi

in_genome "$BCFTOOLS_IMAGE" rm -rf "$C"

finish
