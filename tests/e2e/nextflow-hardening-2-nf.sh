#!/usr/bin/env bash
# The Nextflow pipeline with the docker profile on the fixture, with the tools
# this package changed: manta (inversion conversion, --manta_call_regions),
# pypgx (bundle linked under the task's HOME), telomere_hunter (--cytoband)
# beside clinvar, pharmcat, mosdepth, delly, vcfanno and mito_haplogroup.
#   - every container ran with --network none (no task here needs the network);
#   - no task built a .fai index of the reference in its work directory;
#   - versions.yml has one line per tool, also for the version probes.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

G="$GENOME_DIR"
NF_OUT="${G}/nf-hardening"
WORK="${E2E_WORK}/nf-work-hardening"
SHEET="${CASE_TMP}/samplesheet.csv"
printf 'sample,vcf,vcf_index,bam,bam_index\n%s,%s,%s,%s,%s\n' "$SAMPLE" \
  "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" \
  "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" \
  > "$SHEET"

# Inputs earlier cases leave behind; built here when this case runs alone.
REGIONS="${G}/reference/fixture_regions.bed.gz"
if [ ! -s "${REGIONS}.tbi" ]; then
  sort -k1,1 -k2,2n "${FIXTURE_DIR}/regions.bed" | docker run --rm -i -u "$(id -u):$(id -g)" "$MANTA_IMAGE" \
    bash -c '"$(dirname "$(readlink -f "$(command -v configManta.py)")")/../libexec/bgzip" -c' > "$REGIONS"
  docker run --rm -u "$(id -u):$(id -g)" -v "${G}/reference:/r" "$MANTA_IMAGE" \
    bash -c '"$(dirname "$(readlink -f "$(command -v configManta.py)")")/../libexec/tabix" -f -p bed /r/fixture_regions.bed.gz'
fi
BUNDLE="${G}/reference/pypgx-bundle"
if [ ! -d "$BUNDLE" ]; then
  V="${PYPGX_IMAGE##*:}"
  git clone -q --branch "${V%%--*}" --depth 1 https://github.com/sbslee/pypgx-bundle.git "$BUNDLE"
fi
BANDS="${G}/reference/cytoBand.hg38.txt"
[ -s "$BANDS" ] || bash -c '. "$1/scripts/lib/common.sh" && install_data_file cytoband' _ "$REPO"

# Keep going after a failed task, so one run shows every broken module.
cat > "${CASE_TMP}/e2e.config" <<'NFCONF'
process.errorStrategy = 'ignore'
NFCONF

echo "+ nextflow run main.nf -profile docker (hardening)"
(
  cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
    -c "${CASE_TMP}/e2e.config" \
    -work-dir "$WORK" \
    --input "$SHEET" \
    --reference "${G}/reference/Homo_sapiens_assembly38.fasta" \
    --tools clinvar,pharmcat,mosdepth,delly,manta,vcfanno,pypgx,telomere_hunter,mito_haplogroup \
    --clinvar "${G}/clinvar/clinvar_pathogenic_chr.vcf.gz" \
    --clinvar_index "${G}/clinvar/clinvar_pathogenic_chr.vcf.gz.tbi" \
    --revel "${G}/annotations/revel_grch38.tsv.gz" \
    --revel_index "${G}/annotations/revel_grch38.tsv.gz.tbi" \
    --manta_call_regions "$REGIONS" \
    --pypgx_bundle "$BUNDLE" \
    --cytoband "$BANDS" \
    --outdir "$NF_OUT" \
    --max_cpus 4 --max_memory 14.GB
) 2>&1 | tee "$STEP_LOG"
NF_RC=${PIPESTATUS[0]}
check_eq "nextflow run exits 0" "$NF_RC" 0
cp "${CASE_TMP}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
TRACE=$(find "${NF_OUT}/pipeline_info" -name 'trace_*.txt' 2>/dev/null | LC_ALL=C sort | awk 'END {print}')
FAILED=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
                      $c["status"] != "COMPLETED" && $c["status"] != "CACHED" {print $c["name"] " (" $c["status"] ", exit " $c["exit"] ")"}' \
  "${TRACE:-/dev/null}" 2>/dev/null)
echo "tasks that did not complete: ${FAILED:-none}"
check "the trace lists the tasks" test -s "${TRACE:-/dev/null}"
check_eq "tasks that did not complete" "$(grep -c . <<< "$FAILED" || true)" 0
if [ -n "$FAILED" ]; then
  grep -E 'Error executing process|Command exit status|Command error' -A 6 "${CASE_TMP}/.nextflow.log" \
    | grep -vE 'Pulling|Waiting|Verifying|Download complete|Pull complete|Already exists' | head -80
fi

# --- No network ----------------------------------------------------------------------
TASKS=0 NET=0
while IFS= read -r run; do
  TASKS=$((TASKS + 1))
  grep -m1 'docker run' "$run" | grep -q -- '--network none' && NET=$((NET + 1))
done < <(find "$WORK" -name .command.run)
check_ge "tasks run" "$TASKS" 10
check_eq "tasks whose container ran with --network none" "$NET" "$TASKS"

# --- No .fai built in a task ---------------------------------------------------------------
# The reference index is staged as a link next to the FASTA; a regular .fai
# file in a work directory is one a tool built because the index was missing.
BUILT=$(find "$WORK" -type f -name '*.fai' | sed "s|^${WORK}/||")
echo "index files built in work directories: ${BUILT:-none}"
check_eq ".fai files built in work directories" "$(grep -c . <<< "$BUILT" || true)" 0

# --- Manta -------------------------------------------------------------------------------
R="${NF_OUT}/${SAMPLE}"
MV="${R}/manta/results/variants"
check "MANTA publishes the converted diploidSV.vcf.gz" vcf_ok "nf-hardening/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz"
check "MANTA publishes Manta's own diploidSV.raw.vcf.gz" vcf_ok "nf-hardening/${SAMPLE}/manta/results/variants/diploidSV.raw.vcf.gz"
MANTA_DIR=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $c["name"] ~ /:MANTA / {print $c["hash"]}' "${TRACE:-/dev/null}")
MANTA_OUT=$(cat "$WORK"/"${MANTA_DIR:-none}"*/.command.out 2>/dev/null)
check "the MANTA task logs the inversion conversion" \
  has 'Inversion conversion: [0-9]+ inversion breakend records in, [0-9]+ SVTYPE=INV records out' "$MANTA_OUT"
check "MANTA was configured with the call regions" grep -q -- '--callRegions fixture_regions.bed.gz' "$WORK"/"${MANTA_DIR:-none}"*/.command.sh
check_eq "inversion breakends left in the converted file" \
  "$(gzip -cd "${MV}/diploidSV.vcf.gz" 2>/dev/null | awk -F'\t' '!/^#/ && ($5 ~ /^\[/ || $5 ~ /\]$/) {
       m = $5; sub(/^[^][]*[][]/, "", m); sub(/:.*/, "", m); if (m == $1) n++ } END {print n + 0}')" 0

# --- pypgx: the bundle under the task's HOME -------------------------------------------------
CYP2D6=$(awk -F'\t' '$1 == "CYP2D6" {print $2; exit}' "${R}/pypgx/${SAMPLE}_pypgx_summary.tsv" 2>/dev/null)
echo "PYPGX CYP2D6: ${CYP2D6:-none}"
check "PYPGX returns a CYP2D6 result" test -n "${CYP2D6:-}"
check "the CYP2D6 result is not FAILED or N/A" lacks '^(FAILED|N/A)$' "${CYP2D6:-FAILED}"

# --- TelomereHunter with the GRCh38 bands ------------------------------------------------------
TH_DIR=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $c["name"] ~ /:TELOMERE_HUNTER / {print $c["hash"]}' "${TRACE:-/dev/null}")
check "TELOMERE_HUNTER passes the bands with -b" grep -q -- '-b cytoBand.hg38.txt' "$WORK"/"${TH_DIR:-none}"*/.command.sh
check "no warning about missing bands" lacks 'cytoband is not set' "$(cat "$STEP_LOG")"
SUM=$(find "${R}/telomere" -name '*_summary.tsv' 2>/dev/null | awk 'NR == 1')
check "TelomereHunter's summary has a tel_content column" has 'tel_content' "$(head -n 1 "${SUM:-/dev/null}" 2>/dev/null)"

# --- versions.yml -------------------------------------------------------------------------------
VERS="${NF_OUT}/pipeline_info/software_versions.yml"
cat "$VERS" 2>/dev/null
check_eq "software_versions.yml lines without a key" "$(grep -vc ':' "$VERS" 2>/dev/null || true)" 0
check "TelomereHunter's own version is one value" grep -qE '^ +telomerehunter_reported: [^ ]+$' "$VERS"
check "haplogrep3's own version is one value" grep -qE '^ +haplogrep3_reported: [^ ]+$' "$VERS"

finish
