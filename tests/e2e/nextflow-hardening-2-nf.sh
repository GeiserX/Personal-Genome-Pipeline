#!/usr/bin/env bash
# The Nextflow pipeline with the docker profile on the fixture, with the tools
# this package changed: manta (inversion conversion, --manta_call_regions),
# pypgx (bundle linked under the task's HOME), telomere_hunter (--cytoband)
# beside clinvar, pharmcat, mosdepth, delly, vcfanno and mito_haplogroup,
# and survivor_merge and y_haplogroup, so their process bodies run for real.
#   - every container ran with --network none (no task here needs the network);
#   - no task built a .fai index of the reference in its work directory;
#   - versions.yml has one line per tool, also for the version probes;
#   - SURVIVOR's consensus of the Manta and Delly calls has records, each
#     with SUPP 2; Yleaf read the chrY markers and wrote one prediction row.
# HG002 is male, so the samplesheet says so; indexcov reads the slices as
# female (case 34), hence --sex_check warn. Then the bash steps 22 and 37 run
# on the same calls and BAM in a GENOME_DIR of their own (step 37 reads the
# sex from a male indexcov .ped there), and the E2E workflow compares
# sv_merged/ and y_haplogroup/ with scripts/ci/parity-diff.sh.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

G="$GENOME_DIR"
NF_OUT="${G}/nf-hardening"
WORK="${E2E_WORK}/nf-work-hardening"
SHEET="${CASE_TMP}/samplesheet.csv"
printf 'sample,sex,vcf,vcf_index,bam,bam_index\n%s,male,%s,%s,%s,%s\n' "$SAMPLE" \
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
# Delly with the exclude map step 19 uses, so both sides call the same SVs.
EXCL="${G}/reference/delly_human.hg38.excl.tsv"
[ -s "$EXCL" ] || bash -c '. "$1/scripts/lib/common.sh" && install_data_file delly_exclude' _ "$REPO"
YDATA="${G}/reference/yleaf-${YLEAF_DATA_VERSION}/data"
check "Yleaf's marker tables install" "${REPO}/scripts/setup.sh" --yleaf-data "$G"

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
    --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" \
    --tools clinvar,pharmcat,mosdepth,delly,manta,vcfanno,pypgx,telomere_hunter,mito_haplogroup,survivor_merge,y_haplogroup \
    --clinvar "${G}/clinvar/clinvar_pathogenic_chr.vcf.gz" \
    --clinvar_index "${G}/clinvar/clinvar_pathogenic_chr.vcf.gz.tbi" \
    --revel "${G}/annotations/revel_grch38.tsv.gz" \
    --revel_index "${G}/annotations/revel_grch38.tsv.gz.tbi" \
    --manta_call_regions "$REGIONS" \
    --pypgx_bundle "$BUNDLE" \
    --cytoband "$BANDS" \
    --delly_exclude "$EXCL" \
    --yleaf_data "$YDATA" \
    --sex_check warn \
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

# --- SURVIVOR: the consensus of Manta and Delly ------------------------------------------
CONS="nf-hardening/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf.gz"
check "SURVIVOR_SORT publishes a readable consensus VCF" vcf_ok "$CONS"
check "the consensus VCF is indexed" nonempty "${CONS}.tbi"
bcf query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SUPP\t%INFO/SUPP_VEC\n' "$CONS" 2>/dev/null > "${CASE_TMP}/consensus.tsv"
echo "Consensus records (CHROM POS SVTYPE SUPP SUPP_VEC):"
cat "${CASE_TMP}/consensus.tsv"
check_ge "consensus records" "$(grep -c . "${CASE_TMP}/consensus.tsv" || true)" 1
check_eq "consensus records without SUPP 2 and SUPP_VEC 11 (manta, delly)" \
  "$(awk -F'\t' '$4 != 2 || $5 != "11"' "${CASE_TMP}/consensus.tsv" | grep -c . || true)" 0
check_eq "caller columns in the consensus (named manta, delly)" \
  "$(bcf query -l "$CONS" 2>/dev/null | paste -sd, -)" "manta,delly"

# --- Yleaf ----------------------------------------------------------------------------------
YT="${R}/y_haplogroup/${SAMPLE}_y_haplogroup.txt"
cat "$YT" 2>/dev/null
check_eq "Yleaf's table has a header and one sample row" "$(grep -c . "$YT" 2>/dev/null || true)" 2
check_ge "Y markers Yleaf read (the stub writes 0)" \
  "$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} NR == 2 {print $c["Valid_markers"]}' "$YT" 2>/dev/null)" 1
check "the declared sex goes on with a warning" has "indexcov infers .* --sex_check warn is set: going on with male" "$(cat "$STEP_LOG" "${CASE_TMP}/.nextflow.log" 2>/dev/null)"

# --- The bash steps 22 and 37 on the same calls and BAM, for parity-diff --------------------
# A GENOME_DIR of hard links: the reference, Yleaf's tables, step 04's Manta
# calls, step 19's Delly calls, the BAM and a male indexcov .ped (step 37
# reads the sex from it; the sample's own .ped says female).
B="${E2E_WORK}/genome-hardening"
rm -rf "$B"
mkdir -p "${B}/reference" "${B}/${SAMPLE}/aligned" "${B}/${SAMPLE}/indexcov" \
  "${B}/${SAMPLE}/manta/results/variants" "${B}/${SAMPLE}/delly"
ln -fL "${G}/reference/GRCh38_no_alt_analysis_set.fasta" "${G}/reference/GRCh38_no_alt_analysis_set.fasta.fai" "${B}/reference/"
cp -al "${G}/reference/yleaf-${YLEAF_DATA_VERSION}" "${B}/reference/"
ln -fL "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" "${B}/${SAMPLE}/aligned/"
ln -fL "${G}/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz" "${G}/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz.tbi" \
  "${B}/${SAMPLE}/manta/results/variants/"
ln -fL "${G}/${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz" "${G}/${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz.tbi" "${B}/${SAMPLE}/delly/"
printf '#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\tCNchrX\tCNchrY\n%s\t%s\t-9\t-9\t1\t-9\t1.0\t1.0\n' \
  "$SAMPLE" "$SAMPLE" > "${B}/${SAMPLE}/indexcov/indexcov-indexcov.ped"
for step in 22-survivor-merge.sh 37-y-haplogroup.sh; do
  echo "+ GENOME_DIR=${B} scripts/${step} ${SAMPLE}"
  GENOME_DIR="$B" "${REPO}/scripts/${step}" "$SAMPLE" 2>&1 | tee "$STEP_LOG"
  STEP_RC=${PIPESTATUS[0]}
  check_step_exit "$step"
done
check "step 22 wrote the bash consensus" test -s "${B}/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf.gz"
check "step 37 wrote the bash Y table" test -s "${B}/${SAMPLE}/y_haplogroup/${SAMPLE}_y_haplogroup.txt"

# --- versions.yml -------------------------------------------------------------------------------
VERS="${NF_OUT}/pipeline_info/software_versions.yml"
cat "$VERS" 2>/dev/null
check_eq "software_versions.yml lines without a key" "$(grep -vc ':' "$VERS" 2>/dev/null || true)" 0
check "TelomereHunter's own version is one value" grep -qE '^ +telomerehunter_reported: [^ ]+$' "$VERS"
check "haplogrep3's own version is one value" grep -qE '^ +haplogrep3_reported: [^ ]+$' "$VERS"

finish
