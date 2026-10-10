#!/usr/bin/env bash
# run-all.sh: the whole pipeline for one sample: validate-setup.sh, then main.nf with -resume (a rerun redoes only what changed).
# Usage: GENOME_DIR=/data ./scripts/run-all.sh <sample> <male|female> [nextflow options]
#   Options after the sex go to `nextflow run` as given (--sex_check warn), except -bg. --max_cpus (THREADS) and --max_memory
#   cap each task (default: the host's); they do not limit how much runs at once: Nextflow fills the machine's CPUs and RAM.
#   A full run on 8 CPUs can take more than a day: start it inside tmux or screen, or with nohup. A closed terminal stops
#   it, and -resume then restarts the task that was running (DeepVariant is one task).
# Needs Docker, bash 4.4+, Java 17+ and Nextflow (NEXTFLOW_VERSION in versions.env). Results: GENOME_DIR/<sample>/.
# Input, first match: aligned/<sample>_sorted.bam with .bai (plus vcf/<sample>.vcf.gz with .tbi: not called
#   again, and the gVCF beside it, vcf/<sample>.g.vcf.gz with .tbi, goes to PharmCAT and PRS); fastq/<sample>_R1.fastq.gz
#   and _R2; the VCF alone. Kept in <sample>/nextflow/samplesheet.csv while its files exist and its BAM is this call's.
# Switches: SKIP_VALIDATION=true; THREADS=N (--max_cpus); SKIP_TRIM=true (--skip_trim);
#   INTERVALS="chr20 chr22" (--intervals); ALIGN_DIR=dir (the BAM from <sample>/dir/); TOOLS=a,b
#   (only these --tools names; cyrius, parascopy and y_haplogroup run only when named); KIR=true (KIR genes with HLA);
#   REF_FASTA, EH_CATALOG, MANTA_CALL_REGIONS, PARASCOPY_POPULATION, ANCESTRY_PANEL (--ancestry_ref; none: raw
#   scores) and PGSC_CALC_DIR (--pgsc_calc) as for the single steps; the panel and pgsc_calc are passed when installed.
# Run as scripts after the pipeline: GRIDSS=true (04b), IMPUTATION=true (14),
#   SOMATIC=true (29), EXTRA_CALLERS=gatk,freebayes,strelka2,octopus (03a-03d), BENCHMARK=true (needs
#   EXTRA_CALLERS or a second caller VCF); then the HTML report (24) and the text report, also after a failed pipeline.
#   A step without its data or BAM is skipped. A report-only tool whose task fails is skipped, and its step is marked
#   failed (logs/run_status.tsv). Results are published by hard link when work/ and GENOME_DIR share a filesystem,
#   you passed no -w, -c or --publish_dir_mode, and you run as root or Linux's fs.protected_hardlinks is 0 (then a user
#   may link the tasks' root-owned files).
set -euo pipefail
case "${1:-}" in -h|--help) sed -n '2,/^set -euo/p' "$0" | sed -e '$d' -e 's/^# \{0,1\}//'; exit 0 ;; esac
SAMPLE=${1:-} SEX=${2:-}
[ -n "$SAMPLE" ] && [[ "$SEX" =~ ^(male|female)$ ]] || { { echo "Usage: $0 <sample> <male|female> [nextflow options]"
  [ -z "$SEX" ] || echo "ERROR: sex must be 'male' or 'female', got '${SEX}'."; } >&2; exit 2; }
((BASH_VERSINFO[0] * 100 + BASH_VERSINFO[1] >= 404)) || { echo "ERROR: run-all.sh needs bash 4.4 or later (this is ${BASH_VERSION}); on macOS: brew install bash" >&2; exit 2; }
shift 2; SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd) USER_THREADS=${THREADS:-}
# With -bg nextflow returns at once, and the steps would be recorded and the reports written before the run ends.
for a in "$@"; do [ "$a" != -bg ] || { echo "ERROR: run-all.sh does not take -bg: it waits for the pipeline to record each step and write the reports. Run it inside tmux or screen, or with nohup." >&2; exit 2; }; done
export GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "${SCRIPT_DIR}/lib/common.sh"
validate_sample "$SAMPLE"
G=$GENOME_DIR S="${GENOME_DIR}/${SAMPLE}" LOG_DIR="${GENOME_DIR}/${SAMPLE}/logs" T=${TOOLS:-} C=${EXTRA_CALLERS:-}
on() { [[ "${!1:-}" =~ ^(true|1)$ ]]; }
java_major=$(java -version 2>&1 | awk -F'"' '/version "/ {split($2, v, "."); print (v[1] == 1 ? v[2] : v[1]); exit}') || true
if ! command -v nextflow >/dev/null 2>&1 || [[ ! "${java_major:-}" =~ ^[0-9]+$ ]] || [ "$java_major" -lt 17 ]; then
  printf '%s\n' "ERROR: run-all.sh starts the Nextflow pipeline and needs Java 17 or later (found: ${java_major:-none}) and Nextflow ($(command -v nextflow || echo 'not on PATH'))." \
    "  Install the release CI validates:  curl -s https://get.nextflow.io | NXF_VER=${NEXTFLOW_VERSION} bash && sudo mv nextflow /usr/local/bin/" \
    "  Or run the steps one by one: ./scripts/<step>.sh ${SAMPLE} (docs/getting-started.md)." >&2; exit 2; fi
on SKIP_VALIDATION || "${SCRIPT_DIR}/validate-setup.sh" "$SAMPLE" \
  || { echo "ERROR: setup validation failed; fix the items above. To bypass: SKIP_VALIDATION=true $0 ${SAMPLE} ${SEX}" >&2; exit 1; }
mkdir -p "${S}/nextflow" "$LOG_DIR"
GENOME_DIR="$G" bash "${PGP_ROOT}/bin/write_manifest.sh" "$SAMPLE" run-all.sh "$SEX" || echo "WARNING: could not write ${S}/run_manifest.tsv"
STATUS="${LOG_DIR}/run_status.tsv"  # bin/collect_summary.py marks results older than this run stale
printf '# run-all.sh: when this run started and how each step ended\nmeta\tstarted_epoch\t%s\nmeta\tstarted_utc\t%s\nmeta\tdeclared_sex\t%s\n' \
  "$(date +%s)" "$(date -u '+%Y-%m-%d %H:%M:%S UTC')" "$SEX" > "$STATUS"
NF=() SEL=() RUNS=() KNOWN="" OPTIN=()  # every step of a default run: it runs, or it is skipped with the reason
declare -A TOOL_OF=()                   # step number -> its --tools name, for the steps that run
need() { local f; for f in "$@"; do [ -e "$f" ] || { echo "data not installed: ${f#"$G"/}"; return; }; done; }
optin() { [[ ",${T// /}," == *",$1,"* ]] && return 1; echo "opt-in: add $1 to TOOLS"; }  # prints why it is skipped
plan() {  # plan "STEP Label" TOOL [REASON]: 0 when TOOL runs
  local r=${3:-} t=${TOOLS:-}; KNOWN+=" $2"
  [[ -n "$NOBAM" || " 16 16b 10 20 21 04 19 15 22 08 09 09b 18 05 28 32 35 37 " != *" ${1%% *} "* ]] || r="no BAM"
  [ -z "$t" ] || [[ ",${t// /}," == *",$2,"* ]] || r=${r:-not in TOOLS}
  if [ -z "$r" ]; then SEL+=("$2") RUNS+=("${1%% *}") TOOL_OF[${1%% *}]=$2; printf '  %-28s runs\n' "$1"; return 0; fi
  printf '  %-28s skipped    (%s)\n' "$1" "$r"; printf 'step\t%s\tskipped (%s)\n' "${1%% *}" "$r" >> "$STATUS"; return 1
}
arg() { [ ! -e "$2" ] || NF+=("$1" "$2"); }  # arg --param FILE: pass FILE when it exists
CV="${G}/clinvar/clinvar_pathogenic_chr.vcf.gz" A="${G}/annotations" CAT=${EH_CATALOG:-${G}/reference/expansionhunter_variant_catalog.json}
# Step 09 reads the catalog inside its image; the pipeline takes it as a file.
[ -n "${EH_CATALOG:-}" ] || [ -s "$CAT" ] || { run_in "$EXPANSIONHUNTER_IMAGE" cat /usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json \
  > "${CAT}.tmp" && [ -s "${CAT}.tmp" ] && mv "${CAT}.tmp" "$CAT"; } || rm -f "${CAT}.tmp"
VEPW=$(need "${G}/vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt") SCORES=()
for s in whole_genome_SNVs.tsv.gz:cadd_snv gnomad.genomes.r4.0.indel.tsv.gz:cadd_indel spliceai_scores.raw.snv.hg38.vcf.gz:spliceai_snv \
         spliceai_scores.raw.indel.hg38.vcf.gz:spliceai_indel revel_grch38.tsv.gz:revel AlphaMissense_hg38.tsv.gz:alphamissense; do
  f="${A}/${s%%:*}"; [ -e "${f}.tbi" ] || f=${f/.raw./.masked.}
  [ ! -e "${f}.tbi" ] || SCORES+=("--${s#*:}" "$f" "--${s#*:}_index" "${f}.tbi")
done
# The samplesheet: the last run's while every file it names exists and its BAM is this call's (ALIGN_DIR).
SHEET="${S}/nextflow/samplesheet.csv" row=""
B="${S}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam" V="${S}/vcf/${SAMPLE}.vcf.gz" F="${S}/fastq/${SAMPLE}_R"
[ ! -f "$SHEET" ] || row=$(awk -F, 'NR == 2 {print $2","$3","$4","$5","$6","$7}' "$SHEET")
IFS=, read -r -a cols <<< "$row"; for f in "${cols[@]}"; do [ -z "$f" ] || [ -e "$f" ] || row=""; done
[ "${cols[2]:-}" = "$B" ] || [ -z "${cols[2]:-}${ALIGN_DIR:-}" ] || row=""
[ -n "$row" ] || if [ -f "$B" ] && [ -f "${B}.bai" ] && [ -f "$V" ] && [ -f "${V}.tbi" ]; then row=",,${B},${B}.bai,${V},${V}.tbi"
elif [ -f "$B" ] && [ -f "${B}.bai" ]; then row=",,${B},${B}.bai,,"
elif [ -f "${F}1.fastq.gz" ] && [ -f "${F}2.fastq.gz" ]; then row="${F}1.fastq.gz,${F}2.fastq.gz,,,,"
elif [ -f "$V" ] && [ -f "${V}.tbi" ]; then row=",,,,${V},${V}.tbi"
else echo "ERROR: no input for ${SAMPLE}: no ${B}, no ${F}1/2.fastq.gz and no ${V}, each with its index." >&2; exit 1; fi
IFS=, read -r f1 _ fb _ fv _ <<< "$row"; NOBAM=${f1}${fb} GV="${S}/vcf/${SAMPLE}.g.vcf.gz" GCOL="" GROW=""
# A BAM+VCF row is not called again: the gVCF beside the VCF, when there is one with its index, gives PharmCAT and PRS the reference calls.
if [ -n "$fb" ] && [ -n "$fv" ]; then
  if [ -f "$GV" ] && [ -f "${GV}.tbi" ]; then GCOL=",gvcf,gvcf_index" GROW=",${GV},${GV}.tbi"; echo "NOTE: starting from the existing VCF; PharmCAT and PRS read the gVCF beside it (${GV})."
  else echo "NOTE: starting from the existing VCF, with no gVCF with its index beside it (${GV}): PharmCAT and PRS read its variant sites only. To call again with a gVCF, remove the VCF and its index."; fi
fi
printf 'sample,fastq_1,fastq_2,bam,bam_index,vcf,vcf_index,sex%s\n%s,%s,%s%s\n' "$GCOL" "$SAMPLE" "$row" "$SEX" "$GROW" > "$SHEET"
echo "[Input] ${SHEET}: ${row}"
if [ -n "$NOBAM" ]; then RUNS=(16); else printf 'step\t16\tskipped (no BAM)\n' >> "$STATUS"; fi
echo "[Steps] ${SAMPLE} (${SEX})"
for p in "07 PharmCAT:pharmcat" "27 CPIC lookup:cpic" "11 ROH:roh" "12 Mito haplogroup:mito_haplogroup" "16b mosdepth:mosdepth" \
         "10 TelomereHunter:telomere_hunter" "20 Mito variants (Mutect2):mito_variants" "04 Manta:manta" \
         "19 Delly:delly" "15 duphold:duphold" "22 SV consensus merge:survivor_merge" "28 MultiQC:multiqc"; do plan "${p%:*}" "${p##*:}" || true; done
plan "06 ClinVar screen" clinvar "$(need "$CV" "${CV}.tbi")" && NF+=(--clinvar "$CV" --clinvar_index "${CV}.tbi")
H=$(data_file hla_dat || true) GT=$(data_file gencode_genes || true)
plan "08 HLA typing (T1K)" hla_typing "$(need "$H" "$GT")" && NF+=(--hla_dat "$H" --hla_genes "$GT") \
  && if on KIR; then K="${G}/kir/IPD-KIR_${KIR_DB_RELEASE}/kir.dat"; [ -s "$K" ] && NF+=(--kir true --kir_dat "$K") \
    || echo "  KIR typing skipped: data not installed: ${K#"$G"/} (setup.sh --kir-data)"; fi
plan "09 ExpansionHunter" expansion_hunter "$(need "$CAT")" && NF+=(--expansion_catalog "$CAT")
plan "09b Stranger" stranger "$(need "$CAT")" || true
plan "13 VEP" vep "$VEPW" && NF+=(--vep_cache "${G}/vep_cache") && arg --gnomad_constraint "${A}/gnomad_v4.1_constraint.tsv"
NOSCORE=$([ "${#SCORES[@]}" -gt 0 ] || echo 'data not installed: no score file in annotations/')
plan "30 vcfanno" vcfanno "${VEPW:+needs VEP, }${VEPW:-$NOSCORE}" && NF+=("${SCORES[@]}")
plan "23 Clinical filter" clinical_filter "${VEPW:+needs VEP, ${VEPW}}" || true
plan "31 slivar" slivar "${VEPW:+needs VEP, ${VEPW}}" || true
plan "17 CPSR" cpsr "$(need "${G}/vep_cache/homo_sapiens/${PCGR_VEP_CACHE_RELEASE}_GRCh38/info.txt" "${G}/pcgr_data/${PCGR_DATA_BUNDLE}/data")" \
  && NF+=(--pcgr_data "${G}/pcgr_data/${PCGR_DATA_BUNDLE}" --vep_cache_cpsr "${G}/vep_cache")
plan "18 CNVpytor" cnvpytor "$(need "${G}/reference/cnvpytor/gc_hg38.pytor")" && NF+=(--cnvpytor_resources "${G}/reference/cnvpytor")
plan "05 AnnotSV" annotsv "$(need "${G}/annotsv_annotations/Annotations_Human/Genes/GRCh38")" && NF+=(--annotsv_annotations "${G}/annotsv_annotations")
plan "32 pypgx" pypgx "$(need "${G}/reference/pypgx-bundle")" && NF+=(--pypgx_bundle "${G}/reference/pypgx-bundle")
CY="${G}/tools/cyrius-${CYRIUS_VERSION}" PS="${G}/reference/parascopy-${PARASCOPY_DATA_VERSION}"
CYSTAMP="python=${PYTHON_IMAGE} lock=$(_digest sha256 "${SCRIPT_DIR}/cyrius-constraints.txt")"
plan "21 Cyrius CYP2D6" cyrius "$(optin cyrius || { [ "$(cat "${CY}/INSTALLED" 2>/dev/null)" = "$CYSTAMP" ] \
  || echo "data not installed: ${CY#"$G"/} for this version (setup.sh --cyrius)"; })" && NF+=(--cyrius_install "$CY")
plan "35 Parascopy SMN1/SMN2" parascopy "$(optin parascopy || need "${PS}/homology_table/GRCh38.bed.gz")" && NF+=(--parascopy_data "$PS") \
  && NF+=(--parascopy_population "${PARASCOPY_POPULATION:-EUR}")
plan "25 PRS" prs "$(need "$(compgen -G "${G}/prs_scores/*.txt.gz" | head -n 1 || echo "${G}/prs_scores/<PGS id>.txt.gz")")" && NF+=(--pgs_scoring "${G}/prs_scores") \
  && { PC=${PGSC_CALC_DIR:-${G}/tools/pgsc_calc-${PGSC_CALC_VERSION}}; [ ! -f "${PC}/main.nf" ] || NF+=(--pgsc_calc "$PC"); }
# Step 26 is pgsc_calc's projection inside the PRS run, with its panel and the site list setup.sh writes beside it.
PANEL=${ANCESTRY_PANEL:-${G}/reference/pgsc_calc/${PGSC_PANEL}.tar.zst}
plan "26 Ancestry (pgsc_calc)" ancestry "$(if [[ " ${SEL[*]} " != *" prs "* ]]; then echo 'needs PRS'; elif [ "$PANEL" = none ]; then echo 'ANCESTRY_PANEL=none'
  else need "$PANEL" "${PANEL%.tar.zst}_GRCh38_sites.tsv"; fi)" && NF+=(--ancestry_ref "$PANEL")
YD="${G}/reference/yleaf-${YLEAF_DATA_VERSION}/data"
plan "37 Y haplogroup (Yleaf)" y_haplogroup "$(optin y_haplogroup || need "${YD}/hg38/new_positions.txt")" && NF+=(--yleaf_data "$YD")
for t in ${T//,/ }; do [[ "${KNOWN} " == *" ${t} "* ]] || { echo "ERROR: unknown step '${t}' in TOOLS. Known:${KNOWN}" >&2; exit 2; }; done
arg --cytoband "$(data_file cytoband || true)"; arg --delly_exclude "$(data_file delly_exclude || true)"; arg --manta_call_regions "${MANTA_CALL_REGIONS:-}"
[ -z "${INTERVALS:-}" ] || NF+=(--intervals "$INTERVALS"); [[ " $* " == *" --max_cpus"* ]] || NF+=(--max_cpus "${USER_THREADS:-$(getconf _NPROCESSORS_ONLN)}")
M=$(awk '/^MemTotal:/ {print int($2 / 1048576)}' /proc/meminfo 2>/dev/null || echo $(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1073741824 )))
[[ " $* " == *" --max_memory"* ]] || [ "${M:-0}" -lt 1 ] || NF+=(--max_memory "${M}.GB")
# Publish by hard link when it is sure to work: the BAM and every other output are then stored once, not twice.
# Nextflow does not fall back to a copy when a link fails; the run stops. The tasks run as root (the docker profile)
# and their files can stay root's, and Linux lets another user link a file it does not own only when
# fs.protected_hardlinks is 0. So: root or protected_hardlinks 0, and a probe link from work/ to GENOME_DIR that works.
# The probe only sees the default work/: with -w, -work-dir, NXF_WORK, or a -c/-config that may set workDir, keep the copy.
if [[ " $* " != *" --publish_dir_mode"* && " $* " != *" -w "* && " $* " != *" -work-dir"* && " $* " != *" -c "* \
   && " $* " != *" -config "* ]] && [ -z "${NXF_WORK:-}" ] \
   && { [ "$(id -u)" -eq 0 ] || [ "$(cat /proc/sys/fs/protected_hardlinks 2>/dev/null)" = 0 ]; }; then
  LP="${S}/nextflow/work/.pgp-link-probe.$$"
  mkdir -p "${S}/nextflow/work" && : > "$LP" && ln "$LP" "${S}/.pgp-link-probe.$$" 2>/dev/null && ln "$LP" "${G}/.pgp-link-probe.$$" 2>/dev/null \
    && NF+=(--publish_dir_mode link)
  rm -f "$LP" "${S}/.pgp-link-probe.$$" "${G}/.pgp-link-probe.$$"
fi
if on SKIP_TRIM; then NF+=(--skip_trim true); fi
for v in GRIDSS:04b-gridss IMPUTATION:14-imputation-prep SOMATIC:29-mutect2-somatic; do if on "${v%%:*}"; then OPTIN+=("${v#*:}.sh"); fi; done
for c in ${C//,/ }; do OPTIN+=("$(cd "$SCRIPT_DIR" && compgen -G "03[a-d]-${c}*.sh")") || { echo "ERROR: unknown caller '${c}' in EXTRA_CALLERS (gatk, freebayes, strelka2, octopus)" >&2; exit 2; }; done
if ! on BENCHMARK; then :; elif [ -n "$C" ] || compgen -G "${S}/vcf_*/${SAMPLE}.vcf.gz" >/dev/null || [ -f "${S}/vcf_strelka2/results/variants/variants.vcf.gz" ]; then OPTIN+=(benchmark-variants.sh)
else echo "  benchmark-variants skipped (only one caller VCF: set EXTRA_CALLERS, or run a 03a-03d script first)"; fi
[ -z "${MAX_JOBS:-}" ] || echo "NOTE: MAX_JOBS is no longer read. --max_cpus (THREADS) and --max_memory cap each task; Nextflow fills the machine's CPUs and RAM."
[ -z "${ANCESTRY:-}" ] || echo "NOTE: ANCESTRY is no longer read: step 26 runs in the pipeline whenever the ancestry panel is installed (setup.sh --ancestry-panel)."
export NXF_VER=${NXF_VER:-$NEXTFLOW_VERSION}
echo "[Nextflow ${NXF_VER}] launch directory ${S}/nextflow: .nextflow.log, and every task's files under work/"
FT="${S}/nextflow/failed_tasks.tsv"; rm -f "$FT"  # main.nf writes it when a task failed
rc=0; (cd "${S}/nextflow" && nextflow run "${PGP_ROOT}/main.nf" -profile docker -resume --input "$SHEET" --reference "$REF_FASTA" \
  --outdir "$G" --tools "$(IFS=,; echo "${SEL[*]}")" "${NF[@]}" "$@") || rc=$?
# Each failed task's step is 'failed'. A report-only tool's failure is ignored by the pipeline (nextflow.config), which then
# goes on; any other failure stops it, and then every step it started that is not failed is 'not finished' (a report
# section whose file this run wrote still shows as current). Stranger reads ExpansionHunter's calls, so it does not run
# when ExpansionHunter failed.
declare -A BAD=(); failed=0
while IFS=$'\t' read -r p st; do
  [ "$p" != process ] && [ -n "$p" ] || continue
  p=${p##*:}; case "$p" in Y_*) t=y_haplogroup ;; *) t=${p,,} ;; esac
  [[ " ${SEL[*]} " == *" ${t} "* ]] || { echo "  ${p}: ${st}" >&2; continue; }
  BAD[$t]=1 failed=1 rest=$(IFS=,; echo "${SEL[*]}"); rest=",${rest},"; rest=${rest/,${t},/,}; rest=${rest#,}
  echo "  ${t} (${p}) failed$([ "$st" != ignored ] || echo ', skipped'): see ${S}/nextflow/.nextflow.log. To run without it: TOOLS=${rest%,}" >&2
done < <(cat "$FT" 2>/dev/null)
for t in "${RUNS[@]}"; do
  r=ok; [ "$rc" -eq 0 ] || r="not finished"; [ -z "${BAD[${TOOL_OF[$t]:-none}]:-}" ] || r=failed
  [ "${TOOL_OF[$t]:-}" != stranger ] || [ -z "${BAD[expansion_hunter]:-}" ] || r="skipped (expansion_hunter failed)"
  printf 'step\t%s\t%s\n' "$t" "$r" >> "$STATUS"
done
# The reports are written after a failed pipeline too: they show what finished and mark older results stale.
[ "$rc" -eq 0 ] || { echo "ERROR: the pipeline failed (exit ${rc}): see above and ${S}/nextflow/.nextflow.log. Fix the cause and run the same command; -resume reruns only what did not finish. Writing the reports of what finished." >&2; OPTIN=(); }
for sc in "${OPTIN[@]}" 24-html-report.sh generate-report.sh; do
  [ "$sc" != 24-html-report.sh ] || GENOME_DIR="$G" bash "${PGP_ROOT}/bin/write_manifest.sh" "$SAMPLE" run-all.sh "$SEX" || true
  r=ok; echo "[${sc%.sh}] log: ${LOG_DIR}/${sc%.sh}.log"; bash "${SCRIPT_DIR}/${sc}" "$SAMPLE" > "${LOG_DIR}/${sc%.sh}.log" 2>&1 || { r=failed failed=1; echo "  FAILED"; }
  if [[ "$sc" == [0-9]*-*.sh && "$sc" != 24-* ]]; then printf 'step\t%s\t%s\n' "${sc%%-*}" "$r" >> "$STATUS"; fi
done
[ "$rc" -eq 0 ] || exit "$rc"
echo "Done: results in ${S}/, the HTML report ${S}/${SAMPLE}_report.html, the text report ${S}/${SAMPLE}_report.txt"
exit "$failed"
