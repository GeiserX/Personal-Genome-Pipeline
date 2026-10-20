#!/usr/bin/env bash
# The Nextflow modules CLINICAL_FILTER, SLIVAR_PRIORITIZE and SLIVAR on the
# MANEQ sample of case 1 (clinical-filter-mane-and-quality-1-tiers), with the
# pinned bcftools and slivar images. The pipeline only runs these modules after
# VEP, whose cache is too large here, so a small workflow written below
# includes the three modules and feeds them the VEP-annotated VCF directly.
# bin/ is copied next to it: Nextflow puts the project's bin/ on the task PATH,
# where the modules find bin/csq_rank.awk. The checks are the ones case 1 makes
# on steps 23 and 31, so script and module are held to the same records.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

S=MANEQ
SRC="${GENOME_DIR}/${S}/vep/${S}_vep.vcf"
if [ ! -s "$SRC" ]; then
  fail "${SRC} is missing: case clinical-filter-mane-and-quality-1-tiers writes it"
  finish
fi
O=nf-cfmq
rm -rf "${GENOME_DIR:?}/${O}"
mkdir -p "${GENOME_DIR}/${O}"
bcf view -Oz -o "${O}/${S}_vep.vcf.gz" "${S}/vep/${S}_vep.vcf"
bcf index -t "${O}/${S}_vep.vcf.gz"
check "the MANEQ VCF is compressed and indexed" nonempty "${O}/${S}_vep.vcf.gz.tbi"

P="${CASE_TMP}/project"
rm -rf "$P"
mkdir -p "$P"
cp -R "${REPO}/bin" "${P}/bin"
cat > "${P}/main.nf" <<NF
include { CLINICAL_FILTER } from '${REPO}/modules/local/clinical_filter/main'
include { SLIVAR_PRIORITIZE; SLIVAR } from '${REPO}/modules/local/slivar/main'

workflow {
    ch = channel.of([[id: '${S}'], file('${GENOME_DIR}/${O}/${S}_vep.vcf.gz'), file('${GENOME_DIR}/${O}/${S}_vep.vcf.gz.tbi')])
    CLINICAL_FILTER(ch)
    SLIVAR_PRIORITIZE(ch, [])
    SLIVAR(SLIVAR_PRIORITIZE.out.vcf.join(SLIVAR_PRIORITIZE.out.ped))
}
NF
cat > "${P}/nextflow.config" <<NFCONF
params.outdir = '${GENOME_DIR}/${O}/results'
params.publish_dir_mode = 'copy'
docker.enabled = true
docker.runOptions = '--user $(id -u):$(id -g) --network none'
NFCONF

echo "+ nextflow run (CLINICAL_FILTER, SLIVAR_PRIORITIZE, SLIVAR on ${S})"
(
  cd "$CASE_TMP" && nextflow run "${P}/main.nf" -ansi-log false \
    -c "${REPO}/conf/containers.config" -work-dir "${CASE_TMP}/work"
) 2>&1 | tee "$STEP_LOG"
NF_RC=${PIPESTATUS[0]}
check_eq "nextflow run exits 0" "$NF_RC" 0
[ "$NF_RC" = 0 ] || grep -E 'Error executing process|Command exit status|Command error' -A 8 "${CASE_TMP}/.nextflow.log" | head -60

positions() { bcf query -f '%POS\n' "$1" 2>/dev/null | paste -sd, -; }
cell() {
  awk -F'\t' -v p="$2" -v k="$3" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    $c["POS"] == p && (k in c) {print $c[k]; exit}' "$1" 2>/dev/null
}

# --- CLINICAL_FILTER -----------------------------------------------------------
C="${O}/results/${S}/clinical/${S}"
SUM="${GENOME_DIR}/${C}_clinical_summary.tsv"
cat "$SUM" 2>/dev/null
check_eq "CLINICAL_FILTER: the clinical VCF holds every tier, not the common 4000" \
  "$(positions "${C}_clinical.vcf.gz")" 1000,2000,2500,2600,3000,5000,6000,6100
check_eq "CLINICAL_FILTER: the summary has 21 columns" "$(awk -F'\t' 'NR == 1 {print NF}' "$SUM" 2>/dev/null)" 21
check_eq "CLINICAL_FILTER: 1000 (HIGH off MANE) keeps its MANE consequence" \
  "$(cell "$SUM" 1000 IMPACT) $(cell "$SUM" 1000 Consequence) $(cell "$SUM" 1000 TRANSCRIPT)" "MODIFIER intron_variant MANE"
check_eq "CLINICAL_FILTER: and shows the worst one beside it" \
  "$(cell "$SUM" 1000 WORST_IMPACT) $(cell "$SUM" 1000 WORST_Consequence)" "HIGH stop_gained"
check_eq "CLINICAL_FILTER: 5000 is ranked on its canonical transcript" \
  "$(cell "$SUM" 5000 TRANSCRIPT) $(cell "$SUM" 5000 IMPACT)" "CANONICAL MODERATE"
check_eq "CLINICAL_FILTER: 3000's frequency tier" "$(cell "$SUM" 3000 FREQ_TIER)" unknown
check_eq "CLINICAL_FILTER: 2000 carries its GQ and AD" "$(cell "$SUM" 2000 GQ) $(cell "$SUM" 2000 AD)" "12 3,2"
check_eq "CLINICAL_FILTER: the GQ 12 heterozygote is a weak call" \
  "$(cell "$SUM" 2000 CALL_CHECK) $(cell "$SUM" 2000 CALL_WHY)" "weak GQ12"
check_eq "CLINICAL_FILTER: the GQ 11 1/1 on AD 12,5 is weak by GQ and alt reads" \
  "$(cell "$SUM" 2500 CALL_CHECK) $(cell "$SUM" 2500 CALL_WHY)" "weak GQ11,AB0.29"
check_eq "CLINICAL_FILTER: the GQ 50 1/1 on AD 0,30 is not" "$(cell "$SUM" 2600 CALL_CHECK)" ok
check "CLINICAL_FILTER: the clinical VCF names the floor it used" \
  has '^##PGP_CALL_FLOOR=GQ=30,DP=10,HET_AB=0.2-0.8,HOM_AB=0.8$' "$(bcf view -h "${C}_clinical.vcf.gz" 2>/dev/null)"

# --- SLIVAR_PRIORITIZE and SLIVAR ------------------------------------------------
L="${O}/results/${S}/slivar/${S}"
SUM31="${GENOME_DIR}/${L}_slivar_summary.tsv"
cat "$SUM31" 2>/dev/null
check_eq "SLIVAR_PRIORITIZE: prioritized keeps 3000 for the compound-het search, not 4000" \
  "$(positions "${L}_prioritized.vcf.gz")" 1000,2000,2500,2600,3000,5000,6000,6100
check_eq "SLIVAR_PRIORITIZE: the summary has 15 columns" "$(awk -F'\t' 'NR == 1 {print NF}' "$SUM31" 2>/dev/null)" 15
check_eq "SLIVAR_PRIORITIZE: 3000's frequency tier" "$(cell "$SUM31" 3000 FREQ_TIER)" unknown
check_eq "SLIVAR_PRIORITIZE: 6000 carries its gnomAD AF" \
  "$(cell "$SUM31" 6000 MAX_AF) $(cell "$SUM31" 6000 gnomADe_AF) $(cell "$SUM31" 6000 gnomADg_AF) $(cell "$SUM31" 6000 FREQ_TIER)" \
  "0.0002 0.0001 0.0002 rare"
check_eq "SLIVAR_PRIORITIZE: 1000's MANE consequence beside its worst one" \
  "$(cell "$SUM31" 1000 Consequence) $(cell "$SUM31" 1000 MANE_Consequence)" "stop_gained intron_variant"
check_eq "SLIVAR_PRIORITIZE: no MANE consequence for 5000 (no MANE transcript)" \
  "$(cell "$SUM31" 5000 TRANSCRIPT) $(cell "$SUM31" 5000 MANE_Consequence)" "CANONICAL ."
CH_POS=$(positions "${L}_compound_hets.vcf.gz")
echo "SLIVAR compound-het positions: ${CH_POS:-none}"
check_eq "SLIVAR: the GENEH pair is a compound-het candidate" \
  "$(tr ',' '\n' <<< "$CH_POS" | grep -xE '6000|6100' | paste -sd, -)" 6000,6100

finish
