#!/usr/bin/env bash
# The Nextflow side of step 08's second HLA caller: --hla_qc on a VCF+BAM row
# of HG002 with hla_typing and pharmcat. HLA_QC_READS and HLA_QC run for real;
# the concordance table they publish is the one the bash step wrote for the
# same BAM (case hla-second-caller-1), and PGX_CONSENSUS names both callers
# for what it passes to PharmCAT, which reports HLA-A and HLA-B as outside
# calls. With html_report and --kir in the same run, the Nextflow report
# prints the OptiType agreement per gene (the HLA section of its summary.json
# equals the bash step's, case hla-second-caller-1) and the KIR result, as
# step 08's KIR file says (on the fixture: too few reads).
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

G="$GENOME_DIR"
NF_OUT="${G}/nf-hla-qc"
WORK="${E2E_WORK}/nf-work-hla-qc"
SHEET="${CASE_TMP}/samplesheet.csv"
printf 'sample,sex,vcf,vcf_index,bam,bam_index\n%s,male,%s,%s,%s,%s\n' "$SAMPLE" \
  "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" \
  "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" \
  > "$SHEET"
HLA=$(bash -c '. "$1/scripts/lib/common.sh" && install_data_file hla_dat >/dev/null && data_file hla_dat' _ "$REPO")
GENES=$(bash -c '. "$1/scripts/lib/common.sh" && install_data_file gencode_genes >/dev/null && data_file gencode_genes' _ "$REPO")
check "hla.dat is installed" test -s "$HLA"
check "the GENCODE gene lines are installed" test -s "$GENES"
"${REPO}/scripts/setup.sh" --kir-data "$G" >/dev/null
KIR_DAT="${G}/kir/IPD-KIR_${KIR_DB_RELEASE}/kir.dat"
check "IPD-KIR is installed" test -s "$KIR_DAT"

echo "+ nextflow run main.nf -profile docker (hla_qc)"
(
  cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
    -work-dir "$WORK" \
    --input "$SHEET" \
    --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" \
    --tools hla_typing,pharmcat,html_report \
    --hla_qc true \
    --kir true --kir_dat "$KIR_DAT" \
    --hla_dat "$HLA" \
    --hla_genes "$GENES" \
    --sex_check warn \
    --outdir "$NF_OUT" \
    --max_cpus 4 --max_memory 14.GB
) 2>&1 | tee "$STEP_LOG"
NF_RC=${PIPESTATUS[0]}
check_eq "nextflow run exits 0" "$NF_RC" 0
cp "${CASE_TMP}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
[ "$NF_RC" -eq 0 ] || grep -E 'Error executing process|Command exit status|Command error' -A 12 "${CASE_TMP}/.nextflow.log" | head -80

R="${NF_OUT}/${SAMPLE}"
OT="${R}/hla_optitype/${SAMPLE}_optitype_result.tsv"
CONC="${R}/hla_optitype/${SAMPLE}_hla_concordance.tsv"
cat "$OT" "$CONC" 2>/dev/null
check "HLA_QC published OptiType's result" test -s "$OT"
check "HLA_QC published the concordance table" test -s "$CONC"
check "the Nextflow concordance table is the bash step's (case hla-second-caller-1)" \
  cmp -s "$CONC" "${E2E_WORK}/hla-qc-bash-concordance.tsv"
diff "${E2E_WORK}/hla-qc-bash-concordance.tsv" "$CONC" 2>/dev/null || true
CONS="${R}/pgx_consensus/${SAMPLE}_pgx_consensus.tsv"
cat "$CONS" 2>/dev/null
for g in HLA-A HLA-B; do
  check "PGX_CONSENSUS passes ${g}, naming both callers" \
    grep -Eq "^${g}"$'\t[^\t]+\tyes\tT1K and OptiType \\(step 08\\), which agree at two fields\t' "$CONS"
done
source_of() {
  python3 -c 'import json, sys; g = json.load(open(sys.argv[1]))["genes"][sys.argv[2]]; print(g.get("callSource", ""))' "$1" "$2" 2>/dev/null
}
for g in HLA-A HLA-B; do
  check_eq "PHARMCAT reports ${g} as an outside call" "$(source_of "${R}/pharmcat/${SAMPLE}.report.json" "$g")" OUTSIDE
done

# The Nextflow report: the OptiType agreement per gene and the KIR result.
HTML="${R}/${SAMPLE}_report.html"
KG="${R}/kir/${SAMPLE}_kir_genotype.tsv"
check "HTML_REPORT wrote the report" test -s "$HTML"
check "KIR_TYPING published the KIR genotype file" test -s "$KG"
hla_of() {
  python3 -c 'import json, sys; d = json.load(open(sys.argv[1]))["sections"]["hla"]["data"]; print(json.dumps(d.get(sys.argv[2]), sort_keys=True))' "$1" "$2" 2>/dev/null
}
check "the Nextflow summary has the OptiType rows" has '"agree": "yes"' "$(hla_of "${R}/summary.json" optitype)"
check_eq "and they are the bash report's (case hla-second-caller-1)" "$(hla_of "${R}/summary.json" optitype)" \
  "$(hla_of "${E2E_WORK}/hla-qc-bash-summary.json" optitype)"
# nf_report_has HTML: every concordance row and the KIR result are in the report.
nf_report_has() {
  local g t o a w n=0 ok=0
  while IFS=$'\t' read -r g t o a _; do
    case $a in yes) w=agree ;; no) w=disagree ;; *) w=$a ;; esac
    n=$((n + 1))
    grep -qF "${g}: T1K vs OptiType" "$1" && grep -qF "${w}: T1K ${t}, OptiType ${o}" "$1" && ok=$((ok + 1))
  done < <(tail -n +2 "$CONC")
  [ "$n" -eq 3 ] && [ "$ok" -eq 3 ] || return 1
  if grep -q $'^KIR' "$KG"; then
    while IFS=$'\t' read -r g _; do grep -qF "${g}*" "$1" || return 1; done < <(awk -F'\t' '/^KIR/ && $3 != "."' "$KG")
  else
    grep -qF "$(sed -n '1s/^# *//p' "$KG")" "$1" || return 1
  fi
  grep -qF "No step interprets them" "$1"
}
not_nf_report_has() { ! nf_report_has "$@"; }
grep -A 3 'KIR genes\|OptiType check\|T1K vs OptiType' "$HTML" | head -n 30
check "the Nextflow report prints each gene's agreement, both calls, and the KIR result" nf_report_has "$HTML"
# Red control: the same summary through a report without the OptiType and KIR blocks.
BR="${CASE_TMP}/bin_no_blocks"
rm -rf "$BR"
cp -R "${REPO}/bin" "$BR"
sed -i 's/if hla.get("optitype") is not None:/if False:/; s/if hla.get("kir") is not None:/if False:/' "${BR}/render_report.py"
check "the control's mutations applied" test "$(grep -c '^    if False:$' "${BR}/render_report.py")" -eq 2
python3 -I "${BR}/render_report.py" --summary "${R}/summary.json" -o "${CASE_TMP}/ctl.html" >/dev/null
check "control: the check above rejects that report" not_nf_report_has "${CASE_TMP}/ctl.html"
rm -rf "${NF_OUT:?}" "${WORK:?}"

finish
