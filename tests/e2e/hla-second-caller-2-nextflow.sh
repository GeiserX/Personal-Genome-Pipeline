#!/usr/bin/env bash
# The Nextflow side of step 08's second HLA caller: --hla_qc on a VCF+BAM row
# of HG002 with hla_typing and pharmcat. HLA_QC_READS and HLA_QC run for real;
# the concordance table they publish is the one the bash step wrote for the
# same BAM (case hla-second-caller-1), and PGX_CONSENSUS names both callers
# for what it passes to PharmCAT, which reports HLA-A and HLA-B as outside
# calls.
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

echo "+ nextflow run main.nf -profile docker (hla_qc)"
(
  cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
    -work-dir "$WORK" \
    --input "$SHEET" \
    --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" \
    --tools hla_typing,pharmcat \
    --hla_qc true \
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
rm -rf "${NF_OUT:?}" "${WORK:?}"

finish
