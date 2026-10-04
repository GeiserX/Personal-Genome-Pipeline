#!/usr/bin/env bash
# The Nextflow pipeline from the fixture's FASTQ pair, real containers: FASTP,
# MINIMAP2_INDEX (the reference is given from a folder without the .sr.mmi
# case 20 built, so the run builds its own), ALIGN_MINIMAP2, ALIGN_MARKDUP,
# INDEXCOV, DEEPVARIANT (male: haploid chrX/chrY outside the PARs), then the
# default tools plus clinvar and hla_typing. Not cyrius: it normalises depth
# over bins on every autosome and stops on the first contig the BAM lacks,
# and a BAM aligned from the fixture reads has only the fixture's slices
# (case 37 gives step 21 a separate BAM for that reason). Sample HG002P, the name case
# nextflow-from-fastq-1-bash ran the bash steps under; the E2E workflow
# compares the two with scripts/ci/parity-diff.sh.
#
# --sex_check warn: on the fixture's slices indexcov reads HG002 (male) as
# female (CNchrX and CNchrY near 2, see case 34), so the declared sex would
# stop the run. Case nextflow-from-fastq-3-sexcheck checks that stop.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

P=HG002P
G="$GENOME_DIR"
NF_OUT="${G}/nf-fastq"
WORK="${E2E_WORK}/nf-work-fastq"
SHEET="${CASE_TMP}/samplesheet.csv"
printf 'sample,fastq_1,fastq_2,sex\n%s,%s,%s,male\n' "$P" \
  "${G}/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz" "${G}/${SAMPLE}/fastq/${SAMPLE}_R2.fastq.gz" > "$SHEET"

# The reference without its minimap2 index beside it (hard links: a symlink
# could point outside the folders Nextflow mounts into the containers).
REFDIR="${E2E_WORK}/nf-fastq-ref"
rm -rf "$REFDIR"
mkdir -p "$REFDIR"
for e in fasta fasta.fai dict; do
  ln -f "${G}/reference/GRCh38_no_alt_analysis_set.${e}" "${REFDIR}/GRCh38_no_alt_analysis_set.${e}"
done

# The HLA data case bash-step-08 installed (installed here when it did not)
HLA=$(bash -c '. "$1/scripts/lib/common.sh" && install_data_file hla_dat >/dev/null && data_file hla_dat' _ "$REPO")
GENES=$(bash -c '. "$1/scripts/lib/common.sh" && install_data_file gencode_genes >/dev/null && data_file gencode_genes' _ "$REPO")
check "hla.dat is installed" test -s "$HLA"
check "the GENCODE gene lines are installed" test -s "$GENES"
PGS="${E2E_WORK}/parity-pgs"
check "case nextflow-from-fastq-1-bash wrote the score files" test -s "${PGS}/PGS000018.txt.gz"
INTERVALS=$(awk '{printf "%s%s:%d-%d", (NR > 1 ? " " : ""), $1, $2 + 1, $3}' "${FIXTURE_DIR}/regions.bed")

# Keep going after a failed task, so one run shows every broken module.
cat > "${CASE_TMP}/e2e.config" <<'NFCONF'
process.errorStrategy = 'ignore'
NFCONF

TOOLS=pharmcat,cpic,vcfanno,roh,prs,mito_haplogroup,telomere_hunter,mosdepth,mito_variants,html_report,multiqc,clinvar,hla_typing
echo "+ nextflow run main.nf -profile docker (from FASTQ)"
(
  cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
    -c "${CASE_TMP}/e2e.config" \
    -work-dir "$WORK" \
    --input "$SHEET" \
    --reference "${REFDIR}/GRCh38_no_alt_analysis_set.fasta" \
    --intervals "$INTERVALS" \
    --sex_check warn \
    --tools "$TOOLS" \
    --clinvar "${G}/clinvar/clinvar_pathogenic_chr.vcf.gz" \
    --clinvar_index "${G}/clinvar/clinvar_pathogenic_chr.vcf.gz.tbi" \
    --revel "${G}/annotations/revel_grch38.tsv.gz" \
    --revel_index "${G}/annotations/revel_grch38.tsv.gz.tbi" \
    --pgs_scoring "$PGS" \
    --hla_dat "$HLA" \
    --hla_genes "$GENES" \
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
# tasks NAME: how many tasks of process NAME completed
tasks() {
  awk -F'\t' -v p="$1" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    $c["status"] == "COMPLETED" { n = $c["name"]; sub(/ \(.*/, "", n); sub(/.*:/, "", n); if (n == p) k++ }
    END {print k + 0}' "${TRACE:-/dev/null}" 2>/dev/null
}
for p in FASTP MINIMAP2_INDEX ALIGN_MINIMAP2 ALIGN_MARKDUP INDEXCOV DEEPVARIANT T1K_BUILD HLA_TYPING; do
  check_eq "${p} tasks completed" "$(tasks "$p")" 1
done

R="nf-fastq/${P}"
# --- alignment ---------------------------------------------------------------------
BAM="${R}/aligned/${P}_sorted.bam"
check "BAM passes samtools quickcheck" sam quickcheck -v "$BAM"
check "BAM index exists" nonempty "${BAM}.bai"
HDR=$(sam view -H "$BAM" 2>/dev/null)
TAB=$'\t'
check "@RG SM equals the sample name (${P})" has "^@RG.*${TAB}SM:${P}(${TAB}|\$)" "$HDR"
check "the BAM header records samtools markdup" has '^@PG.*ID:samtools.*markdup' "$HDR"
check_ge "reads flagged as duplicates" "$(sam flagstat "$BAM" 2>/dev/null | awk '/ duplicates$/ {print $1; exit}')" 1
check "fastp's report is published" nonempty "${R}/fastq_trimmed/${P}_fastp.json"
check_eq "trimmed reads published (they stay in the work directory)" \
  "$(find "${G}/${R}" -name '*.fastq.gz' 2>/dev/null | wc -l | tr -d ' ')" 0

# --- sex check ---------------------------------------------------------------------------
check "INDEXCOV wrote the sex check" nonempty "${R}/${P}_sex_check.tsv"
check "the run warns that declared and inferred sex differ, with both" \
  has "says sex male, but indexcov infers (female|unknown) from the BAM index \\(CNchrX=" "$(cat "$STEP_LOG")"

# --- calling -------------------------------------------------------------------------------
VCF="${R}/vcf/${P}.vcf.gz"
GVCF="${R}/vcf/${P}.g.vcf.gz"
check "VCF is readable" vcf_ok "$VCF"
check "VCF index exists" nonempty "${VCF}.tbi"
check "gVCF is readable" vcf_ok "$GVCF"
check_ge "gVCF reference blocks (records with END)" "$(vcf_count -i 'INFO/END>0' "$GVCF")" 10
check_eq "VCF sample column" "$(bcf query -l "$VCF" 2>/dev/null)" "$P"
check_ge "PASS records on the chr20 slice" "$(vcf_count -f PASS -r chr20:10000000-10500000 "$VCF")" 200
GT=$(bcf query -r "$(planted chrom):$(planted pos)" -f '[%GT]\n' "$VCF" 2>/dev/null | awk 'NR == 1')
check "planted SNV is called non-reference (GT ${GT:-none})" has '1' "${GT:-}"
NONPAR="chrX:73700001-74000000,chrY:2700001-3000000"
check_ge "calls with an ALT allele on the non-PAR slices" "$(vcf_count -i 'GT="alt"' -r "$NONPAR" "$VCF")" 5
check_eq "heterozygous calls on non-PAR chrX and chrY (male: haploid)" "$(vcf_count -g het -r "$NONPAR" "$VCF")" 0
check_ge "heterozygous calls on the PAR1 slice (diploid there)" "$(vcf_count -g het -r chrX:1000001-1200000 "$VCF")" 1

# --- the tools on the called VCF and the aligned BAM ------------------------------
GENE=$(planted gene)
check_ge "CLINVAR_SCREEN names the planted hit's gene (${GENE})" "$(sample_side_hits "${G}/${R}/clinvar" "$GENE")" 1
JSON="${G}/${R}/pharmcat/${P}.report.json"
check_ge "PHARMCAT genes called" "$( [ -s "$JSON" ] && pharmcat_called "$JSON" | grep -c . || echo 0)" 1
check_eq "PRS scores in the summary" "$(awk 'NR > 1' "${G}/${R}/prs/${P}_prs_summary.tsv" 2>/dev/null | wc -l | tr -d ' ')" 9
check_eq "PRS input (hom-ref sites from the gVCF)" "$(awk -F'\t' 'NR == 2 {print $NF}' "${G}/${R}/prs/${P}_prs_summary.tsv" 2>/dev/null)" gvcf
check "MOSDEPTH summary has a chr20 row" grep -q '^chr20' "${G}/${R}/coverage/${P}.mosdepth.summary.txt"
check_ge "ROH output lines" "$(cat "${G}/${R}"/roh/*_roh.txt 2>/dev/null | grep -vc '^#' || true)" 1
GT_HLA="${G}/${R}/hla/${P}_hla_genotype.tsv"
head -n 6 "$GT_HLA" 2>/dev/null
for g in HLA-A HLA-B HLA-C; do
  check_ge "${g} rows with a called allele" "$(awk -F'\t' -v g="$g" '$1 == g && $3 ~ /^HLA-/' "$GT_HLA" 2>/dev/null | wc -l | tr -d ' ')" 1
done
check "the HTML report is written" nonempty "${R}/${P}_report.html"
check "MultiQC report is written" nonempty "nf-fastq/multiqc/multiqc_report.html"

finish
