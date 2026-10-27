#!/usr/bin/env bash
# Step 52 (MUC1 VNTR frameshifts, VNtyper) on real reads, and its Nextflow
# mirror:
#   - the fixture's own BAM (case 20's, a slice of HG002 with chr1 but no read
#     over MUC1): VNtyper has nothing to read, and the step says "not
#     assessed" with that reason, never "negative";
#   - two BAMs of VNtyper's own public test data (tests/smoke/
#     vntyper_testdata.py reads them out of its Zenodo archive, sha256
#     pinned): their reads over MUC1 carried onto the fixture's header and
#     added to case 20's BAM, as a whole-genome BAM would hold them. The
#     sample with a frameshift reads positive, with the variant, and the one
#     without reads negative (the red control); both give the depth over the
#     VNTR and the sensitivity line of its band;
#   - nextflow run main.nf --tools vntyper,html_report on the positive sample
#     publishes the bash step's table byte for byte, and the HTML report has
#     the card with the result.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR" T="$CASE_TMP"
BASE="${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
MUC1="chr1:155184000-155194000"
dk() { docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$1" "${@:2}"; }
retry() { local i; for i in 1 2 3; do "$@" && return 0; [ "$i" -lt 3 ] && { echo "  try ${i} of 3 failed; again in $((i * 30)) s"; sleep $((i * 30)); }; done; return 1; }
# metric NAME KEY: one value of the sample's table.
metric() { awk -F'\t' -v k="$2" '$1 == k {print $2}' "${G}/$1/vntyper/$1_vntyper.tsv" 2>/dev/null; }
step() { run_step 52-vntyper.sh "$1"; check_step_exit 52-vntyper.sh; }
sample() {
  rm -rf "${G:?}/$1"; mkdir -p "${G}/$1/aligned"
  cp "$2" "${G}/$1/aligned/$1_sorted.bam"; cp "$2.bai" "${G}/$1/aligned/$1_sorted.bam.bai"
}

# --- the fixture: chr1 is there, MUC1 has no read ------------------------------------
check "case 20's BAM is there" test -s "${BASE}.bai"
step "$SAMPLE"
cat "${G}/${SAMPLE}/vntyper/${SAMPLE}_vntyper.tsv" 2>/dev/null
check_eq "the fixture: not assessed" "$(metric "$SAMPLE" status)" "not assessed"
check "the fixture: the message says the BAM has no read over MUC1" has 'no read in the MUC1 region' "$(metric "$SAMPLE" message)"

# --- VNtyper's positive and negative samples on the fixture's header -------------------
check "VNtyper's two test BAMs download and match their sha256" \
  retry python3 -I "${REPO}/tests/smoke/vntyper_testdata.py" "$T" positive.bam negative.bam
H="${SAMPLE}muc"
for s in positive negative; do
  dk "$SAMTOOLS_IMAGE" samtools index "${T}/${s}.bam"
  # The reads over the region VNtyper slices (mates elsewhere dropped), under
  # case 20's header: the pipeline checks a BAM's header against the
  # reference, and the fixture's chr1 is GRCh38's.
  { dk "$SAMTOOLS_IMAGE" samtools view -H "$BASE" | grep -v '^@PG'
    dk "$SAMTOOLS_IMAGE" samtools view -H "${T}/${s}.bam" | grep '^@RG'
    dk "$SAMTOOLS_IMAGE" samtools view "${T}/${s}.bam" "$MUC1" | awk -F'\t' '$7 == "=" || $7 == "*"'
  } > "${T}/${s}.sam"
  dk "$SAMTOOLS_IMAGE" samtools sort -o "${T}/${s}_muc1.bam" "${T}/${s}.sam"
  check_ge "${s}: reads over MUC1" "$(dk "$SAMTOOLS_IMAGE" samtools view -c "${T}/${s}_muc1.bam" 2>/dev/null || echo 0)" 5000
  dk "$SAMTOOLS_IMAGE" samtools merge -f -o "${T}/${s}_all.bam" "$BASE" "${T}/${s}_muc1.bam"
  dk "$SAMTOOLS_IMAGE" samtools index "${T}/${s}_all.bam"
  sample "${H}${s}" "${T}/${s}_all.bam"
  step "${H}${s}"
  cat "${G}/${H}${s}/vntyper/${H}${s}_vntyper.tsv" 2>/dev/null
  check "${s}: the sensitivity line gives the depth over the VNTR" \
    has '^Mean depth over the MUC1 VNTR: [0-9]+x\. At that depth .* VNtyper found the known variant in [0-9]+ of [0-9]+ runs' \
    "$(metric "${H}${s}" sensitivity)"
  check "${s}: VNtyper's coverage check passes (100x or more)" has '^PASS$' "$(metric "${H}${s}" coverage_qc)"
  echo "vntyper: ${s} test sample on the fixture: $(metric "${H}${s}" status), $(metric "${H}${s}" variant)," \
    "Kestrel $(metric "${H}${s}" confidence), $(metric "${H}${s}" vntr_mean_depth)x over the VNTR" >> "$E2E_NOTES"
done
P="${H}positive" N="${H}negative"
check_eq "positive: VNtyper calls the frameshift" "$(metric "$P" status)" positive
check "positive: with high precision" has '^High_Precision' "$(metric "$P" confidence)"
check "positive: and a variant name" has '[0-9]+_?[0-9]*(del|dup|ins)' "$(metric "$P" variant)"
check "positive: the 100x-or-more band (160 of 161 runs)" has '160 of 161 runs' "$(metric "$P" sensitivity)"
check_eq "negative (red control): no call" "$(metric "$N" status)" negative
check_eq "negative: Kestrel says Negative" "$(metric "$N" confidence)" Negative
check_eq "negative: no variant" "$(metric "$N" variant)" "."
check "VNtyper's own report is kept" test -s "${G}/${P}/vntyper/run/summary_report.html"

# --- Nextflow against the bash step ------------------------------------------------------
command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
check "case 21's VCF is there" test -s "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi"
mkdir -p "${G}/${P}/vcf"
cp "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${G}/${P}/vcf/${P}.vcf.gz"
cp "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" "${G}/${P}/vcf/${P}.vcf.gz.tbi"
SHEET="${T}/samplesheet.csv"
printf 'sample,sex,vcf,vcf_index,bam,bam_index\n%s,male,%s,%s,%s,%s\n' "$P" \
  "${G}/${P}/vcf/${P}.vcf.gz" "${G}/${P}/vcf/${P}.vcf.gz.tbi" \
  "${G}/${P}/aligned/${P}_sorted.bam" "${G}/${P}/aligned/${P}_sorted.bam.bai" > "$SHEET"
NF_OUT="${G}/nf-vntyper"
echo "+ nextflow run main.nf -profile docker --tools vntyper,html_report"
( cd "$T" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -work-dir "${E2E_WORK}/nf-work-vntyper" \
    --input "$SHEET" --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" --tools vntyper,html_report \
    --sex_check warn --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB ) 2>&1 | tee "$STEP_LOG"
NF_RC=${PIPESTATUS[0]}
check_eq "nextflow run exits 0" "$NF_RC" 0
cp "${T}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
a="${G}/${P}/vntyper/${P}_vntyper.tsv" b="${NF_OUT}/${P}/vntyper/${P}_vntyper.tsv"
check "Nextflow wrote ${P}_vntyper.tsv" test -s "$b"
check "Nextflow's ${P}_vntyper.tsv is the bash step's" cmp -s "$a" "$b"
cmp -s "$a" "$b" || diff "$a" "$b" | head -20
check "control: the comparison sees a difference" bash -c '! cmp -s "$1" "$2"' _ \
  "${G}/${N}/vntyper/${N}_vntyper.tsv" "$b"
HTML=$(find "$NF_OUT" -name "${P}_report.html" | head -n 1)
check "the HTML report has the MUC1 VNTR card" has "<h2>MUC1 VNTR Frameshifts \(VNtyper\)</h2>" "$(cat "$HTML" 2>/dev/null)"
check "and the card shows the positive result" has 'badge-red">positive<' "$(cat "$HTML" 2>/dev/null)"

rm -rf "${G:?}/${N}" "${T:?}"/*.bam "${T:?}"/*.bai "${T:?}"/*.sam
finish
