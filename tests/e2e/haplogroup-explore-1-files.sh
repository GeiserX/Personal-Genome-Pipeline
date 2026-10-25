#!/usr/bin/env bash
# Steps 12 and 37 write the files the services that explore a haplogroup
# read (docs/12-mito-haplogroup.md and docs/37-y-haplogroup.md, "Explore your
# branch further"), and both reports link the called branch:
#   1. step 12 on HG002 (step 20's Mutect2 calls, case 36): mito/<s>_chrM.fasta
#      is one record of about 16.6 kb, A/C/G/T/N only, and carries the ALT of a
#      kept call. Control: the same check on the reference chrM fails.
#   2. step 37 on a male copy of HG002 (as in sv-mito-telomere-steps-3):
#      y_haplogroup/<s>_chrY.bam (+ .bai) holds chrY reads and no other contig,
#      and is smaller than the BAM it came from. Control: the same contig
#      check on the full BAM fails. A female copy gets no chrY BAM.
#   3. the text report of HG002 links its mtDNA branch and no Y branch (HG002's
#      Y step is skipped on the slices); the male copy's report links the Y
#      branch when Yleaf called one, and none for 'insufficient markers'.
#   4. the Nextflow processes on the same inputs (mito_haplogroup,
#      y_haplogroup) publish the same two files.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
M="${G}/${SAMPLE}/mito"
REF="${G}/reference/GRCh38_no_alt_analysis_set.fasta"

# --- 1. the mtDNA FASTA ------------------------------------------------------------------
check "step 20's Mutect2 calls are there (case 36)" test -s "${M}/${SAMPLE}_chrM_filtered.vcf.gz"
run_step 12-mito-haplogroup.sh "$SAMPLE"
check_step_exit 12-mito-haplogroup.sh
FA="${M}/${SAMPLE}_chrM.fasta"
check "the log names the FASTA" has "^mtDNA sequence for YFull or FamilyTreeDNA Discover: .*${SAMPLE}_chrM.fasta" "$(cat "$STEP_LOG")"
check_eq "FASTA records" "$(grep -c '^>' "$FA" 2>/dev/null || true)" 1
check_eq "the record is named" "$(head -n 1 "$FA" 2>/dev/null)" ">${SAMPLE}_chrM"
SEQLEN=$(grep -v '^>' "$FA" 2>/dev/null | tr -d '\n' | wc -c | tr -d ' ')
check_ge "sequence length (rCRS is 16569)" "$SEQLEN" 16500
check "sequence length under 16700 (${SEQLEN})" test "${SEQLEN:-99999}" -lt 16700
check_eq "characters other than A, C, G, T, N" "$(grep -v '^>' "$FA" 2>/dev/null | tr -d '\nACGTN' | wc -c | tr -d ' ')" 0
check "the FASTA is small (under 20 kB)" test "$(wc -c < "$FA" 2>/dev/null || echo 999999)" -lt 20000

# The first kept SNV before any kept indel: its position in the FASTA is its
# position on chrM, so the FASTA base there must be the ALT.
bcf view -H -i 'FMT/AF>=0.5' "${SAMPLE}/mito/${SAMPLE}_chrM_for_haplogroup.vcf.gz" 2>/dev/null \
  | awk -F'\t' '{ if (length($4) != 1 || length($5) != 1) exit; print $2 "\t" $4 "\t" $5 }' > "${CASE_TMP}/snvs.tsv"
sam faidx "reference/${REF##*/}" chrM 2>/dev/null > "${CASE_TMP}/ref_chrM.fa"
echo "kept SNVs before the first kept indel: $(grep -c . "${CASE_TMP}/snvs.tsv" || true)"
check_ge "kept SNVs to check" "$(grep -c . "${CASE_TMP}/snvs.tsv" || true)" 1
alt_at() {  # alt_at FASTA [SNVS]: OK when every listed SNV's ALT is the base at its position
  python3 - "$1" "${2:-${CASE_TMP}/snvs.tsv}" <<'PY'
import sys
seq = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">"))
bad = [l for l in open(sys.argv[2]) if l.strip() and seq[int(l.split("\t")[0]) - 1] != l.split("\t")[2].strip()]
print(f"{len(bad)} SNVs whose ALT is not in {sys.argv[1].rsplit('/', 1)[-1]}")
sys.exit(1 if bad or not seq else 0)
PY
}
check "the FASTA carries the ALT of each kept SNV" alt_at "$FA"
check "control: the reference chrM does not" eval '! alt_at "${CASE_TMP}/ref_chrM.fa"'

# --- 2. the chrY BAM -----------------------------------------------------------------------
PED_HEAD=$'#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\tCNchrX\tCNchrY'
copy() {  # copy NAME SEX_CODE: a copy of the sample with an indexcov .ped of that sex
  local t=$1
  mkdir -p "${G}/${t}/aligned" "${G}/${t}/indexcov"
  ln -f "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${G}/${t}/aligned/${t}_sorted.bam"
  ln -f "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" "${G}/${t}/aligned/${t}_sorted.bam.bai"
  printf '%s\n%s\t%s\t-9\t-9\t%s\t-9\t1.0\t1.0\n' "$PED_HEAD" "$t" "$t" "$2" > "${G}/${t}/indexcov/indexcov-indexcov.ped"
}
# other_contigs BAM: reads mapped to a contig other than chrY (idxstats, through the index)
other_contigs() { sam idxstats "$1" 2>/dev/null | awk -F'\t' '$1 != "chrY" { n += $3 } END { print n + 0 }'; }

# Yleaf's marker tables, which the image does not ship (setup.sh --yleaf-data).
check "Yleaf's marker tables install" "${REPO}/scripts/setup.sh" --yleaf-data "$G"
XM="${SAMPLE}xm"
copy "$XM" 1
run_step 37-y-haplogroup.sh "$XM"
check_step_exit 37-y-haplogroup.sh
YB="${XM}/y_haplogroup/${XM}_chrY.bam"
check "the chrY BAM is written" nonempty "$YB"
check "and indexed" nonempty "${YB}.bai"
check "the log names it" has "^Y reads for YFull or a Y-DNA project: .*${XM}_chrY.bam" "$(cat "$STEP_LOG")"
check_ge "reads on chrY" "$(sam idxstats "$YB" 2>/dev/null | awk -F'\t' '$1 == "chrY" { print $3 }')" 1
check_eq "reads on any other contig" "$(other_contigs "$YB")" 0
check_ge "control: the full BAM has reads on other contigs" "$(other_contigs "${XM}/aligned/${XM}_sorted.bam")" 1
check "the chrY BAM is smaller than the full BAM" \
  test "$(wc -c < "${G}/${YB}")" -lt "$(wc -c < "${G}/${XM}/aligned/${XM}_sorted.bam")"
ls -l "${G}/${YB}" "${G}/${XM}/aligned/${XM}_sorted.bam"

XF="${SAMPLE}xf"
copy "$XF" 2
run_step 37-y-haplogroup.sh "$XF"
check_step_exit 37-y-haplogroup.sh
check "a female copy gets no chrY BAM" test ! -e "${G}/${XF}/y_haplogroup/${XF}_chrY.bam"

# --- 3. the links in the reports -------------------------------------------------------------
HG=$(awk -F'\t' 'NR == 2 {gsub(/"/, "", $2); print $2}' "${M}/${SAMPLE}_haplogroup.txt" 2>/dev/null)
HG_Q=$(python3 -c 'import sys, urllib.parse; print(urllib.parse.quote(sys.argv[1], safe=""))' "${HG:-none}")
run_step generate-report.sh "$SAMPLE"
check_step_exit generate-report.sh
TXT=$(cat "${G}/${SAMPLE}/${SAMPLE}_report.txt" 2>/dev/null)
check "HG002's text report links its mtDNA branch on Discover (${HG:-none})" \
  has "^  FamilyTreeDNA Discover \(free, no upload\): https://discover\.familytreedna\.com/mtdna/${HG_Q}\$" "$TXT"
check "and on YFull's MTree" has "^  YFull MTree: https://www\.yfull\.com/mtree/${HG_Q}/\$" "$TXT"
check "and no Y branch (its Y step was skipped)" lacks 'discover\.familytreedna\.com/y-dna/|yfull\.com/tree/' "$TXT"

YT="${G}/${XM}/y_haplogroup/${XM}_y_haplogroup.txt"
check_eq "the male copy has Yleaf's table (a header and one row)" "$(grep -c . "$YT" 2>/dev/null || true)" 2
YHG=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next } NR == 2 { print $c["Hg"] }' "$YT" 2>/dev/null)
run_step generate-report.sh "$XM"
check_step_exit generate-report.sh
TXT=$(cat "${G}/${XM}/${XM}_report.txt" 2>/dev/null)
if [ -n "$YHG" ] && [ "$YHG" != NA ]; then
  B=${YHG%%\**}
  check "the male copy's report links its Y branch (${B})" \
    has "^  YFull YTree: https://www\.yfull\.com/tree/$(python3 -c 'import sys, urllib.parse; print(urllib.parse.quote(sys.argv[1], safe=""))' "$B")/\$" "$TXT"
else
  check "the male copy has no Y call (${YHG:-none}), so no Y link" lacks 'yfull\.com/tree/|y-dna/' "$TXT"
  check "and its report still says insufficient markers" has 'Haplogroup: insufficient markers' "$TXT"
fi

# --- 4. the Nextflow processes -----------------------------------------------------------------
if command -v nextflow >/dev/null; then
  NF_OUT="${G}/nf-explore"
  SHEET="${CASE_TMP}/samplesheet.csv"
  printf 'sample,sex,vcf,vcf_index,bam,bam_index\n%s,male,%s,%s,%s,%s\n' "$SAMPLE" \
    "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" \
    "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" > "$SHEET"
  echo "+ nextflow run main.nf -profile docker (mito_haplogroup, y_haplogroup)"
  (
    cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
      -work-dir "${E2E_WORK}/nf-work-explore" \
      --input "$SHEET" --reference "$REF" \
      --tools mito_haplogroup,y_haplogroup \
      --yleaf_data "${G}/reference/yleaf-${YLEAF_DATA_VERSION}/data" \
      --sex_check warn --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB
  ) 2>&1 | tee "$STEP_LOG"
  check_eq "nextflow run exits 0" "${PIPESTATUS[0]}" 0
  cp "${CASE_TMP}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
  NFA="${NF_OUT}/${SAMPLE}/mito/${SAMPLE}_chrM.fasta"
  check_eq "MITO_EXTRACT_CHRM publishes one FASTA record" "$(grep -c '^>' "$NFA" 2>/dev/null || true)" 1
  check_ge "of mtDNA length" "$(grep -v '^>' "$NFA" 2>/dev/null | tr -d '\n' | wc -c | tr -d ' ')" 16500
  # Without mito_variants its input is the VCF's chrM records: the kept
  # calls are the homozygous ones.
  bcf view -H -i 'GT="AA"' -r chrM "${SAMPLE}/vcf/${SAMPLE}.vcf.gz" 2>/dev/null \
    | awk -F'\t' '{ if (length($4) != 1 || length($5) != 1) exit; print $2 "\t" $4 "\t" $5 }' > "${CASE_TMP}/snvs_vcf.tsv"
  check_ge "homozygous chrM SNVs to check" "$(grep -c . "${CASE_TMP}/snvs_vcf.tsv" || true)" 1
  check "carrying the ALT of each homozygous SNV" alt_at "$NFA" "${CASE_TMP}/snvs_vcf.tsv"
  check "control: the reference chrM does not" eval '! alt_at "${CASE_TMP}/ref_chrM.fa" "${CASE_TMP}/snvs_vcf.tsv"'
  NFY="nf-explore/${SAMPLE}/y_haplogroup/${SAMPLE}_chrY.bam"
  check "Y_PILEUP publishes the chrY BAM" nonempty "$NFY"
  check "and its index" nonempty "${NFY}.bai"
  check_eq "reads on chrY, the same as step 37's" \
    "$(sam idxstats "$NFY" 2>/dev/null | awk -F'\t' '$1 == "chrY" { print $3 }')" \
    "$(sam idxstats "$YB" 2>/dev/null | awk -F'\t' '$1 == "chrY" { print $3 }')"
  check_eq "reads on any other contig (Nextflow)" "$(other_contigs "$NFY")" 0
else
  fail "nextflow is not on PATH"
fi

finish
