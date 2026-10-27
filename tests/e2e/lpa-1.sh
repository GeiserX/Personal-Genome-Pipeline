#!/usr/bin/env bash
# Step 46 (Lp(a), the LPA KIV-2 repeat count) on real reads, and its Nextflow
# mirror:
#   - the fixture's own BAM (case 20's, a slice of HG002 with chr6 but no
#     reads over LPA): the step writes its table and says "too little depth",
#     never a count;
#   - HG002's Illumina reads over LPA (chr6:160.53-160.67 Mb) from the 35x
#     BWA-MEM GRCh38 no-ALT BAM of DeepVariant's case study, carried onto the
#     fixture's header and added to case 20's BAM: a count with its interval, most of the KIV-2 depth
#     from reads of MAPQ 0, and an interval that holds the count. No count for
#     HG002 is published that this case could compare with, so the value goes
#     to the notes and the checks are on its arithmetic:
#   - the same reads with half of those over the KIV-2 block: N + 2 halves
#     (bin/lpa_kiv2.py: N = 14 x M / D - 2);
#   - the same reads at an eighth of the depth: "too little depth";
#   - nextflow run main.nf --tools lpa,html_report on the real reads publishes
#     the same three files as the bash step, byte for byte, and the HTML
#     report has the Lp(a) card.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR" T="$CASE_TMP"
SRC="https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"
BASE="${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
# dk IMAGE CMD...: a pinned image with the work folder at its own path.
dk() { docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$1" "${@:2}"; }
retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && { echo "  try ${i} of 4 failed; again in $((i * 15)) s"; sleep $((i * 15)); }; done; return 1; }
# metric NAME KEY: one value of the sample's table.
metric() { awk -F'\t' -v k="$2" '$1 == k {print $2}' "${G}/$1/lpa/$1_lpa_kiv2.tsv" 2>/dev/null; }
step() { run_step 46-lpa.sh "$1"; check_step_exit 46-lpa.sh; }
# sample NAME BAM: a sample folder with that BAM.
sample() {
  rm -rf "${G:?}/$1"; mkdir -p "${G}/$1/aligned"
  cp "$2" "${G}/$1/aligned/$1_sorted.bam"; cp "$2.bai" "${G}/$1/aligned/$1_sorted.bam.bai"
}

# --- the fixture: chr6 is there, LPA has no reads ------------------------------------
check "case 20's BAM is there" test -s "${BASE}.bai"
step "$SAMPLE"
cat "${G}/${SAMPLE}/lpa/${SAMPLE}_lpa_kiv2.tsv" 2>/dev/null
check_eq "the fixture: too little depth" "$(metric "$SAMPLE" status)" "too little depth"
check_eq "the fixture: no count" "$(metric "$SAMPLE" kiv2_copies)" "."
check "the fixture: the message gives the depth it found" has '^too little depth: [0-9.]+x over the single-copy LPA exons' "$(metric "$SAMPLE" message)"
check_eq "the fixture: a depth row per region" "$(wc -l < "${G}/${SAMPLE}/lpa/${SAMPLE}_lpa_kiv2_depth.tsv" | tr -d ' ')" 29

# --- HG002's reads over LPA, on the fixture's header ---------------------------------
# SAMTOOLS_IMAGE's libcurl has no CA bundle: lend it the host's (htslib reads CURL_CA_BUNDLE).
fetch_reads() {
  rm -f "${T}/win.bam" "${T}/$(basename "$SRC").bai"
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$SAMTOOLS_IMAGE" \
    samtools view -b -o "${T}/win.bam" "$SRC" chr6:160530000-160670000
}
check "HG002's reads over LPA download" retry fetch_reads
# The BAM's header lists every GRCh38 contig; the pipeline checks a BAM's
# header against the reference, and the fixture's has 14 whole contigs, chr6
# among them. So the reads (mates elsewhere set to unpaired positions dropped)
# are written again under the header of case 20's BAM, with their read groups.
{ dk "$SAMTOOLS_IMAGE" samtools view -H "$BASE" | grep -v -e '^@RG' -e '^@PG'
  dk "$SAMTOOLS_IMAGE" samtools view -H "${T}/win.bam" | grep '^@RG'
  dk "$SAMTOOLS_IMAGE" samtools view "${T}/win.bam" | awk -F'\t' '$7 == "=" || $7 == "*" || $7 == "chr6"'
} > "${T}/lpa.sam"
dk "$SAMTOOLS_IMAGE" samtools sort -o "${T}/lpa.bam" "${T}/lpa.sam"
dk "$SAMTOOLS_IMAGE" samtools index "${T}/lpa.bam"
check_ge "HG002 reads over LPA" "$(dk "$SAMTOOLS_IMAGE" samtools view -c "${T}/lpa.bam" 2>/dev/null || echo 0)" 20000
check_eq "the reads carry the fixture's header (chr6 is its 5th contig)" \
  "$(dk "$SAMTOOLS_IMAGE" samtools view -H "${T}/lpa.bam" | awk -F'\t' '$1 == "@SQ" {n++; if ($2 == "SN:chr6") print n}')" 5

# The sample is case 20's BAM with these reads added: the fixture has no read
# over LPA, so the count is the reads' own, and Nextflow's INDEXCOV below gets
# the whole fixture to read.
dk "$SAMTOOLS_IMAGE" samtools merge -f -o "${T}/hg.bam" "$BASE" "${T}/lpa.bam"
dk "$SAMTOOLS_IMAGE" samtools index "${T}/hg.bam"
H="${SAMPLE}lpa"
sample "$H" "${T}/hg.bam"
step "$H"
cat "${G}/${H}/lpa/${H}_lpa_kiv2.tsv" 2>/dev/null
N=$(metric "$H" kiv2_copies) LO=$(metric "$H" ci95_low) HI=$(metric "$H" ci95_high)
check_eq "HG002: a count" "$(metric "$H" status)" estimated
check "HG002: the count (${N}) lies between 10 and 80" awk -v n="$N" 'BEGIN {exit !(n >= 10 && n <= 80)}'
check "HG002: the interval (${LO} to ${HI}) holds the count and is under 20 copies wide" \
  awk -v n="$N" -v l="$LO" -v h="$HI" 'BEGIN {exit !(l < n && n < h && h - l < 20)}'
check "HG002: most of the KIV-2 depth is from reads of MAPQ 0 ($(metric "$H" mapq0_share))" \
  awk -v s="$(metric "$H" mapq0_share)" 'BEGIN {exit !(s >= 0.5)}'
check "HG002: at least 20x over the single-copy exons ($(metric "$H" flank_depth))" \
  awk -v d="$(metric "$H" flank_depth)" 'BEGIN {exit !(d >= 20)}'
echo "lpa: HG002 (35x BWA-MEM) KIV-2 ${N} (95% ${LO} to ${HI}), class $(metric "$H" class);" \
  "KIV-2 exons $(metric "$H" kiv2_depth)x, single-copy exons $(metric "$H" flank_depth)x, MAPQ 0 share $(metric "$H" mapq0_share)" >> "$E2E_NOTES"

# --- half the KIV-2 reads: N + 2 halves ------------------------------------------------
printf 'chr6\t160611561\t160646395\n' > "${T}/block.bed"
dk "$SAMTOOLS_IMAGE" samtools view -b -L "${T}/block.bed" -U "${T}/rest.bam" -o /dev/null "${T}/lpa.bam"
dk "$SAMTOOLS_IMAGE" samtools view -b -L "${T}/block.bed" -s 7.5 -o "${T}/half_block.bam" "${T}/lpa.bam"
dk "$SAMTOOLS_IMAGE" samtools merge -f -o "${T}/half.bam" "${T}/rest.bam" "${T}/half_block.bam"
dk "$SAMTOOLS_IMAGE" samtools index "${T}/half.bam"
V="${SAMPLE}lpah"
sample "$V" "${T}/half.bam"
step "$V"
NH=$(metric "$V" kiv2_copies)
check "half the KIV-2 reads: N + 2 halves (${N} -> ${NH})" \
  awk -v a="$N" -v b="$NH" 'BEGIN {r = (b + 2) / (a + 2); exit !(r > 0.42 && r < 0.58)}'
check_eq "half the KIV-2 reads: the single-copy depth is unchanged" "$(metric "$V" flank_depth)" "$(metric "$H" flank_depth)"
echo "lpa: HG002 with half the KIV-2 reads: ${NH}" >> "$E2E_NOTES"

# --- an eighth of the depth: too little --------------------------------------------------
dk "$SAMTOOLS_IMAGE" samtools view -b -s 42.125 -o "${T}/thin.bam" "${T}/lpa.bam"
dk "$SAMTOOLS_IMAGE" samtools index "${T}/thin.bam"
L="${SAMPLE}lpat"
sample "$L" "${T}/thin.bam"
step "$L"
check_eq "an eighth of the depth: too little depth" "$(metric "$L" status)" "too little depth"
check_eq "an eighth of the depth: no count" "$(metric "$L" kiv2_copies)" "."

# --- Nextflow against the bash step ------------------------------------------------------
command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
check "case 21's VCF is there" test -s "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi"
mkdir -p "${G}/${H}/vcf"
cp "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${G}/${H}/vcf/${H}.vcf.gz"
cp "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" "${G}/${H}/vcf/${H}.vcf.gz.tbi"
SHEET="${T}/samplesheet.csv"
printf 'sample,sex,vcf,vcf_index,bam,bam_index\n%s,male,%s,%s,%s,%s\n' "$H" \
  "${G}/${H}/vcf/${H}.vcf.gz" "${G}/${H}/vcf/${H}.vcf.gz.tbi" \
  "${G}/${H}/aligned/${H}_sorted.bam" "${G}/${H}/aligned/${H}_sorted.bam.bai" > "$SHEET"
NF_OUT="${G}/nf-lpa"
echo "+ nextflow run main.nf -profile docker --tools lpa,html_report"
( cd "$T" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -work-dir "${E2E_WORK}/nf-work-lpa" \
    --input "$SHEET" --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" --tools lpa,html_report \
    --sex_check warn --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB ) 2>&1 | tee "$STEP_LOG"
NF_RC=${PIPESTATUS[0]}
check_eq "nextflow run exits 0" "$NF_RC" 0
cp "${T}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
for f in lpa_kiv2 lpa_kiv2_depth lpa_kiv2_depth_q1; do
  a="${G}/${H}/lpa/${H}_${f}.tsv" b="${NF_OUT}/${H}/lpa/${H}_${f}.tsv"
  check "Nextflow wrote ${H}_${f}.tsv" test -s "$b"
  check "Nextflow's ${H}_${f}.tsv is the bash step's" cmp -s "$a" "$b"
  cmp -s "$a" "$b" || diff "$a" "$b" | head -20
done
check "control: the comparison sees a difference" bash -c '! cmp -s "$1" "$2"' _ \
  "${G}/${V}/lpa/${V}_lpa_kiv2.tsv" "${NF_OUT}/${H}/lpa/${H}_lpa_kiv2.tsv"
HTML=$(find "$NF_OUT" -name "${H}_report.html" | head -n 1)
check "the HTML report has the Lp(a) card with the count" has "<h2>Lp\(a\): LPA KIV-2 Repeat Count</h2>" "$(cat "$HTML" 2>/dev/null)"
check "and the card shows the count" has "KIV-2 repeats \(both chromosomes\): ${N} " "$(cat "$HTML" 2>/dev/null)"

rm -rf "${G:?}/${V}" "${G:?}/${L}" "${T:?}"/*.bam "${T:?}"/*.bai "${T:?}/lpa.sam"
finish
