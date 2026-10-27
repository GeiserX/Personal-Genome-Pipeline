#!/usr/bin/env bash
# Step 43 (archaic ancestry, IBDmix) on a 10 Mb slice of chr22.
#   - without the archaic data, step 43 says so in one line and exits 0;
#   - setup.sh --archaic builds the table for GRCh37 chr22:25-35 Mb only
#     (ARCHAIC_REGIONS): the three archaic chr22 files from MPI-EVA, checked
#     against their pinned md5 lists, and the 1000 Genomes sites, lifted to
#     GRCh38;
#   - HG002 (GIAB's v4.2.1 genotypes on the lifted slice): the step runs, the
#     shares are numbers and the analysed span is most of the slice;
#   - PLANT, HG002 with one chromosome over chr22:28-29 Mb (GRCh38) copied
#     from Vindija 33.19 (0/1 wherever Vindija carries the ALT, no record
#     where it does not): a Neanderthal segment covers 80% of that window or
#     more, and the text report shows the step 43 block;
#   - COMMON, a person built from the table who carries the majority allele
#     at every site and so no allele rarer than 50%: no segment at all, the
#     absence control of the plant;
#   - the Nextflow processes (ARCHAIC_*) on PLANT write the same three files
#     as the bash step; a control shows the comparison sees a difference.
# The time and disk of the install go to the job summary.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR" T="$CASE_TMP"
DATA="${G}/reference/archaic/${ARCHAIC_DATA_VERSION}"
TRUTH="${FIXTURE_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
W_FROM=28000000 W_TO=29000000
dk() { local img=$1; shift; docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$img" "$@"; }
# kv SAMPLE KEY: one value of the sample's summary
kv() { awk -F'\t' -v k="$2" '$1 == k {print $2}' "${G}/$1/archaic/$1_archaic.tsv" 2>/dev/null; }
# covered SAMPLE ANCESTRY FROM TO: bp of FROM-TO under the sample's segments of that ancestry
covered() {
  awk -F'\t' -v a="$2" -v f="$3" -v t="$4" 'NR > 1 && $1 == "chr22" && $5 == a {
      s = ($2 > f ? $2 : f); e = ($3 < t ? $3 : t); if (e > s) n += e - s } END {print n + 0}' \
    "${G}/$1/archaic/$1_archaic_segments.tsv" 2>/dev/null
}
sample() { rm -rf "${G:?}/$1"; mkdir -p "${G}/$1/vcf"; cp "$2" "${G}/$1/vcf/$1.vcf.gz"; dk "$BCFTOOLS_IMAGE" bcftools index -t -f "${G}/$1/vcf/$1.vcf.gz"; }
step() { run_step 43-archaic.sh "$1"; check_step_exit "43-archaic.sh ($1)"; }

GENES=$(bash -c '. "$1/scripts/lib/common.sh" && install_data_file gencode_genes >/dev/null && data_file gencode_genes' _ "$REPO")
check "the GENCODE gene lines are installed" test -s "$GENES"

# --- no data: one line, exit 0 ------------------------------------------------------
ARCHAIC_DIR=none run_step 43-archaic.sh "$SAMPLE"
check_step_exit "43-archaic.sh (no data)"
check_eq "step 43 without data prints one line" "$(grep -c . "$STEP_LOG")" 1
check "that line says it was skipped and how to install the data" has '^Step 43 skipped: .*setup\.sh --archaic' "$(cat "$STEP_LOG")"

# --- the archaic table of the slice ----------------------------------------------------
printf '22\t25000000\t35000000\n' > "${T}/regions37.bed"
echo "+ ARCHAIC_REGIONS=... scripts/setup.sh --archaic"
T0=$(date +%s)
ARCHAIC_REGIONS="${T}/regions37.bed" "${REPO}/scripts/setup.sh" --archaic "$G" 2>&1 | tee "$STEP_LOG"
check_eq "setup.sh --archaic exits 0" "${PIPESTATUS[0]}" 0
T_SETUP=$(( $(date +%s) - T0 ))
check "the table is finished (DONE)" test -s "${DATA}/DONE"
check "only chr22 was built" bash -c '[ "$(cd "$1" && ls chr*.tsv.gz)" = chr22.tsv.gz ]' _ "$DATA"
TABLE="${DATA}/chr22.tsv.gz"
check_ge "table rows on the slice" "$(gzip -dc "$TABLE" 2>/dev/null | awk 'END {print NR - 1}')" 50000
check "every row is on GRCh38 chr22 between 24 and 36 Mb" \
  bash -c 'gzip -dc "$1" | awk -F"\t" "NR > 1 && (\$1 != \"chr22\" || \$2 < 24000000 || \$2 > 36000000) {bad = 1} END {exit bad}"' _ "$TABLE"
check_ge "rows where Vindija carries the ALT inside the plant window" \
  "$(gzip -dc "$TABLE" | awk -F'\t' -v f="$W_FROM" -v t="$W_TO" 'NR > 1 && $2 >= f && $2 < t && ($6 == 1 || $6 == 2) {n++} END {print n + 0}')" 300
printf '#### Step 43 archaic table (GRCh37 chr22:25-35 Mb)\n\n%s s to install, %s on disk, %s rows\n\n' \
  "$T_SETUP" "$(du -sh "$DATA" | cut -f1)" "$(gzip -dc "$TABLE" | awk 'END {print NR - 1}')" >> "$E2E_NOTES"

# --- HG002 ------------------------------------------------------------------------------
dk "$BCFTOOLS_IMAGE" bcftools view -r chr22:24000000-36000000 -Oz -o "${T}/hg002.vcf.gz" "$TRUTH"
H="${SAMPLE}ar"
sample "$H" "${T}/hg002.vcf.gz"
step "$H"
cat "${G}/${H}/archaic/${H}_archaic.txt" 2>/dev/null
check "Neanderthal share is a number" has '^[0-9]+\.[0-9]+$' "$(kv "$H" neanderthal_percent)"
check "Denisovan share is a number" has '^[0-9]+\.[0-9]+$' "$(kv "$H" denisovan_percent)"
check_ge "analysed span (bp) covers most of the slice" "$(kv "$H" analysed_bp)" 8000000
check "the segments file has its header" has '^chrom	start	end	length_bp	ancestry	matched	lod	genes$' \
  "$(head -n 1 "${G}/${H}/archaic/${H}_archaic_segments.tsv" 2>/dev/null)"

# --- PLANT: one chromosome copied from Vindija over 1 Mb -------------------------------------
{ dk "$BCFTOOLS_IMAGE" bcftools view -h "${T}/hg002.vcf.gz"
  dk "$BCFTOOLS_IMAGE" bcftools view -H -t "^chr22:${W_FROM}-$((W_TO - 1))" "${T}/hg002.vcf.gz"
  gzip -dc "$TABLE" | awk -F'\t' -v f="$W_FROM" -v t="$W_TO" 'BEGIN {OFS = "\t"}
    NR > 1 && $2 >= f && $2 < t && ($6 == 1 || $6 == 2) {print "chr22", $2, ".", $3, $4, 50, "PASS", ".", "GT", "0/1"}'
} > "${T}/plant.vcf"
dk "$BCFTOOLS_IMAGE" bcftools sort -Oz -o "${T}/plant.vcf.gz" "${T}/plant.vcf"
P="${SAMPLE}arp"
sample "$P" "${T}/plant.vcf.gz"
step "$P"
cat "${G}/${P}/archaic/${P}_archaic_segments.tsv" 2>/dev/null | head -20
check_ge "a Neanderthal segment covers 80% of the planted window (bp)" "$(covered "$P" Neanderthal "$W_FROM" "$W_TO")" 800000
check "and it names Vindija among the genomes it matches" \
  awk -F'\t' -v f="$W_FROM" -v t="$W_TO" 'NR > 1 && $5 == "Neanderthal" && $2 < t && $3 > f && $6 ~ /Vindija/ {ok = 1} END {exit !ok}' \
  "${G}/${P}/archaic/${P}_archaic_segments.tsv"
check "the segment lists the protein-coding genes it crosses" \
  awk -F'\t' -v f="$W_FROM" -v t="$W_TO" 'NR > 1 && $5 == "Neanderthal" && $2 < t && $3 > f && $8 ~ /[A-Z]/ {ok = 1} END {exit !ok}' \
  "${G}/${P}/archaic/${P}_archaic_segments.tsv"
NP=$(kv "$P" neanderthal_bp) NH=$(kv "$H" neanderthal_bp)
check_ge "the plant adds Neanderthal sequence to HG002's (bp)" "$(( ${NP:-0} - ${NH:-0} ))" 1

run_step generate-report.sh "$P"
check_step_exit "generate-report.sh ($P)"
REPORT_TXT=$(cat "${G}/${P}/${P}_report.txt" 2>/dev/null)
check "the text report has the step 43 block" has 'Archaic Ancestry \(IBDmix, step 43\)' "$REPORT_TXT"
check "with the Neanderthal line of the summary" has "Neanderthal: $(kv "$P" neanderthal_percent)% " "$REPORT_TXT"

# --- COMMON: the majority allele everywhere -----------------------------------------------
{ printf '##fileformat=VCFv4.2\n##contig=<ID=chr22,length=50818468>\n##FILTER=<ID=PASS,Description="All filters passed">\n'
  printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCOMMON\n'
  gzip -dc "$TABLE" | awk -F'\t' 'BEGIN {OFS = "\t"} NR > 1 && $8 >= 0.5 {print "chr22", $2, ".", $3, $4, 50, "PASS", ".", "GT", "1/1"}'
} > "${T}/common.vcf"
dk "$BCFTOOLS_IMAGE" bcftools sort -Oz -o "${T}/common.vcf.gz" "${T}/common.vcf"
C="${SAMPLE}arc"
sample "$C" "${T}/common.vcf.gz"
step "$C"
check_eq "COMMON carries no rare allele: no Neanderthal segment" "$(kv "$C" neanderthal_segments)" 0
check_eq "and no Denisovan segment" "$(kv "$C" denisovan_segments)" 0
check_eq "the segments file is the header alone" "$(awk 'END {print NR}' "${G}/${C}/archaic/${C}_archaic_segments.tsv" 2>/dev/null)" 1

# --- Nextflow against the bash step ----------------------------------------------------------
command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
SHEET="${T}/samplesheet.csv" NF_OUT="${G}/nf-archaic"
printf 'sample,vcf,vcf_index,sex\n%s,%s,%s,male\n' "$P" "${G}/${P}/vcf/${P}.vcf.gz" "${G}/${P}/vcf/${P}.vcf.gz.tbi" > "$SHEET"
echo "+ nextflow run main.nf -profile docker --tools archaic"
( cd "$T" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -work-dir "${E2E_WORK}/nf-work-archaic" \
    --input "$SHEET" --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" --tools archaic \
    --archaic_data "$DATA" --hla_genes "$GENES" --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB ) 2>&1 | tee "$STEP_LOG"
check_eq "nextflow run exits 0" "${PIPESTATUS[0]}" 0
cp "${T}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
for f in archaic.tsv archaic_segments.tsv archaic.txt; do
  a="${G}/${P}/archaic/${P}_${f}" b="${NF_OUT}/${P}/archaic/${P}_${f}"
  check "Nextflow wrote ${P}_${f}" test -s "$b"
  check "Nextflow's ${P}_${f} is the bash step's" cmp -s "$a" "$b"
  cmp -s "$a" "$b" || diff "$a" "$b" | head -20
done
check "control: the comparison sees a difference" bash -c '! cmp -s "$1" "$2"' _ \
  "${G}/${H}/archaic/${H}_archaic_segments.tsv" "${NF_OUT}/${P}/archaic/${P}_archaic_segments.tsv"

finish
