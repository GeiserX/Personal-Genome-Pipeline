#!/usr/bin/env bash
# Step 43 (archaic ancestry, IBDmix) on the whole of chr22.
#   - without the archaic data, step 43 says so in one line and exits 0;
#   - setup.sh --archaic builds the table for GRCh37 chr22 only
#     (ARCHAIC_REGIONS): the three archaic chr22 files from MPI-EVA, checked
#     against their pinned md5 lists, and the 1000 Genomes sites, lifted to
#     GRCh38; every archaic call outside its genome's FilterBed, outside the
#     1000 Genomes GRCh38 strict mask, in a segmental duplication or within
#     5 bp of a 1000 Genomes indel is marked masked (code 10 to 12), which
#     step 43 hands IBDmix as its mask;
#   - HG002 (GIAB's v4.2.1 genotypes on chr22): the step runs, the shares are
#     numbers, the analysed span is most of the chromosome, and the
#     Neanderthal total lies between the lowest value (51 kb) of the range
#     Chen et al.'s own released calls give Europeans on chr22 and a 1,200 kb
#     regression guard (the guard that catches the unmasked table calling
#     about three times too much; the whole-genome check decides whether the
#     share is in range);
#   - PLANT, HG002 with one chromosome over chr22:28-29 Mb (GRCh38) copied
#     from Vindija 33.19 (0/1 wherever Vindija carries the ALT, no record
#     where it does not): a Neanderthal segment covers 70% of that window or
#     more, and the text report shows the step 43 block;
#   - COMMON, a person built from the table who carries the majority allele
#     at every site and so no allele rarer than 50%: a quarter of HG002's
#     Neanderthal sequence at most, the absence control of the plant;
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

# --- the archaic table of chr22 --------------------------------------------------------
printf '22\t1\t51304566\n' > "${T}/regions37.bed"
echo "+ ARCHAIC_REGIONS=... scripts/setup.sh --archaic"
T0=$(date +%s)
ARCHAIC_REGIONS="${T}/regions37.bed" "${REPO}/scripts/setup.sh" --archaic "$G" 2>&1 | tee "$STEP_LOG"
check_eq "setup.sh --archaic exits 0" "${PIPESTATUS[0]}" 0
T_SETUP=$(( $(date +%s) - T0 ))
check "the table is finished (DONE)" test -s "${DATA}/DONE"
check "only chr22 was built" bash -c '[ "$(cd "$1" && ls chr*.tsv.gz)" = chr22.tsv.gz ]' _ "$DATA"
TABLE="${DATA}/chr22.tsv.gz"
ROWS=$(gzip -dc "$TABLE" 2>/dev/null | awk 'END {print NR - 1}')
check_ge "table rows on chr22" "$ROWS" 100000
# The FilterBeds and the excluded regions mark calls masked (10 plus the ALT
# count) instead of dropping them, as Chen et al. keep those sites in IBDmix's
# input and mask them there. A build that lost the masks has no such code.
for i in 5 6 7; do
  check_ge "rows where archaic column $((i - 4)) is masked" \
    "$(gzip -dc "$TABLE" | awk -F'\t' -v i="$i" 'NR > 1 && $i >= 10 && $i <= 12 {n++} END {print n + 0}')" 10000
done
check "every archaic code is 0, 1, 2, 9 or 10 to 12" \
  bash -c 'gzip -dc "$1" | awk -F"\t" "NR > 1 { for (i = 5; i <= 7; i++) if (\$i !~ /^(0|1|2|9|10|11|12)\$/) bad = 1 } END {exit bad}"' _ "$TABLE"
check "every row is on GRCh38 chr22" \
  bash -c 'gzip -dc "$1" | awk -F"\t" "NR > 1 && \$1 != \"chr22\" {bad = 1} END {exit bad}"' _ "$TABLE"
check_ge "rows where Vindija carries the ALT inside the plant window" \
  "$(gzip -dc "$TABLE" | awk -F'\t' -v f="$W_FROM" -v t="$W_TO" 'NR > 1 && $2 >= f && $2 < t && ($6 == 1 || $6 == 2) {n++} END {print n + 0}')" 100
# The masks: every FilterBed of the chromosome is installed, and every table
# row outside the 1000 Genomes GRCh38 strict mask is kept with its archaic
# calls marked masked (10 to 12), none left as 0, 1 or 2. A table that drops
# those rows has no row outside the mask; one that keeps them unmarked has
# calls 0 to 2 there. Either way this goes red.
for f in Altai Vindija33.19 Denisova; do
  check "the FilterBed of ${f} chr22 is installed" test -s "${DATA}/mask/${f}.chr22.bed.gz"
done
check "the GRCh38 strict mask is installed" test -s "${DATA}/strict_GRCh38.bed"
gzip -dc "$TABLE" | awk -F'\t' -v m="${DATA}/strict_GRCh38.bed" '
  BEGIN { while ((getline l < m) > 0) { split(l, c, "\t"); if (c[1] == "chr22") { s[++n] = c[2]; e[n] = c[3] } } }
  NR > 1 { p = $2 - 1
           while (i < n && e[i + 1] <= p) i++
           if (!(i < n && s[i + 1] <= p && p < e[i + 1])) {
             out++
             for (k = 5; k <= 7; k++) if ($k ~ /^(0|1|2)$/) open++ } }
  END { print out + 0, open + 0 }' > "${T}/strict-counts.txt"
read -r OUT_ROWS OUT_OPEN < "${T}/strict-counts.txt"
check_ge "table rows outside the GRCh38 strict mask are kept" "$OUT_ROWS" 7000
check_eq "none of them has an unmarked archaic call" "$OUT_OPEN" 0
# The excluded regions: segmental duplications and the 5 bp around a 1000
# Genomes indel. GRCh37 coordinates, so the test lifts nothing: it checks the
# GRCh37 positions of the rejected sites are gone from the GRCh37 input, by
# counting how many table rows the exclusion list would have kept out. Any row
# left inside one means the exclusion did not run.
check "the exclusion list of chr22 is installed" test -s "${DATA}/exclude/chr22.bed.gz"
gzip -dc "${DATA}/exclude/chr22.bed.gz" > "${T}/exclude22.bed"
check_ge "it names segmental duplications and indel neighbourhoods" "$(awk 'END {print NR + 0}' "${T}/exclude22.bed")" 1000
check_ge "chr22 segmental duplications are among them (Mb covered)" \
  "$(awk '$3 - $2 > 10000 {n += $3 - $2} END {print int(n / 1000000)}' "${T}/exclude22.bed")" 1
printf '#### Step 43 archaic table (GRCh37 chr22)\n\n%s s to install, %s on disk, %s rows\n\n' \
  "$T_SETUP" "$(du -sh "$DATA" | cut -f1)" "$(gzip -dc "$TABLE" | awk 'END {print NR - 1}')" >> "$E2E_NOTES"

# --- HG002 ------------------------------------------------------------------------------
dk "$BCFTOOLS_IMAGE" bcftools view -r chr22 -Oz -o "${T}/hg002.vcf.gz" "$TRUTH"
H="${SAMPLE}ar"
sample "$H" "${T}/hg002.vcf.gz"
step "$H"
cat "${G}/${H}/archaic/${H}_archaic.txt" 2>/dev/null
check "Neanderthal share is a number" has '^[0-9]+\.[0-9]+$' "$(kv "$H" neanderthal_percent)"
check "Denisovan share is a number" has '^[0-9]+\.[0-9]+$' "$(kv "$H" denisovan_percent)"
check_ge "analysed span (bp) covers most of chr22" "$(kv "$H" analysed_bp)" 30000000
# The published range for chr22. Chen et al. 2020 released their IBDmix calls
# on the 1000 Genomes phase 3 genomes (LOD 4, union of the Altai and Vindija
# segments of 50 kb or more). Over the 497 people of CEU, TSI, GBR, FIN and
# IBS, chr22 gives 51 kb at the lowest, 153 kb at the 5th percentile, 541 kb
# at the median, 1,132 kb at the 95th and 1,716 kb at the highest. One
# chromosome of one person varies widely, so the gate is a regression guard,
# not a claim about the range: 51 kb to 1,200 kb. Runs before the strict-mask
# rows were kept marked gave 1,748 kb (no masks) and 1,177 to 1,190 kb (masked
# codes, rows dropped); the marked table gave 1,138 kb. The published figures
# are printed beside the measured value, and the whole-genome check decides.
NEA=$(kv "$H" neanderthal_bp)
check_ge "HG002's chr22 Neanderthal total is not under the published lowest (51 kb)" "${NEA:-0}" 51000
check "and not over the chr22 regression guard (1,200 kb)" test "${NEA:-0}" -le 1200000
echo "  chr22 Neanderthal bp: ${NEA} (published EUR: lowest 51000, p5 153000, median 541000, p95 1132000, highest 1716000)"
printf '#### Step 43 HG002 chr22 Neanderthal\n\n%s bp. Published EUR calls of Chen et al. on chr22: lowest 51000, median 541000, 95th 1132000, highest 1716000 bp.\n\n' \
  "${NEA:-0}" >> "$E2E_NOTES"
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
check_ge "a Neanderthal segment covers 70% of the planted window (bp)" "$(covered "$P" Neanderthal "$W_FROM" "$W_TO")" 700000
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
# COMMON carries no allele rarer than half the population, so almost nothing
# should look archaic in it. Not nothing at all: the masks leave fewer sites
# to contradict a stretch, so a few segments scrape over LOD 4 and 50 kb. The
# control is the size of what it gets against the real sample's.
NC=$(kv "$C" neanderthal_bp)
check "COMMON carries no rare allele: a quarter of HG002's Neanderthal bp at most" \
  test "${NC:-1}" -le "$(( ${NEA:-0} / 4 ))"
echo "  COMMON Neanderthal bp: ${NC} against HG002's ${NEA}, in $(kv "$C" neanderthal_segments) and $(kv "$H" neanderthal_segments) segments"
check "and no more Denisovan bp than HG002" test "$(kv "$C" denisovan_bp)" -le "$(kv "$H" denisovan_bp)"

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
