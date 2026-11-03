#!/usr/bin/env bash
# Step 50 (MEGAnE: Alu, L1 and SVA insertions) on real reads:
#   - HG002's chr20 at 30x: the reads of chr20 in GIAB's 60x novoalign BAM,
#     half the pairs, realigned with minimap2 to chr20 alone. MEGAnE's k-mer
#     set holds 16 bytes per reference base in memory, so the fixture's
#     1.4 Gb reference does not fit a runner and chr20's 64 Mb does. The reads
#     left unmapped are dropped: MEGAnE searches each one against the whole
#     genome, which took 21 minutes of a 23-minute run on this slice;
#   - setup.sh --megane-kmers builds the k-mer set of that reference, and
#     scripts/50-mei.sh calls the sample. One chromosome holds too few
#     insertions for MEGAnE to set its filter thresholds, so the step says so
#     and uses MEGAnE's fixed ones;
#   - its Alu calls against the GIAB v5.0q insertions of HG002 on chr20
#     (RepeatMasker SINE/Alu, interspersed, 250 to 400 bp, a 1 in the
#     genotype, inside the benchmark regions): recall at least 0.70, and
#     precision, the Alu calls inside the benchmark regions within 100 bp of
#     any GIAB insertion, at least 0.80;
#   - the catalog: the BAM has none of chr9, chrX and chr11, so all three
#     known disease insertions are not callable;
#   - AnnotSV on the calls, with the made-up AnnotSV data of
#     tests/fixtures/annotsv/class45.sh plus one made-up gene, FIXCHR20, over
#     all of chr20 (two exons at its ends): every call comes back in
#     FIXCHR20, with its location there, and with AnnotSV's class of an
#     insertion, NA. Red control: the same VCF through bin/mei_report.py
#     without AnnotSV's table fails those checks;
#   - the Nextflow pipeline (--tools mei, with --annotsv_annotations) on the
#     same BAM writes the same calls table as the bash step, and MEI's peak
#     memory from its trace goes to the notes.
# MEI_E2E_DROP_ALU=1 is the red control, run by hand: every read pair within
# 1 kb of a GIAB Alu insertion is removed before the step, and the recall
# check must fail.
. "$(dirname "$0")/lib.sh"
# shellcheck source=../fixtures/annotsv/class45.sh
. "${REPO}/tests/fixtures/annotsv/class45.sh"

G="$GENOME_DIR" T="$CASE_TMP"
SRC="https://giab.s3.amazonaws.com/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam"
TRUTH="https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_stvar"
H="${SAMPLE}mei"
R="${G}/reference/mei"
BAM="${G}/${H}/aligned/${H}_sorted.bam"
dk() { docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$1" "${@:2}"; }
retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && { echo "  try ${i} of 4 failed; again in $((i * 15)) s"; sleep $((i * 15)); }; done; return 1; }

# --- chr20 alone, and its k-mer set ------------------------------------------------------
rm -rf "${R:?}" "${G:?}/reference/megane" "${G:?}/${H}"
mkdir -p "$R" "${G}/${H}/aligned"
# shellcheck disable=SC2016  # $1 and $2 belong to the inner shell
dk "$SAMTOOLS_IMAGE" sh -c 'samtools faidx "$1" chr20 > "$2" && samtools faidx "$2"' _ \
  "${G}/reference/GRCh38_no_alt_analysis_set.fasta" "${R}/chr20.fa"
check_eq "chr20 alone, all of it" "$(cut -f1,2 "${R}/chr20.fa.fai" 2>/dev/null)" $'chr20\t64444167'
export REF_FASTA="${R}/chr20.fa"
echo "+ scripts/setup.sh --megane-kmers"
"${REPO}/scripts/setup.sh" --megane-kmers "$G" 2>&1 | tee "$STEP_LOG"
check_eq "setup.sh --megane-kmers exits 0" "${PIPESTATUS[0]}" 0
MK="${G}/reference/megane/chr20.mk"
check "the k-mer set is where step 50 looks for it" test -s "$MK"
check_ge "repeated 32-mers of chr20" "$(cat "${MK%.mk}.mi" 2>/dev/null || echo 0)" 1000000

# --- HG002's chr20 at 30x, realigned -----------------------------------------------------
# SAMTOOLS_IMAGE's libcurl has no CA bundle: lend it the host's (htslib reads
# CURL_CA_BUNDLE). The index is fetched once, so the eight streams do not each
# write it. A pair is kept when its name hashes into half (-s 42.5); a read is
# taken by the chunk its position falls in, so a read over a chunk boundary
# comes once.
check "the BAM's index downloads" retry curl -sfL -o "${T}/src.bam.bai" "${SRC}.bai"
# shellcheck disable=SC2329  # run through retry
chunk() {
  rm -f "${T}/part$1.bam"
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$SAMTOOLS_IMAGE" \
    bash -c 'set -o pipefail; samtools view -h -s 42.5 "$1##idx##$2" "chr20:$3-$4" | awk -v s="$3" -v e="$4" '"'"'/^@/ || ($4 >= s && $4 <= e)'"'"' | samtools view -b -o "part$5.bam" -' \
    _ "$SRC" "${T}/src.bam.bai" "$2" "$3" "$1"
}
LEN=64444167 N=8
STEP=$(( LEN / N + 1 ))
pids=()
for i in $(seq 0 $((N - 1))); do
  retry chunk "$i" $(( i * STEP + 1 )) $(( (i + 1) * STEP )) & pids+=($!)
done
ok=0
for p in "${pids[@]}"; do wait "$p" && ok=$((ok + 1)); done
check_eq "the eight chr20 streams finish" "$ok" "$N"
dk "$SAMTOOLS_IMAGE" bash -c 'set -o pipefail; samtools cat -o all.bam part*.bam && samtools collate -u -O all.bam coll | samtools fastq -n -1 R1.fq.gz -2 R2.fq.gz -0 /dev/null -s /dev/null -' \
  && rm -f "${T}"/part*.bam "${T}/all.bam"
# shellcheck disable=SC2016  # $1 to $3 belong to the inner bash
dk "$MINIMAP2_SAMTOOLS_IMAGE" bash -c 'set -o pipefail; minimap2 -t 4 -a -x sr -R "@RG\tID:HG002\tSM:HG002\tPL:ILLUMINA\tLB:HG002" "$1" R1.fq.gz R2.fq.gz | samtools view -u -F 4 - | samtools sort -@ 2 -m 1G -T sorttmp -o "$2" - && samtools index "$2"' \
  _ "${R}/chr20.fa" "${T}/all.bam"
rm -f "${T}/R1.fq.gz" "${T}/R2.fq.gz"
check_ge "mapped reads of chr20 at 30x" "$(dk "$SAMTOOLS_IMAGE" samtools view -c "${T}/all.bam" 2>/dev/null || echo 0)" 12000000

# --- the GIAB insertions -----------------------------------------------------------------
check "GIAB v5.0q's SV truth downloads" retry curl -sfL -o "${T}/truth.vcf.gz" "${TRUTH}.vcf.gz"
check "and its benchmark regions" retry curl -sfL -o "${T}/truth.bed" "${TRUTH}.benchmark.bed"
# truth.tsv: every insertion on chr20 (pos, RepeatMasker class, REMAP, length, genotype).
python3 -I - "${T}/truth.vcf.gz" > "${T}/truth.tsv" <<'PY'
import gzip, sys
for line in gzip.open(sys.argv[1], "rt"):
    if line[0] == "#" or not line.startswith("chr20\t"):
        continue
    f = line.split("\t", 10)
    info = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
    if info.get("SVTYPE") == "INS":
        print(f[1], info.get("RM_clsfam", "."), info.get("REMAP", "."), len(f[4]) - len(f[3]), f[9].split(":")[0], sep="\t")
PY
check_ge "GIAB insertions on chr20" "$(wc -l < "${T}/truth.tsv" | tr -d ' ')" 1000

if [ "${MEI_E2E_DROP_ALU:-}" = 1 ]; then
  echo "=== Red control: every pair within 1 kb of a GIAB Alu insertion removed ==="
  awk -F'\t' '$2 == "SINE/Alu" {s = $1 - 1000; print "chr20\t" (s < 0 ? 0 : s) "\t" $1 + 1000}' "${T}/truth.tsv" > "${T}/alu.bed"
  # shellcheck disable=SC2016  # \$1 is awk's
  dk "$SAMTOOLS_IMAGE" bash -c 'set -o pipefail; samtools view -L alu.bed all.bam | cut -f1 | sort -u > drop.txt && samtools view -h all.bam | awk -F "\t" "NR == FNR {d[\$1]; next} /^@/ || !(\$1 in d)" drop.txt - | samtools view -b -o keep.bam - && mv keep.bam all.bam && samtools index all.bam'
fi
mv "${T}/all.bam" "$BAM" && mv "${T}/all.bam.bai" "${BAM}.bai"

# --- AnnotSV's data: the made-up bundle and a gene over all of chr20 ---------------------
# Step 50 reads annotsv_annotations/ in GENOME_DIR; one another case left is
# set aside and put back at the end.
A="${G}/annotsv_annotations"
rm -rf "${T:?}/annotsv_saved"
[ ! -e "$A" ] || mv "$A" "${T}/annotsv_saved"
class45_bundle "$A"
GB="${A}/Annotations_Human/Genes/GRCh38"
for f in RefSeq ENSEMBL; do
  tx=NM_929999; [ "$f" = RefSeq ] || tx=ENST929999
  { printf '20\t0\t64444167\t+\tFIXCHR20\t%s\t500\t64443500\t0,64443167,\t1000,64444167,\n' "$tx"
    cat "${GB}/genes.${f}.sorted.bed"; } | sort -k1,1 -k2,2n > "${T}/genes.bed" && mv "${T}/genes.bed" "${GB}/genes.${f}.sorted.bed"
  printf '%s\t1\n' "$tx" >> "${GB}/transcript_version.${f}.tsv"
done
printf 'NM_929999\n' >> "${A}/mane_select_refseq.txt"

# --- step 50 -----------------------------------------------------------------------------
run_step 50-mei.sh "$H" male
check_step_exit 50-mei.sh
M="${G}/${H}/mei/${H}"
check "the step says MEGAnE used its fixed thresholds on one chromosome" grep -q "used fixed ones" "$STEP_LOG"
check "MEGAnE's VCF is kept" grep -q '^#CHROM' "${M}_megane.vcf"
check_ge "the calls table has insertions" "$(awk 'NR > 1' "${M}_mei.tsv" 2>/dev/null | wc -l | tr -d ' ')" 10
check_eq "the catalog: chr9, chrX and chr11 are not in the BAM, so every row is not callable" \
  "$(awk -F'\t' 'NR > 1 {print $1 ":" $(NF - 2)}' "${M}_mei_sites.tsv" 2>/dev/null | paste -sd, -)" \
  "FKTN:not callable,TAF1:not callable,BBS1:not callable"

# AnnotSV's columns: every call inside FIXCHR20, with a location there, class NA.
annotsv_cols() {  # annotsv_cols TABLE: "<calls>/<in FIXCHR20>/<with a FIXCHR20 location>/<class NA>"
  awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
               {n++; g += ($c["genes"] ~ /(^|,)FIXCHR20(,|$)/); l += ($c["location"] ~ /(^|,)FIXCHR20:/)
                a += ($c["acmg_class"] == "NA")}
               END {print n + 0 "/" g + 0 "/" l + 0 "/" a + 0}' "$1" 2>/dev/null
}
check "the step ran AnnotSV on the insertions" grep -q 'AnnotSV on the insertions' "$STEP_LOG"
check "AnnotSV's table of the insertions is kept" test -s "${M}_mei_annotsv.tsv"
N=$(awk 'NR > 1' "${M}_mei.tsv" 2>/dev/null | wc -l | tr -d ' ')
echo "AnnotSV columns (calls/in FIXCHR20/with its location/class NA): $(annotsv_cols "${M}_mei.tsv")"
check_eq "every call is in FIXCHR20, with its location there and AnnotSV's class NA" \
  "$(annotsv_cols "${M}_mei.tsv")" "${N}/${N}/${N}/${N}"
# Red control: the same VCF without AnnotSV's table.
python3 -I "${REPO}/bin/mei_report.py" --sample "$H" --vcf "${M}_megane.vcf" --sites "${REPO}/assets/mei_sites.tsv" \
  --depth "${M}_mei_depth.tsv" --out-calls "${T}/red_mei.tsv" --out-sites "${T}/red_sites.tsv" > /dev/null
check "red control: without AnnotSV's table the same check fails" \
  test "$(annotsv_cols "${T}/red_mei.tsv")" != "${N}/${N}/${N}/${N}" -a "$(annotsv_cols "${T}/red_mei.tsv")" = "${N}/0/0/0"

# Recall and precision of the Alu calls.
read -r TP NT PP NP < <(python3 -I - "${T}/truth.tsv" "${T}/truth.bed" "${M}_mei.tsv" <<'PY'
import csv, sys
truth = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
bed = [(int(a), int(b)) for c, a, b, *_ in (l.split("\t") for l in open(sys.argv[2])) if c == "chr20"]
inb = lambda p: any(s < p <= e for s, e in bed)
ins = [int(t[0]) for t in truth]
alu = [int(t[0]) for t in truth if t[1] == "SINE/Alu" and t[2] == "interspersed" and 250 <= int(t[3]) <= 400
       and "1" in t[4] and inb(int(t[0]))]
calls = [int(r["pos"]) for r in csv.DictReader(open(sys.argv[3]), delimiter="\t") if r["element"] == "Alu"]
found = sum(any(abs(c - t) <= 100 for c in calls) for t in alu)
inside = [c for c in calls if inb(c)]
right = sum(any(abs(c - t) <= 100 for t in ins) for c in inside)
print(found, len(alu), right, len(inside))
PY
)
echo "Alu recall ${TP}/${NT}, precision ${PP}/${NP}"
check_ge "GIAB Alu insertions to find" "${NT:-0}" 10
check "Alu recall (${TP}/${NT}) is at least 0.70" awk -v a="${TP:-0}" -v b="${NT:-0}" 'BEGIN {exit !(b > 0 && a / b >= 0.70)}'
check "Alu precision (${PP}/${NP}) is at least 0.80" awk -v a="${PP:-0}" -v b="${NP:-0}" 'BEGIN {exit !(b > 0 && a / b >= 0.80)}'
echo "mei: HG002 chr20 30x (minimap2, unmapped reads dropped${MEI_E2E_DROP_ALU:+, Alu reads removed}): Alu recall ${TP}/${NT}," \
  "precision ${PP}/${NP} within 100 bp of GIAB v5.0q; calls per element:" \
  "$(awk -F'\t' 'NR > 1 {n[$4]++} END {for (k in n) printf "%s %d ", k, n[k]}' "${M}_mei.tsv" 2>/dev/null)" >> "$E2E_NOTES"

# --- the Nextflow pipeline on the same BAM ------------------------------------------------
if command -v nextflow > /dev/null; then
  SHEET="${T}/samplesheet.csv" V="${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" NF_OUT="${G}/nf-mei"
  printf 'sample,sex,vcf,vcf_index,bam,bam_index\n%s,male,%s,%s,%s,%s\n' "$H" "$V" "${V}.tbi" "$BAM" "${BAM}.bai" > "$SHEET"
  rm -rf "${NF_OUT:?}"
  echo "+ nextflow run main.nf -profile docker --tools mei --annotsv_annotations"
  ( cd "$T" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -work-dir "${E2E_WORK}/nf-work-mei" \
      --input "$SHEET" --reference "${R}/chr20.fa" --tools mei --megane_kmers "$MK" --annotsv_annotations "$A" \
      --sex_check warn --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB ) 2>&1 | tee "$STEP_LOG"
  check_eq "nextflow run exits 0" "${PIPESTATUS[0]}" 0
  cp "${T}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
  NM="${NF_OUT}/${H}/mei/${H}"
  check "MEI_ANNOTSV's table is published" test -s "${NM}_mei_annotsv.tsv"
  check "Nextflow's calls table is the bash step's" cmp -s "${M}_mei.tsv" "${NM}_mei.tsv"
  cmp -s "${M}_mei.tsv" "${NM}_mei.tsv" || diff "${M}_mei.tsv" "${NM}_mei.tsv" | head -20
  # shellcheck disable=SC2016  # $1 and $2 belong to the inner bash
  check "control: the comparison sees a difference" bash -c '! cmp -s "$1" "$2"' _ "${T}/red_mei.tsv" "${NM}_mei.tsv"
  TRACE=$(find "${NF_OUT}/pipeline_info" -name 'trace_*.txt' 2>/dev/null | LC_ALL=C sort | tail -1)
  PEAK=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
                     $c["process"] == "MEI" {print $c["peak_rss"] " in " $c["realtime"] " on " $c["cpus"] " CPUs"}' "${TRACE:-/dev/null}" 2>/dev/null)
  echo "MEI (MEGAnE) peak_rss: ${PEAK:-not in the trace}"
  echo "mei: Nextflow MEI on HG002 chr20 30x, peak_rss ${PEAK:-not in the trace} (k-mer set $(du -h "$MK" | cut -f1))" >> "$E2E_NOTES"
  rm -rf "${NF_OUT:?}" "${E2E_WORK:?}/nf-work-mei"
else
  fail "nextflow is not on PATH"
fi

# Leave nothing behind: the cases after this one need the disk, and the runner
# ran out of space in them with MEGAnE's output and image still there
# (removing the files alone took the free disk from 7.6 GB to 8.9 GB).
# The image (about 6.4 GB) serves this case alone and is pinned by digest, so
# the e2e image cache never saves it. It goes by id, and the log says if it
# stayed.
echo "disk before the clean-up: $(df -h --output=avail "$E2E_WORK" | tail -1 | tr -d ' ')"
rm -rf "${G:?}/${H}" "${R:?}" "${G:?}/reference/megane" "${T:?}/truth.vcf.gz" "${A:?}"
[ ! -e "${T}/annotsv_saved" ] || mv "${T}/annotsv_saved" "$A"
id=$(docker image inspect -f '{{.Id}}' "$MEGANE_IMAGE" 2>/dev/null) && docker image rm -f "$id" > /dev/null
docker image inspect "$MEGANE_IMAGE" > /dev/null 2>&1 && echo "WARNING: MEGAnE's image is still on the disk"
echo "disk after the clean-up: $(df -h --output=avail "$E2E_WORK" | tail -1 | tr -d ' ')"
unset REF_FASTA
finish
