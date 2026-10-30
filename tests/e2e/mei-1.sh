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
#     known disease insertions are not callable.
# MEI_E2E_DROP_ALU=1 is the red control, run by hand: every read pair within
# 1 kb of a GIAB Alu insertion is removed before the step, and the recall
# check must fail.
. "$(dirname "$0")/lib.sh"

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

rm -rf "${G:?}/${H}/aligned" "${R:?}" "${G:?}/reference/megane" "${T:?}"/*.bam "${T:?}"/*.bai "${T:?}/truth.vcf.gz"
# MEGAnE's image (about 6.4 GB) serves this case alone, and the cases after it
# need the disk: the runner ran out of space in mocha-1 with it still loaded.
# It is pinned by digest, so the e2e image cache never saves it anyway.
docker image rm "$MEGANE_IMAGE" > /dev/null 2>&1 || true
unset REF_FASTA
finish
