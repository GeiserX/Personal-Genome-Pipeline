#!/usr/bin/env bash
# Step 47 (viral DNA in the unmapped reads) on HG002's BAM from step 02, with
# synthetic viral reads added as unmapped reads. Nothing here says anything
# about HG002's own viruses: every expected call comes from the reads the case
# plants. The case:
#   - installs the viral genomes with setup.sh --viral-data (NCBI E-utilities,
#     checked against VIRAL_REF_SHA256) and finds 24 sequences;
#   - measures the human depth of the BAM itself (samtools coverage, chr1 to
#     chr22, over the bases with reads), and plants reads drawn from the HHV-6B
#     genome at half that depth (one integrated copy per cell), 300 reads of
#     CMV (a low level) and 200 reads of phiX (the lab control), with 0.3%
#     substitutions; then the step must call inherited integrated HHV-6B, CMV
#     viral DNA found but no integration, phiX the lab control, and HBV, with
#     no read planted, not found;
#   - plants HHV-6B at a tenth of the human depth instead (a viraemia-like
#     level): the same genome, the same breadth, and the call must not be
#     integration. This is the control that the depth threshold decides;
#   - runs the Nextflow processes (VIRAL_READS, VIRAL_REPORT) on the first
#     sample, whose table and coverages must be the bash step's, byte for byte;
#   - runs generate-report.sh: the text report states the HHV-6 call.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR" T="$CASE_TMP"
BASE="${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
VREF="${G}/reference/viral-${VIRAL_REF_VERSION}/viral_genomes.fa"
dk() { local img=$1; shift; docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$img" "$@"; }

# --- the viral genomes ----------------------------------------------------------------
echo "+ scripts/setup.sh --viral-data"
"${REPO}/scripts/setup.sh" --viral-data "$G" 2>&1 | tail -n 5
check "setup.sh --viral-data wrote the FASTA" test -s "$VREF"
check_eq "the FASTA is the pinned download" "$(sha256sum "$VREF" 2>/dev/null | cut -d' ' -f1)" "$VIRAL_REF_SHA256"
check_eq "24 viral genomes" "$(grep -c '^>' "$VREF" 2>/dev/null)" 24
check "the BAM of step 02 is there" test -s "${BASE}.bai"

# --- the human depth, measured here -------------------------------------------------------
dk "$SAMTOOLS_IMAGE" samtools coverage --ff UNMAP,SECONDARY,QCFAIL "$BASE" > "${T}/human_coverage.tsv"
HUMAN=$(awk -F'\t' '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/ && $5 > 0 {t += $7 * ($3 - $2 + 1); c += $5}
  END {if (c) printf "%.2f", t / c}' "${T}/human_coverage.tsv")
check "a human depth over the autosomes (${HUMAN:-none})" test -n "$HUMAN"
echo "human depth over the autosomal bases with reads: ${HUMAN}"

# plant OUT FRACTION: the step-02 BAM plus unmapped reads from HHV-6B at FRACTION of the human
# depth, 300 CMV reads and 200 phiX reads (150 bp, random strand, 0.3% substitutions, seed fixed).
plant() {
  python3 - "$VREF" "$HUMAN" "$2" > "${T}/spike.sam" <<'PY'
import random, sys
fasta, human, frac = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
seqs, name = {}, None
for line in open(fasta):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:].split()[0]; seqs[name] = []
    elif line:
        seqs[name].append(line.upper())
seqs = {k: "".join(v) for k, v in seqs.items()}
rng, L, comp = random.Random(47), 150, str.maketrans("ACGT", "TGCA")
def reads(acc, n, tag):
    s = seqs[acc]
    for i in range(n):
        p = rng.randrange(0, len(s) - L)
        r = list(s[p:p + L])
        for j in range(L):
            if rng.random() < 0.003:
                r[j] = rng.choice([b for b in "ACGT" if b != r[j]])
        r = "".join(r)
        if rng.random() < 0.5:
            r = r.translate(comp)[::-1]
        print(f"{tag}_{i}\t4\t*\t0\t0\t*\t*\t0\t0\t{r}\t{'I' * L}")
reads("NC_000898.1", round(frac * human * len(seqs["NC_000898.1"]) / L), "spk6b")
reads("NC_006273.2", 300, "spkcmv")
reads("NC_001422.1", 200, "spkphix")
PY
  mkdir -p "${G}/$1/aligned" "${G}/$1/vcf"
  # The planted records go after the BAM's own unmapped reads, at the end of a coordinate-sorted file.
  # shellcheck disable=SC2016
  dk "$SAMTOOLS_IMAGE" bash -euo pipefail -c 'samtools view -h "$1" | cat - "$2" | samtools view -b -o "$3"; samtools index "$3"' \
    _ "$BASE" "${T}/spike.sam" "${G}/$1/aligned/$1_sorted.bam"
  cp "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${G}/$1/vcf/$1.vcf.gz"
  dk "$BCFTOOLS_IMAGE" bcftools index -t -f "${G}/$1/vcf/$1.vcf.gz"
  wc -l < "${T}/spike.sam"
}
# cell NAME ACCESSION COLUMN: one cell of the sample's table.
cell() {
  awk -F'\t' -v a="$2" -v k="$3" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $2 == a {print $c[k]}' \
    "${G}/$1/viral/$1_viral.tsv" 2>/dev/null
}
between() { awk -v x="$1" -v lo="$2" -v hi="$3" 'BEGIN {exit !(x != "" && x + 0 >= lo && x + 0 <= hi)}'; }

# --- HHV-6B at half the human depth -----------------------------------------------------
V="${SAMPLE}vir"
N_PLANTED=$(plant "$V" 0.5)
check_ge "reads planted" "$N_PLANTED" 1000
run_step 47-viral.sh "$V"
check_step_exit 47-viral.sh
cat "${G}/${V}/viral/${V}_viral.tsv" 2>/dev/null
check_eq "a row per viral genome" "$(awk 'END {print NR - 1}' "${G}/${V}/viral/${V}_viral.tsv" 2>/dev/null)" 24
check_eq "the step's human depth is the one measured here" "$(cell "$V" NC_000898.1 human_depth)" "$HUMAN"
N6B=$(grep -c '^spk6b_' "${T}/spike.sam")
check_ge "HHV-6B: the planted reads realigned to it (of ${N6B})" "$(cell "$V" NC_000898.1 reads)" "$((N6B * 95 / 100))"
check "HHV-6B: near the whole genome covered ($(cell "$V" NC_000898.1 breadth))" between "$(cell "$V" NC_000898.1 breadth)" 0.95 1
check "HHV-6B: depth at about half the human depth ($(cell "$V" NC_000898.1 relative_depth))" \
  between "$(cell "$V" NC_000898.1 relative_depth)" 0.40 0.60
check_eq "HHV-6B: inherited integrated HHV-6 likely" "$(cell "$V" NC_000898.1 call)" "inherited integrated HHV-6 likely"
check "HHV-6A is not called integrated" bash -c '[ "$1" != "inherited integrated HHV-6 likely" ]' _ "$(cell "$V" NC_001664.4 call)"
check_ge "CMV: the 300 planted reads" "$(cell "$V" NC_006273.2 reads)" 290
check_eq "CMV: viral DNA found" "$(cell "$V" NC_006273.2 call)" "viral DNA found"
check "CMV: far below the human depth ($(cell "$V" NC_006273.2 copies_per_cell) copies per cell)" \
  between "$(cell "$V" NC_006273.2 copies_per_cell)" 0 0.05
check_eq "phiX: the lab control" "$(cell "$V" NC_001422.1 call)" "lab control (phiX)"
check_eq "HBV, with nothing planted: not found" "$(cell "$V" NC_003977.2 call)" "not found"
for f in viral_reads.bam viral_reads.bam.bai viral_coverage.tsv human_coverage.tsv; do
  check "the step wrote ${V}_${f}" test -s "${G}/${V}/viral/${V}_${f}"
done
check "the step removed its intermediate files" bash -c '! ls "$1"/*_unmapped.fq "$1"/*_viral.sam "$1"/*_viral_depth.tsv 2>/dev/null' \
  _ "${G}/${V}/viral"

# --- the same genome at a tenth of the human depth --------------------------------------------
L10="${SAMPLE}virlow"
plant "$L10" 0.1 >/dev/null
run_step 47-viral.sh "$L10"
check_step_exit 47-viral.sh
check "control: HHV-6B at a tenth of the depth is below 0.30 ($(cell "$L10" NC_000898.1 relative_depth))" \
  between "$(cell "$L10" NC_000898.1 relative_depth)" 0 0.29
check_eq "control: and is viral DNA found, not integration" "$(cell "$L10" NC_000898.1 call)" "viral DNA found"
check "control: over most of the genome still ($(cell "$L10" NC_000898.1 breadth))" between "$(cell "$L10" NC_000898.1 breadth)" 0.5 1

# --- the text report -------------------------------------------------------------------
run_step generate-report.sh "$V"
check_step_exit generate-report.sh
check "the text report states the HHV-6 call" grep -q 'HHV-6: inherited integrated HHV-6B likely' "${G}/${V}/${V}_report.txt"

# --- Nextflow against the bash step --------------------------------------------------------
command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
SHEET="${T}/samplesheet.csv"
printf 'sample,sex,vcf,vcf_index,bam,bam_index\n%s,male,%s,%s,%s,%s\n' "$V" \
  "${G}/${V}/vcf/${V}.vcf.gz" "${G}/${V}/vcf/${V}.vcf.gz.tbi" \
  "${G}/${V}/aligned/${V}_sorted.bam" "${G}/${V}/aligned/${V}_sorted.bam.bai" > "$SHEET"
NF_OUT="${G}/nf-viral"
echo "+ nextflow run main.nf -profile docker --tools viral"
( cd "$T" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -work-dir "${E2E_WORK}/nf-work-viral" \
    --input "$SHEET" --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" --tools viral --viral_ref "$VREF" \
    --sex_check warn --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB ) 2>&1 | tee "$STEP_LOG"
NF_RC=${PIPESTATUS[0]}
check_eq "nextflow run exits 0" "$NF_RC" 0
cp "${T}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
for f in viral viral_coverage human_coverage; do
  a="${G}/${V}/viral/${V}_${f}.tsv" b="${NF_OUT}/${V}/viral/${V}_${f}.tsv"
  check "Nextflow wrote ${V}_${f}.tsv" test -s "$b"
  check "Nextflow's ${V}_${f}.tsv is the bash step's" cmp -s "$a" "$b"
  cmp -s "$a" "$b" || diff "$a" "$b" | head -20
done
check "Nextflow published the viral reads" test -s "${NF_OUT}/${V}/viral/${V}_viral_reads.bam"
check "Nextflow kept the per-base depth out of the results" test ! -e "${NF_OUT}/${V}/viral/${V}_viral_depth.tsv"
check "control: the comparison sees a difference" bash -c '! cmp -s "$1" "$2"' _ \
  "${G}/${L10}/viral/${L10}_viral.tsv" "${NF_OUT}/${V}/viral/${V}_viral.tsv"

finish
