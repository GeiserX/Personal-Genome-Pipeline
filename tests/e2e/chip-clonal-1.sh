#!/usr/bin/env bash
# Step 49 (clonal haematopoiesis) on HG002's real reads over the CHIP genes.
# The fixture's slices hold none of the CHIP genes, so the case installs
# Ensembl's annotation (setup.sh --chip-data), cuts the coding exons of the
# CHIP genes on the fixture reference's contigs (DNMT3A, SF3B1 and IDH1 on
# chr2, TET2 on chr4, ASXL1 on chr20, GNB1 on chr1) and reads HG002's reads
# there from GIAB's 60x GRCh38 BAM, the file the fixture is cut from, on
# GIAB's S3 mirror. Then:
#   - the reads as they are: no call at DNMT3A R882, a depth row for every gene
#     of the table (the genes on contigs the fixture lacks not_on_reference),
#     and a detection limit under 10% at about 60x;
#   - the same reads with DNMT3A R882H (chr2:25234373 C>T) planted in every
#     fifth read over it: a chip row R882H at a VAF near 20%. The reads as they
#     are are the negative control of that check;
#   - planted in every second read: likely_germline, not chip;
#   - both reports show the section and the planted call;
#   - the Nextflow CHIP_* processes on the planted reads merged into step 02's BAM
#     write the bash step's two tables for the same BAM.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR" T="$CASE_TMP"
REF="${G}/reference/GRCh38_no_alt_analysis_set.fasta"
SRC="https://giab.s3.amazonaws.com/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam"
GFF="${G}/reference/ensembl/Homo_sapiens.GRCh38.${ENSEMBL_GFF_RELEASE}.chr.gff3.gz"
R882=25234373
dk() { local img=$1; shift; docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$img" "$@"; }
retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && { echo "  try ${i} of 4 failed; again in $((i * 15)) s"; sleep $((i * 15)); }; done; return 1; }

# --- the annotation and the exons ------------------------------------------------------
echo "+ scripts/setup.sh --chip-data"
"${REPO}/scripts/setup.sh" --chip-data "$G" 2>&1 | tee "$STEP_LOG"
check_eq "setup.sh --chip-data exits 0" "${PIPESTATUS[0]}" 0
check "the annotation is installed" test -s "$GFF"
python3 "${REPO}/bin/chip_report.py" prepare --gff "$GFF" --genes "${REPO}/assets/chip_genes.tsv" --fai "${REF}.fai" \
  --out-gff "${T}/chip.gff3" --out-bed "${T}/chip.bed"
check_eq "genes with exons on the fixture's contigs" "$(cut -f4 "${T}/chip.bed" | sort -u | tr '\n' ' ')" "ASXL1 DNMT3A GNB1 IDH1 SF3B1 TET2 "
mapfile -t WINDOWS < <(awk -F'\t' '{print $1 ":" ($2 - 300) "-" ($3 + 300)}' "${T}/chip.bed")

# --- HG002's reads, under the fixture reference's header ----------------------------
check "GIAB's BAM index downloads" retry curl -fsS -o "${T}/src.bai" "${SRC}.bai"
fetch_reads() {
  rm -f "${T}/win.sam"
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$SAMTOOLS_IMAGE" \
    samtools view -M -X -o "${T}/win.sam" "$SRC" "${T}/src.bai" "${WINDOWS[@]}"
}
check "HG002's reads over the CHIP exons download" retry fetch_reads
check_ge "reads fetched" "$(wc -l < "${T}/win.sam" 2>/dev/null || echo 0)" 20000

# plant N: the reads as SAM under the fixture's @SQ lines, with C>T at R882 in every
# Nth read whose aligned base there is C (0: none). A mate on another contig loses its
# position, so samtools accepts the record under the shorter header.
plant() {
  python3 - "$1" "${REF}.fai" "${T}/win.sam" "$R882" <<'EOF'
import sys
n, fai, sam, pos = int(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4])
contigs = [ln.split("\t")[:2] for ln in open(fai)]
names = {c for c, _ in contigs}
print("@HD\tVN:1.6\tSO:unsorted")
for c, l in contigs:
    print(f"@SQ\tSN:{c}\tLN:{l}")
print("@RG\tID:giab\tSM:HG002\tPL:ILLUMINA")
k = 0
for ln in open(sam):
    f = ln.rstrip("\n").split("\t")
    if f[2] not in names:
        continue
    f = [x for x in f if not x.startswith("RG:Z:")] + ["RG:Z:giab"]
    if f[6] not in ("=", "*") and f[6] not in names:
        f[6], f[7] = "*", "0"
    if n and f[2] == "chr2" and f[5] != "*":
        ref, q, num, hit = int(f[3]), 0, "", None
        for ch in f[5]:
            if ch.isdigit():
                num += ch
                continue
            ln_ = int(num); num = ""
            if ch in "M=X":
                if ref <= pos < ref + ln_:
                    hit = q + pos - ref
                ref += ln_; q += ln_
            elif ch in "IS":
                q += ln_
            elif ch in "DN":
                ref += ln_
        if hit is not None and f[9][hit] == "C":
            if k % n == 0:
                f[9] = f[9][:hit] + "T" + f[9][hit + 1:]
            k += 1
    print("\t".join(f))
EOF
}
# sample NAME N: a sample folder with the reads, planted every Nth read.
sample() {
  rm -rf "${G:?}/$1"; mkdir -p "${G}/$1/aligned"
  plant "$2" > "${T}/$1.sam"
  dk "$SAMTOOLS_IMAGE" samtools sort -o "${G}/$1/aligned/$1_sorted.bam" "${T}/$1.sam" \
    && dk "$SAMTOOLS_IMAGE" samtools index "${G}/$1/aligned/$1_sorted.bam"
}
# row NAME POS: the chip table's row at POS, as status,protein_change,vaf
row() { awk -F'\t' -v p="$2" '$3 == p {print $13 "," $7 "," $9}' "${G}/$1/chip/$1_chip.tsv" 2>/dev/null; }
# depth NAME GENE COLUMN
depth() {
  awk -F'\t' -v g="$2" -v k="$3" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $1 == g {print $c[k]}' \
    "${G}/$1/chip/$1_chip_depth.tsv" 2>/dev/null
}

U="${SAMPLE}chu" P="${SAMPLE}chp" H="${SAMPLE}chh"
check "the reads as they are" sample "$U" 0
check "R882H in every fifth read" sample "$P" 5
check "R882H in every second read" sample "$H" 2
alt_at() { dk "$SAMTOOLS_IMAGE" samtools view "${G}/$1/aligned/$1_sorted.bam" "chr2:${R882}-${R882}" | python3 -c '
import sys
pos = int(sys.argv[1]); a = n = 0
for ln in sys.stdin:
    f = ln.split("\t"); ref, q, num = int(f[3]), 0, ""
    for ch in f[5]:
        if ch.isdigit(): num += ch; continue
        k = int(num); num = ""
        if ch in "M=X":
            if ref <= pos < ref + k:
                n += 1; a += f[9][q + pos - ref] == "T"
            ref += k; q += k
        elif ch in "IS": q += k
        elif ch in "DN": ref += k
print(a, n)' "$R882"; }
read -r A0 N0 < <(alt_at "$U"); read -r A5 N5 < <(alt_at "$P")
check_ge "reads over R882" "${N0:-0}" 30
check_eq "no T at R882 in HG002's own reads" "${A0:-x}" 0
check "about a fifth of the reads carry the planted T (${A5} of ${N5})" test "$((A5 * 6))" -ge "$N5" -a "$((A5 * 4))" -le "$N5"

# --- the bash step ---------------------------------------------------------------------
for s in "$U" "$P" "$H"; do run_step 49-chip.sh "$s"; check_step_exit 49-chip.sh; done
check_eq "control: the reads as they are have no row at R882" "$(row "$U" "$R882")" ""
check_eq "a depth row per gene of the table" "$(awk 'END {print NR - 1}' "${G}/${U}/chip/${U}_chip_depth.tsv" 2>/dev/null)" \
  "$(grep -v '^#' "${REPO}/assets/chip_genes.tsv" | awk 'END {print NR - 1}')"
check_eq "DNMT3A is covered" "$(depth "$U" DNMT3A status)" covered
check_eq "JAK2 (chr9, not in the fixture) is not_on_reference" "$(depth "$U" JAK2 status)" not_on_reference
check "DNMT3A's mean depth is about 60x ($(depth "$U" DNMT3A mean_depth))" \
  awk -v d="$(depth "$U" DNMT3A mean_depth)" 'BEGIN {exit !(d >= 35 && d <= 90)}'
check "DNMT3A's detection limit is under 10% ($(depth "$U" DNMT3A detection_limit_vaf))" \
  awk -v d="$(depth "$U" DNMT3A detection_limit_vaf)" 'BEGIN {exit !(d > 0.02 && d < 0.10)}'
check "the step prints the detection limit and that it is not a diagnosis" \
  has 'Detection limit at this depth.*does not rule one out\. Not a diagnosis\.' "$(cat "$STEP_LOG")"
IFS=, read -r st pc vaf <<< "$(row "$P" "$R882")"
check_eq "planted at a fifth: a chip row" "${st:-none}" chip
check_eq "named R882H on DNMT3A's MANE transcript" "${pc:-none}" R882H
check "its VAF is near 20% (${vaf:-none})" awk -v v="${vaf:-0}" 'BEGIN {exit !(v >= 0.10 && v <= 0.30)}'
check_eq "planted at a half: likely_germline, not chip" "$(row "$H" "$R882" | cut -d, -f1)" likely_germline
check_eq "and the germline-looking call is not counted" \
  "$(awk -F'\t' '$13 == "chip" && $3 == p' p="$R882" "${G}/${H}/chip/${H}_chip.tsv" 2>/dev/null | wc -l | tr -d ' ')" 0

# --- the reports ----------------------------------------------------------------------
python3 "${REPO}/bin/render_report.py" --sample "$P" --sample-dir "${G}/${P}" -o "${T}/P.txt" -o "${T}/P.html" >/dev/null
check "the text report has the CHIP section with R882H" has 'DNMT3A R882H: VAF 0\.[0-9]+' "$(cat "${T}/P.txt" 2>/dev/null)"
check "and the detection limit" has 'Detection limit: a clone of VAF about [0-9]+% or more' "$(cat "${T}/P.txt" 2>/dev/null)"
check "the HTML report has the CHIP card" has 'Clonal Haematopoiesis \(CHIP\)</h2>' "$(cat "${T}/P.html" 2>/dev/null)"
python3 "${REPO}/bin/render_report.py" --sample "$U" --sample-dir "${G}/${U}" -o "${T}/U.txt" >/dev/null
check "control: the unplanted report says no CHIP call" has 'No CHIP call\.' "$(cat "${T}/U.txt" 2>/dev/null)"

# --- Nextflow against the bash step ----------------------------------------------------
command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
# The planted reads merged into the BAM step 02 wrote from the fixture's reads, so the
# pipeline's sex check (INDEXCOV) reads a whole fixture; the bash step runs on the same BAM.
V="${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" BASE="${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" N="${SAMPLE}chn"
check "the VCF of step 03 and the BAM of step 02 are there" test -s "${V}.tbi" -a -s "${BASE}.bai"
rm -rf "${G:?}/${N}"; mkdir -p "${G}/${N}/aligned"
dk "$SAMTOOLS_IMAGE" samtools merge -f -o "${G}/${N}/aligned/${N}_sorted.bam" "$BASE" "${G}/${P}/aligned/${P}_sorted.bam" \
  && dk "$SAMTOOLS_IMAGE" samtools index "${G}/${N}/aligned/${N}_sorted.bam"
run_step 49-chip.sh "$N"; check_step_exit 49-chip.sh
check_eq "the merged BAM keeps the planted chip row" "$(row "$N" "$R882" | cut -d, -f1,2)" "chip,R882H"
SHEET="${T}/samplesheet.csv"
printf 'sample,sex,vcf,vcf_index,bam,bam_index\n%s,male,%s,%s,%s,%s\n' "$N" "$V" "${V}.tbi" \
  "${G}/${N}/aligned/${N}_sorted.bam" "${G}/${N}/aligned/${N}_sorted.bam.bai" > "$SHEET"
NF_OUT="${G}/nf-chip"
echo "+ nextflow run main.nf -profile docker --tools chip"
( cd "$T" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -work-dir "${E2E_WORK}/nf-work-chip" \
    --input "$SHEET" --reference "$REF" --tools chip --chip_gff "$GFF" \
    --sex_check warn --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB ) 2>&1 | tee "$STEP_LOG"
check_eq "nextflow run exits 0" "${PIPESTATUS[0]}" 0
cp "${T}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
for f in chip chip_depth; do
  a="${G}/${N}/chip/${N}_${f}.tsv" b="${NF_OUT}/${N}/chip/${N}_${f}.tsv"
  check "Nextflow wrote ${N}_${f}.tsv" test -s "$b"
  check "Nextflow's ${N}_${f}.tsv is the bash step's" cmp -s "$a" "$b"
  cmp -s "$a" "$b" || diff "$a" "$b" | head -20
done
check "control: the comparison sees a difference" bash -c '! cmp -s "$1" "$2"' _ \
  "${G}/${U}/chip/${U}_chip.tsv" "${NF_OUT}/${N}/chip/${N}_chip.tsv"

finish
