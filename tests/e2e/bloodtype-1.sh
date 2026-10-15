#!/usr/bin/env bash
# Step 41 (blood groups) on HG002, with real reads and GIAB's calls.
# The fixture's slices hold none of the blood group genes, so the case reads
# HG002's reads at every marker site of assets/bloodtype_markers.tsv (500 bp
# each side) and over the RH locus (chr1:25.15-25.45 Mb) from GIAB's 60x
# GRCh38 BAM, the file the fixture is cut from, on GIAB's S3 mirror; and it
# takes HG002's genotypes at those sites from GIAB's v4.2.1 benchmark VCF,
# which the fixture holds whole. What is expected comes from those published
# calls and reads, not from a serological typing of HG002: no c.261 insertion
# in ABO, so both chromosomes carry c.261delG (group O); two RHD copies by
# depth (D+); Fya and Fyb (Fy(a+b+)); no c.578C>T (K-k+); Jkb on both (Jk(a-b+)).
# Then:
#   - planted variants flip Duffy (the GATA box with one Fyb: Fy(a+b-)) and
#     Kell (c.578C>T: K+k+);
#   - the same reads without any RHD read (both copies deleted) give D-, and
#     with half of them (one copy) D+ with one copy;
#   - the VCF without a BAM gives RhD not determined and ABO O at medium;
#   - the Nextflow processes (BLOODTYPE_DEPTH, BLOODTYPE) on a VCF and BAM
#     write the same three files as the bash step on the same two files.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR" T="$CASE_TMP"
MARKERS="${REPO}/assets/bloodtype_markers.tsv"
SRC="https://giab.s3.amazonaws.com/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam"
TRUTH="${FIXTURE_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
# dk IMAGE CMD...: a pinned image with the work folder at its own path.
dk() { local img=$1; shift; docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$img" "$@"; }
# retry CMD...: GIAB's S3 mirror sometimes drops a stream half way; the same call then works.
retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && { echo "  try ${i} of 4 failed; again in $((i * 15)) s"; sleep $((i * 15)); }; done; return 1; }

# --- HG002's reads and calls at the markers ------------------------------------------
mapfile -t WINDOWS < <(awk -F'\t' '!/^#/ && $4 == "site" {print $1 ":" ($2 - 500) "-" ($3 + 500)}' "$MARKERS")
check_ge "marker sites in the table" "${#WINDOWS[@]}" 25
check "GIAB's BAM index downloads" retry curl -fsS -o "${T}/src.bai" "${SRC}.bai"
# SAMTOOLS_IMAGE's libcurl has no CA bundle: lend it the host's (htslib reads CURL_CA_BUNDLE).
fetch_reads() {
  rm -f "${T}/win.bam"
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$SAMTOOLS_IMAGE" \
    samtools view -M -b -X -o "${T}/win.bam" "$SRC" "${T}/src.bai" "${WINDOWS[@]}" chr1:25150000-25450000
}
check "HG002's reads at the markers and over RH download" retry fetch_reads
dk "$SAMTOOLS_IMAGE" samtools index "${T}/win.bam"
check_ge "reads fetched" "$(dk "$SAMTOOLS_IMAGE" samtools view -c "${T}/win.bam" 2>/dev/null || echo 0)" 50000

awk -F'\t' 'BEGIN {OFS = "\t"} !/^#/ && $4 == "site" {print $1, $2 - 50, $3 + 50}' "$MARKERS" | sort -k1,1 -k2,2n > "${T}/sites.bed"
dk "$BCFTOOLS_IMAGE" bcftools view -R "${T}/sites.bed" -Oz -o "${T}/truth.vcf.gz" "$TRUTH"
check_ge "GIAB records near the marker sites" "$(dk "$BCFTOOLS_IMAGE" bcftools view -H "${T}/truth.vcf.gz" 2>/dev/null | wc -l)" 5

# The planted copy: the GATA box (ACKR1 c.-67T>C) and KEL c.578C>T, each on one chromosome.
{ dk "$BCFTOOLS_IMAGE" bcftools view "${T}/truth.vcf.gz"
  printf 'chr1\t159204893\t.\tT\tC\t50\tPASS\t.\tGT\t0/1\nchr7\t142957921\t.\tG\tA\t50\tPASS\t.\tGT\t0/1\n'
} > "${T}/planted.vcf"
dk "$BCFTOOLS_IMAGE" bcftools sort -Oz -o "${T}/planted.vcf.gz" "${T}/planted.vcf"

# The same reads with no read in RHD (both copies deleted), and with half of them (one copy).
printf 'chr1\t25272392\t25330445\n' > "${T}/rhd.bed"
dk "$SAMTOOLS_IMAGE" samtools view -b -L "${T}/rhd.bed" -U "${T}/del.bam" -o /dev/null "${T}/win.bam"
dk "$SAMTOOLS_IMAGE" samtools view -b -L "${T}/rhd.bed" -s 7.5 -o "${T}/half.bam" "${T}/win.bam"
dk "$SAMTOOLS_IMAGE" samtools merge -f -o "${T}/het.bam" "${T}/del.bam" "${T}/half.bam"
for b in del het; do dk "$SAMTOOLS_IMAGE" samtools index "${T}/${b}.bam"; done
check_eq "reads left in RHD after the deletion" \
  "$(dk "$SAMTOOLS_IMAGE" samtools view -c "${T}/del.bam" chr1:25272393-25330445 2>/dev/null)" 0

# sample NAME VCF [BAM]: a sample folder with that VCF and BAM.
sample() {
  rm -rf "${G:?}/$1"; mkdir -p "${G}/$1/vcf" "${G}/$1/aligned"
  cp "$2" "${G}/$1/vcf/$1.vcf.gz"
  [ -z "${3:-}" ] || { cp "$3" "${G}/$1/aligned/$1_sorted.bam"; cp "$3.bai" "${G}/$1/aligned/$1_sorted.bam.bai"; }
}
# cell NAME SYSTEM COLUMN: one cell of the sample's table.
cell() {
  awk -F'\t' -v s="$2" -v k="$3" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $1 == s {print $c[k]}' \
    "${G}/$1/bloodtype/$1_bloodtype.tsv" 2>/dev/null
}
step() { run_step 41-bloodtype.sh "$1"; check_step_exit 41-bloodtype.sh; }

# --- HG002 -------------------------------------------------------------------------
H="${SAMPLE}bt"
sample "$H" "${T}/truth.vcf.gz" "${T}/win.bam"
step "$H"
cat "${G}/${H}/bloodtype/${H}_bloodtype.tsv" 2>/dev/null
check_eq "a row per system" "$(awk 'END {print NR - 1}' "${G}/${H}/bloodtype/${H}_bloodtype.tsv" 2>/dev/null)" 17
check_eq "ABO: GIAB has no c.261 insertion, both alleles O.01" "$(cell "$H" ABO phenotype)/$(cell "$H" ABO genotype)" "O/O.01/O.01"
check_eq "RhD: two RHD copies by depth" "$(cell "$H" RhD phenotype)" "D+"
check "RhD genotype RHD/RHD" has '^RHD/RHD \(RHD copies by read depth (1\.[5-9]|2\.[0-4])' "$(cell "$H" RhD genotype)"
check_eq "Duffy: GIAB's c.125G>A 0/1" "$(cell "$H" Duffy phenotype)" "Fy(a+b+)"
check "Kell: no c.578C>T in GIAB's calls" has '^K-k\+ ' "$(cell "$H" Kell phenotype)"
check_eq "Kidd: GIAB's c.838G>A 1/1" "$(cell "$H" Kidd phenotype)" "Jk(a-b+)"
check_eq "ABO confidence (every site has reads)" "$(cell "$H" ABO confidence)" "high"
check_eq "the evidence file has a row per site and per RhD region" \
  "$(awk 'END {print NR - 1}' "${G}/${H}/bloodtype/${H}_bloodtype_evidence.tsv" 2>/dev/null)" \
  "$(grep -v '^#' "$MARKERS" | awk 'END {print NR - 1}')"

# --- planted -------------------------------------------------------------------------
P="${SAMPLE}btp"
sample "$P" "${T}/planted.vcf.gz" "${T}/win.bam"
step "$P"
check_eq "Duffy flips: the GATA box silences the one Fyb" "$(cell "$P" Duffy phenotype)" "Fy(a+b-)"
check_eq "Duffy genotype" "$(cell "$P" Duffy genotype)" "FY*A/FY*B GATA"
check "Kell flips: c.578C>T 0/1 is K+k+" has '^K\+k\+ ' "$(cell "$P" Kell phenotype)"
check_eq "ABO is unchanged by the planted sites" "$(cell "$P" ABO phenotype)" "O"

# --- RHD deleted -------------------------------------------------------------------
D="${SAMPLE}btd"
sample "$D" "${T}/truth.vcf.gz" "${T}/del.bam"
step "$D"
check_eq "RhD with no RHD read: D-" "$(cell "$D" RhD phenotype)" "D-"
check "RhD genotype: both deleted" has '^RHD deleted/RHD deleted ' "$(cell "$D" RhD genotype)"
E="${SAMPLE}bth"
sample "$E" "${T}/truth.vcf.gz" "${T}/het.bam"
step "$E"
check_eq "RhD with half the RHD reads: D+" "$(cell "$E" RhD phenotype)" "D+"
check "RhD genotype: one copy deleted" has '^RHD/RHD deleted ' "$(cell "$E" RhD genotype)"

# --- no BAM ------------------------------------------------------------------------
V="${SAMPLE}btv"
sample "$V" "${T}/truth.vcf.gz"
step "$V"
check_eq "without a BAM RhD is not determined" "$(cell "$V" RhD phenotype)" "not determined"
check "and the step says why" has 'No BAM .*RhD needs read depth' "$(cat "$STEP_LOG")"
check_eq "without a BAM ABO is O at medium" "$(cell "$V" ABO phenotype)/$(cell "$V" ABO confidence)" "O/medium"

# --- Nextflow against the bash step ------------------------------------------------------
# The same VCF, and the BAM step 02 wrote from the fixture's reads (its header is the
# fixture reference's, which INDEXCOV checks). That BAM has no chr9, chr7, chr17 or
# chr18 contig: the step leaves those rows out of bedcov, which stops on them, and
# the sites there with no record are not determined. Every file must be the bash step's.
command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
N="${SAMPLE}btn" BASE="${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
check "the BAM of step 02 is there" test -s "${BASE}.bai"
sample "$N" "${T}/truth.vcf.gz" "$BASE"
dk "$BCFTOOLS_IMAGE" bcftools index -t -f "${G}/${N}/vcf/${N}.vcf.gz"
step "$N"
check "a contig the BAM lacks (chr9): ABO is not determined" has '^not determined$' "$(cell "$N" ABO phenotype)"
check "and the note says the BAM has no such contig" has 'no depth row \(the BAM has no such contig\)' "$(cell "$N" ABO note)"
SHEET="${T}/samplesheet.csv"
printf 'sample,sex,vcf,vcf_index,bam,bam_index\n%s,male,%s,%s,%s,%s\n' "$N" \
  "${G}/${N}/vcf/${N}.vcf.gz" "${G}/${N}/vcf/${N}.vcf.gz.tbi" \
  "${G}/${N}/aligned/${N}_sorted.bam" "${G}/${N}/aligned/${N}_sorted.bam.bai" > "$SHEET"
NF_OUT="${G}/nf-bloodtype"
echo "+ nextflow run main.nf -profile docker --tools bloodtype"
( cd "$T" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -work-dir "${E2E_WORK}/nf-work-bloodtype" \
    --input "$SHEET" --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" --tools bloodtype \
    --sex_check warn --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB ) 2>&1 | tee "$STEP_LOG"
NF_RC=${PIPESTATUS[0]}
check_eq "nextflow run exits 0" "$NF_RC" 0
cp "${T}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
for f in bloodtype bloodtype_evidence bloodtype_depth; do
  a="${G}/${N}/bloodtype/${N}_${f}.tsv" b="${NF_OUT}/${N}/bloodtype/${N}_${f}.tsv"
  check "Nextflow wrote ${N}_${f}.tsv" test -s "$b"
  check "Nextflow's ${N}_${f}.tsv is the bash step's" cmp -s "$a" "$b"
  cmp -s "$a" "$b" || diff "$a" "$b" | head -20
done
check "control: the comparison sees a difference" bash -c '! cmp -s "$1" "$2"' _ \
  "${G}/${H}/bloodtype/${H}_bloodtype.tsv" "${NF_OUT}/${N}/bloodtype/${N}_bloodtype.tsv"

finish
