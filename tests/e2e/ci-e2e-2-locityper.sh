#!/usr/bin/env bash
# Step 53 (Locityper) and the database it reads, the shipped scripts as a
# user runs them: setup.sh --locityper-data builds the database, then
# scripts/53-locityper.sh genotypes HG002 and writes the table
# bin/locityper_report.py makes. The image's smoke test repeats their
# commands; this case runs the scripts.
#
# The fixture holds reads over the MHC but none over FCGR or KIR, so the
# database is built for the HLA-DQB1 row of assets/locityper_loci.bed
# (LOCITYPER_LOCI). HLA-DRB1, the other MHC row, is left out for time: with
# both, jellyfish took 35 minutes and the pairwise alignments of
# `locityper augment` 23 more (fork run 42301973505). With HLA-DQB1 alone
# augment takes about a minute, but jellyfish, whose time does not depend on
# the loci, still takes 36 (38 minutes in all, fork run 42325306077), so the
# case keeps a limit of 90 minutes instead of the runner's 60 (the
# case-timeout line below). It runs in a GENOME_DIR of its own with a
# reference of chr6 and chr20 cut from the fixture's, so jellyfish
# counts 235 Mb instead of the fixture's 1.8 Gb; the fixture's reads over the
# MHC and chr20:10.0-10.5 Mb are aligned to it, since Locityper needs the
# BAM's contigs to be the reference's. Step 53 measures depth on chr20
# (LOCITYPER_BG_REGION): the fixture has no chr17, where Locityper's own
# background region is.
#
# Red controls: the step stops without a database, without the BAM's index,
# and with a database whose DONE names another reference.
# case-timeout: 5400
. "$(dirname "$0")/lib.sh"

G="${CASE_TMP}/genome"
S="${G}/${SAMPLE}"
rm -rf "$G"
mkdir -p "${G}/reference" "${S}/aligned"
U="$(id -u):$(id -g)"

# --- the reference and the reads -------------------------------------------------
echo "+ a chr6 and chr20 reference cut from the fixture's, and the reads over the MHC and chr20 aligned to it"
in_genome "$SAMTOOLS_IMAGE" samtools faidx reference/GRCh38_no_alt_analysis_set.fasta chr6 chr20 \
  > "${G}/reference/GRCh38_no_alt_analysis_set.fasta"
docker run --rm --user "$U" -v "${G}:/g" "$SAMTOOLS_IMAGE" samtools faidx /g/reference/GRCh38_no_alt_analysis_set.fasta
check_eq "the cut reference has chr6 and chr20" "$(cut -f1 "${G}/reference/GRCh38_no_alt_analysis_set.fasta.fai" 2>/dev/null | paste -sd' ' -)" "chr6 chr20"
docker run --rm --user "$U" -v "${FIXTURE_DIR}:/in:ro" -v "${G}:/g" "$SAMTOOLS_IMAGE" bash -c '
    samtools view -u -F 0x900 /in/HG002_slice.bam chr6:31500001-33100000 chr20:10000001-10500000 \
      | samtools collate -O -u - /g/collate | samtools fastq -n -1 /g/r1.fq.gz -2 /g/r2.fq.gz -0 /dev/null -s /dev/null -'
docker run --rm --user "$U" -v "${G}:/g" "$MINIMAP2_IMAGE" \
    minimap2 -t "${THREADS:-4}" -a -x sr -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" /g/reference/GRCh38_no_alt_analysis_set.fasta /g/r1.fq.gz /g/r2.fq.gz \
  | docker run --rm -i --user "$U" -v "${G}:/g" "$SAMTOOLS_IMAGE" bash -c \
    "samtools sort -o /g/${SAMPLE}/aligned/${SAMPLE}_sorted.bam - && samtools index /g/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
rm -f "${G}/r1.fq.gz" "${G}/r2.fq.gz"
check_ge "reads aligned to the cut reference" \
  "$(docker run --rm -v "${G}:/g:ro" "$SAMTOOLS_IMAGE" samtools view -c -F 0x904 "/g/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" 2>/dev/null)" 10000

LOCI="${CASE_TMP}/loci.bed"
awk -F'\t' '!/^#/ && $4 == "HLA-DQB1"' "${REPO}/assets/locityper_loci.bed" > "$LOCI"
check_eq "the HLA-DQB1 row of assets/locityper_loci.bed" "$(cut -f1,4 "$LOCI" | tr '\t' ' ')" "chr6 HLA-DQB1"

# run53 LOG: step 53 on HG002 in G, its exit code in STEP_RC.
run53() {
  echo "+ scripts/53-locityper.sh ${SAMPLE} (GENOME_DIR=${G})"
  GENOME_DIR="$G" LOCITYPER_BG_REGION=chr20:10000001-10500000 THREADS="${THREADS:-4}" \
    "${REPO}/scripts/53-locityper.sh" "$SAMPLE" > "$1" 2>&1
  STEP_RC=$?
  cat "$1"
}

# --- red control: no database yet ------------------------------------------------
run53 "${CASE_TMP}/nodb.log"
check_eq "red control: step 53 stops (exit 1) before the database is built" "$STEP_RC" 1
check "red control: and says how to build it" grep -q "Locityper's database not found: .*setup.sh --locityper-data" "${CASE_TMP}/nodb.log"

# --- setup.sh --locityper-data ---------------------------------------------------
DB="${G}/reference/locityper-${LOCITYPER_PANEL_VERSION}"
echo "+ LOCITYPER_LOCI=${LOCI} scripts/setup.sh --locityper-data ${G}"
GENOME_DIR="$G" LOCITYPER_LOCI="$LOCI" THREADS="${THREADS:-4}" "${REPO}/scripts/setup.sh" --locityper-data "$G" 2>&1 | tee "$STEP_LOG"
check_eq "setup.sh --locityper-data exits 0" "${PIPESTATUS[0]}" 0
check "and says the database is built" grep -q "\[OK\] Locityper's database (step 53, 1 loci" "$STEP_LOG"
check "the database holds HLA-DQB1 (Locityper's success file)" test -e "${DB}/db/loci/HLA-DQB1/success"
check_eq "its loci.bed is the HLA-DQB1 row" "$(cut -f4 "${DB}/loci.bed" 2>/dev/null | paste -sd' ' -)" "HLA-DQB1"
check "its DONE names the reference's .fai" grep -q "fai=$(sha256sum "${G}/reference/GRCh38_no_alt_analysis_set.fasta.fai" | cut -d' ' -f1) " "${DB}/DONE"
GENOME_DIR="$G" LOCITYPER_LOCI="$LOCI" THREADS="${THREADS:-4}" "${REPO}/scripts/setup.sh" --locityper-data "$G" > "${CASE_TMP}/again.log" 2>&1
check_eq "a second setup.sh --locityper-data exits 0" "$?" 0
check "and builds nothing again" grep -q "already built" "${CASE_TMP}/again.log"

# --- step 53 ---------------------------------------------------------------------
run53 "${CASE_TMP}/step.log"
check_eq "scripts/53-locityper.sh exits 0" "$STEP_RC" 0
TABLE="${S}/locityper/${SAMPLE}_locityper.tsv"
check_eq "the table's columns are those of bin/locityper_report.py" "$(head -n 1 "$TABLE" 2>/dev/null | tr '\t' ' ')" \
  "locus region status genotype haplotype_1 haplotype_2 quality total_reads unexplained_reads weighted_distance warnings filter reason"
check_eq "one row per locus of the database" "$(awk -F'\t' 'NR > 1 {print $1}' "$TABLE" 2>/dev/null | paste -sd' ' -)" "HLA-DQB1"
l=HLA-DQB1
row=$(awk -F'\t' -v l="$l" '$1 == l' "$TABLE" 2>/dev/null)
check_eq "${l}: called" "$(cut -f3 <<< "$row")" called
check "${l}: two panel haplotypes" has '^[^,.][^,]*,[^,.][^,]*$' "$(cut -f4 <<< "$row")"
check "${l}: a numeric quality" has '^[0-9]+(\.[0-9]+)?$' "$(cut -f7 <<< "$row")"
check_ge "${l}: reads used" "$(cut -f8 <<< "$row")" 1
check "${l}: the authors' filter is pass or fail" has '^(pass|fail)$' "$(cut -f12 <<< "$row")"
check "${l}: Locityper's result file" test -s "${S}/locityper/work/gt/loci/${l}/res.json.gz"
check "Locityper's log beside the table" test -s "${S}/locityper/${SAMPLE}_locityper.log"
echo "${CASE_NAME}: $(awk -F'\t' 'NR > 1 {printf "%s %s q%s %s; ", $1, $4, $7, $12}' "$TABLE" 2>/dev/null)" >> "$E2E_NOTES"

# --- red controls: no BAM index, a database of another reference -------------------
mv "${S}/aligned/${SAMPLE}_sorted.bam.bai" "${CASE_TMP}/bai"
run53 "${CASE_TMP}/nobai.log"
check_eq "red control: step 53 stops (exit 1) without the BAM's index" "$STEP_RC" 1
check "red control: and names it" grep -q "File not found: .*${SAMPLE}_sorted.bam.bai" "${CASE_TMP}/nobai.log"
mv "${CASE_TMP}/bai" "${S}/aligned/${SAMPLE}_sorted.bam.bai"
sed -i 's/fai=[0-9a-f]* /fai=0000 /' "${DB}/DONE"
run53 "${CASE_TMP}/stale.log"
check_eq "red control: step 53 stops (exit 1) on a database built from another reference" "$STEP_RC" 1
check "red control: and says so" grep -q "was built from another reference" "${CASE_TMP}/stale.log"

rm -rf "${G:?}"
finish
