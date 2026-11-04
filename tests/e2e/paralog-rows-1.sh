#!/usr/bin/env bash
# Every row of step 35's catalog (assets/paralog_loci.tsv, 37 Parascopy
# models) runs on HG002 reads, and every table reaches the reports:
#   1. the BAM is the fixture's slice plus GIAB's 60x HG002 BAM over every
#      region (main and copies) of every catalog model on the fixture's
#      chromosomes, read over HTTPS and sampled to about 30x as the fixture is;
#      step 35 runs the whole catalog on it, timed, and fails past 30 minutes;
#   2. every row's table holds an integer agCN with a numeric quality, or says
#      not assessed with the reason; a header-only table (what a stub writes)
#      fails that check (red control); at least the GBA row is called;
#   3. step 35 with a catalog of the GBA row alone (no SMN1 row), then both
#      reports: the HBA count and the GBA locus are shown, SMN1 says it has no
#      table; control: without the HBA table the HBA line is gone;
#   4. the pinned Paraphase image embeds data/38/config.yaml of its 4.1.0
#      tag byte for byte, and its top-level regions number what
#      docs/35-paralogs.md and tests/test_paraphase_report.py state.
# The stand-ins of paralog-catalog-1-gba apply (an N-filled chr11 in a copy of
# the reference, background depth from 100 bp windows on the chr20 slice). A
# model with a region on a chromosome the fixture lacks is not assessed. The
# copy numbers on such a slice are not a measurement of HG002.
. "$(dirname "$0")/lib.sh"

echo "+ scripts/setup.sh --parascopy-data"
"${REPO}/scripts/setup.sh" --parascopy-data "$GENOME_DIR"
LOCI_ALL="${REPO}/assets/paralog_loci.tsv"

# --- the stand-in reference (as paralog-catalog-1-gba) ----------------------------
REF="${GENOME_DIR}/reference/paralog_rows_ref.fasta"
python3 - "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta" "$REF" <<'PY'
import shutil, sys
with open(sys.argv[2], "wb") as out, open(sys.argv[1], "rb") as src:
    shutil.copyfileobj(src, out)
    out.write(b">chr11\n")
    left = 135086622            # GRCh38 chr11, as the homology table lists it
    while left:
        n = min(left, 60)
        out.write(b"N" * n + b"\n")
        left -= n
PY
sam faidx "reference/$(basename "$REF")"
check "the stand-in reference has chr11" grep -q $'^chr11\t135086622\t' "${REF}.fai"

# --- 1. HG002 reads over every catalog model ----------------------------------------
P=HG002rows
mkdir -p "${GENOME_DIR}/${P}/aligned"
# Every region of the catalog's models on a chromosome of the fixture's reference, 1 kb padded, BED.
BEDR="${GENOME_DIR}/${P}/model_regions.bed"
python3 - "$LOCI_ALL" "${REPO}/assets/parascopy_models_v1.7.txt" "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta.fai" "$BEDR" <<'PY'
import re, sys
loci, models, fai, out = sys.argv[1:]
with open(loci) as f:
    lines = [x.rstrip("\n").split("\t") for x in f if x.strip() and not x.startswith("#")]
want = {r[lines[0].index("model")][:-3] for r in lines[1:]}
contigs = {x.split("\t")[0] for x in open(fai)}
spans = []
for line in open(models):
    if line.startswith("#") or not line.strip():
        continue
    g, r, _ = line.rstrip("\n").split("\t")
    m = re.fullmatch(r"([^:]+):(\d+)-(\d+)", r)
    if g in want and m.group(1) in contigs:
        spans.append((m.group(1), max(0, int(m.group(2)) - 1001), int(m.group(3)) + 1000))
spans.sort()
merged = []
for c, s, e in spans:
    if merged and merged[-1][0] == c and s <= merged[-1][2]:
        merged[-1][2] = max(merged[-1][2], e)
    else:
        merged.append([c, s, e])
with open(out, "w") as f:
    for c, s, e in merged:
        f.write(f"{c}\t{s}\t{e}\n")
print(f"{len(merged)} regions, {sum(e - s for _, s, e in merged)} bp")
PY
check_ge "regions to read from GIAB" "$(grep -c . "$BEDR" 2>/dev/null || echo 0)" 20
GIAB_BAM="https://giab.s3.amazonaws.com/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam"
sam_https() {  # the image's libcurl has no CA bundle, so lend it the host's
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${GENOME_DIR}:/genome" -w "/genome/${P}" \
    "$SAMTOOLS_IMAGE" samtools "$@"
}
retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && sleep $((i * 15)); done; return 1; }
retry sam_https view -b -M -L model_regions.bed -s 42.5 -o model_reads.bam "$GIAB_BAM"
rm -f "${GENOME_DIR}/${P}/$(basename "$GIAB_BAM").bai"
sam sort -o "${P}/model_reads.sorted.bam" "${P}/model_reads.bam"
cp "${FIXTURE_DIR}/${SAMPLE}_slice.bam" "${GENOME_DIR}/${P}/fixture.bam"
sam merge -f -c -p -o "${P}/aligned/${P}_sorted.bam" "${P}/fixture.bam" "${P}/model_reads.sorted.bam"
rm -f "${GENOME_DIR}/${P}/fixture.bam" "${GENOME_DIR}/${P}/model_reads.bam"
sam index "${P}/aligned/${P}_sorted.bam"
check_ge "${P}: HG002 reads at GBA1" "$(sam view -c -F 0x904 "${P}/aligned/${P}_sorted.bam" chr1:155234452-155244699 2>/dev/null)" 1000
check_ge "${P}: HG002 reads at C4A" "$(sam view -c -F 0x904 "${P}/aligned/${P}_sorted.bam" chr6:31980532-32006219 2>/dev/null)" 1000
BED="${GENOME_DIR}/${P}/parascopy_windows.bed"
awk 'BEGIN { for (s = 10050000; s < 10450000; s += 100) printf "chr20\t%d\t%d\n", s, s + 100 }' > "$BED"

T0=$SECONDS
REF_FASTA="$REF" PARASCOPY_DEPTH_BED="$BED" run_step 35-paralogs.sh "$P"
ELAPSED=$((SECONDS - T0))
check_step_exit 35-paralogs.sh
check "step 35 ran the 37-row catalog in 30 minutes or less (${ELAPSED} s)" test "$ELAPSED" -le 1800

# --- 2. every row: an integer agCN with its quality, or not assessed with the reason --
# table_state TABLE: called, not-assessed, not-resolved or bad. Each row is read on its own: an
# integer agCN with a numeric quality, "." with the filter that says why (not
# resolved), or "not assessed: <reason>"; any other row, or no row, makes it bad.
# shellcheck disable=SC2016  # awk fields
table_state() {
  awk -F'\t' 'NR == 1 {next}
    NF != 11 {b = 1; next}
    $6 ~ /^[0-9]+$/ && $7 ~ /^[0-9.]+$/ {c = 1; next}
    $5 ~ /^not assessed: ./ {n = 1; next}
    $6 == "." && $5 != "" && $5 != "." {r = 1; next}
    {b = 1}
    END {print (b || NR < 2 ? "bad" : c ? "called" : n ? "not-assessed" : r ? "not-resolved" : "bad")}' "$1" 2>/dev/null || echo bad
}
# pass_call TABLE: yes when a row has filter PASS and a quality of 20 or more, the
# rule by which the report calls a locus
# shellcheck disable=SC2016  # awk fields
pass_call() {
  awk -F'\t' 'NR > 1 && $5 == "PASS" && $6 ~ /^[0-9]+$/ && $7 + 0 >= 20 {y = 1} END {print (y ? "yes" : "no")}' "$1" 2>/dev/null || echo no
}
OUTD="${GENOME_DIR}/${P}/paralogs"
TABLES=$(awk -F'\t' '/^#/ || !NF {next} !h {for (i = 1; i <= NF; i++) c[$i] = i; h = 1; next} {print $c["table"]}' "$LOCI_ALL")
check_eq "catalog tables" "$(grep -c . <<< "$TABLES")" 37
called=0 nassessed=0 bad=""
for t in $TABLES; do
  st=$(table_state "${OUTD}/${P}_${t}_copy_number.tsv")
  case "$st" in
    called) called=$((called + 1)) ;;
    not-assessed|not-resolved) nassessed=$((nassessed + 1)) ;;
    *) bad="${bad} ${t}" ;;
  esac
  echo "  ${t}: ${st} $(awk -F'\t' 'NR == 2 {print $5, $6, $7}' "${OUTD}/${P}_${t}_copy_number.tsv" 2>/dev/null)"
done
check_eq "every row has an integer agCN or says not assessed with the reason (bad:${bad:- none})" "${bad:-none}" none
check_eq "the GBA row is called" "$(table_state "${OUTD}/${P}_gba1_copy_number.tsv")" called
check_eq "the GBA call has filter PASS and quality 20 or more" "$(pass_call "${OUTD}/${P}_gba1_copy_number.tsv")" yes
HEAD_ONLY="${CASE_TMP}/header_only.tsv"
head -1 "${OUTD}/${P}_smn_copy_number.tsv" > "$HEAD_ONLY"
check_eq "red control: a header-only table fails the row check" "$(table_state "$HEAD_ONLY")" bad
MIXED="${CASE_TMP}/one_bad_row.tsv"
{ head -2 "${OUTD}/${P}_gba1_copy_number.tsv"; printf 'chr1\t1\t2\tX\t\t\t\t\t\t\t\n'; } > "$MIXED"
check_eq "red control: a called row beside a malformed one fails the row check" "$(table_state "$MIXED")" bad
LOWQ="${CASE_TMP}/lowqual.tsv"
awk -F'\t' -v OFS='\t' 'NR > 1 {$5 = "LowQual"; $7 = 12} 1' "${OUTD}/${P}_gba1_copy_number.tsv" > "$LOWQ"
check_eq "red control: a LowQual quality-12 GBA table is no PASS call" "$(pass_call "$LOWQ")" no
echo "paralog-rows: step 35 on 37 rows in ${ELAPSED} s, ${called} called, ${nassessed} not assessed" >> "$E2E_NOTES"

# The reports read every table: both, on this sample
REP="${CASE_TMP}/report_all"
mkdir -p "$REP"
python3 "${REPO}/bin/render_report.py" --assets "${REPO}/assets" --sample "$P" --sample-dir "${GENOME_DIR}/${P}" \
  --json "${REP}/summary.json" -o "${REP}/report.html" > "${REP}/render.log" 2>&1
loci_in_summary() {
  python3 -c "import json, sys; d = json.load(open(sys.argv[1]))['sections']['parascopy']['data']; print(len(d.get('loci', [])))" \
    "${REP}/summary.json" 2>/dev/null || echo 0
}
check_eq "summary.json has the 36 loci other than SMN1" "$(loci_in_summary)" 36
check "the HTML report lists the GBA locus" grep -q '<td>GBA1, GBAP1</td>' "${REP}/report.html"

# --- 3. a catalog without the SMN1 row ----------------------------------------------
Q=HG002nosmn
mkdir -p "${GENOME_DIR}/${Q}/aligned"
# A copy, not a hard link: the first sample's BAM belongs to the container's user, and the
# kernel refuses a hard link to another user's file (protected_hardlinks).
cp "${GENOME_DIR}/${P}/aligned/${P}_sorted.bam" "${GENOME_DIR}/${Q}/aligned/${Q}_sorted.bam"
cp "${GENOME_DIR}/${P}/aligned/${P}_sorted.bam.bai" "${GENOME_DIR}/${Q}/aligned/${Q}_sorted.bam.bai"
check "the second sample has the BAM" test -s "${GENOME_DIR}/${Q}/aligned/${Q}_sorted.bam"
LOCI_GBA="${GENOME_DIR}/${Q}/paralog_loci.tsv"
awk -F'\t' '/^#/ || !NF {next} !h {print; h = 1; next} $1 == "GBA.gz"' "$LOCI_ALL" > "$LOCI_GBA"
check_eq "the override catalog: header and the GBA row" "$(grep -c . "$LOCI_GBA")" 2
REF_FASTA="$REF" PARASCOPY_DEPTH_BED="$BED" PARALOG_LOCI="$LOCI_GBA" run_step 35-paralogs.sh "$Q"
check_step_exit 35-paralogs.sh
check "no SMN table without the SMN1 row" test ! -e "${GENOME_DIR}/${Q}/paralogs/${Q}_smn_copy_number.tsv"
check "the HBA table is written" test -s "${GENOME_DIR}/${Q}/paralogs/${Q}_hba_copy_number.tsv"
render_q() {
  rm -rf "$1"; mkdir -p "$1"
  python3 "${REPO}/bin/render_report.py" --assets "${REPO}/assets" --sample "$Q" --sample-dir "${GENOME_DIR}/${Q}" \
    --json "$1/summary.json" -o "$1/report.html" > "$1/render.log" 2>&1
}
render_q "${CASE_TMP}/report_nosmn"
H=$(cat "${CASE_TMP}/report_nosmn/report.html" 2>/dev/null)
check "no SMN1 row: the HTML report has the HBA line" has 'HBA1 and HBA2 together \(alpha-globin genes, from depth\)' "$H"
check "no SMN1 row: the HTML report has the GBA locus" has '<td>GBA1, GBAP1</td>' "$H"
check "no SMN1 row: SMN1 says it has no table" has "no SMN1 table \(step 35 wrote no SMN copy-number table\)" "$H"
mv "${GENOME_DIR}/${Q}/paralogs/${Q}_hba_copy_number.tsv" "${CASE_TMP}/hba.tsv"
check "control: the HBA table is set aside" test -s "${CASE_TMP}/hba.tsv"
render_q "${CASE_TMP}/report_nohba"
check "control: without the HBA table the HBA line is gone" \
  lacks 'HBA1 and HBA2 together \(alpha-globin genes, from depth\)' "$(cat "${CASE_TMP}/report_nohba/report.html" 2>/dev/null || echo 'HBA1 and HBA2 together (alpha-globin genes, from depth)')"
mv "${CASE_TMP}/hba.tsv" "${GENOME_DIR}/${Q}/paralogs/${Q}_hba_copy_number.tsv"

# --- 4. Paraphase's region count, from the pinned image -------------------------------
CFG="${CASE_TMP}/paraphase_config.yaml"
retry curl -fsSL -o "$CFG" "https://raw.githubusercontent.com/PacificBiosciences/paraphase/v4.1.0/data/38/config.yaml"
check_eq "Paraphase 4.1.0's data/38/config.yaml (sha256)" "$( (sha256sum "$CFG" 2>/dev/null || shasum -a 256 "$CFG") | cut -d' ' -f1)" \
  a08ffb59c78523b1c8e32a933a021598a939d65bcee81478aa356373386be5d1
BIN="${CASE_TMP}/paraphase.bin"
cid=$(docker create "$PARAPHASE_IMAGE" true)
pbin=$(docker run --rm --entrypoint sh "$PARAPHASE_IMAGE" -c 'command -v paraphase' 2>/dev/null)
docker cp -L "${cid}:${pbin:-/usr/local/bin/paraphase}" "$BIN" >/dev/null 2>&1
docker rm "$cid" >/dev/null 2>&1
embedded() { python3 -c "import sys; sys.exit(0 if open(sys.argv[2], 'rb').read() in open(sys.argv[1], 'rb').read() else 1)" "$1" "$2"; }
check "the pinned image's paraphase embeds that config byte for byte" embedded "$BIN" "$CFG"
printf 'smn1:\n  genes: X\n' > "${CASE_TMP}/other.yaml"
check "control: another config is not found in the binary" sh -c "! python3 -c \"import sys; sys.exit(0 if open(sys.argv[2], 'rb').read() in open(sys.argv[1], 'rb').read() else 1)\" '$BIN' '${CASE_TMP}/other.yaml'"
N=$(python3 -c "import re, sys; print(len(re.findall(r'^[A-Za-z0-9_]+:[ \t]*$', open(sys.argv[1]).read(), re.M)))" "$CFG")
check_eq "regions in the pinned config" "$N" "$(sed -n 's/^PARAPHASE_REGIONS = \([0-9]*\)$/\1/p' "${REPO}/tests/test_paraphase_report.py")"
check "docs/35-paralogs.md states that count" grep -q "It covers ${N} GRCh38 regions" "${REPO}/docs/35-paralogs.md"
rm -f "$BIN" "$REF" "${REF}.fai"

finish
