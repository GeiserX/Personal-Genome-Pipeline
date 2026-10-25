#!/usr/bin/env bash
# Step 45 (Y-STR repeat counts with HipSTR) on HG002, a man, with real reads.
# The fixture's chrY slices hold 11 of the panel's regions, so the case reads
# HG002's reads at every region of assets/ystr_panel.tsv (1 kb each side) from
# GIAB's 60x GRCh38 BAM, the file the fixture is cut from. Then:
#   - the table's regions are repeats in the reference: each region's sequence,
#     read from the fixture FASTA, repeats with the panel's motif length;
#   - the step writes 111 rows in FamilyTreeDNA's order, the multi-copy and
#     absent markers not callable, and most single regions called high;
#   - an independent count from the reads: for every high call, the length of
#     the repeat in the reads that span it with 10 bases on each side
#     (reference length plus the CIGAR's insertions minus its deletions there
#     or in those 10 bases, where an aligner can left-align a repeat's indel)
#     gives the same count in at least 90% of the markers;
#   - a planted change: the same reads with four bases (one AGAT repeat) cut
#     out of DYS393 in every read that spans it give DYS393 one repeat fewer,
#     and leave the other called markers as they were;
#   - a copy whose indexcov row says female is skipped;
#   - both reports carry the section;
#   - the Nextflow process (Y_STR, with html_report) on the fixture's own BAM
#     writes the same table and CSV as the bash step on that BAM.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR" T="$CASE_TMP"
PANEL="${REPO}/assets/ystr_panel.tsv"
REF="${G}/reference/GRCh38_no_alt_analysis_set.fasta"
SRC="https://giab.s3.amazonaws.com/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam"
# dk IMAGE CMD...: a pinned image with the work folder at its own path.
dk() { local img=$1; shift; docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${E2E_WORK}:${E2E_WORK}" -v "${G}:${G}" -w "$T" "$img" "$@"; }
# retry CMD...: GIAB's S3 mirror sometimes drops a stream half way; the same call then works.
retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && { echo "  try ${i} of 4 failed; again in $((i * 15)) s"; sleep $((i * 15)); }; done; return 1; }

# --- the panel against the reference ---------------------------------------------------
python3 "${REPO}/bin/ystr_report.py" bed --panel "$PANEL" --out "${T}/regions.bed"
check_eq "HipSTR regions in the panel" "$(grep -c . "${T}/regions.bed")" 93
awk -F'\t' '{print $1 ":" $2 "-" $3}' "${T}/regions.bed" > "${T}/regions.txt"
dk "$SAMTOOLS_IMAGE" samtools faidx -r "${T}/regions.txt" "$REF" > "${T}/regions.fa"
# periodic REGIONS.BED REGIONS.FA: how many regions repeat with their motif length
# (at least 80% of the bases equal to the base one motif further on).
periodic() {
  python3 -I - "$1" "$2" <<'PY'
import sys
bed = [l.split("\t") for l in open(sys.argv[1])]
seqs, name = {}, None
for line in open(sys.argv[2]):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:]
        seqs[name] = ""
    elif name:
        seqs[name] += line.upper()
ok = 0
for c in bed:
    s, p = seqs.get(f"{c[0]}:{c[1]}-{c[2]}", ""), int(c[3])
    if len(s) == int(c[2]) - int(c[1]) + 1 and len(s) > p:
        same = sum(1 for i in range(len(s) - p) if s[i] == s[i + p])
        ok += same >= 0.8 * (len(s) - p)
print(ok)
PY
}
check_eq "every region is a repeat of its motif length in the reference" "$(periodic "${T}/regions.bed" "${T}/regions.fa")" 93
awk -F'\t' 'BEGIN {OFS = "\t"} {print $1, $2, $3, 1}' "${T}/regions.bed" > "${T}/shifted.bed"
check_eq "control: the same check with a motif length of 1 finds no region" \
  "$(periodic "${T}/shifted.bed" "${T}/regions.fa")" 0

# --- HG002's reads at the regions --------------------------------------------------------
mapfile -t WINDOWS < <(awk -F'\t' '{print $1 ":" ($2 - 1000) "-" ($3 + 1000)}' "${T}/regions.bed")
check "GIAB's BAM index downloads" retry curl -fsS -o "${T}/src.bai" "${SRC}.bai"
# SAMTOOLS_IMAGE's libcurl has no CA bundle: lend it the host's (htslib reads CURL_CA_BUNDLE).
fetch_reads() {
  rm -f "${T}/win.bam"
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$SAMTOOLS_IMAGE" \
    samtools view -M -b -X -o "${T}/win.bam" "$SRC" "${T}/src.bai" "${WINDOWS[@]}"
}
check "HG002's reads at the Y-STR regions download" retry fetch_reads
dk "$SAMTOOLS_IMAGE" samtools index "${T}/win.bam"
check_ge "reads fetched" "$(dk "$SAMTOOLS_IMAGE" samtools view -c "${T}/win.bam" 2>/dev/null || echo 0)" 20000

PED_HEAD=$'#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\tCNchrX\tCNchrY'
# sample NAME BAM SEX_CODE: a sample folder with that BAM and an indexcov row of that sex.
sample() {
  rm -rf "${G:?}/$1"; mkdir -p "${G}/$1/aligned" "${G}/$1/indexcov"
  cp "$2" "${G}/$1/aligned/$1_sorted.bam"; cp "$2.bai" "${G}/$1/aligned/$1_sorted.bam.bai"
  printf '%s\n%s\t%s\t-9\t-9\t%s\t-9\t1.0\t1.0\n' "$PED_HEAD" "$1" "$1" "$3" > "${G}/$1/indexcov/indexcov-indexcov.ped"
}
# col NAME MARKER COLUMN: one cell of the sample's table.
col() {
  awk -F'\t' -v m="$2" -v k="$3" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $3 == m {print $c[k]}' \
    "${G}/$1/y_str/$1_y_str.tsv" 2>/dev/null
}
step() { run_step 45-y-str.sh "$1"; check_step_exit 45-y-str.sh; }

# --- HG002 -------------------------------------------------------------------------------
H="${SAMPLE}ys"
sample "$H" "${T}/win.bam" 1
SECONDS=0
step "$H"
echo "step 45 on HG002's 60x reads at the regions: ${SECONDS} s"
OUT="${G}/${H}/y_str/${H}_y_str.tsv"
cat "$OUT" 2>/dev/null
check_eq "a row per marker" "$(awk 'END {print NR - 1}' "$OUT" 2>/dev/null)" 111
check_eq "rows in FamilyTreeDNA's order (the panel's)" \
  "$(awk -F'\t' 'NR > 1 {print $3}' "$OUT" 2>/dev/null | md5sum)" "$(awk -F'\t' '!/^#/ && $1 != "order" {print $3}' "$PANEL" | md5sum)"
check_eq "the first twelve are Y-12" "$(awk -F'\t' 'NR > 1 && NR <= 13 {print $2}' "$OUT" 2>/dev/null | sort -u)" Y-12
check_eq "multi-copy and absent markers are not callable" \
  "$(awk -F'\t' 'NR > 1 && $5 == "not callable" && $10 ~ /multi-copy|no region/' "$OUT" 2>/dev/null | wc -l | tr -d ' ')" 27
check "DYS464a: not callable from short reads" has '^not callable from short reads: a multi-copy' "$(col "$H" DYS464a note)"
HIGH=$(awk -F'\t' 'NR > 1 && $5 == "high"' "$OUT" 2>/dev/null | wc -l | tr -d ' ')
check_ge "markers called with high confidence at 60x (of 83 callable)" "$HIGH" 60
check "the step prints the tally" has '^Y-STR markers: [0-9]+ high, [0-9]+ low, [0-9]+ no call, [0-9]+ not callable$' "$(cat "$STEP_LOG")"
CSV="${G}/${H}/y_str/${H}_y_str.csv"
check_eq "CSV: header and a row per marker" "$(head -n 1 "$CSV" 2>/dev/null)/$(awk 'END {print NR}' "$CSV" 2>/dev/null)" "marker,value,confidence/112"
check_eq "CSV and table agree on DYS393" "$(awk -F, '$1 == "DYS393" {print $2}' "$CSV" 2>/dev/null)" "$(col "$H" DYS393 count)"
check_eq "DYS393 is called high" "$(col "$H" DYS393 confidence)" high

# The reads' own count: spanning reads (10 bases of flank each side), the region's
# length plus the CIGAR's insertions minus its deletions inside it or its 10-base
# flanks (an aligner left-aligns an indel in a repeat, which can put it a base or
# two before the region's first base, as at DYS426 in HG002's reads); the modal
# length gives the count with the same rounding as the step (reference copies
# rounded half up, plus the length change in whole and partial repeats).
dk "$SAMTOOLS_IMAGE" samtools view -q 20 -F 0xF04 "${T}/win.bam" > "${T}/win.sam"
python3 -I - "$PANEL" "$OUT" "${T}/win.sam" > "${T}/readcount.tsv" <<'PY'
import collections, re, sys
regions = {}
for line in open(sys.argv[1]):
    c = line.rstrip("\n").split("\t")
    if line.startswith("#") or c[0] == "order" or c[3] != "single":
        continue
    where, period, copies, _name = c[4].split("/", 3)
    s, e = (int(x) for x in where.split(":")[1].split("-"))
    regions[c[2]] = (s, e, int(period), int(float(copies) + 0.5))
lens = collections.defaultdict(collections.Counter)
for line in open(sys.argv[3]):
    f = line.split("\t")
    if f[2] != "chrY" or f[5] == "*":
        continue
    pos = int(f[3])
    ops = [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5])]
    end = pos + sum(n for n, op in ops if op in "MDN=X") - 1
    for m, (s, e, p, rc) in regions.items():
        if pos > s - 10 or end < e + 10:
            continue
        ref, diff = pos, 0
        for n, op in ops:
            if op in "M=X":
                ref += n
            elif op in "DN":
                if s - 10 <= ref <= e + 10:
                    diff -= n
                ref += n
            elif op == "I" and s - 10 <= ref <= e + 10:
                diff += n
        lens[m][diff] += 1
for line in open(sys.argv[2]):
    c = line.rstrip("\n").split("\t")
    if c[2] in regions and c[4] == "high" and lens[c[2]]:
        s, e, p, rc = regions[c[2]]
        diff, n = lens[c[2]].most_common(1)[0]
        whole, rest = divmod(diff, p)
        print(f"{c[2]}\t{c[3]}\t{rc + whole}{'.' + str(rest) if rest else ''}\t{n}")
PY
N_CMP=$(grep -c . "${T}/readcount.tsv" || true)
N_SAME=$(awk -F'\t' '$2 == $3' "${T}/readcount.tsv" | wc -l | tr -d ' ')
awk -F'\t' '$2 != $3 {print "  differs: " $0}' "${T}/readcount.tsv"
check_ge "high calls with spanning reads to compare" "$N_CMP" 50
check "the reads' own count agrees in at least 90% of them (${N_SAME} of ${N_CMP})" test $((N_SAME * 10)) -ge $((N_CMP * 9))

# --- planted: one AGAT repeat cut out of DYS393 -------------------------------------------
P="${SAMPLE}ysp"
read -r DS DE < <(awk -F'\t' '$6 == "DYS393" {print $2, $3}' "${T}/regions.bed")
dk "$SAMTOOLS_IMAGE" samtools view -h "${T}/win.bam" > "${T}/all.sam"
# Reads that span DYS393 with a plain match get bases DS..DS+3 removed (CIGAR xM4DyM);
# other reads over the region are dropped, so only the planted allele is left there.
python3 -I - "$DS" "$DE" "${T}/all.sam" > "${T}/planted.sam" <<'PY'
import re, sys
s, e = int(sys.argv[1]), int(sys.argv[2])
cut = 0
for line in open(sys.argv[3]):
    if line.startswith("@"):
        sys.stdout.write(line)
        continue
    f = line.rstrip("\n").split("\t")
    if f[2] != "chrY" or f[5] == "*":
        sys.stdout.write(line)
        continue
    pos = int(f[3])
    ops = re.findall(r"(\d+)([MIDNSHP=X])", f[5])
    end = pos + sum(int(n) for n, op in ops if op in "MDN=X") - 1
    if end < s or pos > e:
        sys.stdout.write(line)
        continue
    m = re.fullmatch(r"(\d+)M", f[5])
    if not m or pos > s - 10 or end < e + 10:
        continue
    off = s - pos
    n = int(m.group(1))
    f[9] = f[9][:off] + f[9][off + 4:]
    f[10] = f[10][:off] + f[10][off + 4:] if f[10] != "*" else "*"
    f[5] = f"{off}M4D{n - off - 4}M"
    f = [x for x in f if not x.startswith(("MD:Z:", "NM:i:"))]
    sys.stdout.write("\t".join(f) + "\n")
    cut += 1
print(f"planted the cut in {cut} reads", file=sys.stderr)
PY
dk "$SAMTOOLS_IMAGE" samtools sort -o "${T}/planted.bam" "${T}/planted.sam"
dk "$SAMTOOLS_IMAGE" samtools index "${T}/planted.bam"
sample "$P" "${T}/planted.bam" 1
step "$P"
C0=$(col "$H" DYS393 count) C1=$(col "$P" DYS393 count)
check "planted: DYS393 has one repeat fewer (${C0} -> ${C1})" test -n "$C0" -a -n "$C1" -a "$C1" = "$((C0 - 1))" 2>/dev/null
# Only the reads that span DYS393 with a plain match keep it (6 at this depth on
# HG002's single Y), fewer than the 10 a high call needs, so the call is low and says so.
check_eq "planted: DYS393 is low on the few reads left" "$(col "$P" DYS393 confidence)" low
check "planted: its note names the read count" has '[0-9]+ reads, fewer than 10' "$(col "$P" DYS393 note)"
# Reads dropped over DYS393 can take their mates' support from a marker nearby, so
# the check is on the counts, not on the confidence.
join -t $'\t' <(awk -F'\t' 'NR > 1 && $5 == "high" && $3 != "DYS393" {print $3 "\t" $4}' "$OUT" | sort) \
               <(awk -F'\t' 'NR > 1 && $4 != "" && $3 != "DYS393" {print $3 "\t" $4}' "${G}/${P}/y_str/${P}_y_str.tsv" 2>/dev/null | sort) \
  > "${T}/planted_pairs.tsv"
check_ge "planted: the other high calls still have a count" "$(grep -c . "${T}/planted_pairs.tsv")" "$((HIGH - 3))"
check_eq "planted: and every one of them is the same count" "$(awk -F'\t' '$2 != $3' "${T}/planted_pairs.tsv" | wc -l | tr -d ' ')" 0

# --- female -----------------------------------------------------------------------------
F="${SAMPLE}ysf"
sample "$F" "${T}/win.bam" 2
step "$F"
check "a female sample is skipped with one line" has '^Step 45 skipped: indexcov infers female from the reads' "$(cat "$STEP_LOG")"
check "and gets no output" test ! -e "${G}/${F}/y_str"

# --- both reports -------------------------------------------------------------------------
python3 "${REPO}/bin/render_report.py" --sample "$H" --sample-dir "${G}/${H}" --json "${T}/summary.json" \
  -o "${T}/report.txt" -o "${T}/report.html" > /dev/null
TXT=$(cat "${T}/report.txt" 2>/dev/null) HTML=$(cat "${T}/report.html" 2>/dev/null)
check "text report: the Y-STR section" has '^## Y-STR Markers \(HipSTR\)$' "$TXT"
check "text report: DYS393 with its count" has "^  Y-12 +DYS393 +${C0} +\\(high\\)$" "$TXT"
check "HTML report: the Y-STR card" has '<h2>Y-STR Markers \(HipSTR\)</h2>' "$HTML"
check "HTML report: the DYS393 row" has "<td>DYS393</td><td>${C0}</td><td>high</td>" "$HTML"

# --- Nextflow against the bash step ----------------------------------------------------------
# The BAM step 02 wrote from the fixture's reads (its header is the fixture reference's,
# which INDEXCOV checks); three of the panel's regions lie in its chrY slices.
command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
N="${SAMPLE}ysn" BASE="${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
check "the BAM of step 02 is there" test -s "${BASE}.bai"
sample "$N" "$BASE" 1
step "$N"
check_eq "the fixture BAM: 111 rows" "$(awk 'END {print NR - 1}' "${G}/${N}/y_str/${N}_y_str.tsv" 2>/dev/null)" 111
check_ge "the fixture BAM: markers with a count (DYS570, DYS388 and DYS442 lie in its slices)" \
  "$(awk -F'\t' 'NR > 1 && $4 != ""' "${G}/${N}/y_str/${N}_y_str.tsv" 2>/dev/null | wc -l | tr -d ' ')" 2
SHEET="${T}/samplesheet.csv"
# With the VCF of step 03 beside the BAM, so the pipeline does not call variants again.
V="${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
check "the VCF of step 03 is there" test -s "${V}.tbi"
printf 'sample,sex,vcf,vcf_index,bam,bam_index\n%s,male,%s,%s,%s,%s\n' "$N" "$V" "${V}.tbi" \
  "${G}/${N}/aligned/${N}_sorted.bam" "${G}/${N}/aligned/${N}_sorted.bam.bai" > "$SHEET"
NF_OUT="${G}/nf-y-str"
echo "+ nextflow run main.nf -profile docker --tools y_str,html_report"
( cd "$T" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -work-dir "${E2E_WORK}/nf-work-y-str" \
    --input "$SHEET" --reference "$REF" --tools y_str,html_report \
    --sex_check warn --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB ) 2>&1 | tee "$STEP_LOG"
NF_RC=${PIPESTATUS[0]}
check_eq "nextflow run exits 0" "$NF_RC" 0
cp "${T}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
for f in tsv csv; do
  a="${G}/${N}/y_str/${N}_y_str.${f}" b="${NF_OUT}/${N}/y_str/${N}_y_str.${f}"
  check "Nextflow wrote ${N}_y_str.${f}" test -s "$b"
  check "Nextflow's ${N}_y_str.${f} is the bash step's" cmp -s "$a" "$b"
  cmp -s "$a" "$b" || diff "$a" "$b" | head -20
done
check "control: the comparison sees a difference" bash -c '! cmp -s "$1" "$2"' _ \
  "$OUT" "${NF_OUT}/${N}/y_str/${N}_y_str.tsv"
check "Nextflow HTML report: the Y-STR card" has '<h2>Y-STR Markers \(HipSTR\)</h2>' "$(cat "${NF_OUT}/${N}/${N}_report.html" 2>/dev/null)"

finish
