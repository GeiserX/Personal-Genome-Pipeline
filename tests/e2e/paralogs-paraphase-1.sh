#!/usr/bin/env bash
# Step 35b (Paraphase, opt-in) on PacBio HiFi reads of HG002, and its refusal
# of a short-read BAM. The fixture has Illumina reads only, and Paraphase is
# written for HiFi reads (its licence allows PacBio data only), so the case
# reads its own HiFi slice:
#   - GIAB's HG002 PacBio Revio BAM (48x, GRCh38-GIABv3, haplotagged), read
#     over HTTPS through its index: the regions of the eight gene families
#     below in Paraphase's own config (data/38/config.yaml of the release the
#     image pins) with 30 kb either side, and its depth windows on chr5, chrX
#     and chrY (data/38/genome_region.bed) so it measures the genome depth as
#     on a whole genome. Methylation tags (MM, ML) are dropped to keep it small;
#   - the reference is the fixture's, plus chr7 and chr15 (PMS2, NCF1, STRC)
#     from UCSC hg38, which is the same sequence as the no-ALT analysis set on
#     those two chromosomes.
# Checks:
#   1. step 35b stops on the fixture's Illumina BAM (no PL:PACBIO read group);
#   2. on the HiFi slice every one of the eight families has its rows, and
#      SMN1, SMN2, PMS2, STRC, NCF1, HBA1, HBA2, the RCCX modules and GBA1
#      have the copy numbers Paraphase calls on this slice (GBA1 and OTOA
#      counted from their gene1 haplotypes, not the family's total);
#   3. a planted -a4.2 deletion: on the reads of haplotype 1 the 4,257 bp
#      between the two homologous HBA "A" boxes (chr16:169,742 and its copy at
#      173,999) are cut out and the reads aligned again. Paraphase must then
#      name the 4.2 kb deletion, report one HBA2 copy fewer than on the
#      unchanged slice, and the same HBA1. A red control swaps the HBA1 and
#      HBA2 labels in that JSON, as swapped regions in the config would, and
#      the HBA2 check must then fail;
#   4. Nextflow on a HiFi row and an Illumina row: the check fails on the
#      Illumina row, the run goes on, and the HiFi row gets the same HBA rows
#      as the bash step and a report card.
. "$(dirname "$0")/lib.sh"

REGIONS=smn1,hba,rccx,GBA,pms2,strc,ncf1,OTOA
PV=${PARAPHASE_IMAGE##*:} && PV=${PV%%--*}
RAW="https://raw.githubusercontent.com/PacificBiosciences/paraphase/v${PV}/data/38"
HIFI="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam"
UCSC=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes

# samtools over HTTPS: the image's libcurl has no CA bundle, so lend it the host's.
sam_https() {
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${GENOME_DIR}:/genome" -w /genome \
    "$SAMTOOLS_IMAGE" samtools "$@"
}
mm2() { in_genome "$MINIMAP2_IMAGE" minimap2 "$@"; }
retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && sleep $((i * 15)); done; return 1; }
tsv() {  # tsv FILE REGION GENE COLUMN: a cell of a step 35b copy-number table
  awk -F'\t' -v r="$2" -v g="$3" -v k="$4" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    $1 == r && $2 == g {print $c[k]; exit}' "$1" 2>/dev/null
}

# --- 1. A short-read BAM stops at the read-group check -----------------------
S=HG002ill
mkdir -p "${GENOME_DIR}/${S}/aligned_longread"
cp "${FIXTURE_DIR}/${SAMPLE}_slice.bam" "${GENOME_DIR}/${S}/aligned_longread/${S}_sorted.bam"
cp "${FIXTURE_DIR}/${SAMPLE}_slice.bam.bai" "${GENOME_DIR}/${S}/aligned_longread/${S}_sorted.bam.bai"
sam addreplacerg -r '@RG\tID:HG002\tSM:HG002\tPL:ILLUMINA' -o "${S}/aligned_longread/rg.bam" "${S}/aligned_longread/${S}_sorted.bam" \
  && mv "${GENOME_DIR}/${S}/aligned_longread/rg.bam" "${GENOME_DIR}/${S}/aligned_longread/${S}_sorted.bam" \
  && sam index "${S}/aligned_longread/${S}_sorted.bam"
run_step 35b-paraphase.sh "$S"
check "step 35b refuses the Illumina BAM (exit ${STEP_RC})" test "$STEP_RC" -ne 0
check "it says why: PL ILLUMINA, not PACBIO" grep -q 'PL ILLUMINA, not PACBIO' "$STEP_LOG"
check "Paraphase did not run on it" test ! -e "${GENOME_DIR}/${S}/paralogs/${S}_paraphase"
rm -rf "${GENOME_DIR:?}/${S}"

# --- The HiFi slice and the reference ----------------------------------------
P=HG002hifi
D="${GENOME_DIR}/${P}/aligned_longread"
mkdir -p "$D"
curl -fsSL --retry 5 -o "${CASE_TMP}/config.yaml" "${RAW}/config.yaml"
curl -fsSL --retry 5 -o "${CASE_TMP}/genome_region.bed" "${RAW}/genome_region.bed"
python3 - "${CASE_TMP}/config.yaml" "${CASE_TMP}/genome_region.bed" "$REGIONS" > "${D}/regions.bed" <<'PY'
import re, sys
cfg, depth, wanted = open(sys.argv[1]).read(), sys.argv[2], set(sys.argv[3].split(","))
iv, block, chrom = [], None, None
for line in cfg.splitlines():
    m = re.match(r"^(\S+):\s*$", line)
    if m:
        block = m.group(1)
        continue
    if block not in wanted:
        continue
    for c, s, e in re.findall(r"(chr[0-9XY]+):(\d+)-(\d+)", line):
        iv.append((c, int(s), int(e)))
        chrom = c
    if "depth_region" in line and chrom:   # [[start, end]] on the family's chromosome
        for s, e in re.findall(r"\[(\d+),\s*(\d+)\]", line):
            iv.append((chrom, int(s), int(e)))
iv = [(c, max(0, s - 30000), e + 30000) for c, s, e in iv]
for line in open(depth):   # positions of 1.6 kb depth windows: 2 kb either side
    f = line.split()
    if len(f) >= 2:
        iv.append((f[0], max(0, int(f[1]) - 2000), int(f[1]) + 2000))
iv.sort()
out = []
for c, s, e in iv:
    if out and out[-1][0] == c and s <= out[-1][2]:
        out[-1][2] = max(out[-1][2], e)
    else:
        out.append([c, s, e])
for c, s, e in out:
    print(f"{c}\t{s}\t{e}")
PY
echo "HiFi slice: $(wc -l < "${D}/regions.bed") regions, $(awk '{n += $3 - $2} END {print n}' "${D}/regions.bed") bp"
check_ge "the slice has the families' regions and the depth windows" "$(wc -l < "${D}/regions.bed" | tr -d ' ')" 20
retry sam_https view -b -M -L "${P}/aligned_longread/regions.bed" --keep-tag RG,NM,SA,rq,HP,PS \
  -o "${P}/aligned_longread/${P}_sorted.bam" "$HIFI"
rm -f "${GENOME_DIR}/$(basename "$HIFI").bai"
sam index "${P}/aligned_longread/${P}_sorted.bam"
check_ge "HiFi reads in the slice" "$(sam view -c "${P}/aligned_longread/${P}_sorted.bam" 2>/dev/null)" 5000

REF="${GENOME_DIR}/reference/paraphase_test_ref.fasta"
cp "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta" "$REF"
for c in chr7 chr15; do
  curl -fsSL --retry 5 "${UCSC}/${c}.fa.gz" | gzip -dc | awk '/^>/ {print; next} {print toupper($0)}' >> "$REF"
done
sam faidx "reference/$(basename "$REF")"
check "the stand-in reference has chr7 and chr15 of GRCh38's lengths" \
  awk -F'\t' '$1 == "chr7" && $2 == 159345973 {a = 1} $1 == "chr15" && $2 == 101991189 {b = 1} END {exit !(a && b)}' "${REF}.fai"

# --- 2. Step 35b on the HiFi slice --------------------------------------------
REF_FASTA="$REF" PARAPHASE_REGIONS="$REGIONS" run_step 35b-paraphase.sh "$P"
check_step_exit 35b-paraphase.sh
CN="${GENOME_DIR}/${P}/paralogs/${P}_paraphase_copy_number.tsv"
cat "$CN" 2>/dev/null
head -c 3000 "${GENOME_DIR}/${P}/paralogs/${P}_paraphase/${P}.paraphase.json" 2>/dev/null; echo
head -n 5 "${GENOME_DIR}/${P}/paralogs/${P}_paraphase_variants.tsv" 2>/dev/null
for r in ${REGIONS//,/ }; do
  check_ge "rows for region ${r}" "$(awk -F'\t' -v r="$r" 'NR > 1 && $1 == r' "$CN" 2>/dev/null | wc -l | tr -d ' ')" 1
done
# The copy numbers Paraphase calls on this HG002 slice with the pinned image,
# the same on every run. They are Paraphase's own calls, not a published truth
# set (none is published per gene for HG002): a change in the slice, the
# image or the parsing that moves one goes red here.
for g in smn1:SMN1:2 smn1:SMN2:2 "smn1:SMN exon 7-8 deletion:0" pms2:PMS2:2 strc:STRC:2 ncf1:NCF1:2 \
         hba:HBA1:2 hba:HBA2:2 "rccx:CYP21A2 (RCCX):4" GBA:GBA1:2; do
  IFS=: read -r reg gene want <<< "$g"
  check_eq "${gene} copies" "$(tsv "$CN" "$reg" "$gene" copy_number)" "$want"
done
check "GBA1 is counted from its gene1 haplotypes, not the family's $(tsv "$CN" GBA GBA1 total_cn)" \
  has "copies of GBA1 alone" "$(tsv "$CN" GBA GBA1 detail)"
OTOA_CN=$(tsv "$CN" OTOA OTOA copy_number) OTOA_TOTAL=$(tsv "$CN" OTOA OTOA total_cn)
check "OTOA counts the gene alone (${OTOA_CN:-none}), fewer than the family's ${OTOA_TOTAL:-none}" \
  test "${OTOA_CN:-x}" -lt "${OTOA_TOTAL:-0}" 2>/dev/null
check_ge "variant rows from Paraphase's VCFs" "$(awk 'NR > 1' "${GENOME_DIR}/${P}/paralogs/${P}_paraphase_variants.tsv" 2>/dev/null | wc -l | tr -d ' ')" 1
HBA2_BASE=$(tsv "$CN" hba HBA2 copy_number) HBA1_BASE=$(tsv "$CN" hba HBA1 copy_number)

# --- 3. A planted 4.2 kb deletion on haplotype 1 -------------------------------
Q=HG002del
QD="${GENOME_DIR}/${Q}/aligned_longread"
mkdir -p "$QD"
DEL_START=169742 DEL_END=173999   # 1-based: bases DEL_START..DEL_END-1 go
# Pass 1: the haplotype-1 reads over the deleted bases lose them: a read that
# spans them is joined across the gap, one that ends or starts inside keeps
# the part outside (at least 1 kb, or it is dropped). They are written as
# FASTQ, their tags in the comment for minimap2 -y, and listed for pass 2.
cat > "${CASE_TMP}/cut.py" <<'PY'
import re, sys
a, b, fq, names = int(sys.argv[1]), int(sys.argv[2]), open(sys.argv[3], "w"), open(sys.argv[4], "w")
spanning = trimmed = dropped = 0
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    tags = f[11:]
    if "HP:i:1" not in tags or int(f[1]) & 0x900:
        continue
    pos = ref = int(f[3])
    qry, qa, qb = 0, None, None   # query offsets of the first base at or after a, and b
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]):
        n = int(n)
        if op in "M=XDN":
            for x, attr in ((a, "qa"), (b, "qb")):
                if (qa if attr == "qa" else qb) is None and ref + n > x:
                    off = qry + max(0, x - ref) if op in "M=X" else qry
                    if attr == "qa":
                        qa = off
                    else:
                        qb = off
            ref += n
            if op in "M=X":
                qry += n
        elif op in "IS":
            qry += n
    if ref <= a or pos >= b:
        continue
    names.write(f[0] + "\n")
    left = f[9][:qa] if pos < a else ""
    right = f[9][qb:] if ref > b else ""
    lq = f[10][:qa] if pos < a else ""
    rq = f[10][qb:] if ref > b else ""
    if left and right:
        spanning += 1
    elif len(left + right) >= 1000:
        trimmed += 1
    else:
        dropped += 1
        continue
    keep = [t for t in tags if t[:2] in ("RG", "rq", "HP", "PS")]
    fq.write(f"@{f[0]}\t" + "\t".join(keep) + f"\n{left + right}\n+\n{lq + rq}\n")
print(f"haplotype-1 reads joined across the deletion: {spanning}, trimmed: {trimmed}, dropped: {dropped}", file=sys.stderr)
print(spanning)
PY
SPAN=$(sam view "${P}/aligned_longread/${P}_sorted.bam" "chr16:$((DEL_START - 50000))-$((DEL_END + 50000))" \
  | python3 "${CASE_TMP}/cut.py" "$DEL_START" "$DEL_END" "${QD}/cut.fq" "${QD}/names.txt")
check_ge "haplotype-1 reads joined across the planted deletion" "${SPAN:-0}" 5
# Pass 2: every other read as it was, and the cut reads aligned again to chr16.
sam view -b -N "^${Q}/aligned_longread/names.txt" -o "${Q}/aligned_longread/kept.bam" "${P}/aligned_longread/${P}_sorted.bam"
sam faidx "reference/$(basename "$REF")" chr16 > "${QD}/chr16.fa"
mm2 -t "${THREADS:-4}" -a -x map-hifi -y "${Q}/aligned_longread/chr16.fa" "${Q}/aligned_longread/cut.fq" > "${QD}/cut.sam"
sam sort -o "${Q}/aligned_longread/cut.bam" "${Q}/aligned_longread/cut.sam"
sam merge -f -o "${Q}/aligned_longread/merged.bam" "${Q}/aligned_longread/kept.bam" "${Q}/aligned_longread/cut.bam"
sam sort -o "${Q}/aligned_longread/${Q}_sorted.bam" "${Q}/aligned_longread/merged.bam"
sam index "${Q}/aligned_longread/${Q}_sorted.bam"
rm -f "${QD}/kept.bam" "${QD}/cut.bam" "${QD}/cut.sam" "${QD}/merged.bam" "${QD}/chr16.fa" "${QD}/chr16.fa.fai"

REF_FASTA="$REF" PARAPHASE_REGIONS=hba run_step 35b-paraphase.sh "$Q"
check_step_exit 35b-paraphase.sh
QCN="${GENOME_DIR}/${Q}/paralogs/${Q}_paraphase_copy_number.tsv"
cat "$QCN" 2>/dev/null
# The named call, not the "no 3.7 kb or 4.2 kb deletion ... found" line.
NAMED='4\.2 kb deletion \(-a4\.2\)'
check "the 4.2 kb deletion is named" has "$NAMED" "$(tsv "$QCN" hba HBA2 detail)"
check "the unchanged slice has no 4.2 kb deletion" lacks "$NAMED" "$(tsv "$CN" hba HBA2 detail)"
check_eq "one HBA2 copy lost (unchanged slice: ${HBA2_BASE:-none})" "$(tsv "$QCN" hba HBA2 copy_number)" \
  "$([[ "${HBA2_BASE:-}" =~ ^[0-9]+$ ]] && echo $((HBA2_BASE - 1)) || echo 'a number')"
check_eq "HBA1 as on the unchanged slice" "$(tsv "$QCN" hba HBA1 copy_number)" "${HBA1_BASE:-none}"
# Red control: the HBA1 and HBA2 regions swapped in Paraphase's config label
# each HBA1 haplotype hba2hap and each HBA2 one hba1hap. The same JSON with the
# two labels swapped must then fail the "one HBA2 copy lost" check above.
SWJ="${CASE_TMP}/swapped.paraphase.json"
sed -e 's/hba1hap/hbaXhap/g' -e 's/hba2hap/hba1hap/g' -e 's/hbaXhap/hba2hap/g' \
  "${GENOME_DIR}/${Q}/paralogs/${Q}_paraphase/${Q}.paraphase.json" > "$SWJ"
python3 "${REPO}/bin/paraphase_report.py" table --sample "$Q" --json "$SWJ" \
  --out-cn "${CASE_TMP}/swapped_cn.tsv" --out-variants "${CASE_TMP}/swapped_variants.tsv" >/dev/null
SW_HBA2=$(tsv "${CASE_TMP}/swapped_cn.tsv" hba HBA2 copy_number)
check "control: with HBA1 and HBA2 swapped the HBA2 check goes red (HBA2 ${SW_HBA2:-none}, the check wants $([[ "${HBA2_BASE:-}" =~ ^[0-9]+$ ]] && echo $((HBA2_BASE - 1))))" \
  test -n "$SW_HBA2" -a "$SW_HBA2" != "$([[ "${HBA2_BASE:-}" =~ ^[0-9]+$ ]] && echo $((HBA2_BASE - 1)) || echo 'a number')"

# --- 4. The same in Nextflow: a HiFi row and an Illumina row ------------------
# PARAPHASE_CHECK fails on the Illumina row and is ignored (report-only), the
# run ends, and the HiFi row gets the table the bash step wrote and a report
# card from it.
if command -v nextflow >/dev/null; then
  NF="${CASE_TMP}/nf"
  rm -rf "$NF" && mkdir -p "$NF"
  V="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
  {
    echo 'sample,vcf,vcf_index,bam,bam_index,sex'
    echo "${P},${V},${V}.tbi,${D}/${P}_sorted.bam,${D}/${P}_sorted.bam.bai,male"
    echo "HG002ill,${V},${V}.tbi,${FIXTURE_DIR}/${SAMPLE}_slice.bam,${FIXTURE_DIR}/${SAMPLE}_slice.bam.bai,male"
  } > "${NF}/samplesheet.csv"
  echo "+ nextflow run main.nf --tools paraphase,html_report --paraphase_regions hba"
  ( cd "$NF" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -work-dir "${NF}/work" \
      --input "${NF}/samplesheet.csv" --reference "$REF" --tools paraphase,html_report \
      --paraphase_regions hba --sex_check warn --outdir "${NF}/out" --max_cpus 4 --max_memory 14.GB ) > "${NF}/run.log" 2>&1
  RC=$?
  grep -vE 'Pulling|Waiting|Verifying|Download complete|Pull complete|Already exists' "${NF}/run.log"
  cp "${NF}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
  NCN="${NF}/out/${P}/paralogs/${P}_paraphase_copy_number.tsv"
  cat "$NCN" 2>/dev/null
  check_eq "Nextflow: the run ends" "$RC" 0
  check "Nextflow: PARAPHASE_CHECK failed on the Illumina row and was ignored" grep -Eq 'PARAPHASE_CHECK[[:space:]]+ignored' "${NF}/failed_tasks.tsv"
  check "Nextflow: no Paraphase output for the Illumina row" test ! -e "${NF}/out/HG002ill/paralogs/HG002ill_paraphase_copy_number.tsv"
  check_eq "Nextflow: HBA2 as the bash step found" "$(tsv "$NCN" hba HBA2 copy_number)" "${HBA2_BASE:-none}"
  check_eq "Nextflow: HBA1 as the bash step found" "$(tsv "$NCN" hba HBA1 copy_number)" "${HBA1_BASE:-none}"
  check "Nextflow: the report has the paralog card with HBA2" grep -q 'HBA2' "${NF}/out/${P}/${P}_report.html"
  # The card counts the variant table too, so HTML_REPORT must be handed it, not only the copy numbers.
  NV=$(awk -F'\t' 'NR > 1 && $1 == "hba"' "${GENOME_DIR}/${P}/paralogs/${P}_paraphase_variants.tsv" 2>/dev/null | wc -l | tr -d ' ')
  check_ge "the bash step has hba variant rows" "$NV" 1
  check "Nextflow: the report counts the ${NV} hba variants the bash step found" \
    grep -q "Variants on phased haplotypes</span><span class=\"value\">${NV}<" "${NF}/out/${P}/${P}_report.html"
  rm -rf "${NF}/work"
else
  fail "nextflow is not on PATH"
fi

rm -f "$REF" "${REF}.fai"
rm -rf "${GENOME_DIR:?}/${Q}" "${GENOME_DIR:?}/${P}/aligned_longread"
finish
