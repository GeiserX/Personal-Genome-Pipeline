#!/usr/bin/env bash
# Step 35 runs Parascopy once per row of the paralog catalog. This case adds a
# GBA row (GBA1 with its pseudogene GBAP1, model GBA.gz of the pinned data) to
# a copy of assets/paralog_loci.tsv and runs the step on HG002 reads:
#   - the genes of assets/parascopy_models_v1.7.txt, each as <gene>.gz, are
#     the names of the list.txt that `setup.sh --parascopy-data` installs, for
#     each of the five populations;
#   - the BAM is the fixture's GIAB slice plus GIAB's 60x HG002 BAM over
#     chr1:155.20-155.26 Mb (GBA1 and GBAP1), read over HTTPS and sampled to
#     about 30x as the fixture is; the GBA table must hold an integer
#     aggregate copy number (agCN) with its quality there;
#   - red control: the fixture's slice alone has no read at GBA1, and its GBA
#     table must say not assessed, with no copy number;
#   - the SMN table is byte for byte what main's commands write (one
#     `parascopy cn-using` of SMN1.gz, the same columns), with a planted
#     one-byte change caught;
#   - the loop of modules/local/parascopy/main.nf, cut from the file, runs on
#     the same BAM in the same image and writes the same two tables.
# The stand-ins of pgx-4-parascopy.sh apply here too (an N-filled chr11 in a
# copy of the reference, background depth from 100 bp windows on the chr20
# slice). The copy numbers on such a slice are not a measurement of HG002; the
# case checks the loop end to end.
. "$(dirname "$0")/lib.sh"

echo "+ scripts/setup.sh --parascopy-data"
"${REPO}/scripts/setup.sh" --parascopy-data "$GENOME_DIR"
DATA="${GENOME_DIR}/reference/parascopy-${PARASCOPY_DATA_VERSION}"
COMMITTED="$(awk -F'\t' '!/^#/ && $3 == "main" {print $1 ".gz"}' "${REPO}/assets/parascopy_models_v1.7.txt" | LC_ALL=C sort)"
for pop in AFR AMR EAS EUR SAS; do
  check "the committed models are the installed ${pop} list.txt" \
    cmp -s <(grep . "${DATA}/models_GRCh38_1KGP/${pop}/list.txt" | LC_ALL=C sort) <(printf '%s\n' "$COMMITTED")
done
check_eq "models on the list" "$(printf '%s\n' "$COMMITTED" | grep -c .)" 167

# --- the stand-in reference (as pgx-4-parascopy.sh) -------------------------------
REF="${GENOME_DIR}/reference/paralog_catalog_ref.fasta"
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

# --- two samples: with and without reads at GBA1 ----------------------------------
GIAB_BAM="https://giab.s3.amazonaws.com/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam"
P=HG002gba
N=HG002nogba
for s in "$P" "$N"; do mkdir -p "${GENOME_DIR}/${s}/aligned"; done
sam_https() {  # the image's libcurl has no CA bundle, so lend it the host's
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${GENOME_DIR}:/genome" -w "/genome/${P}" \
    "$SAMTOOLS_IMAGE" samtools "$@"
}
retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && sleep $((i * 15)); done; return 1; }
retry sam_https view -b -s 42.5 -o gba_reads.bam "$GIAB_BAM" chr1:155200000-155260000
rm -f "${GENOME_DIR}/${P}/$(basename "$GIAB_BAM").bai"
cp "${FIXTURE_DIR}/${SAMPLE}_slice.bam" "${GENOME_DIR}/${N}/aligned/${N}_sorted.bam"
cp "${FIXTURE_DIR}/${SAMPLE}_slice.bam.bai" "${GENOME_DIR}/${N}/aligned/${N}_sorted.bam.bai"
sam merge -f -c -p -o "${P}/aligned/${P}_sorted.bam" "${N}/aligned/${N}_sorted.bam" "${P}/gba_reads.bam"
sam index "${P}/aligned/${P}_sorted.bam"
GBA_REGION=chr1:155234452-155244699
check_ge "${P}: HG002 reads at GBA1" "$(sam view -c -F 0x904 "${P}/aligned/${P}_sorted.bam" "$GBA_REGION" 2>/dev/null)" 1000
check_eq "${N}: no read at GBA1" "$(sam view -c -F 0x904 "${N}/aligned/${N}_sorted.bam" "$GBA_REGION" 2>/dev/null)" 0

LOCI="${GENOME_DIR}/${P}/paralog_loci.tsv"
{ cat "${REPO}/assets/paralog_loci.tsv"
  printf 'GBA.gz\tGBA1,GBAP1\t%s\tagCN of GBA1 and GBAP1 together\tno\tcarrier\tparaphase\tgba\tPMID 35680869\n' "$GBA_REGION"
} > "$LOCI"
BED="${GENOME_DIR}/${P}/parascopy_windows.bed"
awk 'BEGIN { for (s = 10050000; s < 10450000; s += 100) printf "chr20\t%d\t%d\n", s, s + 100 }' > "$BED"

# gba_called TABLE: a row on chr1 over GBA1 with an integer agCN and a numeric quality
# shellcheck disable=SC2016  # awk fields
gba_called() {
  awk -F'\t' 'NR > 1 && $1 == "chr1" && $2 < 155244699 && $3 >= 155234451 && $6 ~ /^[0-9]+$/ && $7 ~ /^[0-9.]+$/ {f = 1}
              END {exit !f}' "$1" 2>/dev/null
}
not_assessed() { awk -F'\t' 'NR > 1 && $5 ~ /^not assessed/ {print $5}' "$1" 2>/dev/null; }

for s in "$P" "$N"; do
  REF_FASTA="$REF" PARASCOPY_DEPTH_BED="$BED" PARALOG_LOCI="$LOCI" run_step 35-paralogs.sh "$s"
  check_step_exit 35-paralogs.sh
  cat "${GENOME_DIR}/${s}/paralogs/${s}_gba_copy_number.tsv" 2>/dev/null
done
PT="${GENOME_DIR}/${P}/paralogs/${P}_gba_copy_number.tsv"
NT="${GENOME_DIR}/${N}/paralogs/${N}_gba_copy_number.tsv"
check "${P}: the GBA table has an integer agCN with its quality over GBA1" gba_called "$PT"
check_eq "${P}: the GBA row is assessed" "$(not_assessed "$PT")" ""
not_called() { ! gba_called "$1"; }
check "red control: the copy-number check fails on the table without reads at GBA1" not_called "$NT"
check_eq "red control: the GBA table says not assessed, with the reason" "$(not_assessed "$NT")" \
  "not assessed: no read over ${GBA_REGION}"
check "${N}: the SMN row still ran" test -s "${GENOME_DIR}/${N}/paralogs/${N}_parascopy/smn/res.samples.bed.gz"
echo "paralog-catalog: GBA agCN on the slice $(awk -F'\t' 'NR > 1 && $1 == "chr1" {printf "%s (q %s) ", $6, $7}' "$PT")" >> "$E2E_NOTES"

# --- the SMN table as main's commands write it ------------------------------------
ST="${GENOME_DIR}/${P}/paralogs/${P}_smn_copy_number.tsv"
O="${GENOME_DIR}/${P}/main_smn"
rm -rf "$O"; mkdir -p "$O"
cat > "${O}/main.sh" <<SH
set -e
parascopy depth -i /genome/${P}/aligned/${P}_sorted.bam::${P} -f /genome/reference/$(basename "$REF") -b /genome/${P}/parascopy_windows.bed --no-gc -o depth -@ 2
parascopy cn-using /genome/reference/parascopy-${PARASCOPY_DATA_VERSION}/models_GRCh38_1KGP/EUR/SMN1.gz -i /genome/${P}/aligned/${P}_sorted.bam::${P} \
  -f /genome/reference/$(basename "$REF") -t /genome/reference/parascopy-${PARASCOPY_DATA_VERSION}/homology_table/GRCh38.bed.gz -d depth -o out -@ 2
{
  printf 'chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions\n'
  gzip -dc out/res.samples.bed.gz | awk -F'\t' -v OFS='\t' '!/^#/ {print \$1, \$2, \$3, \$4, \$6, \$7, \$8, \$9, \$10, \$11, \$13}'
} > smn.tsv
SH
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -v "${GENOME_DIR}:/genome" -w "/genome/${P}/main_smn" \
  "$PARASCOPY_IMAGE" bash main.sh 2>&1 | tail -n 3
check_ge "main's commands wrote SMN rows" "$(awk 'NR > 1' "${O}/smn.tsv" 2>/dev/null | wc -l | tr -d ' ')" 1
check "the SMN table is byte for byte main's" cmp -s "$ST" "${O}/smn.tsv"
check "the comparison sees a planted one-byte change" \
  sh -c "sed '2s/.\$/x/' '${O}/smn.tsv' > '${O}/planted.tsv'; ! cmp -s '$ST' '${O}/planted.tsv'"

# --- the module's loop, cut from main.nf, on the same BAM -------------------------
MOD="${REPO}/modules/local/parascopy/main.nf"
M="${GENOME_DIR}/${P}/module"
rm -rf "$M"; mkdir -p "$M"
SPECS=$(awk -F'\t' '/^#/ || !NF {next} !h {for (i = 1; i <= NF; i++) c[$i] = i; h = 1; next}
  {printf "%s%s|%s|%s", sep, $c["model"], $c["region"], $c["table"]; sep = " "}' "$LOCI")
sed -n '/# Paralog loci (as scripts\/35-paralogs.sh)/,/# end of the paralog loci/p' "$MOD" \
  | sed -e 's/\\\\/\\/g' -e 's/\\\$/$/g' \
        -e "s|\${parascopy_data}|/genome/reference/parascopy-${PARASCOPY_DATA_VERSION}|g" \
        -e "s|\${params.parascopy_population}|EUR|g" \
        -e "s|\${bam}|/genome/${P}/aligned/${P}_sorted.bam|g" \
        -e "s|\${meta.id}|${P}|g" \
        -e "s|\${reference}|/genome/reference/$(basename "$REF")|g" \
        -e "s|\${background}|-b /genome/${P}/parascopy_windows.bed --no-gc|g" \
        -e "s|\${task.cpus}|2|g" \
        -e "s/\${loci}/${SPECS}/g" > "${M}/block.sh"
check_ge "the module's loop was cut from main.nf" "$(wc -l < "${M}/block.sh" | tr -d ' ')" 25
check "the cut loop has no Nextflow variable left" \
  sh -c "! grep -E '[$][{](parascopy_data|params[.]|bam|meta[.]|reference|background|task[.]|loci)' '${M}/block.sh'"
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -v "${GENOME_DIR}:/genome" -w "/genome/${P}/module" \
  "$PARASCOPY_IMAGE" bash -ue block.sh 2>&1 | tail -n 3
for t in smn gba; do
  check "the module's loop writes the same ${t} table as step 35" \
    cmp -s "${GENOME_DIR}/${P}/paralogs/${P}_${t}_copy_number.tsv" "${M}/${P}_${t}_copy_number.tsv"
done
rm -f "$REF" "${REF}.fai"

finish
