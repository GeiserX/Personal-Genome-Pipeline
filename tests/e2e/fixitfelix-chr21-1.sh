#!/usr/bin/env bash
# setup.sh --masked-reference on chr21, and what it changes for KCNE1, CBS and
# CRYAA, three genes GRCh38 holds a false second copy of on chr21's short arm.
# Public data only: GIAB's HG002 and its CMRG v1.00 benchmark. The case:
#   - takes chr21 of the no-ALT analysis set from GIAB's copy (samtools reads
#     only that sequence, checked against its md5, the M5 of chr21), and runs
#     setup.sh --masked-reference on it (MASKED_REF_SHA256 pins chr21's masked
#     FASTA): five stretches of chr21, 771,127 bases, become N, the genes stay,
#     and the file is byte for byte GIAB's maskedGRC_exclusions_v2 chr21;
#   - control: the same run with another MASKED_REF_SHA256 fails and leaves
#     no masked FASTA;
#   - takes HG002's 60x read pairs from GIAB's BAM over both copies of each
#     gene (the five masked stretches and the gene bodies with 5 kb flanks),
#     aligns them with step 02 to the stock chr21 and to the masked chr21 and
#     calls them with step 03 over the three gene bodies;
#   - counts the CMRG v1.00 benchmark variants of each gene (inside its
#     benchmark regions) that a PASS call matches, alleles and genotype, both
#     sides split and left-aligned by bcftools norm: on the masked reference
#     the recall of each gene is 0.9 or more; on the stock reference it is
#     under 0.5, the red control: the same reads, aligner and caller;
#   - runs step 16b's gene pass with the rows of those three genes and U2AF1
#     (no CMRG variant in HG002, so no recall) of assets/callability_genes.tsv:
#     read on the masked reference, not read on the stock one.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
RD="${G}/reference/ffx"           # this case's references and benchmark, under GENOME_DIR
GIAB_REF=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38
STOCK_URL="${GIAB_REF}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz"
GIAB_MASKED_URL="${GIAB_REF}/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz"
BAM_URL=https://giab.s3.amazonaws.com/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam
CMRG=https://giab.s3.amazonaws.com/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00
CMRG_VCF_SHA=80f18efa283d953cfc615d30dd01b1ff80debcc46f52042b942da783acd48da2
CMRG_BED_SHA=dc77803bb474802fd3213d2f3649076dc14675d3b3c0ae3284036be09d2813f7
CHR21_M5=974dc7aec0b755b19f031418fdedf293       # chr21 of GRCh38, upper case, as in a .dict
MASKED_CHR21_SHA=52d5fa9041423d95938a3e8062bc0b2c9e763fed0a84b2694376072a69ed7d1d
# gene, 0-based start, end (GENCODE 50, as in assets/callability_genes.tsv)
GENES="KCNE1:34446687:34512214 CBS:43050899:43076943 CRYAA:43169007:43172810"

rm -rf "$RD" "${G}/ffx-stock" "${G}/ffx-masked"
mkdir -p "$RD"
# SAMTOOLS_IMAGE's libcurl has no CA bundle: lend it the host's (htslib reads CURL_CA_BUNDLE).
dk() {
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${RD}:/w" -w /w "$@"
}
sha() { sha256sum "$1" 2>/dev/null | cut -d' ' -f1; }

# --- chr21, stock and masked ----------------------------------------------------------
dk "$SAMTOOLS_IMAGE" samtools faidx -o chr21.fasta "$STOCK_URL" chr21
dk "$SAMTOOLS_IMAGE" samtools faidx chr21.fasta
M5=$(awk '!/^>/' "${RD}/chr21.fasta" 2>/dev/null | tr -d '\n' | tr 'acgtn' 'ACGTN' | md5sum | cut -d' ' -f1)
check_eq "chr21 of GIAB's no-ALT analysis set has the M5 of GRCh38's chr21" "$M5" "$CHR21_M5"

# setup_mask NAME SHA: setup.sh --masked-reference on RD/NAME/chr21.fasta (a hard
# link of chr21) with MASKED_REF_SHA256=SHA; output in CASE_TMP/NAME.log, exit in MASK_RC.
setup_mask() {
  mkdir -p "${RD}/$1"
  ln -f "${RD}/chr21.fasta" "${RD}/$1/chr21.fasta"
  cp "${RD}/chr21.fasta.fai" "${RD}/$1/chr21.fasta.fai"
  echo "+ REF_FASTA=reference/ffx/$1/chr21.fasta MASKED_REF_SHA256=$2 scripts/setup.sh --masked-reference"
  REF_FASTA="reference/ffx/$1/chr21.fasta" MASKED_REF_SHA256=$2 \
    "${REPO}/scripts/setup.sh" --masked-reference "$G" > "${CASE_TMP}/$1.log" 2>&1
  MASK_RC=$?
  tail -n 6 "${CASE_TMP}/$1.log"
}
setup_mask masked "$MASKED_CHR21_SHA"
MASKED="${RD}/masked/chr21_maskedGRC.fasta"
check_eq "setup.sh --masked-reference exits 0" "$MASK_RC" 0
check "it wrote the masked FASTA, its .fai and its .dict" test -s "$MASKED" -a -s "${MASKED}.fai" -a -s "${RD}/masked/chr21_maskedGRC.dict"
check_eq "the masked FASTA is the pinned one" "$(sha "$MASKED")" "$MASKED_CHR21_SHA"
check "it reports the five chr21 rows" has 'masked 5 of 12 rows \(771127 bases\)' "$(cat "${CASE_TMP}/masked.log")"
check_eq "the .fai lists chr21 at its length" "$(cut -f1,2 "${MASKED}.fai" 2>/dev/null)" "chr21	46709983"
# N over the five stretches, in the masked file and not in the stock one
n_in() { dk "$SAMTOOLS_IMAGE" samtools faidx "$1" chr21:5010001-5166246 chr21:5966594-6161371 chr21:6427259-6580181 \
  chr21:6789086-6934219 chr21:7743701-7865746 | awk '!/^>/' | tr -cd 'Nn' | wc -c | tr -d ' '; }
check_eq "the five stretches are N in the masked chr21" "$(n_in masked/chr21_maskedGRC.fasta)" 771127
STOCK_N=$(n_in chr21.fasta)
check_eq "control: none of them is N in the stock chr21" "$STOCK_N" 0
GENE_SEQ() { dk "$SAMTOOLS_IMAGE" samtools faidx "$1" chr21:34446688-34512214 chr21:43050900-43172810 | md5sum | cut -d' ' -f1; }
check_eq "the gene bodies are the same in both" "$(GENE_SEQ masked/chr21_maskedGRC.fasta)" "$(GENE_SEQ chr21.fasta)"
dk "$SAMTOOLS_IMAGE" samtools faidx -o giab_masked_chr21.fasta "$GIAB_MASKED_URL" chr21
check "the masked chr21 is byte for byte GIAB's maskedGRC_exclusions_v2 chr21" cmp -s "$MASKED" "${RD}/giab_masked_chr21.fasta"

setup_mask gate 0000000000000000000000000000000000000000000000000000000000000000
check "control: another MASKED_REF_SHA256 fails setup.sh (exit ${MASK_RC})" test "$MASK_RC" -ne 0
check "control: it names the sha256 mismatch" has 'is not the file it is pinned for' "$(cat "${CASE_TMP}/gate.log")"
check "control: it leaves no masked FASTA" test ! -e "${RD}/gate/chr21_maskedGRC.fasta" -a ! -e "${RD}/gate/chr21_maskedGRC.fasta.part"

# --- HG002's reads over both copies of the three genes -------------------------------
REGIONS="chr21:5010001-5166246 chr21:5966594-6161371 chr21:6427259-6580181 chr21:6789086-6934219 chr21:7743701-7865746
  chr21:34441688-34517214 chr21:43045900-43177810"
# shellcheck disable=SC2086  # REGIONS is a list of region literals
dk "$SAMTOOLS_IMAGE" bash -euo pipefail -c "samtools view -u -M -F 2316 '${BAM_URL}' ${REGIONS//$'\n'/ } \
  | samtools sort -n -u -T /w/sortn - | samtools fastq -n -1 /w/R1.fastq.gz -2 /w/R2.fastq.gz -0 /dev/null -s /dev/null -"
PAIRS=$(( $(gzip -dc "${RD}/R1.fastq.gz" 2>/dev/null | wc -l) / 4 ))
check_ge "read pairs from GIAB's HG002 BAM" "$PAIRS" 100000
check_eq "R2 holds as many reads" "$(( $(gzip -dc "${RD}/R2.fastq.gz" 2>/dev/null | wc -l) / 4 ))" "$PAIRS"

INTERVALS=""
for g in $GENES; do IFS=: read -r _ s e <<< "$g"; INTERVALS="${INTERVALS:+${INTERVALS} }chr21:$((s + 1))-${e}"; done
export INTERVALS
for side in stock masked; do
  S="ffx-${side}"
  mkdir -p "${G}/${S}/fastq"
  cp "${RD}/R1.fastq.gz" "${G}/${S}/fastq/${S}_R1.fastq.gz"
  cp "${RD}/R2.fastq.gz" "${G}/${S}/fastq/${S}_R2.fastq.gz"
  if [ "$side" = stock ]; then REF=reference/ffx/chr21.fasta; else REF=reference/ffx/masked/chr21_maskedGRC.fasta; fi
  REF_FASTA=$REF run_step 02-alignment.sh "$S"
  check_step_exit "02-alignment.sh (${side})"
  REF_FASTA=$REF run_step 03-deepvariant.sh "$S"
  check_step_exit "03-deepvariant.sh (${side})"
done

# --- recall against CMRG v1.00 --------------------------------------------------------
curl -fsSL --retry 3 --retry-all-errors -o "${RD}/cmrg.vcf.gz" "${CMRG}.vcf.gz"
curl -fsSL --retry 3 --retry-all-errors -o "${RD}/cmrg.bed" "${CMRG}.bed"
check_eq "the CMRG v1.00 VCF is the pinned file" "$(sha "${RD}/cmrg.vcf.gz")" "$CMRG_VCF_SHA"
check_eq "the CMRG v1.00 BED is the pinned file" "$(sha "${RD}/cmrg.bed")" "$CMRG_BED_SHA"
# norm VCF OUT: VCF split and left-aligned on the stock chr21, chr21 only, as text
norm() {
  docker run --rm -u "$(id -u):$(id -g)" -v "${G}:/genome" -w /genome "$BCFTOOLS_IMAGE" \
    bash -c "bcftools view -r chr21 '$1' | bcftools norm -m -any -f reference/ffx/chr21.fasta 2>/dev/null | bcftools view -H" > "$2"
}
dk "$BCFTOOLS_IMAGE" bcftools index -f -t cmrg.vcf.gz
norm reference/ffx/cmrg.vcf.gz "${CASE_TMP}/cmrg.txt"
for side in stock masked; do
  norm "ffx-${side}/vcf/ffx-${side}.vcf.gz" "${CASE_TMP}/${side}.txt"
done
# recall SIDE: "gene matched benchmark recall" lines
recall() {
  python3 - "${CASE_TMP}/cmrg.txt" "${RD}/cmrg.bed" "${CASE_TMP}/$1.txt" "$GENES" <<'PY'
import sys
truth, bed, calls, genes = sys.argv[1:5]
regions = [(int(c[1]), int(c[2])) for c in (l.split("\t") for l in open(bed)) if c[0] == "chr21"]
def alts(gt):
    return sum(1 for a in gt.replace("|", "/").split("/") if a not in ("0", "."))
def rows(path, pass_only):
    out = {}
    for line in open(path):
        c = line.rstrip("\n").split("\t")
        if pass_only and c[6] not in ("PASS", "."):
            continue
        n = alts(c[9].split(":")[0])
        if n:
            out[(int(c[1]), c[3], c[4])] = n
    return out
t, k = rows(truth, False), rows(calls, True)
for g in genes.split():
    name, s, e = g.split(":")
    s, e = int(s), int(e)
    keys = [x for x in t if s < x[0] <= e and any(a < x[0] <= b for a, b in regions)]
    hit = sum(1 for x in keys if k.get(x) == t[x])
    print(name, hit, len(keys), f"{hit / len(keys):.3f}" if keys else "nan")
PY
}
R_MASKED=$(recall masked 2>&1); echo "masked chr21: ${R_MASKED}"
R_STOCK=$(recall stock 2>&1); echo "stock chr21: ${R_STOCK}"
# ok_recall LINES GENE: 0 when GENE's recall is 0.9 or more; low_recall: under 0.5
# shellcheck disable=SC2329  # run through check
ok_recall() { awk -v g="$2" '$1 == g {found = 1; ok = ($4 + 0 >= 0.9)} END {exit !(found && ok)}' <<< "$1"; }
# shellcheck disable=SC2329  # run through check
low_recall() { awk -v g="$2" '$1 == g {found = 1; ok = ($4 + 0 < 0.5)} END {exit !(found && ok)}' <<< "$1"; }
for g in KCNE1:104 CBS:52 CRYAA:9; do
  name=${g%%:*}
  check_eq "${name}: benchmark variants in its CMRG regions" "$(awk -v g="$name" '$1 == g {print $3}' <<< "$R_MASKED")" "${g#*:}"
  check "${name}: recall on the masked chr21 is 0.9 or more ($(awk -v g="$name" '$1 == g {print $2 "/" $3}' <<< "$R_MASKED"))" \
    ok_recall "$R_MASKED" "$name"
  check "control: ${name}: recall on the stock chr21 is under 0.5 ($(awk -v g="$name" '$1 == g {print $2 "/" $3}' <<< "$R_STOCK"))" \
    low_recall "$R_STOCK" "$name"
done

# --- step 16b: the genes are read on the masked reference only -------------------------
CAT="${CASE_TMP}/callability_genes.tsv"
grep -v '^#' "${REPO}/assets/callability_genes.tsv" \
  | awk -F'\t' 'NR == 1 || $1 == "KCNE1" || $1 == "CBS" || $1 == "U2AF1" || $1 == "CRYAA"' > "$CAT"
check_eq "the four genes' rows of assets/callability_genes.tsv" "$(tail -n +2 "$CAT" | cut -f1 | sort | tr '\n' ' ')" "CBS CRYAA KCNE1 U2AF1 "
reading() {
  python3 - "$CAT" "$1" "${REPO}/bin" <<'PY'
import csv, sys
sys.path.insert(0, sys.argv[3])
import collect_summary
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
d = collect_summary.callability(rows, sys.argv[2])
bad = {g["gene"] for g in d["not_read"]}
print(" ".join(f"{r['gene']}={'not_read' if r['gene'] in bad else 'read'}" for r in rows))
PY
}
for side in stock masked; do
  S="ffx-${side}"
  CALLABILITY_GENES="$CAT" run_step 16b-mosdepth.sh "$S"
  check_step_exit "16b-mosdepth.sh (${side})"
done
C_MASKED=$(reading "${G}/ffx-masked/mosdepth/ffx-masked_callability.bed.gz" 2>&1); echo "masked: ${C_MASKED}"
C_STOCK=$(reading "${G}/ffx-stock/mosdepth/ffx-stock_callability.bed.gz" 2>&1); echo "stock: ${C_STOCK}"
check_eq "step 16b reads the four genes on the masked chr21" "$C_MASKED" "KCNE1=read CBS=read U2AF1=read CRYAA=read"
check_eq "control: step 16b reads none of them on the stock chr21" "$C_STOCK" "KCNE1=not_read CBS=not_read U2AF1=not_read CRYAA=not_read"

rm -rf "${RD}" "${G}/ffx-stock" "${G}/ffx-masked"
finish
