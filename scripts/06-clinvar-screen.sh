#!/usr/bin/env bash
# ClinVar Pathogenic Screen — intersect sample VCF with ClinVar pathogenic/LP variants
# Finds known disease-causing variants the person carries
#
# Output: clinvar/${SAMPLE}_clinvar_hits.vcf — the sample's records that match a
# ClinVar Pathogenic/Likely_pathogenic allele and whose genotype carries it,
# annotated with ClinVar's ID, GENEINFO, CLNSIG and CLNREVSTAT. Both reports
# (steps 24 and generate-report.sh) read this file and group the hits by ClinVar
# review stars. clinvar/${SAMPLE}_clinvar_hits.tsv holds the same hits, one row each.
# Two neighbouring records that can be one allele are also matched as their
# combined allele (bin/clinvar_join_adjacent.awk): a matched combined record is
# a hit of its own (INFO/PGP_JOINED_FROM) and its parts are marked PGP_PART_OF.
# clinvar/${SAMPLE}_clinvar_cdna_check.tsv says, for each hit, whether its alt
# reads look like spliced cDNA (bin/clip_check.py; needs the sample's BAM and
# the GENCODE exons setup.sh installs, and says "not checked" without them).
# A hit whose call has a genotype quality below CLINVAR_MIN_GQ (default 30) or
# a depth below CLINVAR_MIN_DP (default 10) is a low-confidence call: it stays
# in the hits, marked, and the reports count it apart from the others. A floor
# other than the default is written into the hits VCF as
# ##PGP_CALL_FLOOR=GQ=..,DP=.., so bin/clinvar_hits.awk and both reports read
# the floor this run used; without the line they use the default.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
VCF_DIR=${VCF_DIR:-vcf}
VCF="${GENOME_DIR}/${SAMPLE}/${VCF_DIR}/${SAMPLE}.vcf.gz"
REF="$REF_FASTA"
# Source file built by setup.sh; the normalised copy beside it is built from it.
CLINVAR="${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz"
CLINVAR_NORM="${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.norm.vcf.gz"
OUTPUT_DIR="${GENOME_DIR}/${SAMPLE}/clinvar"
HITS="${OUTPUT_DIR}/${SAMPLE}_clinvar_hits.vcf"
HITS_TSV="${OUTPUT_DIR}/${SAMPLE}_clinvar_hits.tsv"
CDNA_TSV="${OUTPUT_DIR}/${SAMPLE}_clinvar_cdna_check.tsv"
BAM="${GENOME_DIR}/${SAMPLE}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam"
EXONS_BED=$(data_file gencode_exons) || true
MIN_GQ=${CLINVAR_MIN_GQ:-30}
MIN_DP=${CLINVAR_MIN_DP:-10}
for v in "CLINVAR_MIN_GQ=${MIN_GQ}" "CLINVAR_MIN_DP=${MIN_DP}"; do
  [[ "${v#*=}" =~ ^[0-9]+$ ]] || { echo "ERROR: ${v%%=*} must be a whole number, got '${v#*=}'" >&2; exit 1; }
done

echo "=== ClinVar Pathogenic Screen: ${SAMPLE} ==="

for f in "$VCF" "${VCF}.tbi" "$CLINVAR" "${CLINVAR}.tbi" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

# bcftools_run [--rw DIR] ARGS...: bcftools in its container.
bcftools_run() {
  local -a opt=()
  if [ "$1" = "--rw" ]; then opt=(--rw "$2"); shift 2; fi
  run_in ${opt[@]+"${opt[@]}"} --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" "$@"
}

# Step 0: the contigs both files hold records on, that the reference has too.
# bcftools norm stops at the first record whose contig the reference lacks (a
# vendor scaffold, ClinVar's NT_ contigs), and a record on a contig the other
# file lacks cannot match anyway. No contig in common would give zero hits that
# look clean, so the step stops and says which file is named the other way.
REF_CONTIGS=$(cut -f1 "${REF}.fai" | sort -u)
SAMPLE_CONTIGS=$(bcftools_run bcftools index -s "$(cpath "$VCF")" | awk '$3 > 0 {print $1}' | sort -u)
CLINVAR_CONTIGS=$(bcftools_run bcftools index -s "$(cpath "$CLINVAR")" | awk '$3 > 0 {print $1}' | sort -u)
# ClinVar's side depends on the reference only, so the normalised copy below
# stays valid for every sample.
CLINVAR_KEEP=$(comm -12 <(printf '%s\n' "$CLINVAR_CONTIGS") <(printf '%s\n' "$REF_CONTIGS") | grep . || true)
SAMPLE_KEEP=$(comm -12 <(printf '%s\n' "$SAMPLE_CONTIGS") <(printf '%s\n' "$CLINVAR_KEEP") | grep . || true)
if [ -z "$SAMPLE_KEEP" ]; then
  echo "ERROR: ${VCF} and ${CLINVAR} have no contig name in common." >&2
  echo "  Sample contigs:  $(printf '%s\n' "$SAMPLE_CONTIGS" | head -3 | tr '\n' ' ')..." >&2
  echo "  ClinVar contigs: $(printf '%s\n' "$CLINVAR_CONTIGS" | head -3 | tr '\n' ' ')..." >&2
  if ! grep -q '^chr' <<< "$SAMPLE_CONTIGS"; then
    echo "  The sample VCF is not chr-named: rename its contigs (docs/vcf-first.md)." >&2
  elif ! grep -q '^chr' <<< "$CLINVAR_CONTIGS"; then
    echo "  Rebuild the ClinVar file with chr-prefixed names (scripts/setup.sh does this)." >&2
  fi
  exit 1
fi
# Targets with a constant end: a header without contig lengths cannot shorten them.
CLINVAR_TARGETS="${OUTPUT_DIR}/${SAMPLE}_clinvar_targets.tsv"
SAMPLE_TARGETS="${OUTPUT_DIR}/${SAMPLE}_sample_targets.tsv"
printf '%s\n' "$CLINVAR_KEEP" | awk -v OFS='\t' '{print $1, 1, 2147483647}' > "$CLINVAR_TARGETS"
printf '%s\n' "$SAMPLE_KEEP" | awk -v OFS='\t' '{print $1, 1, 2147483647}' > "$SAMPLE_TARGETS"
SAMPLE_DROPPED=$(comm -23 <(printf '%s\n' "$SAMPLE_CONTIGS") <(printf '%s\n' "$SAMPLE_KEEP") | grep . || true)

# Step 1: normalise ClinVar once (split multiallelics, left-align), rebuilt when the
# source is newer. isec matches on identical REF/ALT, so both sides must be split
# and left-aligned the same way.
if [ ! -f "$CLINVAR_NORM" ] || [ ! -f "${CLINVAR_NORM}.tbi" ] || [ "$CLINVAR" -nt "$CLINVAR_NORM" ]; then
  echo "Normalising ClinVar into ${CLINVAR_NORM} ..."
  rm -f "${CLINVAR_NORM}.part.vcf.gz" "${CLINVAR_NORM}.part.vcf.gz.tbi"
  # The normalised copy is shared by every sample, so clinvar/ is writable here.
  run_in --rw "$(dirname "$CLINVAR_NORM")" --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" \
    sh -c "set -e; bcftools view -T '$(cpath "$CLINVAR_TARGETS")' -Ou '$(cpath "$CLINVAR")' \
      | bcftools norm -m -any -c w -f '$(cpath "$REF")' -Oz -o '$(cpath "${CLINVAR_NORM}.part.vcf.gz")' -"
  bcftools_run --rw "$(dirname "$CLINVAR_NORM")" bcftools index -t "$(cpath "${CLINVAR_NORM}.part.vcf.gz")"
  mv "${CLINVAR_NORM}.part.vcf.gz" "$CLINVAR_NORM"
  mv "${CLINVAR_NORM}.part.vcf.gz.tbi" "${CLINVAR_NORM}.tbi"
fi

# Step 2: filter the sample. Callers that never write PASS (FILTER '.') would lose
# every record under -f PASS, so fall back to '.,PASS' only when no record is PASS.
HAS_PASS=$(run_in --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" \
  sh -c "bcftools view -H -f PASS '$(cpath "$VCF")' | head -n 1 | wc -l")
if [ "$HAS_PASS" -gt 0 ]; then
  FILTER="PASS"
  echo "Filter mode: PASS records only"
else
  FILTER=".,PASS"
  echo "NOTICE: ${VCF} has no PASS record; using records with FILTER '.' or PASS"
fi

if [ -n "$SAMPLE_DROPPED" ]; then
  DROPPED_REGIONS="${OUTPUT_DIR}/${SAMPLE}_dropped_contigs.tsv"
  printf '%s\n' "$SAMPLE_DROPPED" | awk -v OFS='\t' '{print $1, 1, 2147483647}' > "$DROPPED_REGIONS"
  LEFT_OUT=$(run_in --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" \
    sh -c "bcftools view -H -f '${FILTER}' -R '$(cpath "$DROPPED_REGIONS")' '$(cpath "$VCF")' | wc -l" | tr -d ' ')
  echo "NOTICE: ${LEFT_OUT:-0} '${FILTER}' records left out: they lie on $(grep -c . <<< "$SAMPLE_DROPPED") contigs that ClinVar or the reference lacks (first ones: $(head -5 <<< "$SAMPLE_DROPPED" | paste -sd, -))"
  rm -f "$DROPPED_REGIONS"
fi

PASS_VCF="${OUTPUT_DIR}/${SAMPLE}_pass.vcf.gz"
run_in --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" \
  sh -c "set -e; bcftools view -f '${FILTER}' -T '$(cpath "$SAMPLE_TARGETS")' -Ou '$(cpath "$VCF")' \
    | bcftools norm -m -any -c w -f '$(cpath "$REF")' -Oz -o '$(cpath "$PASS_VCF")' -
    bcftools index -f -t '$(cpath "$PASS_VCF")'"

PASS_COUNT=$(bcftools_run bcftools index -n "$(cpath "$PASS_VCF")")
if [ "$PASS_COUNT" -eq 0 ]; then
  echo "ERROR: no records left in ${VCF} after the '${FILTER}' filter; nothing to screen." >&2
  exit 1
fi

# Step 2b: the combined alleles. A caller can write one complex allele (a
# delins) as two records a few bases apart on the same haplotype; matched one by
# one, the parts can hit single-submitter ClinVar records and miss the better
# reviewed record ClinVar holds for the whole allele. bcftools norm -m +any joins
# only records at the same position, so bin/clinvar_join_adjacent.awk builds,
# for each pair of neighbours that can be one allele (the second starts at most
# 10 bases after the first ends, same phase set and GT, or the same unphased
# GT), the allele from the first part's start to the second part's end, with the
# reference bases between the parts (+fill-from-fasta), left-aligned the way
# ClinVar is. The candidates matched below are the split calls plus these.
JOINED="${OUTPUT_DIR}/${SAMPLE}_joined.vcf.gz"
CANDIDATES="${OUTPUT_DIR}/${SAMPLE}_candidates.vcf.gz"
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" \
  bash -o pipefail -c "set -e
    bcftools view '$(cpath "$PASS_VCF")' \
      | awk -v mode=prepare -f /pgp-bin/clinvar_join_adjacent.awk \
      | bcftools +fill-from-fasta --no-version -Ov - -- -c INFO/PGP_JOIN_SEQ -f '$(cpath "$REF")' \
      | awk -v mode=finish -f /pgp-bin/clinvar_join_adjacent.awk \
      | bcftools norm --no-version -c w -d exact -f '$(cpath "$REF")' -Oz -o '$(cpath "$JOINED")' -
    bcftools index -f -t '$(cpath "$JOINED")'
    bcftools concat --no-version -a -Oz -o '$(cpath "$CANDIDATES")' '$(cpath "$PASS_VCF")' '$(cpath "$JOINED")'
    bcftools index -f -t '$(cpath "$CANDIDATES")'"
echo "Combined alleles of neighbouring records: $(bcftools_run bcftools index -n "$(cpath "$JOINED")")"

# Step 3: keep the sample's records whose allele is in ClinVar, then copy ClinVar's
# ID, gene and significance onto them. isec -w1 writes the sample's side, which has
# no GENEINFO/CLNSIG of its own. (annotate -a needs an indexed target, so the
# shared records go to a file first.) Only records whose genotype carries an ALT
# allele are hits: the match is by allele, so after the split above a 0/0 or ./.
# record, the 0/0 half of a multiallelic record and a reference-only ALT '.' row
# would otherwise be listed. The test is a non-zero allele index anywhere in GT,
# not GT="alt", which drops a half call such as ./1 that does carry the ALT allele.
# --no-version keeps command lines out of the header. Last, the tag pass of
# bin/clinvar_join_adjacent.awk marks each part of a matched combined allele
# (PGP_PART_OF), and two hits in one gene that can be one allele but matched no
# combined record (PGP_ONE_ALLELE).
SHARED="${OUTPUT_DIR}/${SAMPLE}_shared.vcf.gz"
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" \
  sh -c "set -e
    bcftools isec --no-version -n=2 -w1 -Oz -o '$(cpath "$SHARED")' '$(cpath "$CANDIDATES")' '$(cpath "$CLINVAR_NORM")'
    bcftools index -f -t '$(cpath "$SHARED")'
    bcftools annotate --no-version -a '$(cpath "$CLINVAR_NORM")' --pair-logic exact \
      -c ID,INFO/GENEINFO,INFO/CLNSIG,INFO/CLNREVSTAT -Ou '$(cpath "$SHARED")' \
      | bcftools view --no-version -i 'GT~\"[1-9]\"' -Ov -o '$(cpath "${HITS}.untagged")'
    awk -v mode=tag -f /pgp-bin/clinvar_join_adjacent.awk '$(cpath "${HITS}.untagged")' '$(cpath "${HITS}.untagged")' \
      > '$(cpath "${HITS}.part")'
    rm -f '$(cpath "${HITS}.untagged")'
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\t%ID\\t%INFO/GENEINFO\\t%INFO/CLNSIG\\t%INFO/CLNREVSTAT\\n' \
      -o '$(cpath "${HITS_TSV}.body")' '$(cpath "${HITS}.part")'"
# A floor other than the default, for bin/clinvar_hits.awk and the reports. At
# the default the hits file stays as it was before the floor existed.
if [ "$MIN_GQ" -ne 30 ] || [ "$MIN_DP" -ne 10 ]; then
  awk -v l="##PGP_CALL_FLOOR=GQ=${MIN_GQ},DP=${MIN_DP}" '/^#CHROM/ {print l} {print}' "${HITS}.part" > "${HITS}.floor"
  mv "${HITS}.floor" "${HITS}.part"
fi
# The same hits as a table, one row per hit.
{
  printf 'chrom\tpos\tref\talt\tgenotype\tclinvar_id\tgeneinfo\tclnsig\tclnrevstat\n'
  cat "${HITS_TSV}.body"
} > "${HITS_TSV}.part"
mv "${HITS}.part" "$HITS"
mv "${HITS_TSV}.part" "$HITS_TSV"
rm -f "$SHARED" "${SHARED}.tbi" "${HITS_TSV}.body" "$CLINVAR_TARGETS" "$SAMPLE_TARGETS" \
  "$JOINED" "${JOINED}.tbi" "$CANDIDATES" "${CANDIDATES}.tbi"

# Step 4: the cDNA read check. Spliced cDNA or expression-vector DNA (a lab
# contaminant) gives reads clipped exactly at an exon edge, with the rest of
# the read on the next exon; a variant it carries is called like any other.
# For each hit, bin/clip_check.py counts the alt and ref reads clipped at an
# exon edge and flags the hit when most alt reads are and few ref reads are.
# It reads only the reads at the hits. pysam is in PYPGX_IMAGE.
if [ ! -f "$BAM" ] || [ ! -f "${BAM}.bai" ]; then
  CDNA_ARGS=(--skip "no BAM (${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam with its .bai)")
elif [ ! -s "$EXONS_BED" ]; then
  CDNA_ARGS=(--skip "no GENCODE exons (scripts/setup.sh installs them)")
else
  CDNA_ARGS=(--bam "$(cpath "$BAM")" --exons "$(cpath "$EXONS_BED")")
fi
run_in --cpus 1 --memory 2g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYPGX_IMAGE" \
  python3 /pgp-bin/clip_check.py --hits "$(cpath "$HITS_TSV")" --out "$(cpath "${CDNA_TSV}.part")" "${CDNA_ARGS[@]}"
mv "${CDNA_TSV}.part" "$CDNA_TSV"

HIT_COUNT=$(grep -c -v '^#' "$HITS" || true)
echo "=== ClinVar screen complete ==="
echo "Hits: ${HITS} (sample records carrying a ClinVar Pathogenic/Likely_pathogenic allele)"
echo "Table: ${HITS_TSV}"
echo "cDNA read check: ${CDNA_TSV}"
echo "Count: ${HIT_COUNT} pathogenic hits"
if [ "$HIT_COUNT" -gt 0 ]; then
  LOW_COUNT=$(awk -f "${PGP_ROOT}/bin/clinvar_hits.awk" "$HITS" | awk -F'\t' '$13 == "yes"' | wc -l | tr -d ' ')
  echo "Low-confidence calls among them (GQ below ${MIN_GQ} or DP below ${MIN_DP}): ${LOW_COUNT}"
fi
# A zero-star submission counts as a hit like an expert-panel one; the review
# status (CLNREVSTAT, copied on above) says how much weight each deserves.
# bin/clinvar_hits.awk maps it to ClinVar's stars, as both reports do, and its
# 10th column says when a hit is a combined allele or part of one, and the
# last when the call is below the floor.
if [ "$HIT_COUNT" -gt 0 ]; then
  echo "By review status (ClinVar stars):"
  awk -f "${PGP_ROOT}/bin/clinvar_hits.awk" "$HITS" | sort -t$'\t' -k1,1nr -k2,2V -k3,3n \
    | awk -F'\t' '{n[$1]++; line[$1] = line[$1] sprintf("    %s %s:%s %s>%s %s %s [%s]%s%s\n", $7, $2, $3, $4, $5, $6, $8, $9, ($10 != "." && $10 != "" ? " (" $10 ")" : ""), ($13 == "yes" ? " LOW-CONFIDENCE CALL (GQ " $12 ", DP " $11 "): check the reads" : ""))}
        END {for (s = 4; s >= 0; s--) if (n[s]) printf "  %d star%s: %d\n%s", s, (s == 1 ? "" : "s"), n[s], line[s]}'
fi
