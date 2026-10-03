#!/usr/bin/env bash
# ClinVar Pathogenic Screen — intersect sample VCF with ClinVar pathogenic/LP variants
# Finds known disease-causing variants the person carries
#
# Output: clinvar/${SAMPLE}_clinvar_hits.vcf — the sample's records that match a
# ClinVar Pathogenic/Likely_pathogenic allele and whose genotype carries it,
# annotated with ClinVar's ID, GENEINFO, CLNSIG and CLNREVSTAT. Both reports
# (steps 24 and generate-report.sh) read this file and group the hits by ClinVar
# review stars. clinvar/${SAMPLE}_clinvar_hits.tsv holds the same hits, one row each.
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

# Step 3: keep the sample's records whose allele is in ClinVar, then copy ClinVar's
# ID, gene and significance onto them. isec -w1 writes the sample's side, which has
# no GENEINFO/CLNSIG of its own. (annotate -a needs an indexed target, so the
# shared records go to a file first.) Only records whose genotype carries an ALT
# allele are hits: the match is by allele, so after the split above a 0/0 or ./.
# record, the 0/0 half of a multiallelic record and a reference-only ALT '.' row
# would otherwise be listed. The test is a non-zero allele index anywhere in GT,
# not GT="alt", which drops a half call such as ./1 that does carry the ALT allele.
# --no-version keeps command lines out of the header.
SHARED="${OUTPUT_DIR}/${SAMPLE}_shared.vcf.gz"
run_in --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" \
  sh -c "set -e
    bcftools isec --no-version -n=2 -w1 -Oz -o '$(cpath "$SHARED")' '$(cpath "$PASS_VCF")' '$(cpath "$CLINVAR_NORM")'
    bcftools index -f -t '$(cpath "$SHARED")'
    bcftools annotate --no-version -a '$(cpath "$CLINVAR_NORM")' --pair-logic exact \
      -c ID,INFO/GENEINFO,INFO/CLNSIG,INFO/CLNREVSTAT -Ou '$(cpath "$SHARED")' \
      | bcftools view --no-version -i 'GT~\"[1-9]\"' -Ov -o '$(cpath "${HITS}.part")'
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\t%ID\\t%INFO/GENEINFO\\t%INFO/CLNSIG\\t%INFO/CLNREVSTAT\\n' \
      -o '$(cpath "${HITS_TSV}.body")' '$(cpath "${HITS}.part")'"
# The same hits as a table, one row per hit.
{
  printf 'chrom\tpos\tref\talt\tgenotype\tclinvar_id\tgeneinfo\tclnsig\tclnrevstat\n'
  cat "${HITS_TSV}.body"
} > "${HITS_TSV}.part"
mv "${HITS}.part" "$HITS"
mv "${HITS_TSV}.part" "$HITS_TSV"
rm -f "$SHARED" "${SHARED}.tbi" "${HITS_TSV}.body" "$CLINVAR_TARGETS" "$SAMPLE_TARGETS"

HIT_COUNT=$(grep -c -v '^#' "$HITS" || true)
echo "=== ClinVar screen complete ==="
echo "Hits: ${HITS} (sample records carrying a ClinVar Pathogenic/Likely_pathogenic allele)"
echo "Table: ${HITS_TSV}"
echo "Count: ${HIT_COUNT} pathogenic hits"
# A zero-star submission counts as a hit like an expert-panel one; the review
# status (CLNREVSTAT, copied on above) says how much weight each deserves.
# bin/clinvar_hits.awk maps it to ClinVar's stars, as both reports do.
if [ "$HIT_COUNT" -gt 0 ]; then
  echo "By review status (ClinVar stars):"
  awk -f "${PGP_ROOT}/bin/clinvar_hits.awk" "$HITS" | sort -t$'\t' -k1,1nr -k2,2V -k3,3n \
    | awk -F'\t' '{n[$1]++; line[$1] = line[$1] sprintf("    %s %s:%s %s>%s %s %s [%s]\n", $7, $2, $3, $4, $5, $6, $8, $9)}
        END {for (s = 4; s >= 0; s--) if (n[s]) printf "  %d star%s: %d\n%s", s, (s == 1 ? "" : "s"), n[s], line[s]}'
fi
