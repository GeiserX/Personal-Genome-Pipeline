#!/usr/bin/env bash
# build-mito-confirmed.sh: write assets/mito_confirmed.tsv, the confirmed mtDNA
# disease variants step 20 names in the report, from ClinVar and the primary
# papers.
#
# The curated part is the CURATED table below: each variant, the disease it
# was first described in, the PubMed ids of the primary papers and a note on
# heteroplasmy. ClinVar gives the rest: the variation id of the same allele,
# the gene and the HGVS name. A variant whose ClinVar record is missing, or is
# no longer classed pathogenic, stops the build, so the list cannot keep a
# variant ClinVar has dropped. MITOMAP is not read: it refuses scripted
# requests and its redistribution terms are not settled (docs/20-mtoolbox.md).
#
# Only stable ClinVar fields go into the table (the variation id, the gene, the
# HGVS name), so any ClinVar release that holds these records writes the same
# file. The review status and the classification of the day come from the
# sample's own annotation table (bin/collect_summary.py mito-annotate).
#
# Usage:
#   scripts/ci/build-mito-confirmed.sh --clinvar VCF           rewrite assets/mito_confirmed.tsv
#   scripts/ci/build-mito-confirmed.sh --clinvar VCF --check   exit 1 unless the committed file is what VCF gives
#   --out FILE                                                 write FILE instead
#
# VCF is ClinVar's GRCh38 VCF with its .tbi (setup.sh installs
# clinvar/clinvar.vcf.gz and clinvar/clinvar_chr.vcf.gz; either works, and so
# does https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz with
# a local bcftools). bcftools runs from PATH, else in BCFTOOLS_IMAGE.
set -euo pipefail
export LC_ALL=C

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
OUT="${ROOT}/assets/mito_confirmed.tsv"
CLINVAR="" CHECK=0
while [ $# -gt 0 ]; do
  case "$1" in
    --clinvar) CLINVAR=${2:?--clinvar needs a VCF}; shift ;;
    --out) OUT=${2:?--out needs a file}; shift ;;
    --check) CHECK=1 ;;
    *) echo "usage: $0 --clinvar VCF [--check] [--out FILE]" >&2; exit 2 ;;
  esac
  shift
done
[ -n "$CLINVAR" ] || { echo "usage: $0 --clinvar VCF [--check] [--out FILE]" >&2; exit 2; }

# kind, position (deletion: the VCF anchor base), ref, alt (deletion: its length),
# disease, PubMed ids (";" between), heteroplasmy note.
CURATED='snv	3243	A	G	MELAS (mitochondrial encephalomyopathy, lactic acidosis and stroke-like episodes); maternally inherited diabetes and deafness	2102678	Usually heteroplasmic; the level differs between tissues and is often lower in blood than in muscle or urine.
snv	8344	A	G	MERRF (myoclonic epilepsy with ragged-red fibres)	2112427	Usually heteroplasmic; the level differs between tissues.
snv	11778	G	A	Leber hereditary optic neuropathy (LHON)	3201231	Usually homoplasmic; not everyone who carries it loses vision.
snv	3460	G	A	Leber hereditary optic neuropathy (LHON)	1928099	Usually homoplasmic; not everyone who carries it loses vision.
snv	14484	T	C	Leber hereditary optic neuropathy (LHON)	1417830	Usually homoplasmic; not everyone who carries it loses vision.
deletion	8469	-	4977	The common 4977 bp deletion, between two 13 bp direct repeats: Kearns-Sayre syndrome, chronic progressive external ophthalmoplegia (CPEO), Pearson syndrome	2711184;2541333;20301382	A single large deletion is heteroplasmic and usually arises new rather than inherited; in Kearns-Sayre syndrome and CPEO it is often found in muscle and not in blood, in Pearson syndrome in blood.'

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

bcf() {
  if command -v bcftools >/dev/null 2>&1; then
    bcftools "$@"
  else
    # shellcheck source=/dev/null
    . "${ROOT}/versions.env"
    local dir
    dir=$(cd "$(dirname "$CLINVAR")" && pwd)
    docker run --rm -u "$(id -u):$(id -g)" -v "${dir}:${dir}:ro" "${BCFTOOLS_IMAGE:?versions.env sets BCFTOOLS_IMAGE}" bcftools "$@"
  fi
}

case "$CLINVAR" in
  http://*|https://*) CV=$CLINVAR ;;
  *) [ -s "$CLINVAR" ] && [ -s "${CLINVAR}.tbi" ] || { echo "ERROR: no ${CLINVAR} with its .tbi" >&2; exit 2; }
     CV="$(cd "$(dirname "$CLINVAR")" && pwd)/$(basename "$CLINVAR")" ;;
esac
# ClinVar's own file names the contig MT, setup.sh's clinvar_chr.vcf.gz chrM.
CHROM=$(bcf index -s "$CV" 2>/dev/null | awk '$1 == "MT" || $1 == "chrM" {print $1; exit}')
[ -n "$CHROM" ] || { echo "ERROR: ${CLINVAR} has no MT or chrM records" >&2; exit 1; }

{
  echo "# assets/mito_confirmed.tsv: confirmed mtDNA disease variants, named in step 20's report section."
  echo "# Written by scripts/ci/build-mito-confirmed.sh from ClinVar (variation id, gene, HGVS) and the primary"
  echo "# papers (disease, pmids, note): edit the script's CURATED table, then rerun it; never edit this file."
  echo "# pos, end: rCRS (chrM) positions of the changed bases; a deletion's are its first and last deleted base."
  printf 'kind\tpos\tend\tref\talt\thgvs\tgene\tdisease\tclinvar_id\tpmids\theteroplasmy_note\n'
} > "${TMP}/out.tsv"

ERR=0
while IFS=$'\t' read -r kind pos ref alt disease pmids note; do
  bcf view -H -r "${CHROM}:${pos}" "$CV" > "${TMP}/rec.vcf"
  # The record of this allele: an SNV by its bases; a deletion by its anchor and length.
  rec=$(awk -F'\t' -v k="$kind" -v p="$pos" -v r="$ref" -v a="$alt" '
    $2 != p {next}
    k == "snv" && $4 == r && $5 == a {print; exit}
    k == "deletion" && length($5) == 1 && length($4) - 1 == a {print; exit}' "${TMP}/rec.vcf")
  if [ -z "$rec" ]; then
    echo "ERROR: ClinVar has no record of ${kind} ${pos} ${ref} ${alt}" >&2; ERR=1; continue
  fi
  info() { awk -F'\t' -v k="$1" '{n = split($8, f, ";"); for (i = 1; i <= n; i++) if (index(f[i], k "=") == 1) {print substr(f[i], length(k) + 2); exit}}' <<< "$rec"; }
  id=$(cut -f3 <<< "$rec")
  sig=$(info CLNSIG)
  case "$sig" in
    Pathogenic|Pathogenic/Likely_pathogenic|Pathogenic\|*) ;;
    *) echo "ERROR: ClinVar ${id} (${kind} ${pos}) is classed '${sig}', not pathogenic" >&2; ERR=1; continue ;;
  esac
  hgvs=$(info CLNHGVS); hgvs=${hgvs#*:}
  gene=$(info GENEINFO | awk -F'|' '{for (i = 1; i <= NF; i++) {split($i, g, ":"); printf "%s%s", (i > 1 ? "," : ""), g[1]}}')
  if [ "$kind" = snv ]; then
    start=$pos end=$pos
  else
    start=$((pos + 1)) end=$((pos + alt)) ref=. alt=.
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$kind" "$start" "$end" "$ref" "$alt" "$hgvs" "${gene:-.}" \
    "$disease" "$id" "$pmids" "$note" >> "${TMP}/out.tsv"
done <<< "$CURATED"
[ "$ERR" = 0 ] || exit 1

if [ "$CHECK" = 1 ]; then
  if cmp -s "${TMP}/out.tsv" "$OUT"; then
    echo "OK: ${OUT#"${ROOT}"/} is what ${CLINVAR} gives ($(($(grep -vc '^#' "$OUT") - 1)) variants)"
    exit 0
  fi
  echo "ERROR: ${OUT#"${ROOT}"/} is not what ${CLINVAR} gives; rerun scripts/ci/build-mito-confirmed.sh --clinvar VCF:" >&2
  diff "$OUT" "${TMP}/out.tsv" >&2 || true
  exit 1
fi
cp "${TMP}/out.tsv" "$OUT"
echo "Wrote ${OUT#"${ROOT}"/}: $(($(grep -vc '^#' "$OUT") - 1)) variants"
