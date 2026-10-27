#!/usr/bin/env bash
# check-rule-alleles.sh: every site of assets/genotype_rules.tsv against dbSNP.
#
# A rule's haplotypes are written as GRCh38 plus-strand bases. TAS2R38 and
# many other genes lie on the minus strand, where the papers name the codon
# bases, the complement; and a ref/alt pair copied the wrong way round reads
# every genotype backwards. Both mistakes give alleles that disagree with
# dbSNP's GRCh38 record. This check compares each site line (rsID, chrom,
# pos, ref, alt) with the dbSNP table at the end of this file, the records
# of the rule sites in dbSNP's GRCh38 VCF of the build CONSUMER_DBSNP_BUILD
# names in versions.env:
#   - the rsID must be in the table, at the same chrom and position;
#   - ref must be dbSNP's REF, and alt one of dbSNP's ALT alleles. A pair that
#     is dbSNP's the wrong way round, or its complement, is named as such.
# A C/G or A/T site keeps both alleles under the complement, so a strand
# mistake in its haplotype alleles cannot show here: the check prints a note
# naming each such site.
#
# Usage:
#   scripts/ci/check-rule-alleles.sh [RULES]         the check, offline (RULES defaults to
#                                                    assets/genotype_rules.tsv)
#   scripts/ci/check-rule-alleles.sh --online [RULES]
#                                                    also: the table below is what dbSNP's
#                                                    VCF holds for those rsIDs, read over
#                                                    HTTPS with bcftools (BCFTOOLS_HTTPS_IMAGE,
#                                                    Docker); the e2e case genotype-rules-1-phase runs it
#   scripts/ci/check-rule-alleles.sh --print-table [RULES]
#                                                    print the table lines of the rule sites
#                                                    from dbSNP (Docker), to put at the end of
#                                                    this file, each after "#: ", when a rule
#                                                    adds a site or the build changes
# Exit 0 when every site matches, 1 when one does not, when RULES has no site
# (an empty check is a failure), or when the table's build is not the one
# versions.env pins. tests/test_traits_report.py runs it on the table and on
# planted copies with a swapped pair, a minus-strand pair and an unknown rsID.
set -euo pipefail
export LC_ALL=C   # one sort order for the table and the dbSNP lines

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
MODE=check
case "${1:-}" in
  --online) MODE=online; shift ;;
  --print-table) MODE=print; shift ;;
  -h|--help) sed -n '2,/^set -euo/p' "$0" | sed '$d; s/^# \{0,1\}//'; exit 0 ;;
esac
RULES=${1:-${ROOT}/assets/genotype_rules.tsv}
[ -f "$RULES" ] || { echo "ERROR: ${RULES} not found" >&2; exit 1; }
# shellcheck source=../../versions.env
. "${ROOT}/versions.env"
: "${CONSUMER_DBSNP_BUILD:?versions.env has no CONSUMER_DBSNP_BUILD}"
DBSNP_URL="https://ftp.ncbi.nih.gov/snp/archive/${CONSUMER_DBSNP_BUILD}/VCF/GCF_000001405.40.gz"

# The rule sites: rsid chrom pos ref alt, one per rsID.
SITES=$(awk -F'\t' '$1 == "site" { print $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 }' "$RULES" | sort -u)
[ -n "$SITES" ] || { echo "ERROR: ${RULES} has no site line: nothing to check" >&2; exit 1; }

# The pinned table: the "#: " lines at the end of this file.
TABLE=$(sed -n 's/^#: //p' "$0")
BUILD=$(sed -n 's/^# build //p' <<< "$TABLE")
if [ "$BUILD" != "$CONSUMER_DBSNP_BUILD" ]; then
  echo "ERROR: the dbSNP table in $0 is of build '${BUILD}', versions.env pins ${CONSUMER_DBSNP_BUILD}." >&2
  echo "  Refresh it: $0 --print-table, and put each line it prints at the end of this file after \"#: \"." >&2
  exit 1
fi

# dbsnp_lines: the table lines of the rule sites, read from dbSNP's VCF over
# HTTPS. The VCF names contigs by RefSeq accession; its header maps them.
dbsnp_lines() {
  local img=${BCFTOOLS_HTTPS_IMAGE:?versions.env has no BCFTOOLS_HTTPS_IMAGE} regions acc
  acc=$(docker run --rm "$img" bcftools view -h "$DBSNP_URL" \
    | sed -n 's/^##contig=<ID=\(NC_0000\([0-9][0-9]\)\.[0-9]*\).*/\2 \1/p' \
    | awk '{ printf "chr%d\t%s\n", $1 + 0, $2 }')
  [ -n "$acc" ] || { echo "ERROR: could not read the contigs of ${DBSNP_URL}" >&2; return 1; }
  regions=$(awk -F'\t' 'NR == FNR { a[$1] = $2; next } ($2 in a) { printf "%s%s:%s", (n++ ? "," : ""), a[$2], $3 }' \
    <(echo "$acc") <(echo "$SITES"))
  [ -n "$regions" ] || { echo "ERROR: no rule site is on a dbSNP contig" >&2; return 1; }
  echo "# build ${CONSUMER_DBSNP_BUILD}"
  docker run --rm "$img" bcftools query -r "$regions" -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' "$DBSNP_URL" \
    | awk -F'\t' 'NR == FNR { chr[$2] = $1; next } { print $1 "\t" chr[$2] "\t" $3 "\t" $4 "\t" $5 }' <(echo "$acc") - \
    | awk -F'\t' 'NR == FNR { want[$1] = 1; next } ($1 in want)' <(echo "$SITES") - \
    | sort -u
}

if [ "$MODE" = print ]; then
  dbsnp_lines
  exit $?
fi

FAILED=0
# The check: each rule site against the table.
OUT=$(awk -F'\t' '
  function comp(b) { return b == "A" ? "T" : b == "T" ? "A" : b == "C" ? "G" : b == "G" ? "C" : b }
  NR == FNR { if ($0 ~ /^#/ || $0 == "") next; chrom[$1] = $2; pos[$1] = $3; ref[$1] = $4; alts[$1] = "," $5 ","; next }
  {
    rs = $1; n++
    if (!(rs in ref)) { print "FAIL\t" rs ": not in the pinned dbSNP table (add it with --print-table)"; next }
    if ($2 != chrom[rs] || $3 != pos[rs]) {
      print "FAIL\t" rs ": the rules place it at " $2 ":" $3 ", dbSNP at " chrom[rs] ":" pos[rs]; next }
    a = alts[rs]
    if ($4 == ref[rs] && index(a, "," $5 ",")) {
      if (comp($4) == $5) print "NOTE\t" rs " " $4 "/" $5 ": the complement of the pair is the pair itself, so this check cannot catch a strand mistake in its haplotype alleles; check them against the codons"
      print "OK\t" rs " " $2 ":" $3 " " $4 ">" $5; next }
    if ($5 == ref[rs] && index(a, "," $4 ",")) {
      print "FAIL\t" rs ": ref and alt are swapped (rules " $4 ">" $5 ", dbSNP REF " ref[rs] ")"; next }
    if (comp($4) == ref[rs] && index(a, "," comp($5) ",")) {
      print "FAIL\t" rs ": the alleles are the minus-strand bases (rules " $4 ">" $5 ", dbSNP " ref[rs] ">" substr(a, 2, length(a) - 2) "): write the plus-strand bases"; next }
    print "FAIL\t" rs ": rules " $4 ">" $5 ", dbSNP REF " ref[rs] " ALT " substr(a, 2, length(a) - 2)
  }
  END { print "SITES\t" n + 0 }' <(echo "$TABLE") <(echo "$SITES"))
while IFS=$'\t' read -r kind msg; do
  case "$kind" in
    OK) echo "ok: ${msg}" ;;
    NOTE) echo "note: ${msg}" ;;
    FAIL) echo "FAIL: ${msg}"; FAILED=1 ;;
    SITES) N=$msg ;;
  esac
done <<< "$OUT"
echo "${N:-0} rule sites checked against dbSNP ${BUILD} (${RULES#"${ROOT}"/})"
[ "${N:-0}" -gt 0 ] || { echo "FAIL: no site was checked"; FAILED=1; }

if [ "$MODE" = online ]; then
  if ! LIVE=$(dbsnp_lines); then
    echo "FAIL: could not read dbSNP ${BUILD} over HTTPS"; FAILED=1
  elif [ "$LIVE" = "$(grep -v '^$' <<< "$TABLE")" ]; then
    echo "ok: the pinned table is what dbSNP ${BUILD} holds for the $(grep -vc '^#' <<< "$LIVE") rule rsIDs"
  else
    echo "FAIL: the pinned table differs from dbSNP ${BUILD} (< pinned, > dbSNP):"
    diff <(grep -v '^$' <<< "$TABLE") <(echo "$LIVE") || true
    FAILED=1
  fi
fi
exit "$FAILED"

# The dbSNP table: one "#: " line per rsID (rsid, chrom, pos, REF, ALT comma
# separated, tab separated), as --print-table writes them, after the build.
#: # build b157
#: rs10246939	chr7	141972804	T	C
#: rs1726866	chr7	141972905	G	A,C,T
#: rs429358	chr19	44908684	T	C
#: rs713598	chr7	141973545	C	A,G,T
#: rs7412	chr19	44908822	C	T
