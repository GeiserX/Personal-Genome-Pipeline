#!/usr/bin/env bash
# Lists every external http(s) address in a src= or href= attribute of the HTML
# reports the run produced (pipeline report, PharmCAT, Nextflow, MultiQC when
# present). src= and <link href=> are fetched when the file is opened; <a href=>
# is a plain link. The list goes to the job summary so the remote assets a
# report pulls in are known. Inline JavaScript that builds a URL is not seen.
. "$(dirname "$0")/lib.sh"

mapfile -t REPORTS < <(find "$GENOME_DIR" -type f -name '*.html' ! -path '*/nf-work/*' | LC_ALL=C sort)
check_ge "HTML reports found" "${#REPORTS[@]}" 2

{
  echo "### External addresses in the generated HTML reports"
  echo
  echo "| Report | Attribute | Address |"
  echo "|---|---|---|"
  for f in "${REPORTS[@]}"; do
    rel="${f#"${GENOME_DIR}"/}"
    grep -oiE '(src|href)[[:space:]]*=[[:space:]]*["'\'']?https?://[^"'\'' >)]+' "$f" 2>/dev/null \
      | sed -E 's/^([A-Za-z]+)[[:space:]]*=[[:space:]]*["'\'']?/\1\t/' \
      | LC_ALL=C sort -u \
      | awk -F'\t' -v r="$rel" '{printf "| %s | %s | %s |\n", r, tolower($1), $2}'
  done
  echo
} > "${CASE_TMP}/assets.md"
cat "${CASE_TMP}/assets.md" | tee -a "$E2E_NOTES"
echo "addresses listed: $(grep -c '^| [^R-]' "${CASE_TMP}/assets.md" || true)"

finish
