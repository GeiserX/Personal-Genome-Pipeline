#!/usr/bin/env bash
# Lists every external http(s) address in a src= or href= attribute of the HTML
# reports the run produced (pipeline report, PharmCAT, Nextflow, MultiQC when
# present). src= and <link href=> are fetched when the file is opened; <a href=>
# is a plain link. The list goes to the job summary so the remote assets a
# report pulls in are known. Inline JavaScript that builds a URL is not seen.
. "$(dirname "$0")/lib.sh"

mapfile -t REPORTS < <(find "$GENOME_DIR" -type f -name '*.html' ! -path '*/nf-work/*' | LC_ALL=C sort)
check_ge "HTML reports found" "${#REPORTS[@]}" 2

# One row per address: how it is referenced and from which reports (a
# directory with many pages, like indexcov's, is named once).
for f in "${REPORTS[@]}"; do
  rel="${f#"${GENOME_DIR}"/}"
  grep -oiE '(src|href)[[:space:]]*=[[:space:]]*["'\'']?https?://[^"'\'' >)]+' "$f" 2>/dev/null \
    | sed -E 's/^([A-Za-z]+)[[:space:]]*=[[:space:]]*["'\'']?/\1\t/' \
    | awk -F'\t' -v r="$rel" 'BEGIN {OFS = "\t"} {print tolower($1), $2, r}'
done > "${CASE_TMP}/assets.tsv"
{
  echo "### External addresses in the generated HTML reports"
  echo
  echo "src= and <link href=> load when the file is opened; <a href=> is a plain link."
  echo
  echo "| Address | Attribute | Found in |"
  echo "|---|---|---|"
  LC_ALL=C sort -u "${CASE_TMP}/assets.tsv" | awk -F'\t' '
    {n = split($3, p, "/"); d = substr($3, 1, length($3) - length(p[n])); k = $2 "\t" $1
     if (!((k, d) in cnt)) {first[k, d] = $3; dirs[k] = dirs[k] SUBSEP d}
     cnt[k, d]++}
    END {for (k in dirs) {m = split(substr(dirs[k], 2), ds, SUBSEP); w = ""
           for (i = 1; i <= m; i++) {e = (cnt[k, ds[i]] > 1) ? ds[i] " (" cnt[k, ds[i]] " files)" : first[k, ds[i]]
                                    w = w (w ? ", " : "") e}
           split(k, a, "\t"); printf "| %s | %s | %s |\n", a[1], a[2], w}}' | LC_ALL=C sort
  echo
} > "${CASE_TMP}/assets.md"
tee -a "$E2E_NOTES" < "${CASE_TMP}/assets.md"
echo "distinct addresses: $(cut -f1,2 "${CASE_TMP}/assets.tsv" | LC_ALL=C sort -u | grep -c . || true)"

finish
