#!/usr/bin/env bash
# check-links.sh: check that every download URL in scripts/ and docs/ answers.
#
# The reference FASTA URL answered 403 and four groups of doc links answered
# 404 while CI stayed green, because nothing ever requested them. This script
# collects the URLs a user is told to download and asks each one:
#   - scripts/**/*.sh (not scripts/ci/): every URL on a line that is not a
#     comment (download calls, URL variables and the "try manually" hints);
#   - docs/**/*.md: every URL inside a fenced code block (the commands a
#     user pastes).
# ${NAME} and $NAME are filled in from versions.env, so a URL built from a pin
# (the PCGR bundle date, the GRIDSS commit) is checked at that pin. A URL that
# still holds a variable, a printf %s or a glob is listed as skipped.
#
# Each URL gets a HEAD request (redirects followed); when that is not 200, a GET
# of one byte (Range: bytes=0-0), since some servers refuse HEAD. 200 and 206
# pass. Anything else fails, and so does no answer at all: a timeout or a DNS
# error is a failure, never a pass.
#
# Usage:
#   scripts/ci/check-links.sh [ROOT]     check the tree at ROOT (default: repo)
#   scripts/ci/check-links.sh --since REF [ROOT]   only URLs that REF did not
#                                        have (a pull request's new links)
#   scripts/ci/check-links.sh --list [ROOT]   print the URLs, request nothing
#   scripts/ci/check-links.sh --self-test
#
# Prints Markdown on stdout. Exit 0 when every URL answers, 1 when one does
# not or when no URL was found (an empty check is a failure).
set -euo pipefail

REPO_ROOT=$(cd "$(dirname "$0")/../.." && pwd)
JOBS=${LINKS_JOBS:-8}

# Prints "URL<TAB>file:line" for every candidate URL under $1, before variable
# expansion.
raw_urls() {
  local root=$1 f
  {
    if [ -d "${root}/scripts" ]; then
      find "${root}/scripts" -name '*.sh' -type f -not -path '*/scripts/ci/*' | sort | while read -r f; do
        awk -v F="${f#"${root}"/}" '!/^[[:space:]]*#/ && /https?:\/\// { print F ":" FNR "\t" $0 }' "$f"
      done
    fi
    if [ -d "${root}/docs" ]; then
      find "${root}/docs" -name '*.md' -type f | sort | while read -r f; do
        awk -v F="${f#"${root}"/}" '
          /^[[:space:]]*(```|~~~)/ { inb = !inb; next }
          inb && /https?:\/\// { print F ":" FNR "\t" $0 }' "$f"
      done
    fi
  } | while IFS=$'\t' read -r loc line; do
    grep -oE "https?://[^[:space:]\"'\`<>()]+" <<<"$line" | while read -r u; do
      printf '%s\t%s\n' "$u" "$loc"
    done
  done
}

# Fills ${NAME} and $NAME from the versions.env under $1, reading stdin, and
# drops trailing punctuation.
expand_vars() {
  local env=${1}/versions.env
  # shellcheck disable=SC2016  # Python source, not shell
  python3 -c '
import re, sys
vals = {}
try:
    for line in open(sys.argv[1]):
        m = re.match(r"^([A-Z][A-Z0-9_]*)=\"([^\"]*)\"", line)
        if m:
            vals[m.group(1)] = m.group(2)
except FileNotFoundError:
    pass
def sub(m):
    name = m.group(1) or m.group(2)
    return vals.get(name, m.group(0))
for line in sys.stdin:
    url, loc = line.rstrip("\n").split("\t", 1)
    url = re.sub(r"\$\{([A-Z][A-Z0-9_]*)\}|\$([A-Z][A-Z0-9_]*)", sub, url)
    # Trailing punctuation, and the brace that closes ${VAR:-https://...}.
    while url and (url[-1] in ".,;:]" or (url[-1] == "}" and url.count("}") > url.count("{"))):
        url = url[:-1]
    print(url + "\t" + loc)
' "$env"
}

# Prints "URL<TAB>locations" (locations joined by ", "), one line per URL.
collect() {
  raw_urls "$1" | expand_vars "$1" | awk -F'\t' '
    { if (!($1 in loc)) { order[++n] = $1; loc[$1] = $2 } else if (index(loc[$1], $2) == 0) { loc[$1] = loc[$1] ", " $2 } }
    END { for (i = 1; i <= n; i++) print order[i] "\t" loc[order[i]] }'
}

# True when URL $1 cannot be requested as written.
unresolved() {
  case "$1" in
    *'$'*|*'%'*|*'*'*|*'{'*|*'}'*) return 0 ;;
  esac
  local host=${1#*://}
  host=${host%%/*}
  case "$host" in
    *.invalid|*.example|example.com|*.example.com|localhost|localhost:*|127.0.0.1*) return 0 ;;
    *.*) ;;
    *) return 0 ;;  # no dot: a placeholder such as https://url/to/your/data
  esac
  return 1
}

# Prints "CODE<TAB>URL"; CODE is 000 when the server never answered.
probe() {
  local url=$1 code retry=${LINKS_RETRY:-2}
  local common=(-s -o /dev/null -w '%{http_code}' -L --connect-timeout 15 --max-time 60
    --retry "$retry" --retry-delay 5 -A 'Personal-Genome-Pipeline link check (curl)')
  code=$(curl "${common[@]}" -I "$url" 2>/dev/null) || true
  if [ "$code" != 200 ]; then
    code=$(curl "${common[@]}" -r 0-0 "$url" 2>/dev/null) || true
  fi
  printf '%s\t%s\n' "${code:-000}" "$url"
}

# check_tree ROOT [SINCE]: check the URLs under ROOT; with SINCE (a git ref of
# ROOT's repository), only the URLs that tree did not already have.
check_tree() {
  local root=$1 since=${2:-} list tmp n=0 nskip=0 nbad=0 nold=0 url loc code
  list=$(collect "$root")
  tmp=$(mktemp -d)
  : > "${tmp}/check"; : > "${tmp}/skip"; : > "${tmp}/old"
  if [ -z "$list" ]; then
    echo "ERROR: no URL found under ${root}; the extraction is broken or the tree is wrong." >&2
    rm -rf "$tmp"
    return 1
  fi
  if [ -n "$since" ]; then
    mkdir "${tmp}/base"
    if ! git -C "$root" archive "$since" 2>/dev/null | tar -x -C "${tmp}/base"; then
      echo "ERROR: could not read the tree of ${since}" >&2
      rm -rf "$tmp"
      return 1
    fi
    collect "${tmp}/base" | cut -f1 > "${tmp}/old"
  fi
  while IFS=$'\t' read -r url loc; do
    [ -n "$url" ] || continue
    if grep -qxF "$url" "${tmp}/old"; then
      nold=$((nold + 1))
    elif unresolved "$url"; then
      printf '%s\t%s\n' "$url" "$loc" >> "${tmp}/skip"
    else
      printf '%s\n' "$url" >> "${tmp}/check"
    fi
  done <<<"$list"
  n=$(grep -c . "${tmp}/check" || true)
  nskip=$(grep -c . "${tmp}/skip" || true)

  if [ -n "$since" ] && [ "$n" -eq 0 ]; then
    echo "No new or changed download URL since ${since} (${nold} unchanged, not requested)."
    rm -rf "$tmp"
    return 0
  fi
  if [ "$n" -eq 0 ]; then
    echo "ERROR: no URL to check under ${root}; the extraction is broken or the tree is wrong." >&2
    rm -rf "$tmp"
    return 1
  fi

  export -f probe
  export LINKS_RETRY
  # shellcheck disable=SC2016  # $1 expands in the child bash
  xargs -P "$JOBS" -I{} bash -c 'probe "$1"' _ {} < "${tmp}/check" > "${tmp}/codes"

  # A URL with no result line was never probed: count it as failed.
  while read -r url; do
    grep -qF "$(printf '\t%s' "$url")" "${tmp}/codes" || printf '000\t%s\n' "$url" >> "${tmp}/codes"
  done < "${tmp}/check"

  if [ -n "$since" ]; then
    echo "Checked ${n} download URLs that are new since ${since} (${nold} unchanged, not requested)."
  else
    echo "Checked ${n} download URLs from scripts/ and docs/ code blocks."
  fi
  echo ""
  while IFS=$'\t' read -r code url; do
    case "$code" in
      200|206) ;;
      *)
        nbad=$((nbad + 1))
        loc=$(awk -F'\t' -v u="$url" '$1 == u { print $2; exit }' <<<"$list")
        [ "$nbad" -eq 1 ] && { echo "| HTTP | URL | Where |"; echo "|---|---|---|"; }
        [ "$code" = 000 ] && code="no answer"
        echo "| ${code} | ${url} | ${loc} |"
        ;;
    esac
  done < <(sort -t$'\t' -k2 "${tmp}/codes")
  if [ "$nbad" -eq 0 ]; then
    echo "All ${n} answered 200."
  else
    echo ""
    echo "${nbad} of ${n} did not answer 200."
  fi
  if [ "$nskip" -gt 0 ]; then
    echo ""
    echo "<details><summary>Not checked, because they hold a run-time value: ${nskip}</summary>"
    echo ""
    while IFS=$'\t' read -r url loc; do
      echo "- \`${url}\` (${loc})"
    done < "${tmp}/skip"
    echo ""
    echo "</details>"
  fi
  rm -rf "$tmp"
  [ "$nbad" -eq 0 ]
}

self_test() {
  local rc out fail=0
  d=$(mktemp -d)
  trap 'rm -rf "$d"' EXIT
  local good=https://github.com/GeiserX/Personal-Genome-Pipeline
  local bad=https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/no-such-file-link-check-control
  mkdir -p "${d}/scripts/ci" "${d}/docs"
  echo 'TESTREPO="Personal-Genome-Pipeline"' > "${d}/versions.env"
  cat > "${d}/scripts/a.sh" <<EOF
#!/usr/bin/env bash
# ${bad}/in-a-comment
fetch "https://github.com/GeiserX/\${TESTREPO}" out
wget -O x "${bad}"
URL=\${URL:-https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/README.md}
echo "https://ftp.example.invalid/x/\${UNKNOWN_VAR}.tgz"
EOF
  echo "wget ${bad}/in-scripts-ci" > "${d}/scripts/ci/tool.sh"
  cat > "${d}/docs/a.md" <<EOF
Prose link, not checked: ${bad}/in-prose

\`\`\`bash
curl -O https://10.255.255.1/unreachable.tgz
wget ${good}/blob/main/LICENSE
\`\`\`
EOF
  rc=0; out=$(LINKS_RETRY=0 check_tree "$d" 2>&1) || rc=$?
  echo "$out"
  [ "$rc" -ne 0 ] || { echo "SELF-TEST FAIL: a tree with a dead URL passed"; fail=1; }
  grep -qF "| 404 | ${bad} | scripts/a.sh:4 |" <<<"$out" || { echo "SELF-TEST FAIL: the 404 URL was not reported"; fail=1; }
  grep -qF "| no answer | https://10.255.255.1/unreachable.tgz |" <<<"$out" || { echo "SELF-TEST FAIL: a URL with no answer was not a failure"; fail=1; }
  grep -qF 'Checked 5 download URLs' <<<"$out" || { echo "SELF-TEST FAIL: expected 5 URLs (expanded \${TESTREPO}, default value, 404, no answer, docs)"; fail=1; }
  grep -qF '2 of 5 did not answer 200' <<<"$out" || { echo "SELF-TEST FAIL: a URL that answers was reported as failed"; fail=1; }
  if grep -qF 'TESTREPO' <<<"$out"; then
    echo "SELF-TEST FAIL: \${TESTREPO} was not filled in from versions.env"; fail=1
  fi
  if grep -qF -e 'in-a-comment' -e 'in-prose' -e 'in-scripts-ci' <<<"$out"; then
    echo "SELF-TEST FAIL: a comment, prose or scripts/ci URL was checked"; fail=1
  fi
  grep -qF 'UNKNOWN_VAR' <<<"$out" || { echo "SELF-TEST FAIL: an unresolved URL was not listed as skipped"; fail=1; }

  # --since: only the URL a commit adds is requested, and a tree that adds
  # none passes.
  local g=${d}/since
  mkdir -p "${g}/scripts"
  echo "wget ${good}/blob/main/README.md" > "${g}/scripts/b.sh"
  git -C "$g" init -q
  git -C "$g" add scripts
  git -C "$g" -c user.name=link-check -c user.email=link-check@example.invalid commit -q -m base
  rc=0; out=$(LINKS_RETRY=0 check_tree "$g" HEAD 2>&1) || rc=$?
  if [ "$rc" -ne 0 ] || ! grep -qF 'No new or changed download URL since HEAD (1 unchanged' <<<"$out"; then
    echo "$out"; echo "SELF-TEST FAIL: --since with no new URL did not pass"; fail=1
  fi
  echo "wget ${bad}/added-by-the-change" >> "${g}/scripts/b.sh"
  rc=0; out=$(LINKS_RETRY=0 check_tree "$g" HEAD 2>&1) || rc=$?
  echo "$out"
  if [ "$rc" -eq 0 ] || ! grep -qF "| 404 | ${bad}/added-by-the-change |" <<<"$out" \
     || ! grep -qF 'Checked 1 download URLs that are new since HEAD (1 unchanged' <<<"$out"; then
    echo "SELF-TEST FAIL: --since did not request exactly the added URL"; fail=1
  fi

  # An empty tree is a failure.
  mkdir -p "${d}/empty"
  rc=0; out=$(check_tree "${d}/empty" 2>&1) || rc=$?
  [ "$rc" -ne 0 ] || { echo "SELF-TEST FAIL: an empty tree passed"; fail=1; }

  if [ "$fail" -eq 0 ]; then
    echo "SELF-TEST OK: a 404, no answer and an empty tree fail; comments, prose and scripts/ci are left out; \${VAR} is filled from versions.env; --since requests only added URLs."
  fi
  return "$fail"
}

case "${1:-}" in
  --self-test) self_test ;;
  --list) collect "${2:-$REPO_ROOT}" ;;
  --since) check_tree "${3:-$REPO_ROOT}" "${2:?--since needs a git ref}" ;;
  -h|--help) sed -n '2,29p' "$0" ;;
  *) check_tree "${1:-$REPO_ROOT}" ;;
esac
