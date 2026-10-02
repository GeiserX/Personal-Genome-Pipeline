#!/usr/bin/env bash
# check-personal-data.sh: fail when a file holds something that looks like
# personal data: a home path, a private address, or a genotype.
#
# This repository is public and the pipeline runs on someone's own genome, so
# an output pasted into a doc, a test or a commit is a leak. Every tracked
# text file is scanned for:
#   home        a home path (/home/x, /Users/x, /mnt/user, C:\Users\x)
#   address     a private-range IPv4 address (10/8, 172.16/12, 192.168/16,
#               100.64/10) or a tailnet name (*.ts.net)
#   diplotype   a star-allele diplotype (*1/*4) on a line with a sample label
#               (sample, patient, subject, proband, NA12878, HG002, ...)
#   rsid        an rsID followed by a genotype (rs4244285 AG, rs123: A/G,
#               rs123 (C;T), or a chip export row rs123<TAB>1<TAB>123<TAB>AG)
#   repeat      a repeat-count genotype (17/19 on a line with REPCN, or after
#               an STR gene such as HTT or FMR1)
# Documented public examples (fixture values, the patterns themselves in the
# contributor guide) are listed in ALLOW below, by file and the text on the
# line, so a new line in the same file is still checked.
#
# Usage:
#   scripts/ci/check-personal-data.sh             every tracked text file
#   scripts/ci/check-personal-data.sh FILE...     only these files
#   scripts/ci/check-personal-data.sh --self-test plant one line of each kind
#                                                 and require each to be caught
#
# Your own identifiers (names, sample ids, hostnames) do not belong in this
# public file. Put one extended regular expression per line in a file outside
# the repository and point PERSONAL_DATA_PATTERNS at it; each match is
# reported as kind "own".
#
# As a pre-commit hook, checking the staged files:
#   cat > .git/hooks/pre-commit <<'HOOK'
#   #!/bin/sh
#   export PERSONAL_DATA_PATTERNS="$HOME/.config/pgp/personal-patterns"   # optional
#   git diff --cached --name-only --diff-filter=ACMR -z |
#     xargs -0 sh -c '[ "$#" -eq 0 ] || exec scripts/ci/check-personal-data.sh "$@"' _
#   HOOK
#   chmod +x .git/hooks/pre-commit
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)

# file (exact path or glob)   text on the line that makes it a documented example
# shellcheck disable=SC2016  # the entries are text, nothing expands
ALLOW='
scripts/ci/check-personal-data.sh             *
.beads/*.jsonl                                *
CLAUDE.md                                     for your own home paths, server hostnames, and names
CONTRIBUTING.md                               /mnt/user\|internal-host\|/home/
CONTRIBUTING.md                               (`/mnt/user/`, `/home/username/`, etc.)
modules/local/cyrius/main.nf                  ${prefix}\\t*1/*1\\tPASS
tests/fake-docker/base-run-all-default.sh     sample1\t*1/*1\tPASS
'

scan() {
  python3 - "$ROOT" "$ALLOW" "${PERSONAL_DATA_PATTERNS:-}" "$@" <<'PY'
import fnmatch, os, re, subprocess, sys
root, allow_text, own_file, args = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
allow = []
for l in allow_text.splitlines():
    if l.strip():
        path, _, text = l.strip().partition(" ")
        allow.append((path, text.strip()))

if args:
    files = args
else:
    out = subprocess.run(["git", "-C", root, "ls-files", "-z"], check=True, capture_output=True).stdout
    files = [os.path.join(root, f) for f in out.decode().split("\0") if f]
if not files:
    print("FAIL: no file to scan")
    sys.exit(1)

LABEL = r'(?:(?i:\b(?:sample|patient|subject|proband|participant|individual|donor)s?[0-9]*\b)|\b(?:NA|HG)[0-9]{3,5}\b)'
GT = r'\(?[ACGTDI][/|;]?[ACGTDI]\)?(?![A-Za-z0-9])'
STR_GENES = (r'HTT|FMR1|C9orf72|C9ORF72|ATXN[0-9]+|DMPK|FXN|CNBP|RFC1|PABPN1|CACNA1A|TCF4|AR|'
             r'PPP2R2B|TBP|JPH3|NOP56|ATN1|CSTB|AFF2|FMR2|NOTCH2NLC|DAB1|BEAN1|STARD7|YEATS2')
OCTET = r'(?:25[0-5]|2[0-4][0-9]|1[0-9][0-9]|[1-9]?[0-9])'
RULES = [
    ("home", re.compile(r'(?<![A-Za-z0-9_])(?:/home/[A-Za-z]|/Users/[A-Za-z]|/mnt/user\b|[A-Za-z]:\\\\?Users\\\\?[A-Za-z])')),
    ("address", re.compile(r'(?<![0-9.])(?:10\.' + OCTET + r'|172\.(?:1[6-9]|2[0-9]|3[01])|192\.168|100\.(?:6[4-9]|[7-9][0-9]|1[01][0-9]|12[0-7]))'
                           r'\.' + OCTET + r'\.' + OCTET + r'(?![0-9]|\.[0-9])')),
    ("address", re.compile(r'[A-Za-z0-9-]\.ts\.net\b')),
    ("diplotype", re.compile(r'^(?=.*' + LABEL + r').*\*[0-9]+[A-Za-z0-9.]*/\*[0-9]+')),
    ("rsid", re.compile(r'\brs[0-9]+\b[\s:=,;|`\'"\[\]-]{1,6}' + GT)),
    ("rsid", re.compile(r'\brs[0-9]+\t[0-9XYMT]{1,2}\t[0-9]+\t[ACGTDI-]{1,2}\b')),
    ("repeat", re.compile(r'^(?=.*\bREPCN\b).*(?<![0-9./])[0-9]{1,4}/[0-9]{1,4}(?![0-9/])')),
    ("repeat", re.compile(r'\b(?:' + STR_GENES + r')\b[^/\n]{0,40}?(?<![0-9.])[0-9]{1,4}/[0-9]{1,4}(?![0-9/])')),
]
if own_file:
    for l in open(own_file):
        if l.strip() and not l.startswith("#"):
            RULES.append(("own", re.compile(l.strip())))

bad = scanned = 0
for path in files:
    full = os.path.abspath(path)
    rel = os.path.relpath(full, root) if full.startswith(root + os.sep) else path
    try:
        data = open(path, "rb").read()
    except (IsADirectoryError, FileNotFoundError):
        continue
    if b"\0" in data[:8192]:
        continue  # binary
    scanned += 1
    texts = [t for p, t in allow if fnmatch.fnmatchcase(rel, p)]
    if "*" in texts:
        continue
    for n, line in enumerate(data.decode("utf-8", "replace").split("\n"), 1):
        for kind, rx in RULES:
            m = rx.search(line)
            if m and not any(t in line for t in texts):
                print("FAIL %-9s %s:%d: %s" % (kind, rel, n, line.strip()[:160]))
                bad = 1
print("scanned %d text files" % scanned)
sys.exit(bad)
PY
}

self_test() {
  local tmp fail=0 rc out k
  tmp=$(mktemp -d)
  # shellcheck disable=SC2064
  trap "rm -rf '$tmp'" EXIT
  # One planted line per kind, each in its own file. All values are made up.
  printf '%s\n' 'cp results.html /home/alice/genome/' > "${tmp}/home.md"
  printf '%s\n' 'scp out.vcf.gz root@192.168.7.20:/data/' > "${tmp}/address.md"
  printf '%s\n' 'ssh nas.tail1234.ts.net' > "${tmp}/tailnet.md"
  printf '%s\n' '| sample7 | CYP2D6 | *1/*4 |' > "${tmp}/diplotype.md"
  printf '%s\n' 'rs4244285 AG' > "${tmp}/rsid.md"
  printf 'rs4477212\t1\t82154\tAA\n' > "${tmp}/chip.txt"
  printf '%s\n' 'HTT CAG repeats: 17/19' > "${tmp}/repeat.md"
  printf '%s\n' 'chr4 3074877 . C <STR17> . PASS . GT:REPCN 0/1:17/19' > "${tmp}/repcn.vcf"
  for k in home address tailnet diplotype rsid chip repeat repcn; do
    f=$(ls "${tmp}/${k}".*)
    rc=0; out=$(scan "$f" 2>&1) || rc=$?
    if [ "$rc" -ne 1 ] || ! grep -q '^FAIL' <<<"$out"; then
      echo "self-test: planted ${k} line was not caught (exit ${rc}):"; printf '%s\n' "$out"; fail=1
    else
      echo "self-test: caught: $(grep '^FAIL' <<<"$out")"
    fi
  done
  # Lines that must pass: doc text without a sample, a version, a public URL.
  # shellcheck disable=SC2016  # the backticks are text
  printf '%s\n' '- `*1/*4` -- Intermediate metabolizer' 'samtools 1.21, GATK 4.6.2.0' \
    'https://example.org/home/page' 'CYP2C19*2 is rs4244285 (G>A)' > "${tmp}/clean.md"
  rc=0; out=$(scan "${tmp}/clean.md" 2>&1) || rc=$?
  [ "$rc" -eq 0 ] || { echo "self-test: clean lines were reported:"; printf '%s\n' "$out"; fail=1; }
  # A pattern of your own, from PERSONAL_DATA_PATTERNS.
  printf '%s\n' 'Report for Jane Example' > "${tmp}/own.md"
  printf '%s\n' '\bJane Example\b' > "${tmp}/patterns"
  rc=0; out=$(PERSONAL_DATA_PATTERNS="${tmp}/patterns" scan "${tmp}/own.md" 2>&1) || rc=$?
  if [ "$rc" -ne 1 ] || ! grep -q '^FAIL own' <<<"$out"; then
    echo "self-test: a PERSONAL_DATA_PATTERNS match was not caught:"; printf '%s\n' "$out"; fail=1
  else
    echo "self-test: caught: $(grep '^FAIL' <<<"$out")"
  fi
  [ "$fail" -eq 0 ] && echo "self-test: OK"
  return "$fail"
}

case "${1:-}" in
  --self-test) self_test ;;
  -*) echo "usage: $0 [FILE... | --self-test]" >&2; exit 2 ;;
  *) scan "$@" ;;
esac
