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
# contributor guide) are listed in ALLOW below, by file and text; public or
# invented test data files are listed whole in SKIP. An entry
# covers only its own text: the rest of that line, and every other line of
# the file, is still checked. This script holds the list, so every listed
# text is allowed here as well.
#
# Usage:
#   scripts/ci/check-personal-data.sh             every tracked text file
#   scripts/ci/check-personal-data.sh FILE...     only these files
#   scripts/ci/check-personal-data.sh --cached FILE...
#                                                 the staged content of these
#                                                 files (git index), not the
#                                                 working tree
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
#     xargs -0 sh -c '[ "$#" -eq 0 ] || exec scripts/ci/check-personal-data.sh --cached "$@"' _
#   HOOK
#   chmod +x .git/hooks/pre-commit
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)

# file (exact path or glob)   text on the line that makes it a documented example
# shellcheck disable=SC2016  # the entries are text, nothing expands
ALLOW='
CLAUDE.md                                     /mnt/user\|/Users/\|/home/
scripts/ci/check-personal-data.sh             (/home/x, /Users/x, /mnt/user, C:\Users\x)
scripts/ci/check-personal-data.sh             (*1/*4) on a line with a sample label
scripts/ci/check-personal-data.sh             (rs4244285 AG, rs123: A/G,
scripts/ci/check-personal-data.sh             rs123 (C;T), or a chip export row
scripts/ci/check-personal-data.sh             (17/19 on a line with REPCN
CONTRIBUTING.md                               /mnt/user\|internal-host\|/home/
CONTRIBUTING.md                               (`/mnt/user/`, `/home/username/`, etc.)
modules/local/cyrius/main.nf                  ${prefix}\\t*1/*1\\tPASS
tests/fake-docker/reports-stale-step.sh       sample1\t*1/*1\tPASS
tests/test_collect_summary.py                 S\t*1/*4\tPASS
tests/smoke/commands.tsv                      HTT at 19/45 repeats
tests/smoke/commands.tsv                      ATXN3 at 22/24 repeats
'

# Whole files that are not scanned (glob), with the reason. Keep this list
# short: an entry hides every line of the file.
SKIP='
tests/fixtures/*                public or invented test data; each fixture directory says where it came from
tests/smoke/*.vcf.in            invented tool inputs for the image smoke tests
'

scan() {
  python3 - "$ROOT" "$ALLOW" "$SKIP" "${PERSONAL_DATA_PATTERNS:-}" "${CACHED:-}" "$@" <<'PY'
import fnmatch, os, re, subprocess, sys
root, allow_text, skip_text, own_file, cached, args = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6:]
skip = [l.split()[0] for l in skip_text.splitlines() if l.strip()]
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

# A label may carry a number after a separator: sample7, sample_01, donor-3.
LABEL = r'(?:(?i:(?<![A-Za-z])(?:sample|patient|subject|proband|participant|individual|donor)s?(?:[_-]?[0-9]+)?(?![A-Za-z]))|\b(?:NA|HG)[0-9]{3,5}\b)'
GT = r'\(?[ACGTDI][/|;]?[ACGTDI]\)?(?![A-Za-z0-9])'
STR_GENES = (r'HTT|FMR1|C9orf72|C9ORF72|ATXN[0-9]+|DMPK|FXN|CNBP|RFC1|PABPN1|CACNA1A|TCF4|AR|'
             r'PPP2R2B|TBP|JPH3|NOP56|ATN1|CSTB|AFF2|FMR2|NOTCH2NLC|DAB1|BEAN1|STARD7|YEATS2')
OCTET = r'(?:25[0-5]|2[0-4][0-9]|1[0-9][0-9]|[1-9]?[0-9])'
# (kind, pattern, the line must also match this or None). The pattern
# marks only the leaked value, so an allowed example covers just its span.
RULES = [
    ("home", re.compile(r'(?<![A-Za-z0-9_])(?:/home/[A-Za-z]|/Users/[A-Za-z]|/mnt/use[r]\b|[A-Za-z]:\\\\?Users\\\\?[A-Za-z])'), None),
    ("address", re.compile(r'(?<![0-9.])(?:10\.' + OCTET + r'|172\.(?:1[6-9]|2[0-9]|3[01])|192\.168|100\.(?:6[4-9]|[7-9][0-9]|1[01][0-9]|12[0-7]))'
                           r'\.' + OCTET + r'\.' + OCTET + r'(?![0-9]|\.[0-9])'), None),
    ("address", re.compile(r'[A-Za-z0-9-]\.ts\.net\b'), None),
    ("diplotype", re.compile(r'\*[0-9]+[A-Za-z0-9.]*/\*[0-9]+'), re.compile(LABEL)),
    ("rsid", re.compile(r'\brs[0-9]+\b[\s:=,;|`\'"\[\]-]{1,6}' + GT), None),
    ("rsid", re.compile(r'\brs[0-9]+\t[0-9XYMT]{1,2}\t[0-9]+\t[ACGTDI-]{1,2}\b'), None),
    ("repeat", re.compile(r'(?<![0-9./])[0-9]{1,4}/[0-9]{1,4}(?![0-9/])'), re.compile(r'\bREPCN\b')),
    ("repeat", re.compile(r'\b(?:' + STR_GENES + r')\b[^/\n]{0,40}?(?<![0-9.])[0-9]{1,4}/[0-9]{1,4}(?![0-9/])'), None),
]
if own_file:
    for l in open(own_file):
        if l.strip() and not l.startswith("#"):
            RULES.append(("own", re.compile(l.strip()), None))

bad = scanned = 0
for path in files:
    full = os.path.abspath(path)
    rel = os.path.relpath(full, root) if full.startswith(root + os.sep) else path
    if any(fnmatch.fnmatchcase(rel, g) for g in skip):
        continue
    if cached:
        # The blob in the index: what the commit will hold.
        r = subprocess.run(["git", "-C", root, "show", ":" + rel], capture_output=True)
        if r.returncode != 0:
            continue
        data = r.stdout
    else:
        try:
            data = open(path, "rb").read()
        except (IsADirectoryError, FileNotFoundError):
            continue
    if b"\0" in data[:8192]:
        continue  # binary
    scanned += 1
    if rel == "scripts/ci/check-personal-data.sh":
        texts = [t for _, t in allow]
    else:
        texts = [t for p, t in allow if fnmatch.fnmatchcase(rel, p)]
    for n, line in enumerate(data.decode("utf-8", "replace").split("\n"), 1):
        # Spans of the line that are documented examples.
        spans = [(i.start(), i.end()) for t in texts for i in re.finditer(re.escape(t), line)]
        for kind, rx, cond in RULES:
            if cond is not None and not cond.search(line):
                continue
            for m in rx.finditer(line):
                if any(a <= m.start() and m.end() <= b for a, b in spans):
                    continue  # inside an allowed example; the rest of the line is still checked
                print("FAIL %-9s %s:%d: %s" % (kind, rel, n, line.strip()[:160]))
                bad = 1
                break
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
  # plant drops every "~", which keeps this file itself free of the patterns.
  plant() { printf '%s\n' "${2//\~/}" > "${tmp}/$1"; }
  plant home.md 'cp results.html /home/~alice/genome/'
  plant address.md 'scp out.vcf.gz root@192.~168.7.20:/data/'
  plant tailnet.md 'ssh nas.tail1234.ts.~net'
  plant diplotype.md '| sample7 | CYP2D6 | *1/~*4 |'
  plant label.md 'sample_01 CYP2D6 *1/~*4'
  plant rsid.md 'rs4244285~ AG'
  printf 'rs4477212\t1\t82154\tAA\n' > "${tmp}/chip.txt"
  plant repeat.md 'HTT CAG repeats: 17/~19'
  plant repcn.vcf 'chr4 3074877 . C <STR17> . PASS . GT:REPCN 0/~1:17/~19'
  for k in home address tailnet diplotype label rsid chip repeat repcn; do
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
  # An allowed example covers its own text only: alone it passes, with a home
  # path appended to the same line it is reported.
  mkdir -p "${tmp}/tree"
  # shellcheck disable=SC2016  # the backticks are text
  plant tree/CONTRIBUTING.md 'Personal file paths (`/mnt/user/`, `/home/username/`, etc.)'
  rc=0; out=$(ROOT="${tmp}/tree" scan "${tmp}/tree/CONTRIBUTING.md" 2>&1) || rc=$?
  [ "$rc" -eq 0 ] || { echo "self-test: an allowed example was reported:"; printf '%s\n' "$out"; fail=1; }
  # shellcheck disable=SC2016  # the backticks are text
  plant tree/CONTRIBUTING.md 'Personal file paths (`/mnt/user/`, `/home/username/`, etc.) see /home/~bob/x'
  rc=0; out=$(ROOT="${tmp}/tree" scan "${tmp}/tree/CONTRIBUTING.md" 2>&1) || rc=$?
  if [ "$rc" -ne 1 ] || ! grep -q '^FAIL home' <<<"$out"; then
    echo "self-test: a home path after an allowed example was not caught:"; printf '%s\n' "$out"; fail=1
  else
    echo "self-test: caught: $(grep '^FAIL' <<<"$out")"
  fi
  # --cached reads the index: a staged leak is caught after the working copy
  # drops it.
  git -C "${tmp}/tree" init -q
  plant tree/staged.md 'cp results.html /home/~alice/genome/'
  git -C "${tmp}/tree" add staged.md
  plant tree/staged.md 'cp results.html genome/'
  rc=0; out=$(cd "${tmp}/tree" && ROOT="${tmp}/tree" CACHED=1 scan staged.md 2>&1) || rc=$?
  if [ "$rc" -ne 1 ] || ! grep -q '^FAIL home' <<<"$out"; then
    echo "self-test: --cached did not read the staged content:"; printf '%s\n' "$out"; fail=1
  else
    echo "self-test: caught: $(grep '^FAIL' <<<"$out")"
  fi
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
  --cached) shift; CACHED=1 scan "$@" ;;
  -*) echo "usage: $0 [FILE... | --self-test]" >&2; exit 2 ;;
  *) scan "$@" ;;
esac
