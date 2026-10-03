#!/usr/bin/env bash
# check-images.sh: fail when an image is named anywhere but versions.env.
#
# versions.env is the one list of container images. conf/containers.config
# and docs/versions.md are generated from it. Every other file refers to an
# image through its *_IMAGE variable, so a bump is one line in versions.env
# and nothing is left on the old tag. This script checks three things:
#
#   1. No literal image (registry/name:tag, name@sha256:..., or a Docker Hub
#      official image such as python:3.11) in scripts/, modules/, workflows/,
#      conf/, docs/, tests/, .github/workflows/, main.nf or nextflow.config.
#   2. conf/containers.config is what versions.env produces, and every
#      process has a selector (scripts/ci/gen-containers-config.sh --check).
#   3. docs/versions.md is what versions.env produces
#      (scripts/ci/gen-versions-doc.sh --check).
#
# Usage:
#   scripts/ci/check-images.sh               check this repository
#   scripts/ci/check-images.sh --root DIR    check another tree
#   scripts/ci/check-images.sh --self-test   plant each fault in a copy and
#                                            require this script to report it
set -euo pipefail
export LC_ALL=C

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
MODE=check
while [ $# -gt 0 ]; do
  case "$1" in
    --self-test) MODE=self-test ;;
    --root) ROOT=$(cd "$2" && pwd); shift ;;
    *) echo "usage: $0 [--root DIR | --self-test]" >&2; exit 2 ;;
  esac
  shift
done

# Paths whose files are scanned.
SCAN='scripts modules workflows conf docs tests .github/workflows main.nf nextflow.config'

# Files that may name images, with the reason.
EXEMPT_FILES='
versions.env                          the one list of images
conf/containers.config                generated from versions.env
docs/versions.md                      generated from versions.env
docs/lessons-learned.md               history: names the tags that failed
docs/research/                        dated research notes, kept as written
.github/workflows/container-test.yml  its own matrix, held to versions.env by its sync step
scripts/ci/check-images.sh            this file: the exemptions and the faults its self-test plants
scripts/ci/freshness.py               reads image references; its self-test plants made-up ones
scripts/ci/changed-images.sh          reads image references; its self-test plants made-up ones
tests/fixtures/                       provenance notes name the image that wrote each fixture
'

# Single literals that stay on purpose: file, image, reason. Each line allows
# one occurrence: a change to the image or one more copy of it in the same
# file is reported again.
EXEMPT_PAIRS='
docs/quick-test.md                quay.io/biocontainers/samtools:1.20--h50ea8bc_0   needs CA certificates, which SAMTOOLS_IMAGE lacks
docs/08-hla-typing.md             jiachenzdocker/hla-la@sha256:ecca23de6635aa85e60b4ee39dd4e15341b5febb514e5478f2b2a086f05a447c   HLA-LA is not a pipeline step
docs/stylesheets/theme.css        quay.io/biocontainers/expansionhunter:5.0.0   example in a CSS comment
scripts/cyrius-constraints.txt    python:3.11   the image pip resolved these constraints in
scripts/ci/settle-doubts.sh       jmcdani20/hap.py:v0.3.12   repeats versions.env, to be removed
scripts/ci/settle-doubts.sh       quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_5   repeats versions.env, to be removed
scripts/ci/settle-doubts.sh       hkubal/clair3:v2.0.2   repeats versions.env, to be removed
docs/troubleshooting.md           quay.io/biocontainers/toolname:tag   a placeholder, not an image
'

# scan: print one line per literal image in the tree at $ROOT; exit 1 if any.
scan() {
  local list
  if [ -e "${ROOT}/.git" ]; then
    # shellcheck disable=SC2086  # SCAN splits into paths on purpose
    list=$(git -C "$ROOT" ls-files -- $SCAN)
  else
    # shellcheck disable=SC2086
    list=$(cd "$ROOT" && find $SCAN -type f 2>/dev/null || true)
  fi
  [ -n "$list" ] || { echo "FAIL: found no file to scan under ${ROOT}"; return 1; }
  # The file list goes in as an argument: stdin carries the program.
  python3 - "$ROOT" "$EXEMPT_FILES" "$EXEMPT_PAIRS" "$list" <<'PY'
import collections, os, re, sys
root, exempt_files, exempt_pairs = sys.argv[1], sys.argv[2], sys.argv[3]
files = [f for f in sys.argv[4].split("\n") if f]
if not files:
    print("FAIL: no file to scan")
    sys.exit(1)
skip = [l.split()[0] for l in exempt_files.splitlines() if l.strip()]
# Each line allows one occurrence; list a pair twice to allow two.
pairs = collections.Counter(tuple(l.split()[:2]) for l in exempt_pairs.splitlines() if l.strip())

# Docker Hub official images have no slash: those versions.env uses, plus
# the usual base images.
env = open(os.path.join(root, "versions.env")).read()
bare = set(re.findall(r'^[A-Z0-9_]+_IMAGE="([^"/:@]+)[:@]', env, re.M))
bare |= {"ubuntu", "debian", "alpine", "busybox", "centos", "rockylinux", "fedora",
         "python", "node", "perl", "r-base", "openjdk", "eclipse-temurin", "golang", "rust"}

# Any tag: a version, or a word such as latest or stable.
TAG = r'(?::[A-Za-z0-9_][A-Za-z0-9._-]*|@sha256:[0-9a-f]{64})'
LEAD = r'(?<![A-Za-z0-9_./:@$-])'
SHAPE = re.compile(LEAD + r'([a-z][a-z0-9._-]*/[a-z][a-z0-9._/-]*' + TAG + r')(?![A-Za-z0-9_/-])')
BARE = re.compile(LEAD + r'((?:' + "|".join(map(re.escape, sorted(bare))) + r')' + TAG + r')(?![A-Za-z0-9_/-])')
URL = re.compile(r'[a-z][a-z0-9+.-]*://\S+')
# A path to a file, such as modules/local/vep/main.nf:23, is not an image.
FILE = re.compile(r'\.(nf|sh|md|py|config|ya?ml|json|txt|tsv|csv|vcf|gz|html|css|js|R|r):[0-9]+(-[0-9]+)?$')

bad = 0
for f in files:
    if any(f == s or (s.endswith("/") and f.startswith(s)) for s in skip):
        continue
    try:
        text = open(os.path.join(root, f), encoding="utf-8").read()
    except (UnicodeDecodeError, IsADirectoryError, FileNotFoundError):
        continue
    left = collections.Counter({k: v for k, v in pairs.items() if k[0] == f})
    for n, line in enumerate(text.split("\n"), 1):
        # docker://image:tag (a step image in a workflow) is an image, not a URL.
        line = URL.sub(" ", line.replace("docker://", " "))
        for rx in (SHAPE, BARE):
            for m in rx.finditer(line):
                image = m.group(1).rstrip(".")  # a full stop ends the sentence, not the tag
                if FILE.search(image) or image.startswith("example/"):
                    continue  # example/ is the placeholder namespace of the self-tests
                if left[(f, image)] > 0:
                    left[(f, image)] -= 1
                    continue
                print("FAIL: %s:%d names the image %s; use its *_IMAGE variable from versions.env" % (f, n, image))
                bad = 1
print("scanned %d files" % len(files))
sys.exit(bad)
PY
}

check() {
  local fail=0
  if scan; then echo "OK: no file names an image outside versions.env and the files generated from it."; else fail=1; fi
  "${ROOT}/scripts/ci/gen-containers-config.sh" --check || fail=1
  "${ROOT}/scripts/ci/gen-versions-doc.sh" --check || fail=1
  return "$fail"
}

# self_test: plant each fault in a copy of this tree and require check() on
# the copy to name it; the unchanged copy must pass.
self_test() {
  local tmp fail=0 rc out
  tmp=$(mktemp -d)
  # shellcheck disable=SC2064
  trap "rm -rf '$tmp'" EXIT
  copy() {
    mkdir -p "${tmp}/$1"
    (cd "$ROOT" && git ls-files -z | xargs -0 tar -cf -) | tar -xf - -C "${tmp}/$1"
    cp "${BASH_SOURCE[0]}" "${tmp}/$1/scripts/ci/check-images.sh"   # this script, even before it is committed
  }
  # expect NAME PATTERN: the check of copy NAME exits 1 and prints PATTERN.
  expect() {
    rc=0; out=$("${tmp}/$1/scripts/ci/check-images.sh" 2>&1) || rc=$?
    if [ "$rc" -ne 1 ] || ! grep -qE -- "$2" <<<"$out"; then
      echo "self-test: '$1' exited ${rc} and did not report /$2/:"; printf '%s\n' "$out"; fail=1
    else
      echo "self-test: '$1' caught: $(grep -E -- "$2" <<<"$out" | head -1)"
    fi
  }

  copy clean
  rc=0; out=$("${tmp}/clean/scripts/ci/check-images.sh" 2>&1) || rc=$?
  [ "$rc" -eq 0 ] || { echo "self-test: the unchanged copy failed:"; printf '%s\n' "$out"; fail=1; }

  # A step script that names its image.
  copy script-literal
  echo 'docker run --rm staphb/bcftools:1.21 bcftools --version' >> "${tmp}/script-literal/scripts/06-clinvar-screen.sh"
  expect script-literal '^FAIL: scripts/06-clinvar-screen\.sh:[0-9]+ names the image staphb/bcftools:1\.21;'

  # A test that names its image at the end of a sentence: the full stop is
  # not part of the tag.
  copy test-literal
  echo '# Needs docker: runs staphb/samtools:1.19.' >> "${tmp}/test-literal/tests/test_indexcov_sex.sh"
  expect test-literal '^FAIL: tests/test_indexcov_sex\.sh:[0-9]+ names the image staphb/samtools:1\.19;'

  # A Docker Hub official image in a module, and a floating tag in a workflow.
  copy official
  echo '// container python:3.12' >> "${tmp}/official/modules/local/roh/main.nf"
  expect official '^FAIL: modules/local/roh/main\.nf:[0-9]+ names the image python:3\.12;'
  copy floating
  echo '#   image: staphb/samtools:latest' >> "${tmp}/floating/.github/workflows/lint.yml"
  expect floating '^FAIL: \.github/workflows/lint\.yml:[0-9]+ names the image staphb/samtools:latest;'
  # A word tag, and a step image written as docker://.
  echo '#   image: staphb/samtools:stable' >> "${tmp}/floating/.github/workflows/lint.yml"
  echo '#   - uses: docker://alpine:3.8' >> "${tmp}/floating/.github/workflows/lint.yml"
  expect floating '^FAIL: \.github/workflows/lint\.yml:[0-9]+ names the image staphb/samtools:stable;'
  expect floating '^FAIL: \.github/workflows/lint\.yml:[0-9]+ names the image alpine:3\.8;'

  # A doc that repeats a tag, and an exempt literal that changed.
  copy doc-literal
  echo 'Run docker pull quay.io/biocontainers/mosdepth:0.3.13--h05c3d44_0 first.' >> "${tmp}/doc-literal/docs/quick-test.md"
  expect doc-literal '^FAIL: docs/quick-test\.md:[0-9]+ names the image quay\.io/biocontainers/mosdepth:0\.3\.13--h05c3d44_0;'
  copy exempt-changed
  sed -i.bak 's|samtools:1.20--h50ea8bc_0|samtools:1.21--h50ea8bc_0|' "${tmp}/exempt-changed/docs/quick-test.md"
  expect exempt-changed '^FAIL: docs/quick-test\.md:[0-9]+ names the image quay\.io/biocontainers/samtools:1\.21--h50ea8bc_0;'
  copy exempt-twice
  echo 'Then docker pull quay.io/biocontainers/samtools:1.20--h50ea8bc_0 again.' >> "${tmp}/exempt-twice/docs/quick-test.md"
  expect exempt-twice '^FAIL: docs/quick-test\.md:[0-9]+ names the image quay\.io/biocontainers/samtools:1\.20--h50ea8bc_0;'

  # A module process with no selector in conf/containers.config.
  copy no-selector
  printf '%s\n' 'process NEW_TOOL {' '    script:' '    """' '    true' '    """' '}' \
    >> "${tmp}/no-selector/modules/local/roh/main.nf"
  expect no-selector '^FAIL: process NEW_TOOL has no container'

  # A tag bumped in versions.env without running the generators.
  copy stale
  sed -i.bak -E 's|^(MOSDEPTH_IMAGE="[^"]*:)[^"]*"|\19.99"|' "${tmp}/stale/versions.env"
  expect stale 'conf/containers.config does not match versions.env'
  expect stale 'docs/versions.md is out of date'

  [ "$fail" -eq 0 ] && echo "self-test: OK"
  return "$fail"
}

case "$MODE" in
  check) check ;;
  self-test) self_test ;;
esac
