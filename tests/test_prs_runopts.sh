#!/usr/bin/env bash
# The outer run's docker flags reach pgsc_calc's own containers.
#
# PRS (modules/local/prs/main.nf) runs pgsc_calc, a Nextflow pipeline of its
# own, with a config it writes, images.config. Its docker.runOptions line
# must start with -u UID:GID and --network none (and the sex folder's mount),
# then carry the outer run's docker.runOptions and the PRS process's
# containerOptions, less their user and network flags.
#
#   bash tests/test_prs_runopts.sh [MODULE]
#       No Docker, no Nextflow (guard.yml runs it on every pull request).
#       Runs the task script's "outer docker flags" block, with Groovy's
#       escapes undone, on several outer values and reads the line back from
#       the images.config it wrote. MODULE defaults to the PRS module.
#       Without the block (the module before it existed) the test fails.
#
# What the Groovy side hands the block (the run's docker.runOptions and
# task.containerOptions) is checked on real runs by the e2e case
# tests/e2e/prs-runopts-and-trace-1.sh.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
MODULE=${1:-${REPO}/modules/local/prs/main.nf}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
FAILS=0

# The block between the markers, as bash sees it: Groovy turns \\ into \ and
# \$ into $; any other escape in the block is refused, so it cannot be misread.
if ! python3 - "$MODULE" > "${WORK}/block.sh" <<'PY'
import sys
text = open(sys.argv[1]).read()
begin, end = "# >>> outer docker flags", "# <<< outer docker flags"
if begin not in text or end not in text:
    sys.exit(f"no '{begin}' ... '{end}' block in {sys.argv[1]}")
block = text.split(begin, 1)[1].split(end, 1)[0].split("\n", 1)[1]
out, i = [], 0
while i < len(block):
    c = block[i]
    if c == "\\":
        nxt = block[i + 1] if i + 1 < len(block) else ""
        if nxt not in ("\\", "$"):
            sys.exit(f"unexpected Groovy escape \\{nxt} in the block")
        out.append(nxt)
        i += 2
    else:
        out.append(c)
        i += 1
print("".join(out))
PY
then
  echo "FAIL: the outer docker flags block of ${MODULE} could not be read"
  exit 1
fi

# line OUTER MOUNT: the docker.runOptions line the block writes for those
# values, after the process block the task writes first (a file a glob would find).
line() {
  local d
  d=$(mktemp -d "${WORK}/case.XXXXXX")
  printf 'process {\n}\n' > "${d}/images.config"
  (cd "$d" && OUTER_DOCKER=$1 MOUNT=$2 bash -euo pipefail "${WORK}/block.sh" > /dev/null) || { echo "the block failed"; return; }
  grep '^docker\.runOptions' "${d}/images.config" || echo "no docker.runOptions line"
}

FIXED="-u $(id -u):$(id -g) --network none"
# expect DESCRIPTION OUTER MOUNT WANT
expect() {
  local got
  got=$(line "$2" "$3")
  if [ "$got" = "docker.runOptions = '$4'" ]; then
    echo "ok:   $1"
  else
    echo "FAIL: $1"
    echo "      outer: $2"
    echo "      got:   ${got}"
    echo "      want:  docker.runOptions = '$4'"
    FAILS=$((FAILS + 1))
  fi
}

expect "no outer flags: the line of before" "" "" "$FIXED"
expect "the docker profile's --user root:root and --network none are dropped" \
  "--user root:root --network none" "" "$FIXED"
expect "a low-priority runOptions reaches pgsc_calc after -u and --network none" \
  "--user root:root --cpu-shares 256" "" "${FIXED} --cpu-shares 256"
expect "every user and network form is dropped, the rest kept in order" \
  "-u 0:0 --cpu-shares 256 --user=1:1 --net host -u0 --network=bridge --memory-swappiness 0 --network none --label a=b" "" \
  "${FIXED} --cpu-shares 256 --memory-swappiness 0 --label a=b"
expect "the sex folder's mount stays before the outer flags" \
  "--cpu-shares 256" ' -v "/w/sex:/w/sex:ro"' "${FIXED}"' -v "/w/sex:/w/sex:ro" --cpu-shares 256'
expect "a glob is passed as it is, not expanded to the files of the task folder" \
  "--cpu-shares 256 *" "" "${FIXED} --cpu-shares 256 *"
# A quote and a backslash, escaped for the single-quoted config string.
expect "a quote and a backslash are escaped for the config" \
  "--label 'a=b c' --env X=a\\b" "" "${FIXED} --label \\'a=b c\\' --env X=a\\\\b"

# Control: a block that drops nothing must fail the second case, so the
# comparison above can fail.
# shellcheck disable=SC2016  # $out and $o belong to the block
sed 's/skip=true ;;/out="$out $o" ;;/' "${WORK}/block.sh" > "${WORK}/broken.sh"
mv "${WORK}/block.sh" "${WORK}/good.sh"
mv "${WORK}/broken.sh" "${WORK}/block.sh"
if [ "$(line "--user root:root --network none" "")" = "docker.runOptions = '${FIXED}'" ]; then
  echo "FAIL: control: a block that keeps --user root:root was not caught"
  FAILS=$((FAILS + 1))
else
  echo "ok:   control: a block that keeps --user root:root is caught"
fi

echo "${FAILS} failure(s)"
[ "$FAILS" -eq 0 ]
