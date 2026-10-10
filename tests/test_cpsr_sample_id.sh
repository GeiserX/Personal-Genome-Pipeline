#!/usr/bin/env bash
# CPSR (PCGR 2.3.2) stops on a --sample_id shorter than 3 or longer than 40
# characters (pcgr/pcgr_vars.py SAMPLE_ID_MIN_LENGTH and SAMPLE_ID_MAX_LENGTH,
# checked in pcgr/arg_checker.py). A sample called S1 made the step, and so the
# whole run, fail after DeepVariant. Both entry points now give CPSR an id
# inside 3..40 (bin/cpsr_sample_id) and rename its files back to the sample's:
#   - the CPSR process of modules/local/cpsr/main.nf: its script block and its
#     stub, rendered for each id and run with a fake `cpsr` on PATH that
#     records its argv and writes the files the real one names after
#     --sample_id;
#   - scripts/17-cpsr.sh, with a fake `docker` that does the same.
# Ids of 1, 2, 3, 40 and 41 characters. Each must reach CPSR as 3..40
# characters (unchanged when already inside), and leave
# <id>.cpsr.grch38.html and <id>.cpsr.grch38.classification.tsv.gz, and no
# file under the other name.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

L40=$(printf 'a%.0s' $(seq 40))
L41=$(printf 'b%.0s' $(seq 41))
IDS=(A S1 abc "$L40" "$L41")

# want_id ID: the id CPSR must get (the rule of bin/cpsr_sample_id)
want_id() {
  if [ "${#1}" -lt 3 ]; then echo "${1}_cpsr"; elif [ "${#1}" -gt 40 ]; then echo "${1:0:40}"; else echo "$1"; fi
}

# --- bin/cpsr_sample_id -----------------------------------------------------------
for id in "${IDS[@]}"; do
  got=$(sh "${REPO}/bin/cpsr_sample_id" "$id" 2>&1) || got="(exit $?) ${got}"
  if [ "$got" = "$(want_id "$id")" ]; then pass "cpsr_sample_id: ${#id} characters -> ${got}"
  else fail "cpsr_sample_id ${id} printed '${got}', want '$(want_id "$id")'"; fi
done

# --- a fake cpsr: records its argv, writes what CPSR writes ------------------------
mkdir -p "${WORK}/bin"
cat > "${WORK}/bin/cpsr" <<'FAKE'
#!/usr/bin/env bash
printf '%s\n' "$*" >> "$CPSR_ARGV"
out="" sid=""
while [ $# -gt 0 ]; do
  case "$1" in --output_dir) out=$2; shift ;; --sample_id) sid=$2; shift ;; esac
  shift
done
for e in html classification.tsv.gz conf.yaml json.gz; do : > "${out}/${sid}.cpsr.grch38.${e}"; done
FAKE
chmod +x "${WORK}/bin/cpsr"

# check_argv LABEL ID ARGV_FILE: the --sample_id CPSR got
check_argv() {
  local label=$1 id=$2 sid want
  want=$(want_id "$id")
  sid=$(sed -n 's/.*--sample_id \([^ ]*\).*/\1/p' "$3" 2>/dev/null | tail -n 1)
  if [ -z "$sid" ]; then
    fail "${label} (${#id} characters): cpsr was not run with --sample_id"
  elif [ "${#sid}" -lt 3 ] || [ "${#sid}" -gt 40 ]; then
    fail "${label} (${#id} characters): CPSR got --sample_id '${sid}' (${#sid} characters), outside 3..40"
  elif [ "$sid" != "$want" ]; then
    fail "${label} (${#id} characters): CPSR got --sample_id '${sid}', want '${want}'"
  else
    pass "${label} (${#id} characters): CPSR got --sample_id of ${#sid} characters"
  fi
}
# check_files LABEL ID DIR EXT...: DIR holds <ID>.cpsr.grch38.<EXT> for each EXT,
# and no CPSR file named after another id
check_files() {
  local label=$1 id=$2 dir=$3 e f missing="" other
  shift 3
  for e in "$@"; do [ -e "${dir}/${id}.cpsr.grch38.${e}" ] || missing+=" ${id}.cpsr.grch38.${e}"; done
  if [ -n "$missing" ]; then fail "${label} (${#id} characters): missing${missing}; found: $(cd "$dir" && printf '%s ' *)"
  else pass "${label} (${#id} characters): the files carry the sample's id"; fi
  other=$(cd "$dir" && for f in *.cpsr.*; do [ ! -e "$f" ] || [[ "$f" == "${id}.cpsr."* ]] || echo "$f"; done)
  [ -z "$other" ] || fail "${label} (${#id} characters): CPSR files under another name are left: ${other//$'\n'/ }"
}

# --- the Nextflow process: its script and stub blocks, rendered -------------------
MOD="${REPO}/modules/local/cpsr/main.nf"
# render BLOCK ID: the block's Groovy string as Nextflow would run it for meta.id ID
render() {
  python3 - "$MOD" "$1" "$2" <<'PY'
import re, sys, textwrap
src, block, mid = open(sys.argv[1]).read(), sys.argv[2], sys.argv[3]
m = re.search(r"^[ \t]*" + block + r':[ \t]*\n(.*?)^[ \t]*"""[ \t]*\n(.*?)\n[ \t]*"""[ \t]*$', src, re.S | re.M)
if not m:
    sys.exit(f"no {block}: block with a triple-quoted string in {sys.argv[1]}")
if m.group(1).strip():
    sys.exit(f"{block}: code before the string is not rendered here: {m.group(1).strip()!r}")
values = {"meta.id": mid, "vcf": "sample.vcf.gz", "vep_cache_cpsr": "vep_cache", "pcgr_data": "pcgr_data",
          "task.process": "CPSR", "task.container.replaceFirst(/^[^:@]+[:@]/, '')": "2.3.2"}
s, out, i = m.group(2), [], 0
while i < len(s):
    c = s[i]
    if c == "\\" and i + 1 < len(s):
        out.append(s[i + 1]); i += 2
    elif s.startswith("${", i):
        j = s.index("}", i)
        expr = s[i + 2:j]
        if expr not in values:
            sys.exit(f"{block}: unknown Groovy expression ${{{expr}}}")
        out.append(values[expr]); i = j + 1
    elif c == "$" and i + 1 < len(s) and (s[i + 1].isalpha() or s[i + 1] == "_"):
        sys.exit(f"{block}: unescaped ${s[i + 1:i + 20]!r}: Groovy would read it as a variable")
    else:
        out.append(c); i += 1
print(textwrap.dedent("".join(out)))
PY
}
for block in script stub; do
  for id in "${IDS[@]}"; do
    d="${WORK}/task-${block}-${#id}"
    mkdir -p "$d"
    if ! render "$block" "$id" > "${d}/.command.sh" 2> "${d}/render.err"; then
      fail "CPSR ${block} block: $(cat "${d}/render.err")"; continue
    fi
    : > "${d}/argv"
    if [ "$block" = stub ]; then
      # A stub runs no cpsr: only the files it leaves count
      (cd "$d" && PATH="${REPO}/bin:${PATH}" bash -euo pipefail .command.sh > run.log 2>&1) \
        || fail "CPSR stub (${#id} characters) exited non-zero: $(tail -3 "${d}/run.log")"
      check_files "CPSR stub" "$id" "$d" html classification.tsv.gz
    else
      (cd "$d" && CPSR_ARGV="${d}/argv" PATH="${WORK}/bin:${REPO}/bin:${PATH}" bash -euo pipefail .command.sh > run.log 2>&1) \
        || fail "CPSR script (${#id} characters) exited non-zero: $(tail -3 "${d}/run.log")"
      check_argv "CPSR script" "$id" "${d}/argv"
      check_files "CPSR script" "$id" "$d" html classification.tsv.gz conf.yaml
    fi
  done
done

# --- scripts/17-cpsr.sh, with a fake docker ---------------------------------------
cat > "${WORK}/bin/docker" <<'FAKE'
#!/usr/bin/env bash
printf '%s\n' "$*" >> "$CPSR_ARGV"
out="" sid="" prev=""
for a in "$@"; do
  case "$prev" in
    -v) case "$a" in *:/mnt/outputs) out=${a%:/mnt/outputs} ;; esac ;;
    --sample_id) sid=$a ;;
  esac
  prev=$a
done
[ -n "$out" ] && [ -n "$sid" ] || exit 0
for e in html classification.tsv.gz conf.yaml json.gz; do : > "${out}/${sid}.cpsr.grch38.${e}"; done
FAKE
chmod +x "${WORK}/bin/docker"
GD="${WORK}/genome"
# shellcheck source=../versions.env
. "${REPO}/versions.env"
mkdir -p "${GD}/vep_cache/homo_sapiens/${PCGR_VEP_CACHE_RELEASE}_GRCh38" "${GD}/pcgr_data/${PCGR_DATA_BUNDLE}/data"
echo "species homo_sapiens" > "${GD}/vep_cache/homo_sapiens/${PCGR_VEP_CACHE_RELEASE}_GRCh38/info.txt"
for id in "${IDS[@]}"; do
  mkdir -p "${GD}/${id}/vcf"
  : > "${GD}/${id}/vcf/${id}.vcf.gz"
  : > "${WORK}/argv-17-${#id}"
  if ! CPSR_ARGV="${WORK}/argv-17-${#id}" PATH="${WORK}/bin:${PATH}" GENOME_DIR="$GD" \
       bash "${REPO}/scripts/17-cpsr.sh" "$id" > "${WORK}/17-${#id}.log" 2>&1; then
    fail "scripts/17-cpsr.sh (${#id} characters) exited non-zero: $(tail -3 "${WORK}/17-${#id}.log" | tr '\n' '|')"
  fi
  check_argv "scripts/17-cpsr.sh" "$id" "${WORK}/argv-17-${#id}"
  check_files "scripts/17-cpsr.sh" "$id" "${GD}/${id}/cpsr" html classification.tsv.gz conf.yaml
done

if [ "$FAILS" -ne 0 ]; then
  echo "${FAILS} check(s) failed"
  exit 1
fi
echo "all CPSR sample id checks passed"
