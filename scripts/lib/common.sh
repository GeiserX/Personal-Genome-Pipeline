# shellcheck shell=bash
# common.sh: sourced by every script in scripts/, once SAMPLE and GENOME_DIR
# are read from the arguments and the environment:
#   # shellcheck source=lib/common.sh
#   . "$(dirname "$0")/lib/common.sh"
#   validate_sample "$SAMPLE"
#
# It gives each script the same:
#   versions.env   every image tag and coupled data version, sourced with no
#                  fallback: a missing file or a typo stops the script.
#   umask 077      files the host side writes are readable by their owner only.
#   GENOME_DIR     as exported by the caller (not required here: setup.sh and
#                  validate-setup.sh run without it).
#   REF_FASTA      host path of the GRCh38 FASTA (override from the environment;
#                  it must lie inside GENOME_DIR, the only directory containers
#                  see). REF_FASTA_C is the same file inside a container,
#                  REF_DICT its GATK sequence dictionary.
#   THREADS        CPU budget of one step. Default 8; a script that wants another
#                  default sets THREADS=${THREADS:-N} before sourcing this file.
#   CONTAINER_ENGINE  docker by default.
# and the helpers validate_sample, require_image, cpath, run_in, fetch,
# lock_acquire / lock_release and pipeline_images, described where they are
# defined. Keep it bash 3.2 compatible: macOS runs setup.sh with /bin/bash.

PGP_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
# shellcheck source=../versions.env
. "${PGP_ROOT}/versions.env"

umask 077
CONTAINER_ENGINE=${CONTAINER_ENGINE:-docker}
THREADS=${THREADS:-8}

# validate_sample NAME: exit 2 unless NAME is a plain sample name. Scripts put
# it into container paths and `bash -c` bodies, so only letters, digits, '.',
# '_' and '-' are allowed (the same rule as main.nf), and '.' and '..' are not
# names at all.
validate_sample() {
  local s=${1:-}
  if [[ ! "$s" =~ ^[A-Za-z0-9._-]+$ ]] || [ "$s" = "." ] || [ "$s" = ".." ]; then
    echo "ERROR: invalid sample name '${s}'. Use only letters, digits, '.', '_' and '-' (and not '.' or '..')." >&2
    exit 2
  fi
}

# require_image VAR...: exit 1 unless each VAR is set by versions.env.
require_image() {
  local v
  for v in "$@"; do
    if [ -z "${!v:-}" ]; then
      echo "ERROR: ${v} is not set; add it to versions.env." >&2
      exit 1
    fi
  done
}

# cpath HOST_PATH: the path a container sees for a host path inside GENOME_DIR.
cpath() {
  case "$1" in
    "${GENOME_DIR:?}") printf '/genome' ;;
    "${GENOME_DIR}"/*) printf '/genome%s' "${1#"$GENOME_DIR"}" ;;
    *) echo "ERROR: ${1} is outside GENOME_DIR (${GENOME_DIR}); containers only see GENOME_DIR." >&2
       return 1 ;;
  esac
}

# shellcheck disable=SC2034  # REF_DICT and REF_FASTA_C are read by the scripts
if [ -n "${GENOME_DIR:-}" ]; then
  GENOME_DIR=${GENOME_DIR%/}
  REF_FASTA=${REF_FASTA:-${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta}
  REF_DICT="${REF_FASTA%.*}.dict"
  REF_FASTA_C=$(cpath "$REF_FASTA") || exit 2
fi

# run_in [--net] [--root] [--rw DIR]... [docker run options] IMAGE [COMMAND...]
#
# `docker run --rm` with these defaults:
#   --network none            the analysis steps need no network; a tool that
#                             reaches out at run time fails loudly instead.
#   GENOME_DIR at /genome:ro  references, ClinVar, caches and bundles are
#                             read-only, and only
#   GENOME_DIR/SAMPLE         (when SAMPLE is set) is writable, at /genome/SAMPLE.
#   --user UID:GID, HOME=/tmp outputs belong to the caller, not to root.
# Opt-outs, each named at the call site with the reason:
#   --net      the step downloads something (a database, a pip package).
#   --root     the image cannot run as an unprivileged user.
#   --rw DIR   DIR (inside GENOME_DIR) is writable too, e.g. a shared index.
# Everything after the opt-outs goes to `docker run` unchanged.
run_in() {
  local isolate=true root=false d
  local -a extra=()
  while [ $# -gt 0 ]; do
    case "$1" in
      --net) isolate=false; shift ;;
      --root) root=true; shift ;;
      --rw)
        d=${2:?run_in --rw needs a directory}
        mkdir -p "$d"
        extra+=(-v "${d}:$(cpath "$d")")
        shift 2 ;;
      *) break ;;
    esac
  done
  local -a args=(run --rm)
  if $isolate; then args+=(--network none); fi
  if ! $root; then args+=(--user "$(id -u):$(id -g)" -e HOME=/tmp); fi
  args+=(-v "${GENOME_DIR:?}:/genome:ro")
  if [ -n "${SAMPLE:-}" ]; then
    mkdir -p "${GENOME_DIR}/${SAMPLE}"
    args+=(-v "${GENOME_DIR}/${SAMPLE}:/genome/${SAMPLE}")
  fi
  "$CONTAINER_ENGINE" "${args[@]}" ${extra[@]+"${extra[@]}"} "$@"
}

# _digest md5|sha256|sum FILE: the checksum of FILE (sum = BSD sum, the
# format of Ensembl's CHECKSUMS files; only its first field is printed).
_digest() {
  case "$1" in
    md5)
      if command -v md5sum >/dev/null 2>&1; then md5sum "$2" | awk '{print $1}'; else md5 -q "$2"; fi ;;
    sha256)
      if command -v sha256sum >/dev/null 2>&1; then sha256sum "$2" | awk '{print $1}'; else shasum -a 256 "$2" | awk '{print $1}'; fi ;;
    sum) sum "$2" | awk '{print $1 + 0}' ;;
    *) echo "ERROR: unknown checksum type '${1}'" >&2; return 2 ;;
  esac
}

# _get URL OUT: download URL to OUT (resuming a partial OUT), curl first and
# wget when curl is missing. OUT "-" writes to stdout.
_get() {
  if command -v curl >/dev/null 2>&1; then
    if [ "$2" = "-" ]; then curl -fsSL "$1"; else curl -fL -C - -o "$2" "$1"; fi
  elif command -v wget >/dev/null 2>&1; then
    if [ "$2" = "-" ]; then wget -q -O - "$1"; else wget -c -O "$2" "$1"; fi
  else
    echo "ERROR: neither curl nor wget is installed." >&2
    return 1
  fi
}

# fetch URL DEST [md5|sha256|sum VALUE|CHECKSUM_URL]
#
# Downloads URL to DEST.part (resuming one left by an interrupted run), checks
# it, then renames it to DEST. Nothing ever writes DEST directly, so a DEST
# that exists is whole. The check is the checksum when one is given (VALUE, or
# the file at CHECKSUM_URL: its only line, or the line naming URL's file), else
# `gzip -t` for a .gz file, else non-empty. A failed check deletes DEST.part
# and returns 1, so the next run downloads it again.
fetch() {
  local url=$1 dest=$2 kind=${3:-} want=${4:-} part="${2}.part" name got line i ok=false
  name=$(basename "$url")
  if [ -n "$kind" ]; then
    case "$want" in
      http://*|https://*|ftp://*)
        if ! line=$(_get "$want" - | awk -v n="$name" '
              NF { lines++; first = $1; f = $NF; sub(/.*\//, "", f); if (f == n) hit = $1 }
              END { if (lines == 1) print first; else if (hit != "") print hit }'); then
          line=""
        fi
        if [ -z "$line" ]; then
          echo "ERROR: could not read the checksum of ${name} from ${want}" >&2
          return 1
        fi
        want=$line ;;
    esac
  fi
  mkdir -p "$(dirname "$dest")"
  for i in 1 2 3; do
    if _get "$url" "$part"; then ok=true; break; fi
    echo "  Download attempt ${i}/3 failed: ${url}" >&2
    [ "$i" -lt 3 ] && sleep $((i * 5))
  done
  if ! $ok; then
    echo "ERROR: could not download ${url} (the partial file ${part} is kept and resumed next time)" >&2
    return 1
  fi
  if [ -n "$kind" ]; then
    got=$(_digest "$kind" "$part") || return 1
    if [ "$kind" = sum ]; then want=$((10#${want})); fi
    if [ "$got" != "$want" ]; then
      rm -f "$part"
      echo "ERROR: ${kind} checksum of ${name} is ${got}, expected ${want}. Removed the download; run again to fetch it anew." >&2
      return 1
    fi
  elif [[ "$dest" == *.gz ]]; then
    if ! gzip -t "$part" 2>/dev/null; then
      rm -f "$part"
      echo "ERROR: ${name} is not a complete gzip file. Removed the download; run again to fetch it anew." >&2
      return 1
    fi
  elif [ ! -s "$part" ]; then
    rm -f "$part"
    echo "ERROR: ${name} downloaded empty." >&2
    return 1
  fi
  mv -f "$part" "$dest"
}

# lock_acquire DIR / lock_release DIR: a lock shared by concurrent runs. The
# lock is a directory holding the owner's PID; a lock whose PID is gone (a
# killed run) is stale and taken over.
lock_acquire() {
  local lock=$1 pid
  until mkdir "$lock" 2>/dev/null; do
    pid=$(cat "${lock}/pid" 2>/dev/null || true)
    if [ -n "$pid" ] && ! kill -0 "$pid" 2>/dev/null; then
      echo "  Removing stale lock ${lock} (process ${pid} is gone)."
      rm -rf "$lock"
      continue
    fi
    echo "  Waiting for ${lock} (held by process ${pid:-starting})..."
    sleep 5
  done
  echo "$$" > "${lock}/pid"
}
lock_release() {
  rm -rf "$1"
}

# pipeline_images: the value of every *_IMAGE in versions.env, one per line,
# except those whose line is marked `# optional` (pulled on first use).
pipeline_images() {
  local line name
  while IFS= read -r line || [ -n "$line" ]; do
    case "$line" in
      *"# optional"*) continue ;;
    esac
    if [[ "$line" =~ ^([A-Z0-9_]+_IMAGE)= ]]; then
      name=${BASH_REMATCH[1]}
      printf '%s\n' "${!name}"
    fi
  done < "${PGP_ROOT}/versions.env"
}
