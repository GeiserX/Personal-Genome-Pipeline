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
#   REF_FASTA      host path of the GRCh38 FASTA. Default: the NCBI GRCh38
#                  no-ALT analysis set setup.sh installs. Override from the
#                  environment, absolute or relative to GENOME_DIR; it must lie
#                  inside GENOME_DIR, the only directory containers see.
#                  REF_FASTA_C is the same file inside a container, REF_DICT
#                  its GATK sequence dictionary.
#   THREADS        CPU budget of one step. Default 8; a script that wants another
#                  default sets THREADS=${THREADS:-N} before sourcing this file.
#   CONTAINER_ENGINE  docker by default.
# and the helpers validate_sample, require_image, cpath, run_in, have_output,
# wrote_vcf, atomic_out, fetch, data_file / install_data_file, install_vep_cache,
# lock_acquire / lock_release, pipeline_images, scatter_beds and run_parallel,
# described where they are
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
# names at all. Nor is a name that starts with '-': CPSR, bin/pgx_parse.py and
# bin/collect_summary.py take it as an argument value, and argparse reads '-x'
# as an option. CPSR's 3 to 40 characters are bin/cpsr_sample_id's job.
validate_sample() {
  local s=${1:-}
  if [[ ! "$s" =~ ^[A-Za-z0-9._-]+$ ]] || [ "$s" = "." ] || [ "$s" = ".." ]; then
    echo "ERROR: invalid sample name '${s}'. Use only letters, digits, '.', '_' and '-' (and not '.' or '..')." >&2
    exit 2
  fi
  if [[ "$s" == -* ]]; then
    echo "ERROR: invalid sample name '${s}': it starts with '-', which tools read as an option. Start it with a letter or digit." >&2
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
  # A new file name on purpose: a reference with ALT contigs left from an
  # older version is never picked up by accident (docs/realignment.md).
  REF_FASTA=${REF_FASTA:-${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta}
  case "$REF_FASTA" in
    /*) ;;
    *) REF_FASTA="${GENOME_DIR}/${REF_FASTA}" ;;
  esac
  # GATK names the dictionary after the FASTA without its extensions:
  # ref.fasta and ref.fa.gz both give ref.dict.
  REF_DICT="${REF_FASTA%.gz}"
  REF_DICT="${REF_DICT%.*}.dict"
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
#   --root     the image cannot run as an unprivileged user: the container runs
#              as root (--user 0:0), also when the image names another user.
#   --rw DIR   DIR (inside GENOME_DIR) is writable too, e.g. a shared index.
# Everything after the opt-outs goes to `docker run` unchanged.
# When an unprivileged run fails and the sample directory holds a path the
# caller does not own (versions before this library ran every container as
# root), run_in says how to take the directory back.
run_in() {
  local isolate=true root=false d c rc=0 other
  local -a extra=()
  while [ $# -gt 0 ]; do
    case "$1" in
      --net) isolate=false; shift ;;
      --root) root=true; shift ;;
      --rw)
        d=${2:?run_in --rw needs a directory}
        c=$(cpath "$d") || exit 2
        mkdir -p "$d"
        extra+=(-v "${d}:${c}")
        shift 2 ;;
      *) break ;;
    esac
  done
  local -a args=(run --rm)
  if $isolate; then args+=(--network none); fi
  # --root is explicit: the container runs as root (--user 0:0) even when
  # the image names another user. It is for an image that has not been shown
  # to run unprivileged, for example PCGR in scripts/17-cpsr.sh.
  if $root; then
    args+=(--user 0:0)
  else
    args+=(--user "$(id -u):$(id -g)" -e HOME=/tmp)
  fi
  args+=(-v "${GENOME_DIR:?}:/genome:ro")
  if [ -n "${SAMPLE:-}" ]; then
    mkdir -p "${GENOME_DIR}/${SAMPLE}"
    args+=(-v "${GENOME_DIR}/${SAMPLE}:/genome/${SAMPLE}")
  fi
  "$CONTAINER_ENGINE" "${args[@]}" ${extra[@]+"${extra[@]}"} "$@" || rc=$?
  if [ "$rc" -ne 0 ] && ! $root && [ -n "${SAMPLE:-}" ]; then
    # A root-owned *.lock (left by a Nextflow CRAM_ARCHIVE task under
    # --user root) is opened read-only and blocks nothing: not worth the hint.
    other=$(find "${GENOME_DIR}/${SAMPLE}" ! -user "$(id -u)" ! -name '*.lock' -print 2>/dev/null | head -n 1)
    if [ -n "$other" ]; then
      echo "NOTE: ${other} is not yours. Older versions of this pipeline ran every container as root;" >&2
      echo "  steps now run as you and cannot write over such files. If the error above is" >&2
      echo "  'Permission denied', take the sample directory back and run the step again:" >&2
      echo "    sudo chown -R \"$(id -u):$(id -g)\" \"${GENOME_DIR}/${SAMPLE}\"" >&2
    fi
  fi
  return "$rc"
}

# have_output FILE...: true when every FILE is a finished output, the test a
# step uses before it skips work. Every FILE must be non-empty, and by name:
#   .vcf                  starts with ##fileformat=VCF and has a #CHROM line;
#   .vcf.gz, .bcf, .bam   ends with the BGZF end-of-file block, which htslib
#                         and htsjdk write only when they close the file, and
#                         starts with the format's own header.
# So an empty file left by a failed redirect, or one cut short by a killed
# run, is redone instead of trusted.
have_output() {
  local f first
  for f in "$@"; do
    [ -s "$f" ] || return 1
    case "$f" in
      *.vcf)
        first=$(head -n 1 "$f")
        [[ "$first" == "##fileformat=VCF"* ]] || return 1
        grep -q '^#CHROM' "$f" || return 1 ;;
      *.vcf.gz|*.bcf|*.bam)
        _bgzf_complete "$f" || return 1
        first=$(gzip -cd "$f" 2>/dev/null | head -c 16 | tr -d '\000') || true
        case "$f" in
          *.vcf.gz) [[ "$first" == "##fileformat=VCF"* ]] || return 1 ;;
          *.bcf) [[ "$first" == BCF* ]] || return 1 ;;
          *.bam) [[ "$first" == BAM* ]] || return 1 ;;
        esac ;;
    esac
  done
}

# wrote_vcf FILE: true when FILE is non-empty and starts with a VCF header,
# plain or gzip-compressed. The check after a tool that exited 0: it catches
# a run that wrote nothing where the tool said it succeeded.
wrote_vcf() {
  local first
  [ -s "$1" ] || return 1
  first=$(gzip -cdf "$1" 2>/dev/null | head -c 16) || true
  [[ "$first" == "##fileformat=VCF"* ]]
}

# vcf_header_samples: the sample names of the VCF header on stdin, one per
# line: the #CHROM columns after FORMAT. Nothing for a sites-only VCF. The
# pipeline analyses one sample per run (VCF_PRECHECK applies the same rule).
vcf_header_samples() {
  awk -F'\t' '/^#CHROM/ { for (i = 10; i <= NF; i++) print $i }'
}

# _bgzf_complete FILE: true when FILE ends with the 28-byte BGZF EOF block.
_bgzf_complete() {
  local last
  last=$(tail -c 28 "$1" | od -An -v -tx1 | tr -d ' \n')
  [ "$last" = "1f8b08040000000000ff0600424302001b0003000000000000000000" ]
}

# atomic_out FILE COMMAND [ARGS...]: run COMMAND with its standard output in
# FILE.tmp, and rename that to FILE only when COMMAND exits 0 and wrote
# something. Otherwise FILE.tmp is removed, FILE is left as it was and the
# exit code is returned, so a failed image pull or a crash never leaves an
# empty FILE behind for the next run to accept.
atomic_out() {
  local out=$1 rc=0
  shift
  rm -f "${out}.tmp"
  "$@" > "${out}.tmp" || rc=$?
  if [ "$rc" -eq 0 ] && [ ! -s "${out}.tmp" ]; then
    echo "ERROR: $1 wrote nothing for ${out}." >&2
    rc=1
  fi
  if [ "$rc" -ne 0 ]; then
    rm -f "${out}.tmp"
    return "$rc"
  fi
  mv -f "${out}.tmp" "$out"
}

# _digest md5|sha256|sum FILE: the checksum of FILE (sum = BSD sum, the
# format of Ensembl's CHECKSUMS files; only its first field is printed).
_digest() {
  case "$1" in
    md5)
      if command -v md5sum >/dev/null 2>&1; then md5sum "$2" | awk '{print $1}'
      elif command -v md5 >/dev/null 2>&1; then md5 -q "$2"
      else openssl md5 -r "$2" | awk '{print $1}'; fi ;;
    sha256)
      if command -v sha256sum >/dev/null 2>&1; then sha256sum "$2" | awk '{print $1}'; else shasum -a 256 "$2" | awk '{print $1}'; fi ;;
    sum) sum "$2" | awk '{print $1 + 0}' ;;
    *) echo "ERROR: unknown checksum type '${1}'" >&2; return 2 ;;
  esac
}

# _get URL OUT: download URL to OUT (resuming a partial OUT), curl first and
# wget when curl is missing. OUT "-" writes to stdout. A server that cannot
# resume (curl exit 33) gets the partial OUT dropped and one fresh download.
_get() {
  local rc=0
  if command -v curl >/dev/null 2>&1; then
    if [ "$2" = "-" ]; then
      curl -fsSL "$1"
    else
      curl -fL -C - -o "$2" "$1" || rc=$?
      if [ "$rc" -eq 33 ]; then
        rm -f "$2"
        curl -fL -o "$2" "$1"
      else
        return "$rc"
      fi
    fi
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
#
# A failed download, of the file or of its CHECKSUM_URL, is tried FETCH_TRIES
# times (default 3), FETCH_WAIT seconds apart (default 5, then 10, ...). A
# larger FETCH_TRIES and FETCH_WAIT wait out a server that answers 404 for
# minutes at a time and then comes back.
fetch() {
  local url=$1 dest=$2 kind=${3:-} want=${4:-} part="${2}.part" name got line i ok=false
  local tries=${FETCH_TRIES:-3} wait=${FETCH_WAIT:-}
  case "$tries" in ''|*[!0-9]*) tries=0 ;; *) tries=$((10#$tries)) ;; esac
  [ "$tries" -gt 0 ] || { echo "ERROR: FETCH_TRIES must be a whole number above 0, got '${FETCH_TRIES}'" >&2; return 1; }
  case "$wait" in *[!0-9]*) echo "ERROR: FETCH_WAIT must be a whole number of seconds, got '${wait}'" >&2; return 1 ;; esac
  name=$(basename "$url")
  if [ -n "$kind" ]; then
    case "$want" in
      http://*|https://*|ftp://*|file://*)
        # A read that worked but names no checksum is not retried.
        for ((i = 1; i <= tries; i++)); do
          if line=$(_get "$want" - | awk -v n="$name" '
                NF { lines++; first = $1; f = $NF; sub(/.*\//, "", f); if (f == n) hit = $1 }
                END { if (lines == 1) print first; else if (hit != "") print hit }'); then
            break
          fi
          line=""
          echo "  Checksum download attempt ${i}/${tries} failed: ${want}" >&2
          [ "$i" -lt "$tries" ] && sleep "${wait:-$((i * 5))}"
        done
        if [ -z "$line" ]; then
          echo "ERROR: could not read the checksum of ${name} from ${want}" >&2
          return 1
        fi
        want=$line ;;
    esac
  fi
  mkdir -p "$(dirname "$dest")"
  for ((i = 1; i <= tries; i++)); do
    if _get "$url" "$part"; then ok=true; break; fi
    echo "  Download attempt ${i}/${tries} failed: ${url}" >&2
    [ "$i" -lt "$tries" ] && sleep "${wait:-$((i * 5))}"
  done
  if ! $ok; then
    echo "ERROR: could not download ${url} (the partial file ${part} is kept and resumed next time)" >&2
    return 1
  fi
  if [ -n "$kind" ]; then
    got=$(_digest "$kind" "$part") || return 1
    if [ "$kind" = sum ]; then
      case "$want" in
        ''|*[!0-9]*) echo "ERROR: '${want}' is not a sum checksum for ${name}" >&2; return 1 ;;
      esac
      want=$((10#${want}))
    fi
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

# --- Pinned data files ----------------------------------------------------------
# Small reference files from a fixed commit or release, each checked before it
# is stored. setup.sh installs them (install_data_file); a step only reads
# them (data_file) and runs no download of its own. HLA_DB_RELEASE can be
# overridden from the environment; step 08 keys its T1K index on it, so a new
# release builds a new index.
HLA_DB_RELEASE=${HLA_DB_RELEASE:-3.65.0}
GENCODE_RELEASE=50
DELLY_EXCL_COMMIT=72a05f3cd23bf55f1edd723b9504a0faf1fba230
# shellcheck disable=SC2034  # read by setup.sh and validate-setup.sh
DATA_FILES="delly_exclude cytoband hla_dat gencode_genes"

# data_file NAME: print the host path of data file NAME and return 0 when it
# is installed, 1 when it is not (the path is printed either way).
#   delly_exclude  Delly's GRCh38 exclude map: telomeres, centromeres and the
#                  extra contigs (step 19)
#   cytoband       UCSC's GRCh38 chromosome bands, chr1-22, X and Y (step 10)
#   hla_dat        hla.dat of IPD-IMGT/HLA release HLA_DB_RELEASE (step 08)
#   gencode_genes  the gene lines of GENCODE's basic annotation, release
#                  GENCODE_RELEASE: T1K takes HLA gene coordinates from it (step 08)
data_file() {
  local out
  case "$1" in
    delly_exclude) out="${GENOME_DIR:?}/reference/delly_human.hg38.excl.tsv" ;;
    cytoband) out="${GENOME_DIR:?}/reference/cytoBand.hg38.txt" ;;
    hla_dat) out="${GENOME_DIR:?}/hla/IPD-IMGT-HLA_${HLA_DB_RELEASE}/hla.dat" ;;
    gencode_genes) out="${GENOME_DIR:?}/reference/gencode.v${GENCODE_RELEASE}.basic.genes.gtf" ;;
    *) echo "ERROR: unknown data file '${1}'" >&2; return 2 ;;
  esac
  printf '%s\n' "$out"
  [ -s "$out" ]
}

# install_data_file NAME: download, check and store data file NAME unless it
# is installed. Progress goes to stderr. Returns 1 when it cannot.
install_data_file() {
  local out url gz md5
  out=$(data_file "$1") && return 0
  [ -n "$out" ] || return 2
  case "$1" in
    delly_exclude)
      fetch "https://raw.githubusercontent.com/dellytools/delly/${DELLY_EXCL_COMMIT}/excludeTemplates/human.hg38.excl.tsv" \
        "$out" sha256 caef8593f82f513694ae41b429680eab8abe39c057a119dae7760715aef006ee >&2 || return 1 ;;
    cytoband)
      gz="${out}.download.gz"
      fetch https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz "$gz" \
        sha256 e514e48b0a12a516fd5231b44241a20a282f290bdd6131cadcddc6276ee26082 >&2 || return 1
      # TelomereHunter reads chr1-22, X and Y and reports every other line
      # (the ALT, random and chrM rows) as invalid.
      _keep_lines "$gz" "$out" '$1 ~ /^chr([0-9]+|X|Y)$/' || return 1 ;;
    hla_dat)
      url="https://raw.githubusercontent.com/ANHIG/IMGTHLA/v${HLA_DB_RELEASE}-alpha"
      # md5checksum.txt lines read "MD5 (hla.dat.zip) = <md5>".
      md5=$(_get "${url}/md5checksum.txt" - 2>/dev/null | awk '$2 == "(hla.dat.zip)" {print $NF}') || md5=""
      if [ -z "$md5" ]; then
        echo "ERROR: no md5 for hla.dat.zip in ${url}/md5checksum.txt (is ${HLA_DB_RELEASE} an IPD-IMGT/HLA release?)" >&2
        return 1
      fi
      fetch "${url}/hla.dat.zip" "${out}.zip" md5 "$md5" >&2 || return 1
      # The T1K image has Perl's unzip module; the host may have no unzip.
      # shellcheck disable=SC2154  # T1K_IMAGE comes from versions.env
      run_in --rw "$(dirname "$out")" "$T1K_IMAGE" perl -e \
        'use IO::Uncompress::Unzip qw(unzip $UnzipError); unzip($ARGV[0] => $ARGV[1]) or die "unzip failed: $UnzipError\n";' \
        "$(cpath "${out}.zip")" "$(cpath "${out}.tmp")" >&2 || { rm -f "${out}.tmp"; return 1; }
      if ! grep -q 'IPD-IMGT/HLA' "${out}.tmp" 2>/dev/null; then
        rm -f "${out}.tmp"
        echo "ERROR: ${out}.zip held no IPD-IMGT/HLA hla.dat" >&2
        return 1
      fi
      mv -f "${out}.tmp" "$out"
      rm -f "${out}.zip" ;;
    gencode_genes)
      url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}/gencode.v${GENCODE_RELEASE}.basic.annotation.gtf.gz"
      gz="${out}.download.gz"
      fetch "$url" "$gz" md5 "$(dirname "$url")/MD5SUMS" >&2 || return 1
      _keep_lines "$gz" "$out" '$3 == "gene"' || return 1 ;;
  esac
  [ -s "$out" ]
}

# _keep_lines GZ OUT AWK_CONDITION: write the tab-separated lines of GZ that
# match AWK_CONDITION to OUT (through OUT.tmp) and remove GZ. Fails, keeping
# no OUT, when no line matches.
_keep_lines() {
  if gzip -cd "$1" | awk -F'\t' "$3" > "${2}.tmp" && [ -s "${2}.tmp" ]; then
    mv -f "${2}.tmp" "$2"
    rm -f "$1"
  else
    rm -f "${2}.tmp" "$1"
    echo "ERROR: no usable lines in $(basename "$1")" >&2
    return 1
  fi
}

# vep_cache_url RELEASE: where Ensembl publishes the indexed GRCh38 VEP cache.
vep_cache_url() {
  printf 'https://ftp.ensembl.org/pub/release-%s/variation/indexed_vep_cache/homo_sapiens_vep_%s_GRCh38.tar.gz' "$1" "$1"
}

# install_vep_cache DIR RELEASE: download the cache (checked against Ensembl's
# CHECKSUMS file), unpack it in a temporary directory and move it to
# DIR/homo_sapiens/RELEASE_GRCh38 only when complete. The tarball is deleted.
install_vep_cache() {
  local dir=$1 rel=$2 url tarball tmp final
  url=$(vep_cache_url "$rel")
  tarball="${dir}/$(basename "$url")"
  final="${dir}/homo_sapiens/${rel}_GRCh38"
  if [ -e "$final" ]; then
    echo "ERROR: ${final} exists but has no info.txt (an interrupted extraction?)." >&2
    echo "  Remove it and run again." >&2
    return 1
  fi
  mkdir -p "${dir}/homo_sapiens"
  fetch "$url" "$tarball" sum "$(dirname "$url")/CHECKSUMS" || return 1
  echo "  Extracting $(basename "$tarball")..."
  tmp=$(mktemp -d "${dir}/.extract.XXXXXX")
  if ! tar xzf "$tarball" -C "$tmp" || [ ! -f "${tmp}/homo_sapiens/${rel}_GRCh38/info.txt" ]; then
    rm -rf "$tmp"
    echo "ERROR: could not unpack ${tarball} into a ${rel}_GRCh38 cache." >&2
    return 1
  fi
  mv "${tmp}/homo_sapiens/${rel}_GRCh38" "$final"
  rm -rf "$tmp" "$tarball"
}

# lock_acquire DIR / lock_release DIR: a lock shared by concurrent runs. The
# lock is a directory holding the owner's PID. It is stale when that PID is
# gone (a killed run), or when it has no PID a minute after it was made (a run
# killed between mkdir and writing the PID). A stale lock is renamed before it
# is removed, so of two waiters that judged it stale only one takes it over.
lock_acquire() {
  local lock=$1 pid stale
  until mkdir "$lock" 2>/dev/null; do
    pid=$(cat "${lock}/pid" 2>/dev/null || true)
    stale=false
    if [ -n "$pid" ]; then
      kill -0 "$pid" 2>/dev/null || stale=true
    elif [ -n "$(find "$lock" -maxdepth 0 -mmin +1 2>/dev/null)" ]; then
      stale=true
    fi
    if $stale; then
      if mv "$lock" "${lock}.stale.$$" 2>/dev/null; then
        if [ -n "$pid" ]; then
          echo "  Removing stale lock ${lock} (process ${pid} is gone)."
        else
          echo "  Removing stale lock ${lock} (no owner after a minute)."
        fi
        rm -rf "${lock}.stale.$$"
      fi
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

# scatter_beds DIR [INTERVALS]: split the calling of a single-process caller
# into units, one BED file each (DIR/001.bed, ...), and print their paths in
# reference order. Without INTERVALS: chr1-22, chrX, chrY and chrM of
# REF_FASTA's .fai one unit each, and every other contig in one last unit.
# With INTERVALS (space-separated contigs or contig:start-end regions, 1-based
# inclusive): one unit per region. SCATTER=false puts every region in one unit,
# so the caller runs once over all of them.
scatter_beds() {
  local dir=$1 intervals=${2:-} fai="${REF_FASTA}.fai" r
  rm -rf "$dir"
  mkdir -p "$dir"
  [ -s "$fai" ] || { echo "ERROR: ${fai} not found" >&2; return 1; }
  {
    if [ -n "$intervals" ]; then
      for r in $intervals; do printf 'R\t%s\n' "$r"; done
    fi
  } | awk -F'\t' -v d="$dir" -v all="${SCATTER:-true}" -v given="${intervals:+1}" '
    FNR == NR { len[$1] = $2; order[++n] = $1; rank[$1] = n; next }
    function unit(name) { if (all == "false") return sprintf("%s/001.bed", d); return sprintf("%s/%03d.bed", d, ++u) }
    { r = $2; gsub(/,/, "", r); c = r; s = 0; e = ""
      if (match(r, /:[0-9]+-[0-9]+$/)) { c = substr(r, 1, RSTART - 1); split(substr(r, RSTART + 1), p, "-"); s = p[1] - 1; e = p[2] }
      if (!(c in len)) { print "ERROR: contig " c " (INTERVALS) is not in the reference" > "/dev/stderr"; bad = 1; exit 1 }
      if (e == "") e = len[c]
      k++; rk[k] = rank[c]; st[k] = s + 0; reg[k] = c "\t" s "\t" e }
    END {
      if (bad) exit 1
      if (given) {
        # Reference order, then start, whatever the order of INTERVALS: the
        # units are joined in this order.
        for (i = 2; i <= k; i++)
          for (j = i; j > 1 && (rk[j - 1] > rk[j] || (rk[j - 1] == rk[j] && st[j - 1] > st[j])); j--) {
            t = rk[j]; rk[j] = rk[j - 1]; rk[j - 1] = t
            t = st[j]; st[j] = st[j - 1]; st[j - 1] = t
            t = reg[j]; reg[j] = reg[j - 1]; reg[j - 1] = t
          }
        for (i = 1; i <= k; i++) { f = unit(); print reg[i] >> f; close(f) }
        exit 0
      }
      for (i = 1; i <= n; i++) if (order[i] ~ /^chr([0-9]+|X|Y|M)$/) { f = unit(); print order[i] "\t0\t" len[order[i]] >> f; close(f) }
      rest = (all == "false") ? sprintf("%s/001.bed", d) : sprintf("%s/%03d.bed", d, u + 1)
      for (i = 1; i <= n; i++) if (order[i] !~ /^chr([0-9]+|X|Y|M)$/) print order[i] "\t0\t" len[order[i]] >> rest
    }' "$fai" - || return 1
  ls "$dir"/*.bed
}

# run_parallel MAX FUNC ARG...: run `FUNC ARG` for every ARG in the
# background, at most MAX at a time, and return 1 when any of them failed,
# after all have ended. Needs bash 4.3 (wait -n); the step scripts need 4.4.
run_parallel() {
  local max=$1 fn=$2 a running=0 rc=0
  shift 2
  for a in "$@"; do
    "$fn" "$a" &
    running=$((running + 1))
    if [ "$running" -ge "$max" ]; then
      wait -n || rc=1
      running=$((running - 1))
    fi
  done
  while [ "$running" -gt 0 ]; do
    wait -n || rc=1
    running=$((running - 1))
  done
  return "$rc"
}
