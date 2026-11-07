#!/usr/bin/env bash
# e2e-run.sh — run every end-to-end case in tests/e2e/ on the HG002 fixture.
#
# Usage: scripts/ci/e2e-run.sh [pattern]
#        scripts/ci/e2e-run.sh --self-test-pull | --image-groups | --pull-images
#        scripts/ci/e2e-run.sh --pull-image REF
#   pattern  optional shell glob; runs only the cases whose file name matches
#            (e.g. '2*' or '*cyrius*'). Later cases read earlier cases' outputs,
#            so a partial run is for debugging only.
#   --self-test-pull  checks the pull retry of step 3 against a fake docker
#            (seconds, no network) and the image cache groups, and exits.
#   --image-groups  prints "GROUP IMAGE" for every image the cases use (the
#            e2e workflow's image cache entries) and exits.
#   --pull-images  pulls those images, as step 3 below does, and exits.
#   --pull-image REF  pulls one image with the same retry (the docker shim's
#            pulls of the lazy images at their first use), and exits.
#   --data-key NAME  prints the Actions cache key of data entry NAME (hla, pgs,
#            catalogs).
#   --data-check NAME  exits 0 when entry NAME in E2E_DATA_CACHE is whole and
#            matches its checksums (the e2e workflow saves only such an entry).
#   --data-fill NAME  fills entry NAME in E2E_DATA_CACHE as a case's run does,
#            within E2E_DATA_BUDGET, and exits (124 when the budget ran out).
#
# Env: E2E_WORK  work area (default ${RUNNER_TEMP:-/tmp}/e2e-work), about 25 GB
#      GH_TOKEN  token for `gh release download` (the job's GITHUB_TOKEN)
#      E2E_CASES  a file of case names (without .sh): run only those, as a
#                full run (pre-pull, all-cases-ran marker). The e2e workflow
#                gives each of its parallel jobs its list from e2e-shards.sh.
#      E2E_JOB   the job name written beside each result (default local)
#      E2E_MIN_FREE_GB  fail before the first case when less is free under
#                E2E_WORK (default 0)
#      E2E_BASE_OUT  after the cases, pack every file they wrote or changed
#                (under E2E_WORK and the checkout) into this .tar.zst, and
#                write "E2E_WORK CHECKOUT" to E2E_BASE_OUT.from
#      E2E_BASE_IN   after the layout, unpack this .tar.zst from E2E_BASE_OUT,
#                which must have been packed with the same E2E_WORK and
#                checkout, then delete it to give the cases its disk space
#      E2E_TRACE=true  after each case, record how long it took and which
#                earlier case wrote each file it read (atime and ctime, with
#                the file systems remounted strictatime; needs sudo), in
#                $E2E_WORK/report/trace.tsv for `e2e-shards.sh times`
#      E2E_DATA_CACHE  a folder of downloads kept between runs (see "Downloads
#                kept between runs" below); E2E_DATA_KEEP=true keeps it whole
#      E2E_DATA_MIRROR  a folder URL holding copies of the pinned downloads,
#                tried before their sources
#      E2E_DATA_BUDGET  the most seconds one data cache fill may take
#                (default 2400); a fill still running then fails its case
#      E2E_MIRROR_BUDGET  the most seconds one file may take on the mirror
#                before its source is tried (default 900, or half of
#                E2E_DATA_BUDGET when that is less)
#
# Every run writes $E2E_WORK/report/results.tsv ("CASE RESULT SECONDS JOB",
# tab-separated) and report/extra.md (failed checks and notes), which the
# e2e workflow merges into one table across its jobs.
#
# What it does:
#   1. downloads the release named in tests/fixtures/VERSION and checks
#      SHA256SUMS; when a pull request bumps VERSION, waits up to 75 minutes
#      for the build-fixture job to publish it;
#   2. lays out GENOME_DIR the way setup.sh and step 13 would leave it;
#   3. pulls every image the cases use that is not on the machine yet,
#      trying a pull again only when the registry answered with a 5xx (the
#      five largest wait for the first case that runs them), and checks the
#      free disk against E2E_MIN_FREE_GB;
#   4. runs tests/e2e/*.sh in name order (C locale: numbered cases first, then
#      the <package-key>-*.sh cases later packages add), or those of E2E_CASES,
#      and keeps going after a failure, so one run lists every broken step;
#   5. prints a case / result / time / log table to the job summary;
#   6. with E2E_BASE_OUT, packs what the cases wrote; exits 1 if any case failed.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
PATTERN=${1:-*}

export REPO
export SAMPLE=HG002
export THREADS=4
export E2E_WORK="${E2E_WORK:-${RUNNER_TEMP:-/tmp}/e2e-work}"
export FIXTURE_DIR="${E2E_WORK}/fixture"
export GENOME_DIR="${E2E_WORK}/genome"
export E2E_NOTES="${E2E_WORK}/notes.md"
LOG_DIR="${E2E_WORK}/logs"
CASE_TIMEOUT=${CASE_TIMEOUT:-3600}
GH_REPO=${GITHUB_REPOSITORY:-GeiserX/OwnGenome}
TAG=$(tr -d '[:space:]' < "${REPO}/tests/fixtures/VERSION")
SUMMARY=${GITHUB_STEP_SUMMARY:-/dev/null}

# The docker shim clamps --cpus to this machine's CPU count.
export PATH="${REPO}/tests/e2e/bin:${PATH}"

# --- Image pulls ----------------------------------------------------------------
# A `docker run` whose image is missing pulls it with one manifest request, and
# a registry 5xx on that request ends the step: in run 37980248545 Docker Hub
# answered 500 for hap.py and the GIAB case failed; the rerun passed. So the
# images are pulled before the cases, and a pull that fails with a registry
# 5xx is tried again, up to PULL_TRIES tries in all. Any other error (a wrong
# tag, a rate limit) fails at once. The tools' own `docker run` is never retried.
PULL_TRIES=${PULL_TRIES:-4}
PULL_WAIT=${PULL_WAIT:-20}   # seconds before the second try; doubles after each
REGISTRY_5XX='(HTTP status|status code):? 5[0-9][0-9]|50[0-4] (Internal Server Error|Bad Gateway|Service Unavailable|Gateway Time-?out)'
# pull_image IMAGE: docker pull with the retry above.
pull_image() {
  local image=$1 try=1 wait=$PULL_WAIT out
  while :; do
    if out=$(docker pull -q "$image" 2>&1); then echo "pulled ${image}"; return 0; fi
    echo "$out" >&2
    if ! grep -qE "$REGISTRY_5XX" <<< "$out"; then
      echo "ERROR: pulling ${image} failed, and not with a registry 5xx: not tried again." >&2
      return 1
    fi
    if [ "$try" -ge "$PULL_TRIES" ]; then
      echo "ERROR: pulling ${image}: a registry 5xx on each of ${try} tries." >&2
      return 1
    fi
    echo "Registry 5xx pulling ${image} (try ${try} of ${PULL_TRIES}); trying again in ${wait}s." >&2
    sleep "$wait"
    try=$((try + 1)) wait=$((wait * 2))
  done
}
# Images no case uses: those of the steps in docs/testing.md, "What no e2e case
# runs", plus BWA_IMAGE (the classic index GRIDSS needs, step 04b) and
# SAMTOOLS_HTTPS_IMAGE (docs/quick-test.md). Pulling them would only cost time
# and cache space. PICARD_IMAGE is not one of them: consumer-export-1 runs
# chip-to-vcf.sh, whose liftover is Picard's.
UNUSED_IMAGES=" SAMTOOLS_HTTPS_IMAGE BWA_IMAGE STRELKA_IMAGE OCTOPUS_IMAGE CLAIR3_IMAGE TIDDIT_IMAGE SNIFFLES_IMAGE GRIDSS_IMAGE DUPHOLD_IMAGE CNVPYTOR_IMAGE ANNOTSV_IMAGE VEP_IMAGE PCGR_IMAGE "
# case_images: the NAME_IMAGE="ref" lines of versions.env and of the cases,
# without the images above, once each.
case_images() {
  grep -hoE '^[A-Z0-9_]+_IMAGE="[^"]+"' "${REPO}/versions.env" "${REPO}"/tests/e2e/*.sh \
    | while IFS='=' read -r name ref; do
        [[ "$UNUSED_IMAGES" == *" ${name} "* ]] || echo "${ref//\"/}"
      done | awk '!seen[$0]++'
}

# Images pulled at first use, not before the cases: the five largest (1.9 to
# 6.0 GB each compressed; DeepVariant's, GATK's and VNtyper's alone take
# 17.6 GB unpacked), which no image cache entry holds (pinned by digest or over
# the group budget) and which only some legs run. Pulled up front, every leg
# carried all five. A `docker run` or `docker create` of one that is missing
# pulls it first, with the retry above (tests/e2e/bin/docker reads
# E2E_LAZY_IMAGES).
LAZY_IMAGES=" DEEPVARIANT_IMAGE PHARMCAT_IMAGE GATK_IMAGE MEGANE_IMAGE VNTYPER_IMAGE "
# lazy_images: the references of LAZY_IMAGES, from versions.env.
lazy_images() {
  grep -hoE '^[A-Z0-9_]+_IMAGE="[^"]+"' "${REPO}/versions.env" \
    | while IFS='=' read -r name ref; do
        [[ "$LAZY_IMAGES" == *" ${name} "* ]] && echo "${ref//\"/}"
      done
  return 0
}

# --- Image cache groups -------------------------------------------------------
# The e2e workflow caches the images of case_images in IMAGE_GROUPS cache
# entries, so no single entry, and not their sum, reaches the repository's
# 10 GB cache limit. An image's group comes from its repository name without
# the tag or digest, so a version bump moves no image to another group and
# changes the cache key of one group only.
IMAGE_GROUPS=4
# image_group REF: the group (0 to IMAGE_GROUPS-1) of an image reference.
image_group() {
  local repo=${1%@*}
  [[ "${repo##*/}" == *:* ]] && repo=${repo%:*}
  echo $(( $(printf '%s' "$repo" | cksum | cut -d' ' -f1) % IMAGE_GROUPS ))
}
# image_groups: one "GROUP REF" line per image of case_images.
image_groups() {
  local ref
  case_images | while IFS= read -r ref; do echo "$(image_group "$ref") ${ref}"; done
}

if [ "${1:-}" = --image-groups ]; then
  image_groups
  exit 0
fi

# --- Data cache keys and checks -------------------------------------------------
# The e2e workflow asks here for the Actions cache key of the entries hla,
# pgs and catalogs (`--data-key NAME`) and whether an entry is whole before it saves it
# (`--data-check NAME`), so the key and the check are written once. The other
# entries (mocha, archaic) are keyed in e2e.yml on their versions.env pins.
DATA_CACHE=${E2E_DATA_CACHE:-}
HLA_CASE="${REPO}/tests/e2e/bash-step-08-hla.sh"
PGS_CASE="${REPO}/tests/e2e/prs-scores-and-run-options-1-strand-sex.sh"
# PGS_MD5: "ID MD5" of the score files the e2e-data-1 release holds, the md5
# the PGS Catalog served beside each on 2026-11-06 (its .md5 file).
# E2E_PGS_MD5 replaces the list (tests/test_e2e_data_cache.sh).
PGS_MD5=${E2E_PGS_MD5-"PGS002257 bb1b2b11a499c5930e8c835647e48296
PGS000065 90ac7aabbbde9bdaa681dbdb66c0ebd8
PGS000064 43b39c4450f423f0772bbd8779dd640f
PGS000066 b34f31b8b3a7dbbbb58882a70e14649d"}
# CATALOGS: the builders of the committed tables two cases rebuild and
# compare, as "FOLDER COMMAND...": each lists its pinned downloads with
# --list ("NAME URL SHA256"), and entry catalogs holds them in FOLDER under
# NAME, where the builder's cache reads them (callability-1-catalog,
# cn-intervals-1-catalog). GENCODE's GTFs come from
# EBI's FTP server, which refused every connection on 2026-11-06 (run
# 42259432681).
CATALOGS="callability bash ${REPO}/scripts/ci/build-callability.sh
cn-intervals python3 ${REPO}/scripts/ci/build-cn-intervals.py"
# catalog_files: "FOLDER NAME URL SHA256" of every download of CATALOGS.
catalog_files() {
  local dir interp builder out
  while read -r dir interp builder; do
    out=$("$interp" "$builder" --list) && [ -n "$out" ] || { echo "ERROR: ${builder} --list gave nothing" >&2; return 1; }
    sed "s|^|${dir} |" <<< "$out"
  done <<< "$CATALOGS"
}
# hla_releases: the IPD-IMGT/HLA releases bash-step-08-hla installs (the
# release before the pin, then the pin), from the case and scripts/lib/common.sh.
hla_releases() {
  local other pin
  other=$(sed -n 's/^OTHER=//p' "$HLA_CASE")
  pin=$(sed -n 's/^HLA_DB_RELEASE=\${HLA_DB_RELEASE:-\(.*\)}$/\1/p' "${REPO}/scripts/lib/common.sh")
  [ -n "$other" ] && [ -n "$pin" ] || { echo "ERROR: no OTHER= line in ${HLA_CASE} or no HLA_DB_RELEASE default in scripts/lib/common.sh" >&2; return 1; }
  echo "$other $pin"
}
# data_key NAME: the Actions cache key of entry NAME. hla: the two releases,
# GENCODE's release and a hash of the code that installs them (the "Pinned
# data files" part of scripts/lib/common.sh), since the entry holds what that
# code writes. pgs: a hash of the score ids and the server they come from.
data_key() {
  local r g h
  case "$1" in
    hla)
      r=$(hla_releases) || return 1
      g=$(sed -n 's/^GENCODE_RELEASE=//p' "${REPO}/scripts/lib/common.sh")
      h=$(awk '/^# --- Pinned data files/,/^_keep_lines\(\)/' "${REPO}/scripts/lib/common.sh" | sha256sum | cut -c1-16)
      echo "e2e-data-hla-${r// /-}-gencode${g}-${h}" ;;
    pgs)
      h=$(grep -E '^(NEW_IDS|PGS_BASE)=' "$PGS_CASE" | sha256sum | cut -c1-16)
      echo "e2e-data-pgs-${h}" ;;
    catalogs)
      h=$(catalog_files | sha256sum | cut -c1-16) && [ "${#h}" -eq 16 ] || return 1
      echo "e2e-data-catalogs-${h}" ;;
    *) echo "ERROR: no data cache key for ${1}" >&2; return 1 ;;
  esac
}
# data_ok NAME: 0 when entry NAME in E2E_DATA_CACHE is whole and checks out.
# hla: every file its SHA256SUMS lists, which data_fetch writes after the
# downloads passed their md5 checks, matches. pgs: each score of the case is
# there with the md5 file the PGS Catalog serves beside it, and matches it.
data_ok() {
  local d="${DATA_CACHE:?}/$1" id pin sum
  case "$1" in
    hla)
      [ -s "${d}/SHA256SUMS" ] && [ "$(wc -l < "${d}/SHA256SUMS")" -eq 4 ] \
        && (cd "$d" && sha256sum --quiet -c SHA256SUMS) ;;
    pgs)
      # A score pinned in PGS_MD5 must match its pin as well as its .md5.
      for id in $(sed -n 's/^NEW_IDS="\(.*\)"$/\1/p' "$PGS_CASE"); do
        pin=$(awk -v i="$id" '$1 == i {print $2}' <<< "$PGS_MD5")
        sum=$(md5sum "${d}/${id}.txt.gz" 2> /dev/null | cut -d' ' -f1)
        [ -s "${d}/${id}.txt.gz" ] && [ -s "${d}/${id}.txt.gz.md5" ] \
          && [ "$sum" = "$(awk '{print $1; exit}' "${d}/${id}.txt.gz.md5")" ] \
          && { [ -z "$pin" ] || [ "$sum" = "$pin" ]; } \
          || return 1
      done ;;
    catalogs)
      # Every download of every builder, matching its pin.
      local list dir name url
      list=$(catalog_files) || return 1
      while read -r dir name url pin; do
        [ -s "${d}/${dir}/${name}" ] && [ "$(sha256sum "${d}/${dir}/${name}" | cut -d' ' -f1)" = "$pin" ] || return 1
      done <<< "$list" ;;
    *) return 1 ;;
  esac
}
# --- Downloads kept between runs ----------------------------------------------
# The two largest downloads of the cases come from EBI's 1000 Genomes server,
# which is slow and stalls at times: MoChA's chr20 panel (mocha-1, about
# 655 MB, with its manifest) and the 1000 Genomes phase 3 sites VCF
# (archaic-ancestry-1, about 1.5 GB). With E2E_DATA_CACHE set (the e2e
# workflow restores it from the Actions cache), the runner gives each of
# these cases its files before it starts, at the paths setup.sh reads them
# from, so setup.sh finds them and downloads nothing. A file missing from
# E2E_DATA_CACHE is downloaded there first, checked against its pin. The
# files are hard links where they can be (setup.sh deletes its copy, the
# cache keeps its own). Without E2E_DATA_KEEP=true (the workflow sets it on
# a run that saves the cache) a case's cached files are deleted after it,
# to give the disk back.
# E2E_DATA_MIRROR, when set, is a folder URL (a GitHub release) holding copies
# of these pinned files under their own names. Each file is tried there once
# first, then at its source. EBI's 1000 Genomes hosts swing between a few KB/s
# and MB/s from one request to the next, and a leg with no Actions cache could
# spend its whole time limit on one download. The pin checks both copies the
# same way, so the mirror can only save time, never change what a case reads.
# kill_tree PID: PID and every process under it (data_fetch's subshells and
# their curl).
# The children are listed first and the parent is stopped before them: a
# parent still running when its child dies would start the next command
# (mfetch's source download after the mirror try), which the list missed.
kill_tree() {
  local c kids
  kids=$(pgrep -P "$1" 2> /dev/null || true)
  kill -TERM "$1" 2> /dev/null || true
  for c in $kids; do kill_tree "$c"; done
}
# within SECS CMD...: runs CMD (a function too) and stops it with everything
# under it after SECS seconds; returns CMD's status, or 124 when it was stopped.
within() {
  local secs=$1 pid wd rc=0 over
  shift
  over=$(mktemp)
  "$@" & pid=$!
  ( sleep "$secs"; echo over > "$over"; kill_tree "$pid" ) & wd=$!
  wait "$pid" || rc=$?
  kill_tree "$wd"
  wait "$wd" 2> /dev/null || true
  if [ -s "$over" ]; then rm -f "$over"; return 124; fi
  rm -f "$over"
  return "$rc"
}
# mfetch URL DEST KIND SUM: fetch from the mirror, then from URL. Called
# inside data_fetch, after it sources scripts/lib/common.sh for fetch. The
# mirror try gets E2E_MIRROR_BUDGET seconds (default 900, or half of
# E2E_DATA_BUDGET when that is less), so a mirror that trickles leaves the
# source its turn within the fill's budget.
mfetch() {
  local secs=${E2E_MIRROR_BUDGET:-} rc=0
  if [ -n "${E2E_DATA_MIRROR:-}" ]; then
    if [ -z "$secs" ]; then
      secs=$(( ${E2E_DATA_BUDGET:-2400} / 2 ))
      [ "$secs" -le 900 ] || secs=900
      [ "$secs" -ge 1 ] || secs=1
    fi
    FETCH_TRIES=1 within "$secs" fetch "${E2E_DATA_MIRROR%/}/${1##*/}" "$2" "$3" "$4" 2> /dev/null || rc=$?
    if [ "$rc" -eq 0 ]; then
      echo "Data cache: ${1##*/} from the mirror"
      return 0
    fi
    if [ "$rc" -eq 124 ]; then
      echo "Data cache: ${1##*/} took more than ${secs} s on the mirror; downloading it from its source"
    else
      echo "Data cache: ${1##*/} not on the mirror; downloading it from its source"
    fi
  fi
  fetch "$@"
}
# data_fetch NAME: fills E2E_DATA_CACHE/NAME with what is missing, in a
# subshell that sources scripts/lib/common.sh for fetch and the pins.
data_fetch() (
  d="${DATA_CACHE}/$1"
  # shellcheck source=../lib/common.sh
  . "${REPO}/scripts/lib/common.sh" || exit 1
  # More and longer-spaced tries than a user's run: EBI and GitHub's raw host
  # stayed down for minutes at a time on 2026-11-04 (curl exit 7 and 28).
  FETCH_TRIES=${E2E_FETCH_TRIES:-6} FETCH_WAIT=${E2E_FETCH_WAIT:-30}
  case "$1" in
    mocha)
      m="${d}/manifest.txt"
      { [ -s "$m" ] && [ "$(_digest sha256 "$m")" = "$MOCHA_PANEL_MANIFEST_SHA256" ]; } \
        || mfetch "${MOCHA_PANEL_URL}/20220804_manifest.txt" "$m" sha256 "$MOCHA_PANEL_MANIFEST_SHA256" || exit 1
      # mocha-1 cuts the panel to the fixture's chr20 slice.
      f=1kGP_high_coverage_Illumina.chr20.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
      for e in "" .tbi; do
        md5=$(awk -v n="${f}${e}" '$1 == n {print $2}' "$m")
        [ -n "$md5" ] || { echo "ERROR: ${f}${e} is not in the panel manifest" >&2; exit 1; }
        [ -s "${d}/raw/${f}${e}" ] && [ "$(_digest md5 "${d}/raw/${f}${e}")" = "$md5" ] && continue
        mfetch "${MOCHA_PANEL_URL}/${f}${e}" "${d}/raw/${f}${e}" md5 "$md5" || exit 1
      done ;;
    archaic)
      url=$(sed -n 's/^ARCHAIC_1KG_SITES_URL=//p' "${REPO}/scripts/setup.sh")
      [ -n "$url" ] || { echo "ERROR: no ARCHAIC_1KG_SITES_URL line in scripts/setup.sh" >&2; exit 1; }
      s="${d}/1kg_sites.vcf.gz"
      [ -s "$s" ] && [ "$(_digest md5 "$s")" = "$ARCHAIC_1KG_SITES_MD5" ] && exit 0
      mfetch "$url" "$s" md5 "$ARCHAIC_1KG_SITES_MD5" || exit 1 ;;
    hla)
      # The files install_data_file writes, made by install_data_file itself
      # with the entry as GENOME_DIR: hla.dat of both releases (md5 from
      # IPD-IMGT/HLA's md5checksum.txt) and GENCODE's genes and exons (md5
      # from its MD5SUMS). Their sha256 go to SHA256SUMS, which data_ok checks.
      data_ok hla && exit 0
      rm -rf "$d" && mkdir -p "$d" || exit 1
      read -r other pin <<< "$(hla_releases)"
      [ -n "$pin" ] || exit 1
      # Each install is tried again as a whole: install_data_file reads
      # hla.dat.zip's md5 from md5checksum.txt in one request, outside fetch's
      # tries. What an earlier try installed is kept.
      ( GENOME_DIR=$d
        unset SAMPLE
        for job in "$other hla_dat" "$pin hla_dat" "$pin gencode_genes" "$pin gencode_exons"; do
          read -r rel name <<< "$job"
          for ((t = 1; ; t++)); do
            HLA_DB_RELEASE=$rel install_data_file "$name" && break
            [ "$t" -lt 3 ] || exit 1
            echo "Data cache: installing ${name} (${rel}) failed, try ${t} of 3; again in ${FETCH_WAIT} s" >&2
            sleep "$FETCH_WAIT"
          done
        done ) || exit 1
      (cd "$d" && find . -type f ! -name 'SHA256SUMS*' ! -name '*.part' | LC_ALL=C sort | xargs sha256sum) > "${d}/SHA256SUMS.tmp" \
        && mv "${d}/SHA256SUMS.tmp" "${d}/SHA256SUMS" && data_ok hla || exit 1 ;;
    pgs)
      # The four real PGS Catalog files of prs-scores-and-run-options-1, each
      # with the .md5 file the PGS Catalog serves beside it. A score pinned
      # in PGS_MD5 is tried on the mirror first, checked against the pin, and
      # its .md5 file is written from the pin; any other score comes from
      # the PGS Catalog with its served .md5, as before.
      base=$(sed -n 's/^PGS_BASE=//p' "$PGS_CASE")
      [ -n "$base" ] || { echo "ERROR: no PGS_BASE= line in ${PGS_CASE}" >&2; exit 1; }
      for id in $(sed -n 's/^NEW_IDS="\(.*\)"$/\1/p' "$PGS_CASE"); do
        f="${d}/${id}.txt.gz" url="${base}/${id}/ScoringFiles/Harmonized/${id}_hmPOS_GRCh38.txt.gz"
        pin=$(awk -v i="$id" '$1 == i {print $2}' <<< "$PGS_MD5")
        [ -s "$f" ] && [ -s "${f}.md5" ] && [ "$(_digest md5 "$f")" = "$(awk '{print $1; exit}' "${f}.md5")" ] \
          && { [ -z "$pin" ] || [ "$(_digest md5 "$f")" = "$pin" ]; } && continue
        rm -f "$f" "${f}.md5"
        if [ -n "$pin" ]; then
          mfetch "$url" "$f" md5 "$pin" || exit 1
          echo "${pin}  ${url##*/}" > "${f}.md5"
        else
          fetch "${url}.md5" "${f}.md5" || exit 1
          fetch "$url" "$f" md5 "$(awk '{print $1; exit}' "${f}.md5")" || exit 1
        fi
      done
      data_ok pgs || exit 1 ;;
    catalogs)
      list=$(catalog_files) || exit 1
      while read -r dir name url sum; do
        f="${d}/${dir}/${name}"
        [ -s "$f" ] && [ "$(_digest sha256 "$f")" = "$sum" ] && continue
        rm -f "$f"
        mfetch "$url" "$f" sha256 "$sum" || exit 1
      done <<< "$list"
      data_ok catalogs || exit 1 ;;
    *) echo "ERROR: no data cache entry ${1}" >&2; exit 1 ;;
  esac
)
# data_fill NAME: data_fetch NAME within E2E_DATA_BUDGET seconds (default
# 2400). A fill that is still running then is stopped, its partial files
# kept for the next try, and data_fill returns 124 with a message that says
# so. Each try of fetch has its own stall limit, but a source that trickles
# without stalling is never cut by it: on 2026-11-05 EBI served the 1000
# Genomes sites VCF slowly enough that two shard legs reached their 240-minute
# limit on that one download (fork run 41957608395).
data_fill() {
  local budget=${E2E_DATA_BUDGET:-2400} rc=0
  case "$budget" in ''|*[!0-9]*) echo "ERROR: E2E_DATA_BUDGET must be a whole number of seconds, got '${budget}'" >&2; return 1 ;; esac
  within "$budget" data_fetch "$1" || rc=$?
  if [ "$rc" -eq 124 ]; then
    echo "ERROR: Data cache: filling ${1} took more than E2E_DATA_BUDGET=${budget} s and was stopped; a source this slow would eat the job's time limit" >&2
  fi
  return "$rc"
}
case "${1:-}" in
  --data-key) data_key "${2:?--data-key NAME}"; exit ;;
  --data-check) data_ok "${2:?--data-check NAME}"; exit ;;
  --data-fill) data_fill "${2:?--data-fill NAME}"; exit ;;
esac

# pull_all: pull every image of case_images that is not here yet, but the
# lazy ones.
pull_all() {
  local image pulled=0 lazy
  mapfile -t IMAGES < <(case_images)
  [ "${#IMAGES[@]}" -gt 0 ] || { echo "ERROR: no NAME_IMAGE line found in versions.env" >&2; return 1; }
  lazy=$(lazy_images)
  for image in "${IMAGES[@]}"; do
    grep -qxF "$image" <<< "$lazy" && continue
    docker image inspect "$image" > /dev/null 2>&1 && continue
    pull_image "$image" || return 1
    pulled=$((pulled + 1))
  done
  echo "${#IMAGES[@]} images: $(( ${#IMAGES[@]} - pulled - $(wc -l <<< "$lazy") )) already here, ${pulled} pulled, $(wc -l <<< "$lazy") left to their first use"
}

# The e2e workflow's shard legs pull while the base leg runs.
if [ "${1:-}" = --pull-images ]; then
  pull_all
  exit
fi
# --pull-image REF: one pull with the retry (the docker shim's lazy pulls).
if [ "${1:-}" = --pull-image ]; then
  pull_image "${2:?--pull-image needs an image}"
  exit
fi

# --- Free disk ------------------------------------------------------------------
# The legs ran out of disk pulling images on 2026-11-05 (run 41938436657). Before
# the cases a leg needs E2E_MIN_FREE_GB free under E2E_WORK: room for what its
# cases write and for the lazy images they pull.
# disk_check NEED_GB: prints the free space; fails when it is under NEED_GB.
disk_check() {
  local free
  free=$(( $(df --output=avail -k "$E2E_WORK" | tail -1) / 1048576 ))
  echo "Free disk before the cases: ${free} GB under ${E2E_WORK} (the leg needs ${1} GB)"
  if [ "$free" -lt "$1" ]; then
    echo "ERROR: ${free} GB free under ${E2E_WORK}, under the ${1} GB the cases need." >&2
    return 1
  fi
}

if [ "${1:-}" = --self-test-pull ]; then
  T=$(mktemp -d)
  trap 'rm -rf "$T"' EXIT
  # A fake docker that only pulls: each image name has one behaviour, and every
  # pull is counted in $T/<image>.count.
  cat > "${T}/docker" <<'FAKE'
#!/usr/bin/env bash
[ "$1" = pull ] || { echo "fake docker: only pull is expected, got: $*" >&2; exit 2; }
img=${*: -1}
c="${FAKE_DIR}/${img//[\/:]/_}.count"
echo x >> "$c"
n=$(wc -l < "$c")
case "$img" in
  flaky500)
    [ "$n" -ge 2 ] && exit 0
    echo 'Error response from daemon: Head "https://registry-1.docker.io/v2/library/flaky500/manifests/latest": received unexpected HTTP status: 500 Internal Server Error' >&2 ;;
  down503)
    echo 'Error response from daemon: Head "https://registry-1.docker.io/v2/library/down503/manifests/latest": received unexpected HTTP status: 503 Service Unavailable' >&2 ;;
  notag)
    echo 'Error response from daemon: manifest for notag not found: manifest unknown: manifest unknown' >&2 ;;
  ratelimit)
    echo 'Error response from daemon: toomanyrequests: You have reached your unauthenticated pull rate limit. https://www.docker.com/increase-rate-limit' >&2 ;;
esac
exit 1
FAKE
  chmod +x "${T}/docker"
  export FAKE_DIR="$T"
  PATH="${T}:${PATH}" PULL_WAIT=0
  st_failed=0
  # expect DESC pass|fail PULLS IMAGE
  expect() {
    local rc=0 got n
    pull_image "$4" > /dev/null 2> "${T}/err" || rc=$?
    got=$([ "$rc" -eq 0 ] && echo pass || echo fail)
    n=$({ cat "${T}/${4//[\/:]/_}.count" 2> /dev/null || true; } | wc -l | tr -d ' ')
    if [ "$got" = "$2" ] && [ "$n" = "$3" ]; then
      echo "ok: $1"
    else
      echo "FAIL: $1 (${got} after ${n} pulls; want ${2} after ${3})"
      sed 's/^/  /' "${T}/err"
      st_failed=1
    fi
  }
  expect "a registry 500 on the first pull is tried again and passes" pass 2 flaky500
  expect "a registry 503 on every pull gives up after ${PULL_TRIES} tries" fail "$PULL_TRIES" down503
  expect "a missing tag fails at once" fail 1 notag
  expect "a rate limit (429) fails at once" fail 1 ratelimit
  LIST=$(case_images)
  # shellcheck source=../../versions.env
  HAPPY=$(. "${REPO}/versions.env" && echo "$HAPPY_IMAGE") VEP=$(. "${REPO}/versions.env" && echo "$VEP_IMAGE")
  if grep -qxF "$HAPPY" <<< "$LIST"; then echo "ok: the pull list holds hap.py's image"; else echo "FAIL: the pull list lacks ${HAPPY}"; st_failed=1; fi
  if grep -qxF "$VEP" <<< "$LIST"; then echo "FAIL: the pull list holds ${VEP}, which no case runs"; st_failed=1; else echo "ok: the pull list leaves VEP out"; fi
  # Cache groups: a version bump keeps an image in its group, and every image
  # of the pull list is in exactly one group.
  for pair in "example/tool:2.31--h118bc1c_0 example/tool:2.32--h0_0" \
              "example/pinned:3.15.0@sha256:16e7e5fc example/pinned:3.16.0@sha256:0123abcd" \
              "registry.example:5000/tool:1.0 registry.example:5000/tool:2.0"; do
    read -r a b <<< "$pair"
    if [ "$(image_group "$a")" = "$(image_group "$b")" ]; then echo "ok: ${a} and ${b} share a cache group"
    else echo "FAIL: ${a} and ${b} fall in different cache groups"; st_failed=1; fi
  done
  GROUPS_LIST=$(image_groups)
  if [ "$(cut -d' ' -f2 <<< "$GROUPS_LIST")" = "$LIST" ] && ! grep -qvE "^[0-9]+ " <<< "$GROUPS_LIST" \
     && [ -z "$(awk -v n="$IMAGE_GROUPS" '$1 >= n' <<< "$GROUPS_LIST")" ]; then
    echo "ok: each of the $(wc -l <<< "$LIST" | tr -d ' ') listed images is in one cache group of ${IMAGE_GROUPS}"
  else
    echo "FAIL: the cache groups do not hold the pull list once each"; st_failed=1
  fi
  # Lazy images: each one a case image, five in all.
  LAZY=$(lazy_images)
  if [ "$(wc -l <<< "$LAZY")" -eq 5 ] && ! grep -qvxF -f <(echo "$LIST") <<< "$LAZY"; then
    echo "ok: the 5 lazy images are case images"
  else
    echo "FAIL: the lazy images are not 5 case images: ${LAZY//$'\n'/ }"; st_failed=1
  fi
  # disk_check against a fake df: 10 GB free fails a need of 11 and passes 10.
  # shellcheck disable=SC2016  # the fake df expands FAKE_AVAIL_K when it runs
  printf '#!/usr/bin/env bash\nprintf "Avail\\n%%s\\n" "$FAKE_AVAIL_K"\n' > "${T}/df"
  chmod +x "${T}/df"
  E2E_WORK=$T
  if FAKE_AVAIL_K=$((10 * 1048576)) disk_check 11 > /dev/null 2>&1; then
    echo "FAIL: 10 GB free passed a need of 11 GB"; st_failed=1
  else echo "ok: 10 GB free fails a need of 11 GB"; fi
  if FAKE_AVAIL_K=$((10 * 1048576)) disk_check 10 > /dev/null 2>&1; then echo "ok: 10 GB free passes a need of 10 GB"
  else echo "FAIL: 10 GB free failed a need of 10 GB"; st_failed=1; fi
  # The docker shim: a run of a missing lazy image pulls it once, then runs;
  # a run of any other image only runs.
  mkdir -p "${T}/real"
  cat > "${T}/real/docker" <<'FAKE'
#!/usr/bin/env bash
echo "$*" >> "${FAKE_DIR}/calls"
case "$1 ${2:-}" in
  "image inspect") awk -v i="${*: -1}" '$1 == "pull" && $NF == i {f = 1} END {exit !f}' "${FAKE_DIR}/calls" ;;
  *) exit 0 ;;
esac
FAKE
  chmod +x "${T}/real/docker"
  : > "${T}/calls"
  PATH="${REPO}/tests/e2e/bin:${T}/real:${PATH}" E2E_LAZY_IMAGES="lazyone lazytwo" \
    docker run --rm -v /a:/a lazytwo tool --x > /dev/null 2>&1
  PATH="${REPO}/tests/e2e/bin:${T}/real:${PATH}" E2E_LAZY_IMAGES="lazyone lazytwo" \
    docker run --rm lazytwo tool > /dev/null 2>&1
  PATH="${REPO}/tests/e2e/bin:${T}/real:${PATH}" E2E_LAZY_IMAGES="lazyone lazytwo" \
    docker run --rm plainthree tool > /dev/null 2>&1
  if [ "$(grep -c '^pull ' "${T}/calls")" = 1 ] && grep -qx 'pull -q lazytwo' "${T}/calls" \
     && [ "$(grep -c '^run ' "${T}/calls")" = 3 ]; then
    echo "ok: the docker shim pulls a missing lazy image once, before its run"
  else
    echo "FAIL: the docker shim's lazy pulls:"; sed 's/^/  /' "${T}/calls"; st_failed=1
  fi
  exit "$st_failed"
fi

REPORT_DIR="${E2E_WORK}/report"
JOB=${E2E_JOB:-local}
TRACE=${E2E_TRACE:-false}
CASE_LIST=${E2E_CASES:-}
if [ -n "$CASE_LIST" ]; then
  [ -f "$CASE_LIST" ] || { echo "ERROR: E2E_CASES names no file: ${CASE_LIST}" >&2; exit 1; }
  [ "$PATTERN" = "*" ] || { echo "ERROR: give E2E_CASES or a pattern, not both" >&2; exit 1; }
fi
mkdir -p "$E2E_WORK" "$FIXTURE_DIR" "$LOG_DIR" "$REPORT_DIR"
: > "$E2E_NOTES"
: > "${REPORT_DIR}/results.tsv"
: > "${REPORT_DIR}/extra.md"

fixture_complete() {
  [ -f "${FIXTURE_DIR}/SHA256SUMS" ] && (cd "$FIXTURE_DIR" && sha256sum --quiet -c SHA256SUMS)
}

# --- 1. Fixture ---------------------------------------------------------------
echo "=== Fixture ${TAG} ==="
deadline=$(( $(date +%s) + 75 * 60 ))
until fixture_complete; do
  if gh release view "$TAG" -R "$GH_REPO" >/dev/null 2>&1; then
    # Inside the condition, so a failed download (assets still uploading, a
    # network error) waits for the next try instead of ending the run.
    if gh release download "$TAG" -R "$GH_REPO" -D "$FIXTURE_DIR" --clobber && fixture_complete; then
      break
    fi
    echo "Release ${TAG} could not be downloaded, is incomplete or its checksums do not match."
  else
    echo "Release ${TAG} does not exist yet; the build-fixture job publishes it."
  fi
  if [ "$(date +%s)" -ge "$deadline" ]; then
    echo "ERROR: fixture ${TAG} not usable after 75 minutes. Run the E2E workflow with job=build-fixture." >&2
    exit 1
  fi
  sleep 60
done
(cd "$FIXTURE_DIR" && sha256sum -c SHA256SUMS)

# --- 2. GENOME_DIR layout -----------------------------------------------------
# The reference under the name scripts/lib/common.sh reads by default.
echo "=== Layout ${GENOME_DIR} ==="
mkdir -p "${GENOME_DIR}/reference" "${GENOME_DIR}/clinvar" "${GENOME_DIR}/annotations" \
  "${GENOME_DIR}/${SAMPLE}/fastq" "${GENOME_DIR}/${SAMPLE}/vep"
REF="${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set"
if [ ! -s "${REF}.fasta" ]; then
  gzip -dc "${FIXTURE_DIR}/fixture_ref.fa.gz" > "${REF}.fasta.tmp"
  mv "${REF}.fasta.tmp" "${REF}.fasta"
fi
cp "${FIXTURE_DIR}/fixture_ref.fa.gz.fai" "${REF}.fasta.fai"
cp "${FIXTURE_DIR}/fixture_ref.dict" "${REF}.dict"
for f in clinvar.vcf.gz clinvar_chr.vcf.gz clinvar_pathogenic_chr.vcf.gz; do
  cp "${FIXTURE_DIR}/${f}" "${FIXTURE_DIR}/${f}.tbi" "${GENOME_DIR}/clinvar/"
done
# The synthetic score file stands in for REVEL (step 30 looks for this name).
cp "${FIXTURE_DIR}/revel_synthetic.tsv.gz" "${GENOME_DIR}/annotations/revel_grch38.tsv.gz"
cp "${FIXTURE_DIR}/revel_synthetic.tsv.gz.tbi" "${GENOME_DIR}/annotations/revel_grch38.tsv.gz.tbi"
# Reads for step 02, and the VEP output step 13 would have written.
cp "${FIXTURE_DIR}/${SAMPLE}_R1.fastq.gz" "${FIXTURE_DIR}/${SAMPLE}_R2.fastq.gz" "${GENOME_DIR}/${SAMPLE}/fastq/"
cp "${FIXTURE_DIR}/${SAMPLE}_vep.vcf" "${GENOME_DIR}/${SAMPLE}/vep/"
df -h "$E2E_WORK"
# Files whose ctime is newer than this were written by the cases (or unpacked
# from the base job), not by the layout above.
touch "${E2E_WORK}/.layout-done"

# The files the cases of the base job wrote, at the paths they were written to.
# base-from names the work area and checkout of the packing run; anywhere else
# the paths in the archive would land beside the files the cases read.
if [ -n "${E2E_BASE_IN:-}" ]; then
  echo "=== Unpacking the base cases' outputs ==="
  # Checked before anything is written: an archive from another work area
  # would put its files where no case reads them.
  from=$(zstd -dc "$E2E_BASE_IN" | tar -P -xOf - "${E2E_WORK}/base-from" 2> /dev/null || true)
  if [ "$from" != "${E2E_WORK} ${REPO}" ]; then
    echo "ERROR: ${E2E_BASE_IN} was not packed with E2E_WORK=${E2E_WORK} and the checkout at ${REPO} (it says: ${from:-nothing})" >&2
    exit 1
  fi
  zstd -dc "$E2E_BASE_IN" | sudo tar -P -xpf -
  echo "Unpacked $(zstd -dc "$E2E_BASE_IN" | tar -tf - | wc -l) entries"
  # The archive is about 4 GB on the disk the cases write to. Shard 1 of run
  # 41524201317 filled that disk while the archive was still there.
  rm -f "$E2E_BASE_IN"
fi

# --- 3. Images ----------------------------------------------------------------
# Only a full run pre-pulls. A case reaches most images through the steps it
# runs, so the images of a few cases cannot be told apart from the rest; a
# partial run lets its cases pull their own, without the retry.
if [ "$PATTERN" = "*" ]; then
  echo "=== Pulling the images the cases use ==="
  pull_all
  E2E_LAZY_IMAGES=$(lazy_images | paste -sd' ' -)
  export E2E_LAZY_IMAGES
else
  echo "=== Partial run ('${PATTERN}'): no pre-pull, each case pulls what it uses ==="
fi
docker system df
disk_check "${E2E_MIN_FREE_GB:-0}"
# The least free space while the cases run, sampled every 10 s, for the report.
LOW_FILE="${REPORT_DIR}/disk-low"
df --output=avail -k "$E2E_WORK" | tail -1 > "$LOW_FILE"
( while sleep 10; do
    a=$(df --output=avail -k "$E2E_WORK" | tail -1)
    [ "$a" -lt "$(cat "$LOW_FILE")" ] && echo "$a" > "$LOW_FILE"
  done ) &
LOW_PID=$!
trap 'kill "$LOW_PID" 2> /dev/null || true' EXIT

# data_name CASE: the E2E_DATA_CACHE entry CASE reads, and where setup.sh
# wants it ("-": the case reads the entry itself), then "offline" when the
# case runs with the network blocked once its files are in place.
#   bash-step-08-hla: hla.dat of two IPD-IMGT/HLA releases (raw.githubusercontent.com,
#     which served it under 1 KB/s at times) and GENCODE's genes and exons
#     (ftp.ebi.ac.uk), at the paths data_file names, so install_data_file
#     finds them and downloads nothing.
#   prs-scores-and-run-options-1: the four real PGS Catalog score files
#     (ftp.ebi.ac.uk); the case copies them from E2E_DATA_CACHE/pgs.
#   callability-1-catalog, cn-intervals-1-catalog: the pinned downloads of their builders (CATALOGS above); each case
#     points its builder's cache at E2E_DATA_CACHE/catalogs/<its folder>.
data_name() {
  # shellcheck source=../../versions.env
  case "$1" in
    mocha-1.sh) echo "mocha ${GENOME_DIR}/reference/mocha-$(. "${REPO}/versions.env" && echo "$MOCHA_PANEL_VERSION")" ;;
    archaic-ancestry-1.sh) echo "archaic ${GENOME_DIR}/reference/archaic/$(. "${REPO}/versions.env" && echo "$ARCHAIC_DATA_VERSION")/download" ;;
    bash-step-08-hla.sh) echo "hla ${GENOME_DIR} offline" ;;
    prs-scores-and-run-options-1-strand-sex.sh) echo "pgs -" ;;
    callability-1-catalog.sh|cn-intervals-1-catalog.sh) echo "catalogs -" ;;
  esac
}
# Every HTTP(S) request of a process run under it goes to a closed port.
OFFLINE=(env http_proxy=http://127.0.0.1:9 https_proxy=http://127.0.0.1:9 HTTP_PROXY=http://127.0.0.1:9
         HTTPS_PROXY=http://127.0.0.1:9 no_proxy= NO_PROXY= FETCH_TRIES=1 FETCH_WAIT=0)
# offline_blocked CURL_ERROR: true when curl's error says the request went to
# the closed proxy of OFFLINE, in either wording curl uses.
offline_blocked() {
  grep -qE '127\.0\.0\.1( port |:)9\b|via 127\.0\.0\.1\b' <<< "$1"
}
# data_seed CASE: before CASE, its files from E2E_DATA_CACHE. A failure here
# only warns and costs the download: setup.sh in the case fetches what is
# not there. A fill that runs past E2E_DATA_BUDGET sets CASE_DATA_FAILED
# instead: the case does not run (its own download would be as slow) and
# counts as failed. Sets CASE_ENV to OFFLINE when the case's entry says "offline"
# and every file is in place, after a control shows the block holds.
data_seed() {
  local n dest flag f t placed=true
  CASE_ENV=() CASE_CONTROL_FAILED=false CASE_DATA_FAILED=false
  [ -n "$DATA_CACHE" ] || return 0
  read -r n dest flag <<< "$(data_name "$1")"
  [ -n "$n" ] || return 0
  if [ -n "$(find "${DATA_CACHE}/${n}" -type f ! -name '*.part' 2> /dev/null | head -n 1)" ]; then
    echo "Data cache: ${n} restored from the Actions cache"
  else
    echo "Data cache: ${n} not in the Actions cache; downloading it"
  fi
  local rc=0
  data_fill "$n" || rc=$?
  if [ "$rc" -eq 124 ]; then
    echo "::error::Data cache: filling ${n} took more than E2E_DATA_BUDGET=${E2E_DATA_BUDGET:-2400} s; ${1%.sh} does not run and counts as failed"
    CASE_DATA_FAILED=true
    return 0
  fi
  [ "$rc" -eq 0 ] || { echo "::warning::Data cache: could not fill ${n}; ${1%.sh} downloads it itself"; return 0; }
  [ "$dest" = - ] && { echo "Data cache: ${n} is whole; ${1%.sh} reads it from ${DATA_CACHE}/${n}"; return 0; }
  while IFS= read -r f; do
    t="${dest}/${f#./}"
    mkdir -p "${t%/*}" \
      && { ln -f "${DATA_CACHE}/${n}/${f#./}" "$t" 2> /dev/null || cp -f "${DATA_CACHE}/${n}/${f#./}" "$t"; } \
      || { echo "::warning::Data cache: could not place ${f#./}; ${1%.sh} downloads it itself"; placed=false; }
  done < <(cd "${DATA_CACHE}/${n}" && find . -type f ! -name '*.part' ! -name SHA256SUMS)
  echo "Data cache: ${n}: $(cd "${DATA_CACHE}/${n}" && find . -type f ! -name '*.part' ! -name SHA256SUMS | wc -l | tr -d ' ') files in ${dest}"
  [ "$flag" = offline ] && $placed || return 0
  # Control: under OFFLINE a request to a host that answers must fail, or
  # the offline run would prove nothing.
  # Read from curl's own error, so a host that is down cannot pass for it:
  # the error must name the closed proxy. curl words it two ways: older
  # releases "Failed to connect to 127.0.0.1 port 9", newer ones "Failed to
  # connect to <host> port 443 via 127.0.0.1" (the 26.04 runners, 2026-11-05).
  # (A variable, not a pipe: under pipefail curl's own failure would fail it.)
  local out
  out=$("${OFFLINE[@]}" curl -fsS --max-time 20 -o /dev/null https://raw.githubusercontent.com/ 2>&1 || true)
  if ! offline_blocked "$out"; then
    echo "$out"
    echo "::error::Data cache: a request under the network block did not fail at the closed port; ${1%.sh} runs with the network and counts as failed"
    CASE_CONTROL_FAILED=true
    return 0
  fi
  echo "Data cache: control passed (a request under the block failed); ${1%.sh} runs with every HTTP(S) request sent to a closed port"
  CASE_ENV=("${OFFLINE[@]}")
}
# data_drop CASE: after CASE, its cache files go, unless this run saves them.
data_drop() {
  local n dest
  [ -n "$DATA_CACHE" ] && [ "${E2E_DATA_KEEP:-false}" != true ] || return 0
  read -r n dest <<< "$(data_name "$1")"
  [ -z "$n" ] || rm -rf "${DATA_CACHE:?}/${n}"
}

# --- Trace -------------------------------------------------------------------
# With E2E_TRACE=true each case's reads are told apart by atime, which only
# moves on every read under strictatime, and its writes by ctime, which no
# copy can carry over from the source. trace_scan prints "CTIME ATIME PATH" for
# every file under the work area and the checkout, without the fixture, the
# report and .git. The logs stay in: a case may read an earlier case's log.
TRACE_DIR="${REPORT_DIR}/trace-state"
trace_scan() {
  sudo find "$E2E_WORK" "$REPO" \( -path "$FIXTURE_DIR" -o -path "$REPORT_DIR" -o -path "${REPO}/.git" \) \
    -prune -o -type f -printf '%C@\t%A@\t%p\n' 2> /dev/null || true
}
if [ "$TRACE" = true ]; then
  mkdir -p "$TRACE_DIR"
  for d in "$E2E_WORK" "$REPO"; do sudo mount -o remount,strictatime "$(findmnt -n -o TARGET -T "$d")"; done
  # Positive control: a second read must move atime too (relatime moves it on
  # the first read after a write only), or the trace would miss readers.
  for d in "$E2E_WORK" "$REPO"; do
    probe="${d}/.atime-probe"
    echo x > "$probe"; sleep 1.1; cat "$probe" > /dev/null
    before=$(stat -c %x "$probe"); sleep 1.1; cat "$probe" > /dev/null
    if [ "$(stat -c %x "$probe")" = "$before" ]; then
      echo "ERROR: E2E_TRACE: a second read did not move atime under ${d}" >&2; exit 1
    fi
    rm -f "$probe"
  done
  # Every file there before the first case was written by the layout or the base job.
  trace_scan | awk -F'\t' '{print $3 "\tsetup"}' > "${TRACE_DIR}/writers.tsv"
  : > "${REPORT_DIR}/trace.tsv"
fi
# trace_case NAME START SECONDS: one trace.tsv line for the case that just ran,
# "NAME SECONDS WRITERS READS", WRITERS being the comma-separated cases (or
# "setup") that wrote the files it read, READS how many such files it read.
trace_case() {
  trace_scan > "${TRACE_DIR}/scan.tsv"
  awk -F'\t' -v s="$2" -v k="$1" -v secs="$3" -v out="${TRACE_DIR}/writers.new" '
    FILENAME == ARGV[1] { w[$1] = $2; next }
    {
      p = $3
      if (p in w) {
        if ($2 + 0 >= s + 0) { reads++; if (!(w[p] in seen)) { seen[w[p]] = 1; list = (list == "" ? w[p] : list "," w[p]) } }
        if ($1 + 0 >= s + 0) w[p] = k
      } else w[p] = k
      print p "\t" w[p] > out
    }
    END { printf "%s\t%s\t%s\t%d\n", k, secs, (list == "" ? "-" : list), reads }' \
    "${TRACE_DIR}/writers.tsv" "${TRACE_DIR}/scan.tsv" >> "${REPORT_DIR}/trace.tsv"
  mv "${TRACE_DIR}/writers.new" "${TRACE_DIR}/writers.tsv"
}

# --- 4. Cases -----------------------------------------------------------------
mapfile -t CASES < <(cd "${REPO}/tests/e2e" && find . -maxdepth 1 -type f -name '*.sh' ! -name lib.sh -printf '%f\n' | LC_ALL=C sort)
if [ -n "$CASE_LIST" ]; then
  unknown=$(grep -vxF -f <(printf '%s\n' "${CASES[@]%.sh}") "$CASE_LIST" || true)
  [ -z "$unknown" ] || { echo "ERROR: E2E_CASES lists names that are no case: ${unknown//$'\n'/ }" >&2; exit 1; }
fi
names=() results=() times=()
failed=0
for c in "${CASES[@]}"; do
  name="${c%.sh}"
  if [ -n "$CASE_LIST" ]; then
    grep -qxF "$name" "$CASE_LIST" || continue
  else
    # shellcheck disable=SC2053  # PATTERN is a glob on purpose
    [[ "$c" == $PATTERN ]] || continue
  fi
  log="${LOG_DIR}/${name}.log"
  # A case may raise its own limit with a line "# case-timeout: SECONDS"
  # (never below CASE_TIMEOUT).
  ctimeout=$(sed -n 's/^# case-timeout: \([0-9][0-9]*\)$/\1/p' "${REPO}/tests/e2e/${c}" | head -n 1)
  [ "${ctimeout:-0}" -gt "$CASE_TIMEOUT" ] || ctimeout=$CASE_TIMEOUT
  start=$(date +%s)
  tstart=$(date +%s.%N)
  echo "::group::${name}"
  data_seed "$c"
  if $CASE_DATA_FAILED; then
    echo "${name}: not run: its data fill ran past E2E_DATA_BUDGET" | tee "$log"
    rc=-2
  else
    set +e
    timeout -k 30 "$ctimeout" "${CASE_ENV[@]}" bash "${REPO}/tests/e2e/${c}" 2>&1 | tee "$log"
    rc=${PIPESTATUS[0]}
    set -e
  fi
  data_drop "$c"
  # A lazy image goes when its case ends, so a leg holds one at a time; a
  # later case that runs it pulls it again (DeepVariant's and PharmCAT's images
  # took 15 GB together in the base leg of run 41960912641).
  for i in ${E2E_LAZY_IMAGES:-}; do
    docker image inspect "$i" > /dev/null 2>&1 && { docker image rm -f "$i" > /dev/null || true; }
  done
  echo "::endgroup::"
  # A case whose network block could not be shown to hold proves nothing offline.
  [ "$rc" -eq 0 ] && $CASE_CONTROL_FAILED && rc=-1
  # Nor does a case that exited 0 without one [PASS] line and without saying
  # why it had nothing to check (NO_CHECKS_REASON, finish in lib.sh).
  [ "$rc" -eq 0 ] && ! grep -q '^\[PASS\] ' "$log" && ! grep -q "^${name}: no checks run: " "$log" && rc=-3
  case "$rc" in
    0)   result=pass ;;
    -1)  result="FAIL (network-block control)"; failed=$((failed + 1)) ;;
    -3)  result="FAIL (no check ran)"; failed=$((failed + 1)) ;;
    -2)  result="FAIL (data fill over ${E2E_DATA_BUDGET:-2400}s)"; failed=$((failed + 1)) ;;
    124) result="FAIL (timeout ${ctimeout}s)"; failed=$((failed + 1)) ;;
    *)   result="FAIL (exit ${rc})"; failed=$((failed + 1)) ;;
  esac
  secs=$(( $(date +%s) - start ))
  names+=("$name"); results+=("$result"); times+=("${secs}s")
  printf '%s\t%s\t%s\t%s\n' "$name" "$result" "$secs" "$JOB" >> "${REPORT_DIR}/results.tsv"
  [ "$TRACE" = true ] && trace_case "$name" "$tstart" "$secs"
  echo "${name}: ${result}"
done
if [ "${#names[@]}" -eq 0 ]; then
  echo "ERROR: no case matches '${PATTERN}'${CASE_LIST:+ or is listed in ${CASE_LIST}}" >&2
  exit 1
fi
# Tells the workflow every case ran, so the pulled images are complete enough to cache.
[ "$PATTERN" = "*" ] && touch "${E2E_WORK}/all-cases-ran"
kill "$LOW_PID" 2> /dev/null || true
echo "Disk (${JOB}): $(( $(cat "$LOW_FILE") / 1048576 )) GB free at the lowest during the cases, $(( $(df --output=avail -k "$E2E_WORK" | tail -1) / 1048576 )) GB after them" | tee -a "$SUMMARY"
docker system df

# --- 5. Report ----------------------------------------------------------------
{
  echo "### E2E on fixture ${TAG}: ${failed} of ${#names[@]} cases failed"
  echo
  echo "| Case | Result | Time | Log |"
  echo "|---|---|---|---|"
  for i in "${!names[@]}"; do
    echo "| ${names[$i]} | ${results[$i]} | ${times[$i]} | logs/${names[$i]}.log |"
  done
  echo
  echo "Logs are in the e2e-logs artifacts of this run."
  {
    for i in "${!names[@]}"; do
      [ "${results[$i]}" = pass ] && continue
      echo
      echo "#### ${names[$i]}"
      echo '```'
      grep -E '^\[FAIL\]' "${LOG_DIR}/${names[$i]}.log" || tail -n 15 "${LOG_DIR}/${names[$i]}.log"
      echo '```'
    done
    if [ -s "$E2E_NOTES" ]; then
      echo
      cat "$E2E_NOTES"
    fi
  } | tee "${REPORT_DIR}/extra.md"
} | tee -a "$SUMMARY"

# --- 6. Pack for the shards ---------------------------------------------------
# Everything the cases wrote or changed, at its absolute path, owners, modes
# and times (ns) kept, the base cases' logs too (a case may read one).
# Directories go in as entries alone (--no-recursion), so the reference and
# the fixture copies of the layout stay out.
if [ -n "${E2E_BASE_OUT:-}" ] && [ "$failed" -eq 0 ]; then
  echo "${E2E_WORK} ${REPO}" > "${E2E_WORK}/base-from"
  # Beside the archive too, so a shard learns the work area before it unpacks.
  cp "${E2E_WORK}/base-from" "${E2E_BASE_OUT}.from"
  sudo find "$E2E_WORK" "$REPO" \( -path "$FIXTURE_DIR" -o -path "$REPORT_DIR" \
      -o -path "${REPO}/.git" -o -path "$E2E_NOTES" -o -path "${E2E_WORK}/all-cases-ran" \) -prune \
      -o -cnewer "${E2E_WORK}/.layout-done" -print0 \
    | sudo tar -P --format=posix --no-recursion --null -T - -cf - | zstd -T0 -3 -q -o "$E2E_BASE_OUT"
  echo "Packed the base outputs: $(du -h "$E2E_BASE_OUT" | cut -f1), $(zstd -dc "$E2E_BASE_OUT" | tar -tf - | wc -l) entries. The largest:"
  zstd -dc "$E2E_BASE_OUT" | tar -tvf - | sort -k3,3n | tail -n 5
fi

[ "$failed" -eq 0 ]
