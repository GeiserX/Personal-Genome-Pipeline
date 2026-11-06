#!/usr/bin/env bash
# check-release.sh: the checks a release must pass before it is published.
#
#   check-release.sh tag TAG [--notes-out FILE]
#       TAG (vMAJOR.MINOR.PATCH) must equal "v" + the manifest version that
#       nextflow.config declares at TAG, so `nextflow run -r TAG` reports the
#       version the release is named after. A major or minor release
#       (PATCH = 0) must carry its hand-written notes, docs/releases/TAG.md,
#       at TAG. With --notes-out the release body is written to FILE: the
#       notes file when TAG has one, else (a patch release only) the commit
#       subjects since the previous tag.
#
#   check-release.sh checks RUNS.TSV [--ignore NAME]...
#       RUNS.TSV is the tagged commit's check runs, one per line:
#       name<TAB>status<TAB>conclusion (as `gh api .../commits/SHA/check-runs`
#       gives them). Every run must have finished, and finished as success,
#       skipped or neutral; Guard's ci-ok, which runs on every push to main,
#       must be among them. NAME (the release job itself) is left out.
#
#   check-release.sh declared
#       The working tree: the version nextflow.config declares, and every
#       docs/releases/*.md. Each notes file is named vX.Y.Z.md, opens with
#       "# OwnGenome X.Y.Z", is linked from CHANGELOG.md, and a declared
#       major or minor version has one. Guard runs this through
#       tests/test_release.sh, so a bump to X.Y.0 cannot merge without its
#       notes.
#
# --repo DIR runs against another git work tree (the tests plant one); the
# default is the checkout this script is in.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
if [ "${1:-}" = --repo ]; then ROOT=${2:?--repo needs a directory}; shift 2; fi
MODE=${1:-}
shift || true

die() { echo "ERROR: $*" >&2; exit 1; }
semver='^v([0-9]+)\.([0-9]+)\.([0-9]+)$'

# manifest_version: the version = '...' line inside the manifest { } block of
# nextflow.config on stdin. Exactly one, or an error.
manifest_version() {
  local v
  v=$(awk '
    /^[[:space:]]*manifest[[:space:]]*\{/ {inside = 1; next}
    inside && /^[[:space:]]*\}/ {inside = 0}
    inside && /^[[:space:]]*version[[:space:]]*=/ {
      s = $0; sub(/^[^=]*=[[:space:]]*/, "", s); sub(/[[:space:]]*(\/\/.*)?$/, "", s)
      gsub(/^["\047]|["\047]$/, "", s); print s
    }')
  [ -n "$v" ] || die "nextflow.config declares no version in its manifest block"
  [ "$(printf '%s\n' "$v" | wc -l)" -eq 1 ] || die "nextflow.config declares more than one manifest version"
  printf '%s\n' "$v"
}

# check_notes NAME TEXT: NAME is vX.Y.Z.md; TEXT is its content.
check_notes() {
  local name=$1 text=$2 tag v
  tag=${name%.md}
  [[ "$name" == *.md && "$tag" =~ $semver ]] || die "docs/releases/${name}: a notes file is named vMAJOR.MINOR.PATCH.md"
  v=${tag#v}
  [ "$(printf '%s\n' "$text" | sed -n '/[^[:space:]]/{p;q;}')" = "# OwnGenome ${v}" ] \
    || die "docs/releases/${name}: the first line must be '# OwnGenome ${v}'"
  [ "$(printf '%s\n' "$text" | grep -c '[^[:space:]]')" -ge 2 ] \
    || die "docs/releases/${name}: holds no notes, only a heading"
}

case "$MODE" in
  tag)
    TAG=${1:?Usage: $0 tag TAG [--notes-out FILE]}
    shift
    OUT=""
    if [ "${1:-}" = --notes-out ]; then OUT=${2:?--notes-out needs a file}; fi
    [[ "$TAG" =~ $semver ]] || die "the tag is not of the form vMAJOR.MINOR.PATCH"
    patch=${BASH_REMATCH[3]}
    git -C "$ROOT" rev-parse -q --verify "refs/tags/${TAG}" >/dev/null || die "tag ${TAG} does not exist in this repository"
    declared=$(git -C "$ROOT" show "${TAG}:nextflow.config" | manifest_version)
    [ "v${declared}" = "$TAG" ] \
      || die "tag ${TAG} does not match the version nextflow.config declares at that commit (${declared}): bump manifest.version in a pull request, then tag the merged commit"
    notes="docs/releases/${TAG}.md"
    if text=$(git -C "$ROOT" show "${TAG}:${notes}" 2>/dev/null); then
      check_notes "${TAG}.md" "$text"
      [ -z "$OUT" ] || printf '%s\n' "$text" > "$OUT"
      echo "OK: ${TAG} matches manifest.version, notes from ${notes}"
    elif [ "$patch" = 0 ]; then
      die "${TAG} is a major or minor release and has no ${notes} at that commit: write its notes first"
    else
      if [ -n "$OUT" ]; then
        prev=$(git -C "$ROOT" tag --sort=-version:refname | grep -E "$semver" | grep -vx "$TAG" \
          | while read -r t; do git -C "$ROOT" merge-base --is-ancestor "$t" "$TAG" 2>/dev/null && { echo "$t"; break; }; done || true)
        {
          echo "## What's Changed"
          echo ""
          git -C "$ROOT" log "${prev:+${prev}..}${TAG}" --pretty=format:"- %s" --no-merges
          echo ""
          echo ""
          echo "All image versions are listed in versions.env at this tag."
        } > "$OUT"
      fi
      echo "OK: ${TAG} matches manifest.version, a patch release: notes from the commit subjects"
    fi
    ;;
  checks)
    RUNS=${1:?Usage: $0 checks RUNS.TSV [--ignore NAME]...}
    shift
    ignore=()
    while [ $# -gt 0 ]; do
      [ "$1" = --ignore ] || die "unknown option: $1"
      ignore+=("${2:?--ignore needs a name}")
      shift 2
    done
    [ -f "$RUNS" ] || die "no check-run list at ${RUNS}"
    n=0 bad=0 cios=0
    while IFS=$'\t' read -r name status conclusion || [ -n "$name" ]; do
      [ -n "$name" ] || continue
      skip=0
      for i in "${ignore[@]+"${ignore[@]}"}"; do [ "$name" != "$i" ] || skip=1; done
      [ "$skip" = 0 ] || continue
      n=$((n + 1))
      if [ "$status" != completed ]; then
        echo "NOT FINISHED: ${name} (${status})"; bad=1
      else
        case "$conclusion" in
          success|skipped|neutral) [ "$name" != ci-ok ] || [ "$conclusion" != success ] || cios=1 ;;
          *) echo "NOT GREEN: ${name} (${conclusion:-no conclusion})"; bad=1 ;;
        esac
      fi
    done < "$RUNS"
    [ "$n" -gt 0 ] || die "the tagged commit has no check runs: tag a commit that main's CI ran on"
    [ "$cios" = 1 ] || die "Guard's ci-ok did not succeed on the tagged commit: tag a commit that main's CI ran on"
    [ "$bad" = 0 ] || die "the tagged commit's checks are not all green: wait for them or fix main, then run the release again"
    echo "OK: all ${n} check runs of the tagged commit finished green"
    ;;
  declared)
    declared=$(manifest_version < "${ROOT}/nextflow.config")
    [[ "v${declared}" =~ $semver ]] || die "nextflow.config declares version '${declared}', which is not MAJOR.MINOR.PATCH"
    patch=${BASH_REMATCH[3]}
    changelog=""
    [ ! -f "${ROOT}/CHANGELOG.md" ] || changelog=$(cat "${ROOT}/CHANGELOG.md")
    shopt -s nullglob
    count=0
    for f in "${ROOT}"/docs/releases/*; do
      name=$(basename "$f")
      check_notes "$name" "$(cat "$f")"
      grep -qF "docs/releases/${name}" <<< "$changelog" || die "CHANGELOG.md does not link docs/releases/${name}"
      count=$((count + 1))
    done
    if [ "$patch" = 0 ] && [ ! -f "${ROOT}/docs/releases/v${declared}.md" ]; then
      die "nextflow.config declares ${declared}, a major or minor release, and docs/releases/v${declared}.md does not exist: write its notes in the same pull request"
    fi
    echo "OK: declared version ${declared}, ${count} notes file(s) well formed"
    ;;
  *)
    die "usage: $0 [--repo DIR] {tag TAG [--notes-out FILE] | checks RUNS.TSV [--ignore NAME]... | declared}"
    ;;
esac
