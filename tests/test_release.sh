#!/usr/bin/env bash
# scripts/ci/check-release.sh on planted trees, then on the real one.
#   tag:      a tag that differs from manifest.version fails (the state main
#             was in before the release check: 0.8.2 declared, v1.0.0 pushed);
#             a minor tag without docs/releases/<tag>.md fails; with it, the
#             notes file is the release body; a patch tag without notes gets
#             the commit subjects; a missing tag fails.
#   checks:   a run still going, a failed or cancelled run, no ci-ok, or no
#             runs at all fail; the ignored release job itself does not count.
#   declared: a declared X.Y.0 without notes, notes with the wrong heading,
#             a wrong file name or no CHANGELOG link fail.
set -uo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
CHECK="${REPO}/scripts/ci/check-release.sh"
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
fails=0
ok()  { echo "PASS $*"; }
bad() { echo "FAIL $*"; fails=$((fails + 1)); }
# expect RC WHAT ARGS...: run the check, compare its exit code (0 or non-0).
expect() {
  local want=$1 what=$2 rc=0
  shift 2
  out=$(bash "$CHECK" "$@" 2>&1) || rc=$?
  if { [ "$want" = 0 ] && [ "$rc" = 0 ]; } || { [ "$want" != 0 ] && [ "$rc" != 0 ]; }; then
    ok "$what"
  else
    bad "$what (exit ${rc}): ${out}"
  fi
}

g() { git -C "$T" -c user.name=t -c user.email=t@example.invalid -c tag.gpgSign=false -c commit.gpgSign=false -c core.hooksPath=/dev/null "$@"; }
config() {
  cat > "${T}/nextflow.config" <<NF
params {
    cyrius_install = null   // tools/cyrius-<version>
    version        = 'not this one'
}
manifest {
    name            = 'OwnGenome'
    nextflowVersion = '!>=25.10.4'
    version         = '$1'
}
NF
}
notes() {
  mkdir -p "${T}/docs/releases"
  printf '# OwnGenome %s\n\nWhat you get.\n\n- one thing\n' "$1" > "${T}/docs/releases/v$1.md"
}
g init -q
config 0.8.2
g add -A && g commit -qm "feat: the last commit before the cut"
g tag v0.8.2
g commit -q --allow-empty -m "fix: a change after 0.8.2"
g tag v1.0.0

# --- tag ---
expect 1 "tag v1.0.0 on a commit that declares 0.8.2 is refused" --repo "$T" tag v1.0.0
grep -q "does not match the version" <<< "$out" || bad "the mismatch message names the cause: ${out}"
expect 1 "a tag that does not exist is refused" --repo "$T" tag v9.9.9
expect 1 "a tag not of the form vX.Y.Z is refused" --repo "$T" tag 1.0.0

config 1.0.0
g add -A && g commit -qm "chore: declare 1.0.0"
g tag -f v1.0.0 >/dev/null
expect 1 "minor tag v1.0.0 without docs/releases/v1.0.0.md is refused" --repo "$T" tag v1.0.0
grep -q "major or minor release" <<< "$out" || bad "the missing-notes message names the cause: ${out}"

notes 1.0.0
g add -A && g commit -qm "docs: the 1.0.0 notes"
g tag -f v1.0.0 >/dev/null
expect 0 "v1.0.0 with matching manifest.version and its notes passes" --repo "$T" tag v1.0.0 --notes-out "${T}/body.md"
if cmp -s "${T}/body.md" "${T}/docs/releases/v1.0.0.md"; then ok "the release body is the notes file"; else bad "the release body is not the notes file"; fi

: > "${T}/docs/releases/v1.0.0.md"; printf '# OwnGenome 1.0.0\n' > "${T}/docs/releases/v1.0.0.md"
g add -A && g commit -qm "docs: empty the notes"
g tag -f v1.0.0 >/dev/null
expect 1 "notes that hold only a heading are refused" --repo "$T" tag v1.0.0

config 1.0.1
g add -A && g commit -qm "fix: a patch after 1.0.0"
g tag v1.0.1
expect 0 "patch tag v1.0.1 without notes passes" --repo "$T" tag v1.0.1 --notes-out "${T}/body.md"
if grep -qx -- "- fix: a patch after 1.0.0" "${T}/body.md"; then ok "a patch release body lists the commits since the previous tag"
else bad "patch body: $(cat "${T}/body.md")"; fi
if grep -q "a change after 0.8.2" "${T}/body.md"; then bad "the patch body reaches past the previous tag"; else ok "the patch body stops at the previous tag"; fi

# --- checks ---
tab=$'\t'
runs() { printf '%s\n' "$@" > "${T}/runs.tsv"; }
runs "ci-ok${tab}completed${tab}success" "shellcheck${tab}completed${tab}success" "deploy${tab}completed${tab}skipped" \
     "create-release${tab}in_progress${tab}"
expect 0 "all checks green, the release job itself left out" checks "${T}/runs.tsv" --ignore create-release
expect 1 "the release job counts when it is not left out" checks "${T}/runs.tsv"
runs "ci-ok${tab}completed${tab}success" "nextflow (pinned)${tab}in_progress${tab}"
expect 1 "a check still running is refused" checks "${T}/runs.tsv"
runs "ci-ok${tab}completed${tab}success" "unit-tests${tab}completed${tab}failure"
expect 1 "a failed check is refused" checks "${T}/runs.tsv"
runs "ci-ok${tab}completed${tab}success" "e2e${tab}completed${tab}cancelled"
expect 1 "a cancelled check is refused" checks "${T}/runs.tsv"
runs "shellcheck${tab}completed${tab}success"
expect 1 "no ci-ok among the checks is refused" checks "${T}/runs.tsv"
runs "ci-ok${tab}completed${tab}failure"
expect 1 "a failed ci-ok is refused" checks "${T}/runs.tsv"
: > "${T}/runs.tsv"
expect 1 "a commit with no check runs is refused" checks "${T}/runs.tsv"

# --- declared ---
D="${T}/declared"
mkdir -p "${D}/docs/releases"
cat > "${D}/nextflow.config" <<'NF'
manifest {
    version = '1.0.0'
}
NF
echo "# Changelog" > "${D}/CHANGELOG.md"
expect 1 "declared 1.0.0 without its notes is refused" --repo "$D" declared
printf '# OwnGenome 1.0.0\n\n- a thing\n' > "${D}/docs/releases/v1.0.0.md"
expect 1 "notes that CHANGELOG.md does not link are refused" --repo "$D" declared
echo "- 1.0.0: docs/releases/v1.0.0.md" >> "${D}/CHANGELOG.md"
expect 0 "declared 1.0.0 with linked notes passes" --repo "$D" declared
printf '# Release 1.0.0\n\n- a thing\n' > "${D}/docs/releases/v1.0.0.md"
expect 1 "notes whose heading is not '# OwnGenome 1.0.0' are refused" --repo "$D" declared
printf '# OwnGenome 1.0.0\n\n- a thing\n' > "${D}/docs/releases/v1.0.0.md"
printf '# OwnGenome 1.1\n\n- a thing\n' > "${D}/docs/releases/v1.1.md"
expect 1 "a notes file not named vX.Y.Z.md is refused" --repo "$D" declared
rm "${D}/docs/releases/v1.1.md"
sed -i.bak "s/'1.0.0'/'1.0.2'/" "${D}/nextflow.config"
expect 0 "a declared patch version needs no notes of its own" --repo "$D" declared
printf 'manifest {\n    name = "x"\n}\n' > "${D}/nextflow.config"
expect 1 "a manifest with no version is refused" --repo "$D" declared

# --- the real tree ---
expect 0 "the real tree: declared version and every notes file" declared
echo "$out"

[ "$fails" -eq 0 ] || { echo "${fails} release check case(s) failed"; exit 1; }
echo "All release check cases passed"
