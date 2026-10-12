#!/usr/bin/env bash
# versions.env pins a few images by digest (name:tag@sha256:...). Docker gives
# an image pulled that way no tag, lists it as dangling, and a dangling-image
# cleanup (`docker image prune`) deletes it. setup.sh gives each one a second
# name, pinned_tag in scripts/lib/common.sh, and this test checks that name:
#   - every digest-pinned image of versions.env gets one, a tag-only image none;
#   - it starts with localhost/pgp-pinned/ (no public registry owns localhost)
#     and is a valid Docker reference;
#   - it is the same on every call and carries the digest: one tag per digest;
#   - tag_pinned_image is safe to repeat, and a failed `docker tag` is reported.
# No Docker: a fake engine on CONTAINER_ENGINE logs the calls. The real pull,
# tag and `docker images --filter dangling=true` run in Container Test.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

cat > "${WORK}/engine" <<'FAKE'
#!/usr/bin/env bash
printf '%s\n' "$*" >> "$ENGINE_LOG"
[ "${1:-}" = tag ] && [ -n "${TAG_FAILS:-}" ] && exit 1
exit 0
FAKE
chmod +x "${WORK}/engine"
export ENGINE_LOG="${WORK}/engine.log"
: > "$ENGINE_LOG"

# shellcheck source=../scripts/lib/common.sh
. "${REPO}/scripts/lib/common.sh"
CONTAINER_ENGINE="${WORK}/engine"

D1=$(printf 'a%.0s' {1..64})
D2=$(printf 'b%.0s' {1..64})

# A Docker reference: lowercase path components joined by '/', each one
# alphanumerics with '.', '_', '__' or '-' between them, then a tag of up to
# 128 characters. The whole name stays under 256 characters.
COMP='[a-z0-9]+(([._]|__|-+)[a-z0-9]+)*'
REF_RE="^${COMP}(/${COMP})*:[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}$"

# 1. Every digest-pinned image of versions.env gets a tag; a tag-only one none.
pinned=0
while IFS= read -r img; do
  tag=$(pinned_tag "$img")
  case "$img" in
    *@sha256:*)
      pinned=$((pinned + 1))
      hex=${img##*@sha256:}
      [ -n "$tag" ] || { fail "no tag for digest-pinned ${img}"; continue; }
      [[ "$tag" == localhost/pgp-pinned/* ]] || fail "tag of ${img} is outside localhost/pgp-pinned/: ${tag}"
      [[ "$tag" == *":sha256-${hex}" ]] || fail "tag of ${img} does not carry its digest: ${tag}"
      [[ "$tag" =~ $REF_RE ]] && [ "${#tag}" -le 255 ] || fail "tag of ${img} is not a valid Docker reference: ${tag}"
      [ "$(pinned_tag "$img")" = "$tag" ] || fail "tag of ${img} changes from one call to the next"
      ;;
    *)
      [ -z "$tag" ] || fail "tag-only image ${img} got a pinned tag: ${tag}"
      ;;
  esac
done < <(grep -E '^[A-Z0-9_]+_IMAGE=' "${REPO}/versions.env" | sed -E 's/^[^=]+="([^"]*)".*/\1/')
[ "$pinned" -ge 1 ] || fail "found no digest-pinned image in versions.env; the test reads nothing"
pass "${pinned} digest-pinned image(s) of versions.env get a localhost/pgp-pinned/ tag; tag-only images get none"

# 2. The derivation itself.
check() {  # check IMAGE WANT
  local got
  got=$(pinned_tag "$1")
  if [ "$got" = "$2" ]; then pass "$1 -> ${2:-(no tag)}"; else fail "$1 -> '${got}', expected '${2}'"; fi
}
check "python:3.14@sha256:${D1}"                 "localhost/pgp-pinned/python:sha256-${D1}"
check "python@sha256:${D1}"                      "localhost/pgp-pinned/python:sha256-${D1}"
check "ghcr.io/org/tool:2-beta@sha256:${D2}"     "localhost/pgp-pinned/ghcr.io/org/tool:sha256-${D2}"
check "registry.example:5000/org/tool@sha256:${D1}"     "localhost/pgp-pinned/registry.example-5000/org/tool:sha256-${D1}"
check "registry.example:5000/org/tool:1.0@sha256:${D1}" "localhost/pgp-pinned/registry.example-5000/org/tool:sha256-${D1}"
check "python:3.14"                              ""
check "python:3.14@sha256:abc"                   ""
check ""                                         ""
[ "$(pinned_tag "python@sha256:${D1}")" != "$(pinned_tag "python@sha256:${D2}")" ] \
  && pass "two digests of one image get two tags" || fail "two digests of one image share a tag"
[[ "$(pinned_tag "registry.example:5000/org/tool@sha256:${D1}")" =~ $REF_RE ]] \
  && pass "a registry port leaves a valid reference" || fail "a registry port leaves an invalid reference"

# 3. tag_pinned_image: one `tag` call per digest-pinned image, the same on a
# second run; none for a tag-only image; a failed `docker tag` is reported.
for run in 1 2; do
  : > "$ENGINE_LOG"
  while IFS= read -r img; do
    tag_pinned_image "$img" || fail "run ${run}: tag_pinned_image ${img} failed"
  done < <(pipeline_images)
  cp "$ENGINE_LOG" "${WORK}/run${run}.log"
done
want=$(pipeline_images | grep -c '@sha256:' || true)
got=$(grep -c '^tag ' "${WORK}/run1.log" || true)
[ "$got" -eq "$want" ] && [ "$want" -ge 1 ] \
  && pass "setup's list: ${got} tag call(s) for ${want} digest-pinned image(s)" \
  || fail "setup's list: ${got} tag call(s) for ${want} digest-pinned image(s)"
if cmp -s "${WORK}/run1.log" "${WORK}/run2.log"; then pass "a second run makes the same calls"; else fail "a second run differs"; fi
while read -r _ src dst; do
  [ "$dst" = "$(pinned_tag "$src")" ] || fail "tagged ${src} as ${dst}"
done < "${WORK}/run1.log"

: > "$ENGINE_LOG"
tag_pinned_image "python:3.14" || fail "tag_pinned_image failed on a tag-only image"
[ ! -s "$ENGINE_LOG" ] && pass "a tag-only image is left alone" || fail "a tag-only image was tagged: $(cat "$ENGINE_LOG")"

if TAG_FAILS=1 tag_pinned_image "python@sha256:${D1}"; then
  fail "a failed docker tag was not reported"
else
  pass "a failed docker tag is reported"
fi

[ "$FAILS" -eq 0 ] || { echo "${FAILS} check(s) failed"; exit 1; }
echo "All pinned-tag checks passed."
