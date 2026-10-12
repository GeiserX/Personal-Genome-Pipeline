#!/usr/bin/env bash
# versions.env pins a few images by digest (name:tag@sha256:...). Docker gives
# an image pulled that way no tag, lists it as dangling, and a dangling-image
# cleanup (`docker image prune`) deletes it. setup.sh gives each one a second
# name, pinned_tag in scripts/lib/common.sh, and this test checks that name:
#   - every digest-pinned image of versions.env gets one, a tag-only image none;
#   - it starts with localhost/pgp-pinned/ (no public registry owns localhost)
#     and is a valid Docker reference;
#   - it is the same on every call and carries the digest: one tag per digest;
#   - tag_pinned_image is safe to repeat, and a failed `docker tag` is reported;
#   - setup.sh --pull-only says how many pinned tags it added (0 on a rerun);
#   - validate-setup.sh says a cleanup would delete an image only when the
#     image is dangling: another local tag on the same image id keeps it.
# No Docker: a fake engine on CONTAINER_ENGINE logs the calls, and a second one
# keeps a table of names and image ids. The real pull, tag and
# `docker images --filter dangling=true` run in Container Test.
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

# 2. The derivation itself. Image names use the example/ placeholder
# namespace, so scripts/ci/check-images.sh does not read them as pins.
P=localhost/pgp-pinned
check() {  # check IMAGE WANT
  local got
  got=$(pinned_tag "$1")
  if [ "$got" = "$2" ]; then pass "$1 -> ${2:-(no tag)}"; else fail "$1 -> '${got}', expected '${2}'"; fi
}
check "example/tool:3.14@sha256:${D1}"           "${P}/example/tool:sha256-${D1}"
check "example/tool@sha256:${D1}"                "${P}/example/tool:sha256-${D1}"
H=ghcr.io  # a registry host, kept apart so the line is not read as a pin
check "${H}/example/tool:2-beta@sha256:${D2}"    "${P}/${H}/example/tool:sha256-${D2}"
check "registry.example:5000/example/tool@sha256:${D1}"     "${P}/registry.example-5000/example/tool:sha256-${D1}"
check "registry.example:5000/example/tool:1.0@sha256:${D1}" "${P}/registry.example-5000/example/tool:sha256-${D1}"
check "example/tool:3.14"                        ""
check "example/tool:3.14@sha256:abc"             ""
check ""                                         ""
if [ "$(pinned_tag "example/tool@sha256:${D1}")" != "$(pinned_tag "example/tool@sha256:${D2}")" ]; then
  pass "two digests of one image get two tags"
else
  fail "two digests of one image share a tag"
fi
if [[ "$(pinned_tag "registry.example:5000/example/tool@sha256:${D1}")" =~ $REF_RE ]]; then
  pass "a registry port leaves a valid reference"
else
  fail "a registry port leaves an invalid reference"
fi

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
if [ "$got" -eq "$want" ] && [ "$want" -ge 1 ]; then
  pass "setup's list: ${got} tag call(s) for ${want} digest-pinned image(s)"
else
  fail "setup's list: ${got} tag call(s) for ${want} digest-pinned image(s)"
fi
if cmp -s "${WORK}/run1.log" "${WORK}/run2.log"; then pass "a second run makes the same calls"; else fail "a second run differs"; fi
while read -r _ src dst; do
  [ "$dst" = "$(pinned_tag "$src")" ] || fail "tagged ${src} as ${dst}"
done < "${WORK}/run1.log"

: > "$ENGINE_LOG"
tag_pinned_image "example/tool:3.14" || fail "tag_pinned_image failed on a tag-only image"
if [ -s "$ENGINE_LOG" ]; then fail "a tag-only image was tagged: $(cat "$ENGINE_LOG")"; else pass "a tag-only image is left alone"; fi

if TAG_FAILS=1 tag_pinned_image "example/tool@sha256:${D1}"; then
  fail "a failed docker tag was not reported"
else
  pass "a failed docker tag is reported"
fi

# 4. setup.sh --pull-only and validate-setup.sh against an engine that keeps
# state: one "NAME ID" line per local name in $ENGINE_STATE. A pull adds the
# digest reference, `tag` adds a name, `rmi` removes one; an image id is
# dangling when none of its names is a tag (every name has '@').
cat > "${WORK}/stateful" <<'FAKE'
#!/usr/bin/env bash
S=$ENGINE_STATE
id_of() { awk -v r="$1" '$1 == r { print $2; f = 1 } END { exit !f }' "$S"; }
put() { awk -v r="$1" '$1 != r' "$S" > "${S}.new"; echo "$1 $2" >> "${S}.new"; mv "${S}.new" "$S"; }
case "${1:-}" in
  --version) echo "Docker version 0.0.0-stateful" ;;
  info) case " $* " in *" --format "*) echo 17179869184 ;; *) echo ok ;; esac ;;
  pull) put "$2" "sha256:$(printf '%s' "$2" | cksum | cut -d' ' -f1)" ;;
  tag) id=$(id_of "$2") || exit 1; put "$3" "$id" ;;
  rmi) awk -v r="$2" '$1 != r' "$S" > "${S}.new"; mv "${S}.new" "$S" ;;
  image)
    [ "${2:-}" = inspect ] || exit 0
    shift 2
    fmt=""; [ "${1:-}" = --format ] && { fmt=$2; shift 2; }
    id=$(id_of "$1") || { echo "Error: No such image: $1" >&2; exit 1; }
    [ "$fmt" = '{{.Id}}' ] && echo "$id"
    exit 0 ;;
  images)
    case " $* " in
      *" dangling=true "*)
        awk '{ if (!($2 in t)) t[$2] = 0; if ($1 !~ /@/) t[$2] = 1 } END { for (i in t) if (!t[i]) print i }' "$S" ;;
    esac ;;
esac
exit 0
FAKE
chmod +x "${WORK}/stateful"
export ENGINE_STATE="${WORK}/state"
: > "$ENGINE_STATE"
STATEFUL="${WORK}/stateful"

npinned=$(pipeline_images | grep -c '@sha256:' || true)
first=$(pipeline_images | grep '@sha256:' | head -1)
first_tag=$(pinned_tag "$first")

pull_only() { CONTAINER_ENGINE=$STATEFUL bash "${REPO}/scripts/setup.sh" --pull-only > "${WORK}/setup.log" 2>&1; }
summary_says() {  # summary_says N: the summary line reports N pinned tags added
  if grep -qE "^\[OK\] Docker images: .*, $1 pgp-pinned tag\(s\) added\.$" "${WORK}/setup.log"; then
    pass "setup.sh --pull-only reports $1 pgp-pinned tag(s) added ($2)"
  else
    fail "setup.sh --pull-only should report $1 pgp-pinned tag(s) added ($2); it printed: $(grep '^\[OK\] Docker images' "${WORK}/setup.log" || tail -3 "${WORK}/setup.log")"
  fi
}
validate() { CONTAINER_ENGINE=$STATEFUL bash "${REPO}/scripts/validate-setup.sh" > "${WORK}/validate.log" 2>&1 || true; }
# warned: validate-setup.sh lists IMAGE under the "would delete" warning.
warned() {
  awk -v i="$1" '/digest-pinned image\(s\) have no local tag/ { w = 1; next } /Tag them with/ { w = 0 }
                 w && NF == 1 && $1 == i { f = 1 } END { exit !f }' "${WORK}/validate.log"
}
# noted: validate-setup.sh lists IMAGE as missing its pinned tag only.
noted() {
  awk -v i="$1" '/lack their localhost\/pgp-pinned\/ tag/ { w = 1; next } /Add it with/ { w = 0 }
                 w && NF == 1 && $1 == i { f = 1 } END { exit !f }' "${WORK}/validate.log"
}

pull_only || fail "setup.sh --pull-only failed: $(tail -5 "${WORK}/setup.log")"
summary_says "$npinned" "first run, every image pulled"
[ "$(awk -v r="$first_tag" '$1 == r' "$ENGINE_STATE")" ] || fail "setup.sh did not tag ${first} as ${first_tag}"
pull_only || fail "setup.sh --pull-only failed on a rerun"
summary_says 0 "rerun, every tag already there"

validate
if warned "$first" || noted "$first"; then fail "validate-setup.sh reports ${first} although it has its pinned tag"; else pass "no report for an image with its pinned tag"; fi

# The pinned tag removed: the image has no tag left, so it is dangling.
CONTAINER_ENGINE=$STATEFUL bash -c '"$CONTAINER_ENGINE" rmi "$1"' _ "$first_tag"
validate
if warned "$first"; then pass "an image with no tag at all gets the would-delete warning"; else fail "validate-setup.sh did not warn that a cleanup would delete ${first} (no tag at all): $(grep -A3 'Docker Images' "${WORK}/validate.log" | tr '\n' ' ')"; fi

# Another tag on the same image id: not dangling, so no deletion claim.
CONTAINER_ENGINE=$STATEFUL bash -c '"$CONTAINER_ENGINE" tag "$1" localhost/pgp-test/foreign:keep' _ "$first"
validate
if warned "$first"; then
  fail "validate-setup.sh claims a cleanup would delete ${first}, but another local tag keeps it"
else
  pass "an image kept by another tag gets no would-delete warning"
fi
if noted "$first"; then pass "an image kept by another tag is noted as missing its pinned tag"; else fail "validate-setup.sh did not note that ${first} lacks its pinned tag"; fi
pull_only || fail "setup.sh --pull-only failed after the tag was removed"
summary_says 1 "one pinned tag missing"

[ "$FAILS" -eq 0 ] || { echo "${FAILS} check(s) failed"; exit 1; }
echo "All pinned-tag checks passed."
