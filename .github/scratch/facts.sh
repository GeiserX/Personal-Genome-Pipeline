#!/usr/bin/env bash
# Image-level facts behind beads pgp-9ms.4, .5 and .6. Temporary.
# shellcheck source=lib.sh
. "$(dirname "$0")/lib.sh"
prelude
. "$NEW/versions.env"

echo "== bgzip/tabix in ${BCFTOOLS_IMAGE} =="
OUT=$(docker run --rm "$BCFTOOLS_IMAGE" sh -c 'for b in bcftools bgzip tabix; do p=$(command -v $b) && echo "$b: $p" || echo "$b: absent"; done' 2>&1)
echo "$OUT"
check "bcftools image has bcftools" "$(grep '^bcftools:' <<< "$OUT")" "on PATH" "$(grep -q '^bcftools: /' <<< "$OUT" && echo 1 || echo 0)"
check "bcftools image has no bgzip" "$(grep '^bgzip:' <<< "$OUT")" "absent" "$(grep -q '^bgzip: absent' <<< "$OUT" && echo 1 || echo 0)"
check "bcftools image has no tabix" "$(grep '^tabix:' <<< "$OUT")" "absent" "$(grep -q '^tabix: absent' <<< "$OUT" && echo 1 || echo 0)"

echo "== haplogrep3 image =="
docker image inspect "$HAPLOGREP3_IMAGE" >/dev/null 2>&1 || docker pull -q "$HAPLOGREP3_IMAGE" >/dev/null
echo "Entrypoint/Cmd: $(docker image inspect "$HAPLOGREP3_IMAGE" --format '{{json .Config.Entrypoint}} {{json .Config.Cmd}}')"
OUT=$(docker run --rm "$HAPLOGREP3_IMAGE" sh -c 'for b in haplogrep3 classify; do p=$(command -v $b) && echo "$b: $p" || echo "$b: absent"; done' 2>&1)
echo "$OUT"
check "haplogrep3 on PATH" "$(grep '^haplogrep3:' <<< "$OUT")" "on PATH" "$(grep -q '^haplogrep3: /' <<< "$OUT" && echo 1 || echo 0)"
check "classify is not a command" "$(grep '^classify:' <<< "$OUT")" "absent" "$(grep -q '^classify: absent' <<< "$OUT" && echo 1 || echo 0)"

echo "== slivar image =="
docker pull -q "$SLIVAR_IMAGE" >/dev/null
OUT=$(docker run --rm "$SLIVAR_IMAGE" slivar compound-hets --help 2>&1 || true)
echo "$OUT"
check "slivar compound-hets help" "$(grep -m1 -i 'allow-non-trios' <<< "$OUT" || echo '<none>')" "mentions --allow-non-trios" "$(grep -q 'allow-non-trios' <<< "$OUT" && echo 1 || echo 0)"

echo "== Cyrius 1.1.1 in ${PYTHON_IMAGE}: resolved versions and console scripts =="
echo "python image digest: $(docker image inspect "$PYTHON_IMAGE" --format '{{index .RepoDigests 0}}' 2>/dev/null || docker pull -q "$PYTHON_IMAGE")"
OUT=$(docker run --rm "$PYTHON_IMAGE" bash -c '
  pip install --no-cache-dir --disable-pip-version-check -q "cyrius==1.1.1" >/dev/null 2>&1 &&
  python --version && echo "--- pip freeze ---" && pip freeze --all --exclude pip --exclude setuptools --exclude wheel &&
  echo "--- console scripts ---" &&
  for b in cyrius star_caller; do p=$(command -v $b) && echo "$b: $p" || echo "$b: absent"; done' 2>&1)
echo "$OUT"
check "cyrius console script" "$(grep '^cyrius:' <<< "$OUT")" "on PATH" "$(grep -q '^cyrius: /' <<< "$OUT" && echo 1 || echo 0)"
check "no star_caller command" "$(grep '^star_caller:' <<< "$OUT")" "absent" "$(grep -q '^star_caller: absent' <<< "$OUT" && echo 1 || echo 0)"

finish
