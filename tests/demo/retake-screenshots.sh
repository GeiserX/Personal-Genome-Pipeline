#!/usr/bin/env bash
# retake-screenshots.sh — rebuild DEMO-001 and retake the three docs pictures
# Usage: tests/demo/retake-screenshots.sh <work_dir> [<pictures_dir>]
#
#   1. tests/demo/make_demo_sample.py writes the invented sample DEMO-001 in
#      <work_dir>/genome (fixed seed; no real sample is read);
#   2. steps 36, 27, 24 and 28 run on it unchanged, with the images of
#      versions.env;
#   3. `make_demo_sample.py check` stops here if a report section is not ok,
#      so no picture shows a card that says "Not run";
#   4. tests/demo/render_pictures.py draws demo-html-report.png,
#      demo-cpic-report.png and demo-multiqc-report.png into <pictures_dir>
#      (default <work_dir>/pictures), with pinned playwright and pillow in a
#      venv (<work_dir>/venv, or $VENV);
#   5. tests/demo/check_pictures.py runs its self-test, then checks the
#      pictures: image chunks only and, on macOS, OCR finds DEMO-001 and no
#      private term.
#
# Env: DEMO_RECORDS     variant records of the VCF (default 4700000)
#      DEMO_TERMS_FILE  your private terms for the OCR check, one regular
#                       expression per line; keep it outside the repository
#      VENV             the venv to use (created when missing)
#      PLAYWRIGHT_WITH_DEPS  set (Linux) to let Playwright install the
#                       libraries Chromium needs, with apt through sudo
#
# Docker must see <work_dir> and this checkout. With Colima on macOS both
# must sit in a folder Colima mounts (the home folder and /private/tmp by
# default): a folder it does not mount looks empty inside the container. The
# script checks this before it starts.
#
# It writes nothing in docs/: look at the pictures, then copy them into
# docs/images/ (docs/testing.md, "Retaking the report pictures").
set -euo pipefail

WORK=${1:?Usage: $0 <work_dir> [<pictures_dir>]}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
SAMPLE="DEMO-001"
DEMO="${REPO}/tests/demo"
mkdir -p "$WORK"
WORK=$(cd "$WORK" && pwd)
OUT=${2:-${WORK}/pictures}
export GENOME_DIR="${WORK}/genome"
# shellcheck source=../../versions.env
. "${REPO}/versions.env"
ENGINE=${CONTAINER_ENGINE:-docker}

# Only a folder this script made is deleted.
if [ -e "$GENOME_DIR" ] && [ ! -f "${GENOME_DIR}/.demo-genome" ]; then
  echo "ERROR: ${GENOME_DIR} exists and was not written by this script; give another <work_dir>." >&2
  exit 1
fi
rm -rf "$GENOME_DIR"
mkdir -p "$GENOME_DIR"
: > "${GENOME_DIR}/.demo-genome"

echo "=== Docker sees ${GENOME_DIR} and ${REPO}/bin?"
if ! "$ENGINE" run --rm --network none -v "${GENOME_DIR}:/probe:ro" -v "${REPO}/bin:/pgp-bin:ro" "$PYTHON_IMAGE" \
    test -f /probe/.demo-genome -a -f /pgp-bin/collect_summary.py; then
  echo "ERROR: inside a container ${GENOME_DIR} or ${REPO}/bin is empty: Docker does not mount that folder." >&2
  echo "  Put the checkout and <work_dir> under a folder it mounts (Colima: the home folder or /private/tmp)." >&2
  exit 1
fi

echo "=== 1. DEMO-001"
python3 "${DEMO}/make_demo_sample.py" --genome-dir "$GENOME_DIR" --records "${DEMO_RECORDS:-4700000}"

echo "=== 2. Steps 36, 27, 24 and 28"
bash "${REPO}/scripts/36-pgx-consensus.sh" "$SAMPLE"
bash "${REPO}/scripts/27-cpic-lookup.sh" "$SAMPLE"
bash "${REPO}/scripts/24-html-report.sh" "$SAMPLE"
bash "${REPO}/scripts/28-multiqc.sh" "$SAMPLE"

echo "=== 3. Every section of the summary is ok?"
python3 "${DEMO}/make_demo_sample.py" check --genome-dir "$GENOME_DIR"

echo "=== 4. Pictures"
VENV=${VENV:-${WORK}/venv}
[ -x "${VENV}/bin/python" ] || python3 -m venv "$VENV"
PKGS=(playwright==1.63.0 pillow==12.3.0)
[ "$(uname -s)" != Darwin ] || PKGS+=(pyobjc-framework-Vision==12.2.2)
"${VENV}/bin/python" -m pip install -q "${PKGS[@]}"
export PLAYWRIGHT_BROWSERS_PATH=${PLAYWRIGHT_BROWSERS_PATH:-${VENV}/browsers}
# On a fresh Linux box Chromium also needs system libraries: --with-deps
# installs them with apt (through sudo).
if [ -n "${PLAYWRIGHT_WITH_DEPS:-}" ]; then
  "${VENV}/bin/python" -m playwright install --with-deps chromium
else
  "${VENV}/bin/python" -m playwright install chromium
fi
"${VENV}/bin/python" "${DEMO}/render_pictures.py" --genome-dir "$GENOME_DIR" --out "$OUT"

echo "=== 5. Checks"
ARGS=()
[ "$(uname -s)" != Darwin ] || ARGS+=(--ocr)
[ -z "${DEMO_TERMS_FILE:-}" ] || ARGS+=(--terms-file "$DEMO_TERMS_FILE")
"${VENV}/bin/python" "${DEMO}/check_pictures.py" --self-test ${ARGS[@]+"${ARGS[@]}"}
"${VENV}/bin/python" "${DEMO}/check_pictures.py" ${ARGS[@]+"${ARGS[@]}"} "${OUT}"/demo-*.png
[ "$(uname -s)" = Darwin ] || echo "NOTE: no OCR check here (it needs macOS); run check_pictures.py --ocr on a Mac before committing."

echo ""
echo "Pictures in ${OUT}:"
"${VENV}/bin/python" - "${OUT}"/demo-*.png <<'PY'
import sys
from PIL import Image
for p in sys.argv[1:]:
    with Image.open(p) as im:
        print(f"  {p}  {im.width} x {im.height}")
PY
echo "Look at them, then: cp ${OUT}/demo-*.png ${REPO}/docs/images/"
