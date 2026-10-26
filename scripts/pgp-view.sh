#!/usr/bin/env bash
# pgp-view.sh — Open a sample's reads and variants in a local genome viewer
# Usage: ./scripts/pgp-view.sh <sample_name> [locus]
#   locus: chr1, chr1:12345 or chr1:12000-13000 (default: the first contig)
#
# Serves an igv.js page on http://127.0.0.1:PGP_VIEW_PORT (default 8765) with
# the reference (REF_FASTA), the sample's BAM or CRAM, its small-variant VCF
# and its SV VCF, and prints the address to open. The report's positions link
# to the same address, so with the viewer running a click on a ClinVar hit
# opens the reads at it. The server (assets/viewer/serve.py) listens on
# 127.0.0.1 only and its pages may not call any other server, so the data never
# leaves this machine; docs/viewer.md says how to reach it from a remote box
# over ssh. Needs python3 on the host. The first run downloads the pinned
# igv.js release (IGV_JS_VERSION below) into GENOME_DIR/tools/ and checks its
# sha256; later runs need no network.
#
# Files, each overridable from the environment:
#   READS   aligned/<sample>_sorted.bam, else aligned/<sample>_sorted.cram
#           (ALIGN_DIR as in the steps); "none" leaves the reads out
#   VCF     vcf/<sample>.vcf.gz (VCF_DIR as in the steps)
#   SV_VCF  sv_merged/<sample>_sv_consensus.vcf.gz (step 22), else Manta's
#           manta/results/variants/diploidSV.vcf.gz (step 4)
# A file that is not there is left out with a note; the reference is required.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> [locus]}
LOCUS=${2:-}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"

# igv.js, MIT licence, https://github.com/igvteam/igv.js. Bump both lines
# together: the sha256 is of dist/igv.min.js of that npm release.
IGV_JS_VERSION=3.8.10
IGV_JS_SHA256=87b991ae2a720978172e269dd5127c495e56df5af31144dff19c83124e2f5744
IGV_JS_URL="https://cdn.jsdelivr.net/npm/igv@${IGV_JS_VERSION}/dist/igv.min.js"

PORT=${PGP_VIEW_PORT:-8765}
S="${GENOME_DIR}/${SAMPLE}"
ALN="${S}/${ALIGN_DIR:-aligned}"

if [ ! -d "$S" ]; then
  echo "ERROR: Sample directory not found: ${S}" >&2
  exit 1
fi
if ! command -v python3 >/dev/null 2>&1; then
  echo "ERROR: python3 is needed on this computer to serve the viewer." >&2
  exit 1
fi
case "$PORT" in
  ''|*[!0-9]*) echo "ERROR: PGP_VIEW_PORT must be a port number, got '${PORT}'" >&2; exit 2 ;;
esac

first_file() {
  local f
  for f in "$@"; do
    if [ -f "$f" ]; then echo "$f"; return 0; fi
  done
  return 1
}

READS=${READS:-$(first_file "${ALN}/${SAMPLE}_sorted.bam" "${ALN}/${SAMPLE}_sorted.cram" || true)}
VCF=${VCF:-$(first_file "${S}/${VCF_DIR:-vcf}/${SAMPLE}.vcf.gz" || true)}
SV_VCF=${SV_VCF:-$(first_file "${S}/sv_merged/${SAMPLE}_sv_consensus.vcf.gz" \
  "${S}/manta/results/variants/diploidSV.vcf.gz" || true)}
[ "$READS" != none ] || READS=""

ARGS=(--sample "$SAMPLE" --fasta "$REF_FASTA" --port "$PORT")
[ -n "$LOCUS" ] && ARGS+=(--locus "$LOCUS")
for pair in "reads:${READS}" "vcf:${VCF}" "sv-vcf:${SV_VCF}"; do
  name=${pair%%:*} f=${pair#*:}
  if [ -z "$f" ]; then
    echo "Note: no ${name} file for ${SAMPLE}; the viewer shows it without that track."
  elif [ ! -f "$f" ]; then
    echo "ERROR: ${name} file not found: ${f}" >&2
    exit 1
  else
    ARGS+=("--${name}" "$f")
  fi
done

IGV_JS="${GENOME_DIR}/tools/igv.js-${IGV_JS_VERSION}/igv.min.js"
if [ ! -f "$IGV_JS" ]; then
  echo "Downloading igv.js ${IGV_JS_VERSION} (once)..."
  fetch "$IGV_JS_URL" "$IGV_JS" sha256 "$IGV_JS_SHA256"
fi

echo "============================================"
echo "  Genome viewer: ${SAMPLE}"
echo "  Reference: ${REF_FASTA}"
echo "  Reads:     ${READS:-none}"
echo "  VCF:       ${VCF:-none}"
echo "  SV VCF:    ${SV_VCF:-none}"
echo "============================================"
exec python3 "${PGP_ROOT}/assets/viewer/serve.py" "${ARGS[@]}" --igv-js "$IGV_JS"
