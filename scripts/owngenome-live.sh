#!/usr/bin/env bash
# owngenome-live.sh: start the live site of the runs under GENOME_DIR (bin/owngenome_live.py, docs/live-site.md)
# Usage: GENOME_DIR=/data ./scripts/owngenome-live.sh [sample] [extra owngenome_live.py flags]
#        GENOME_DIR=/data ./scripts/owngenome-live.sh hash-password
#
# With a sample, it checks the name and prints the address of that run's page.
#
# hash-password asks for a password and writes its hash (never the password) to
# OWNGENOME_LIVE_PASSWORD_FILE, mode 600. The site then asks for that password.
# Without a hash file the site listens on loopback only and asks for nothing.
#
# Settings, each from the environment:
#   OWNGENOME_LIVE_BIND           browser address (default 127.0.0.1; any other needs the hash file)
#   OWNGENOME_LIVE_PORT           browser port (default 8787)
#   OWNGENOME_LIVE_INGEST_PORT    port Nextflow and run-all.sh post to, on 127.0.0.1 (default 8788)
#   OWNGENOME_LIVE_STATE_DIR      where messages and the access log go (default GENOME_DIR/.owngenome-live)
#   OWNGENOME_LIVE_PASSWORD_FILE  the hash file (default OWNGENOME_LIVE_STATE_DIR/password)
#   OWNGENOME_LIVE_ALLOWED_HOSTS  space-separated host names the browser may use besides loopback
#   OWNGENOME_LIVE_VIEWER         "off" leaves the genome viewer out (default on)
#
# The genome viewer (bin/owngenome_live_viewer.py) shows each run's reads at
# /runs/<run>/viewer/, behind the same password: the reference is REF_FASTA, as
# in every step, and the first start downloads the pinned igv.js
# (IGV_JS_VERSION in versions.env) into GENOME_DIR/tools/, checked by sha256,
# as owngenome-view.sh does. Without the reference or igv.js the site starts
# without the viewer and says why.
# Then start a run with OWNGENOME_LIVE_URL=http://127.0.0.1:8788 ./scripts/run-all.sh <sample> <sex>.
# Needs python3 (3.8 or later) on the host; nothing is installed.
set -euo pipefail

GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
STATE=${OWNGENOME_LIVE_STATE_DIR:-${GENOME_DIR}/.owngenome-live}
PWF=${OWNGENOME_LIVE_PASSWORD_FILE:-${STATE}/password}
command -v python3 >/dev/null 2>&1 || { echo "ERROR: owngenome-live.sh needs python3 on the host" >&2; exit 2; }

if [ "${1:-}" = hash-password ]; then
  mkdir -p "$(dirname "$PWF")"
  h=$(python3 "${PGP_ROOT}/bin/owngenome_live.py" hash-password)
  # One command per line: under set -e a failure inside an && list would not stop the script.
  (umask 077; printf '%s\n' "$h" > "${PWF}.tmp")
  chmod 600 "${PWF}.tmp"
  mv "${PWF}.tmp" "$PWF"
  echo "Password hash written: ${PWF} (mode 600). Start the site with: GENOME_DIR=${GENOME_DIR} $0"
  exit 0
fi

OPEN=""
if [ $# -gt 0 ] && [[ "$1" != --* ]]; then validate_sample "$1"; OPEN=$1; shift; fi

args=(--root "$GENOME_DIR" --state-dir "$STATE" --bind "${OWNGENOME_LIVE_BIND:-127.0.0.1}" --port "${OWNGENOME_LIVE_PORT:-8787}"
  --ingest-port "${OWNGENOME_LIVE_INGEST_PORT:-8788}")
if [ -s "$PWF" ]; then args+=(--password-file "$PWF"); fi
for h in ${OWNGENOME_LIVE_ALLOWED_HOSTS:-}; do args+=(--allowed-host "$h"); done
IGV_JS="${GENOME_DIR}/tools/igv.js-${IGV_JS_VERSION}/igv.min.js"
if [ "${OWNGENOME_LIVE_VIEWER:-on}" = off ]; then
  :
elif [ ! -f "${REF_FASTA}.fai" ]; then
  echo "Note: no reference index ${REF_FASTA}.fai; the site starts without the genome viewer." >&2
else
  if [ ! -f "$IGV_JS" ]; then
    echo "Downloading igv.js ${IGV_JS_VERSION} for the genome viewer (once)..."
    fetch "https://cdn.jsdelivr.net/npm/igv@${IGV_JS_VERSION}/dist/igv.min.js" "$IGV_JS" sha256 "$IGV_JS_SHA256" \
      || echo "Note: igv.js could not be downloaded; the site starts without the genome viewer." >&2
  fi
  if [ -f "$IGV_JS" ]; then args+=(--viewer-fasta "$REF_FASTA" --viewer-igv-js "$IGV_JS"); fi
fi
scheme=http
for a in "$@"; do [ "$a" != --tls-cert ] && [[ "$a" != --tls-cert=* ]] || scheme=https; done
[ -z "$OPEN" ] || echo "The run's page: ${scheme}://${OWNGENOME_LIVE_BIND:-127.0.0.1}:${OWNGENOME_LIVE_PORT:-8787}/runs/${OPEN}/"
exec python3 "${PGP_ROOT}/bin/owngenome_live.py" "${args[@]}" "$@"
