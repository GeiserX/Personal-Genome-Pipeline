#!/usr/bin/env bash
# pgp-live.sh: start the live site of the runs under GENOME_DIR (bin/pgp_live.py, docs/live-site.md)
# Usage: GENOME_DIR=/data ./scripts/pgp-live.sh [sample] [extra pgp_live.py flags]
#        GENOME_DIR=/data ./scripts/pgp-live.sh hash-password
#
# With a sample, it checks the name and prints the address of that run's page.
#
# hash-password asks for a password and writes its hash (never the password) to
# PGP_LIVE_PASSWORD_FILE, mode 600. The site then asks for that password.
# Without a hash file the site listens on loopback only and asks for nothing.
#
# Settings, each from the environment:
#   PGP_LIVE_BIND           browser address (default 127.0.0.1; any other needs the hash file)
#   PGP_LIVE_PORT           browser port (default 8787)
#   PGP_LIVE_INGEST_PORT    port Nextflow and run-all.sh post to, on 127.0.0.1 (default 8788)
#   PGP_LIVE_STATE_DIR      where messages and the access log go (default GENOME_DIR/.pgp-live)
#   PGP_LIVE_PASSWORD_FILE  the hash file (default PGP_LIVE_STATE_DIR/password)
#   PGP_LIVE_ALLOWED_HOSTS  space-separated host names the browser may use besides loopback
# Then start a run with PGP_LIVE_URL=http://127.0.0.1:8788 ./scripts/run-all.sh <sample> <sex>.
# Needs python3 (3.8 or later) on the host; nothing is installed.
set -euo pipefail

GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
STATE=${PGP_LIVE_STATE_DIR:-${GENOME_DIR}/.pgp-live}
PWF=${PGP_LIVE_PASSWORD_FILE:-${STATE}/password}
command -v python3 >/dev/null 2>&1 || { echo "ERROR: pgp-live.sh needs python3 on the host" >&2; exit 2; }

if [ "${1:-}" = hash-password ]; then
  mkdir -p "$(dirname "$PWF")"
  h=$(python3 "${PGP_ROOT}/bin/pgp_live.py" hash-password)
  (umask 077 && printf '%s\n' "$h" > "${PWF}.tmp") && chmod 600 "${PWF}.tmp" && mv "${PWF}.tmp" "$PWF"
  echo "Password hash written: ${PWF} (mode 600). Start the site with: GENOME_DIR=${GENOME_DIR} $0"
  exit 0
fi

OPEN=""
if [ $# -gt 0 ] && [[ "$1" != --* ]]; then validate_sample "$1"; OPEN=$1; shift; fi

args=(--root "$GENOME_DIR" --state-dir "$STATE" --bind "${PGP_LIVE_BIND:-127.0.0.1}" --port "${PGP_LIVE_PORT:-8787}"
  --ingest-port "${PGP_LIVE_INGEST_PORT:-8788}")
if [ -s "$PWF" ]; then args+=(--password-file "$PWF"); fi
for h in ${PGP_LIVE_ALLOWED_HOSTS:-}; do args+=(--allowed-host "$h"); done
[ -z "$OPEN" ] || echo "The run's page: http://${PGP_LIVE_BIND:-127.0.0.1}:${PGP_LIVE_PORT:-8787}/runs/${OPEN}/"
exec python3 "${PGP_ROOT}/bin/pgp_live.py" "${args[@]}" "$@"
