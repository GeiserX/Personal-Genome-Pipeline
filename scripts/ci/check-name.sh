#!/usr/bin/env bash
# check-name.sh: fail when a tracked file still uses the project's old name.
#
# The project was called Personal-Genome-Pipeline until November 2026 and is
# now OwnGenome. This gate keeps the old name, and the old command and
# environment-variable names, from coming back through a branch that started
# before the rename. It checks:
#   name      any spelling of the old name (Personal-Genome-Pipeline,
#             Personal Genome Pipeline, personal_genome_pipeline, any case)
#   command   the old command and variable names: pgp-view.sh, pgp-lookup.sh,
#             pgp-live.sh, pgp_live.py, PGP_LIVE_*,
#             PGP_VIEW_PORT, PGP_DISK_WARN_PERCENT, PGP_MEM_WARN_PERCENT
# Lines that keep the old name on purpose are listed in ALLOW below, by file
# and exact line text. Internal identifiers that keep the PGP prefix (the
# PGP_* VCF tags, PGP_ROOT, localhost/pgp-pinned/, /pgp-* mounts, bead ids)
# are not checked.
#
# Usage:
#   scripts/ci/check-name.sh               every tracked file
#   scripts/ci/check-name.sh --self-test   plant one line of each kind in a
#                                          scratch copy and require a failure
set -euo pipefail

NAME_RE='personal[-_ ]genome[-_ ]pipeline'
CMD_RE='pgp-view\.sh|pgp-lookup\.sh|pgp-live\.sh|pgp_live|PGP_LIVE_|PGP_VIEW_PORT|PGP_DISK_WARN_PERCENT|PGP_MEM_WARN_PERCENT'

# file<TAB>exact line text
ALLOW=$(cat <<'EOF'
README.md	OwnGenome was called Personal-Genome-Pipeline until November 2026.
docs/index.md	OwnGenome was called Personal-Genome-Pipeline until November 2026.
EOF
)

# Ordinary words, not the name: "A personal genome pipeline focused on ...".
ALLOW_WORDS_FILE="docs/research/nf-core-pipeline-requirements.md"
ALLOW_WORDS_TEXT="A personal genome pipeline focused on"

scan() {
  # $1 = directory to scan (a git work tree). Prints file:line:text per hit.
  local dir=$1 hits=0 file line text
  while IFS= read -r hit; do
    file=${hit%%:*}; rest=${hit#*:}; line=${rest%%:*}; text=${rest#*:}
    if printf '%s\n' "$ALLOW" | grep -qxF "${file}	${text}"; then continue; fi
    if [ "$file" = "$ALLOW_WORDS_FILE" ] && [[ "$text" == *"$ALLOW_WORDS_TEXT"* ]]; then
      lower=$(printf '%s' "$text" | tr '[:upper:]' '[:lower:]')
      # only that phrase may carry the words; any other spelling on the line fails
      if [ "$(printf '%s\n' "$lower" | grep -oE "$NAME_RE" | wc -l)" -eq 1 ]; then continue; fi
    fi
    echo "OLD NAME: ${file}:${line}: ${text}"
    hits=$((hits + 1))
  done < <(git -C "$dir" grep -n -I -i -E "$NAME_RE" -- ':!scripts/ci/check-name.sh' || true)
  while IFS= read -r hit; do
    echo "OLD COMMAND: ${hit}"
    hits=$((hits + 1))
  done < <(git -C "$dir" grep -n -I -E "$CMD_RE" -- ':!scripts/ci/check-name.sh' || true)
  return $((hits > 0))
}

if [ "${1:-}" = "--self-test" ]; then
  root=$(git rev-parse --show-toplevel)
  scratch=$(mktemp -d)
  trap 'rm -rf "$scratch"' EXIT
  # a copy of the committed tree, as its own repository
  git -C "$root" archive HEAD | tar -x -C "$scratch"
  git -C "$scratch" init -q && git -C "$scratch" add -A
  fail=0
  if ! out=$(scan "$scratch"); then echo "self-test: the clean copy already fails:"; echo "$out"; exit 1; fi
  for plant in "See the Personal Genome Pipeline docs." "Run ./scripts/pgp-view.sh first."; do
    cp "$scratch/docs/index.md" "$scratch/docs/index.md.orig"
    echo "$plant" >> "$scratch/docs/index.md"
    git -C "$scratch" add docs/index.md
    if scan "$scratch" >/dev/null; then echo "self-test FAIL: not caught: $plant"; fail=1
    else echo "self-test ok: caught: $plant"; fi
    mv "$scratch/docs/index.md.orig" "$scratch/docs/index.md"
    git -C "$scratch" add docs/index.md
  done
  exit $fail
fi

if scan "$(git rev-parse --show-toplevel)"; then
  echo "OK: no old project name"
else
  echo "Rename these to OwnGenome (owngenome-view.sh, owngenome-lookup.sh, owngenome-live.sh, owngenome_live.py, OWNGENOME_*), or list a line that keeps the old name on purpose in ALLOW in scripts/ci/check-name.sh." >&2
  exit 1
fi
