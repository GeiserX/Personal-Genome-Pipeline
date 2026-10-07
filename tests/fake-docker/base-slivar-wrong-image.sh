#!/usr/bin/env bash
# Step 31 when the slivar image cannot run (a wrong image name): the step
# must exit non-zero, not report "no compound heterozygote candidates".
# The hook fails every run of the slivar image like a missing image does,
# answers `bcftools +split-vep -l` with a CSQ field list, and creates the
# files the other bcftools calls write.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
mkdir -p "${GENOME_DIR}/sample1/vep"
printf '##fileformat=VCFv4.2\n' | gzip -c > "${GENOME_DIR}/sample1/vep/sample1_vep.vcf.gz"
: > "${GENOME_DIR}/sample1/vep/sample1_vep.vcf.gz.tbi"
use_output_hook

cat > "${CASE_WORK}/slivar-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
case "$1" in
  *slivar*)
    echo "docker: Error response from daemon: manifest for $1 not found: manifest unknown (fake)" >&2
    exit 125 ;;
esac
case " $* " in
  *" +split-vep -l "*) printf 'Allele\nConsequence\nIMPACT\nSYMBOL\nCLIN_SIG\n' ;;
esac
exec "${CASE_WORK:?}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/slivar-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/slivar-hook"

run_rc slivar "${SCRIPTS}/31-slivar.sh" sample1
output_lacks slivar 'unbound variable'
docker_log_has '^run image=[^ ]*slivar' "step 31 never ran the slivar image"
[ "$RC" -ne 0 ] || fail "step 31 exited 0 although slivar could not run"
output_lacks slivar 'No compound heterozygote candidates found'
