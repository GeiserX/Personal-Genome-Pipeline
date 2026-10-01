#!/usr/bin/env bash
# 04-manta.sh run twice on the same sample. The hook below acts like Manta
# 1.6.0: configManta.py refuses a runDir that already holds runWorkflow.py,
# and runWorkflow.py writes results/variants/. The second run must exit 0
# and say the work is already done.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1

cat > "${CASE_WORK}/manta-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
shift   # the image
host() { printf '%s%s' "${GENOME_DIR:?}" "${1#/genome}"; }
case "${1##*/}" in
  configManta.py)
    run_dir=""
    while [ $# -gt 1 ]; do
      if [ "$1" = "--runDir" ]; then run_dir=$2; fi
      shift
    done
    d=$(host "$run_dir")
    if [ -e "${d}/runWorkflow.py" ]; then
      echo "fake Manta: run directory ${run_dir} already holds runWorkflow.py; Manta 1.6.0 refuses to configure it again" >&2
      exit 1
    fi
    mkdir -p "$d"
    printf '#!/usr/bin/env python\n' > "${d}/runWorkflow.py"
    ;;
  runWorkflow.py)
    d="$(dirname "$(host "$1")")/results/variants"
    mkdir -p "$d"
    for f in diploidSV candidateSV candidateSmallIndels; do
      printf '##fileformat=VCFv4.1\n' | gzip -c > "${d}/${f}.vcf.gz"
      : > "${d}/${f}.vcf.gz.tbi"
    done
    ;;
esac
HOOK
chmod +x "${CASE_WORK}/manta-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/manta-hook"

run_expect 0 manta-1 "${SCRIPTS}/04-manta.sh" sample1
[ -f "${GENOME_DIR}/sample1/manta/results/variants/diploidSV.vcf.gz.tbi" ] \
  || fail "the first run left no diploidSV.vcf.gz.tbi (the hook did not run)"
run_expect 0 manta-2 "${SCRIPTS}/04-manta.sh" sample1
output_has manta-2 '[Aa]lready done'
