#!/usr/bin/env bash
# The README says bash 4.0 or newer. Before 4.4, "${ARR[@]}" of an empty
# array is an unbound variable under `set -u`, so a step that expands an
# optional, empty argument list that way dies before its tool runs. This case
# runs the steps that build such lists (03d, 04a, 10, 17, 19) under bash 4.2
# in a container, against the fake docker of the fake-docker suite.
. "$(dirname "$0")/lib.sh"

BASH42_IMAGE="bash:4.2.53-alpine3.22"
B="${CASE_TMP}/b42"
rm -rf "$B"
mkdir -p "$B"
cat > "${B}/inner.sh" <<'INNER'
set -uo pipefail
echo "bash ${BASH_VERSION}"
case "$BASH_VERSION" in 4.2.*) ;; *) echo "FAIL: not bash 4.2"; exit 1 ;; esac
export REPO_ROOT=/repo CASE_WORK=/work FAKE_DOCKER_LOG=/work/docker.log
export PATH="/repo/scripts/ci/fake-docker:${PATH}" GENOME_DIR=/work/genome
G=$GENOME_DIR
mkdir -p "$G/reference" "$G/s1/aligned" "$G/vep_cache" "$G/pcgr_data"
for f in reference/GRCh38_no_alt_analysis_set.fasta reference/GRCh38_no_alt_analysis_set.fasta.fai \
         s1/aligned/s1_sorted.bam s1/aligned/s1_sorted.bam.bai s1/vcf/s1.vcf.gz; do
  mkdir -p "$(dirname "$G/$f")"; echo placeholder > "$G/$f"
done
# Both BWA indexes: step 04a then runs TIDDIT with no extra argument (and so
# did the version before, which looked for the BWA-MEM2 one).
for ext in amb ann bwt pac sa bwt.2bit.64; do echo placeholder > "$G/reference/GRCh38_no_alt_analysis_set.fasta.$ext"; done
. /repo/versions.env
mkdir -p "$G/vep_cache/homo_sapiens/${PCGR_VEP_CACHE_RELEASE}_GRCh38" "$G/pcgr_data/${PCGR_DATA_BUNDLE}/data"
echo species > "$G/vep_cache/homo_sapiens/${PCGR_VEP_CACHE_RELEASE}_GRCh38/info.txt"
# TIDDIT writes <prefix>.vcf; bcftools stats answers a count.
cat > /work/hook <<'HOOK'
#!/usr/bin/env bash
args=" ${*:2} "
case "$args" in
  *" tiddit "*)
    o=$(echo "$args" | sed -n 's/.* -o \([^ ]*\) .*/\1/p')
    mkdir -p "$(dirname "/work/genome${o#/genome}")"
    printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > "/work/genome${o#/genome}.vcf" ;;
  *" bcftools stats "*) printf 'SN\t0\tnumber of records:\t0\n' ;;
esac
exit 0
HOOK
chmod +x /work/hook
export FAKE_DOCKER_RUN_HOOK=/work/hook
rc_all=0
for step in 03d-octopus.sh 04a-tiddit.sh 10-telomere-hunter.sh 19-delly.sh 17-cpsr.sh; do
  rc=0
  CPSR_SECONDARY_FINDINGS=false bash "/repo/scripts/${step}" s1 > "/work/${step}.log" 2>&1 || rc=$?
  if [ "$rc" -ne 0 ] || grep -q 'unbound variable' "/work/${step}.log"; then
    echo "FAIL: ${step} under bash 4.2 (exit ${rc}):"
    tail -n 5 "/work/${step}.log"
    rc_all=1
  else
    echo "ok: ${step} under bash 4.2"
  fi
done
exit "$rc_all"
INNER

docker run --rm -v "${REPO}:/repo:ro" -v "${B}:/work" "$BASH42_IMAGE" bash /work/inner.sh 2>&1 | tee "${CASE_TMP}/inner.log"
check_eq "every step runs under bash 4.2 (inner exit)" "${PIPESTATUS[0]}" 0
check_eq "steps that ran under bash 4.2" "$(grep -c '^ok: ' "${CASE_TMP}/inner.log" || true)" 5
docker run --rm -v "${B}:/work" "$BASH42_IMAGE" rm -rf /work/genome /work/hook 2>/dev/null || true

finish
