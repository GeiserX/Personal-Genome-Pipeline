#!/usr/bin/env bash
# setup.sh creates the reference's sequence dictionary (Picard and GATK need
# it; chip-to-vcf.sh now checks for it) and tries to install the small pinned
# data files. The fake downloads fail their checksums, which must not stop
# setup: each file gets a [WARN], and validate-setup.sh names what is missing
# and what the step does without it.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
use_output_hook
# CreateSequenceDictionary writes a SAM header; the generic hook would leave an
# empty file, which setup.sh rightly refuses.
cat > "${CASE_WORK}/dict-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
case " ${*:2} " in
  *" CreateSequenceDictionary "*)
    out=$(sed -n 's/.* -O \([^ ]*\) .*/\1/p' <<<" ${*:2} ")
    printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n' > "$(host_path "$out")"
    exit 0 ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/dict-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/dict-hook"

# chip-to-vcf.sh stops on a missing .dict before it converts anything.
mkdir -p "${GENOME_DIR}/reference_hg19" "${GENOME_DIR}/liftover" "${GENOME_DIR}/chip1/raw"
for f in reference_hg19/human_g1k_v37.fasta reference_hg19/human_g1k_v37.fasta.fai liftover/hg19ToHg38.over.chain.gz; do
  printf 'placeholder\n' > "${GENOME_DIR}/${f}"
done
printf 'rs1\t1\t100\tAG\n' > "${GENOME_DIR}/chip1/raw/chip1_raw.txt"
: > "$FAKE_DOCKER_LOG"
run_expect 1 chip-no-dict "${SCRIPTS}/chip-to-vcf.sh" chip1
output_has chip-no-dict 'Required file not found: .*Homo_sapiens_assembly38\.dict'
if grep -q '^run ' "$FAKE_DOCKER_LOG"; then fail "chip-to-vcf.sh started a container without the sequence dictionary"; fi

run_expect 0 setup "${SCRIPTS}/setup.sh" "$GENOME_DIR"
docker_log_has '^run image=[^ ]*gatk.* gatk CreateSequenceDictionary -R /genome/reference/Homo_sapiens_assembly38\.fasta ' \
  "setup.sh did not create the sequence dictionary"
[ -f "${GENOME_DIR}/reference/Homo_sapiens_assembly38.dict" ] || fail "setup.sh left no Homo_sapiens_assembly38.dict"
for what in 'Delly exclude map' 'GRCh38 chromosome bands' 'IPD-IMGT/HLA' 'GENCODE'; do
  output_has setup "\[WARN\] Could not install the ${what}"
done
docker_log_has '^curl https://raw\.githubusercontent\.com/dellytools/delly/[0-9a-f]{40}/excludeTemplates/human\.hg38\.excl\.tsv ' \
  "setup.sh did not fetch Delly's exclude map from a pinned commit"

run_rc validate "${SCRIPTS}/validate-setup.sh" sample1
output_has validate 'Delly exclude map \(step 19 runs without it'
output_has validate 'IPD-IMGT/HLA [0-9.]+ \(step 08 is skipped without it\) not found'
output_has validate 'GRCh38 chromosome bands \(step 10 falls back to hg19 bands\) not found'
