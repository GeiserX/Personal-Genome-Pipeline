#!/usr/bin/env bash
# Step 08 run with KIR=true, then again without it, as after a realignment:
# the second run must not leave the first run's KIR genotype files, or the
# report (collect_summary.kir_genotypes reads kir_t1k/<sample>_kir_genotype.tsv
# whenever it exists) shows them as this run's. A run that stops early, here
# with the HLA database removed, must not leave them either.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
G=$GENOME_DIR
seed_reference "$G"
seed_sample "$G" sample1
use_output_hook
# t1k-build.pl writes an index with coordinates; run-t1k writes nothing, so
# the KIR pass records "too few reads", which is a KIR result the report shows.
cat > "${CASE_WORK}/t1k-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
args=" ${*:2} "
case "$args" in
  *" t1k-build.pl "*)
    o=$(host_path "$(sed -n 's/.* -o \([^ ]*\) .*/\1/p' <<<"$args")")
    p=$(sed -n 's/.* --prefix \([^ ]*\) .*/\1/p' <<<"$args")
    p=${p:-hla}
    mkdir -p "$o"
    printf '>%s*01:01\nACGT\n' "$p" > "${o}/${p}_dna_seq.fa"
    : > "${o}/${p}_dna_coord.fa"
    for g in HLA-A HLA-B HLA-C HLA-DRB1 HLA-DQB1 HLA-DPB1; do
      printf '>%s*01:01 chr6 29942532 29945870 +\nACGT\n' "$g" >> "${o}/${p}_dna_coord.fa"
    done
    exit 0 ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/t1k-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/t1k-hook"

# What setup.sh installs: the HLA database, the GENCODE genes and IPD-KIR.
# shellcheck source=../../scripts/lib/common.sh
HLA_DAT=$(. "${REPO_ROOT}/scripts/lib/common.sh"; data_file hla_dat || true)
# shellcheck source=../../scripts/lib/common.sh
GENES=$(. "${REPO_ROOT}/scripts/lib/common.sh"; data_file gencode_genes || true)
mkdir -p "$(dirname "$HLA_DAT")" "$(dirname "$GENES")"
printf 'ID   HLA00001; SV 4; standard; DNA; HUM; 3503 BP.\nCC   IPD-IMGT/HLA Release Version 3.0.0\n' > "$HLA_DAT"
printf 'chr6\tHAVANA\tgene\t29941260\t29949572\t.\t+\t.\tgene_name "HLA-A";\n' > "$GENES"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"
mkdir -p "${G}/kir/IPD-KIR_${KIR_DB_RELEASE}"
printf 'CC   IPD-KIR Release Version %s\n' "$KIR_DB_RELEASE" > "${G}/kir/IPD-KIR_${KIR_DB_RELEASE}/kir.dat"

K="${G}/sample1/kir_t1k"
# shown: what the report's summary reads for KIR (None = no KIR block).
shown() {
  python3 -c 'import sys; sys.path.insert(0, sys.argv[1]); import collect_summary as c; print(c.kir_genotypes(sys.argv[2], "sample1") is not None)' \
    "${REPO_ROOT}/bin" "${G}/sample1"
}

KIR=true run_expect 0 kir-on "${SCRIPTS}/08-hla-typing.sh" sample1
output_lacks kir-on 'unbound variable'
[ -s "${K}/sample1_kir_genotype.tsv" ] || fail "KIR=true wrote no KIR genotype file (the control)"
[ "$(shown)" = True ] || fail "the report's summary does not read the KIR result of the KIR=true run (the control)"

: > "${K}/sample1_kir_t1k_genotype.tsv"   # T1K's raw table, as a real run leaves it
run_expect 0 kir-off "${SCRIPTS}/08-hla-typing.sh" sample1
for f in sample1_kir_genotype.tsv sample1_kir_t1k_genotype.tsv; do
  [ ! -e "${K}/${f}" ] || fail "a run without KIR=true left the earlier ${f}"
done
[ "$(shown)" = False ] || fail "the report's summary still shows a KIR result after a run without KIR=true"

KIR=true run_expect 0 kir-on-again "${SCRIPTS}/08-hla-typing.sh" sample1
[ -s "${K}/sample1_kir_genotype.tsv" ] || fail "the second KIR=true run wrote no KIR genotype file"
rm -f "$HLA_DAT"
run_expect 0 kir-early-stop "${SCRIPTS}/08-hla-typing.sh" sample1
output_has kir-early-stop 'SKIPPED: IPD-IMGT/HLA'
[ ! -e "${K}/sample1_kir_genotype.tsv" ] || fail "a run that stopped early left the earlier KIR genotype file"
echo "Step 08 removes an earlier KIR result before its KIR gate, so a run without KIR=true shows none."
