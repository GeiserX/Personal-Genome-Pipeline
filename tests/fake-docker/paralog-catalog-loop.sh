#!/usr/bin/env bash
# Step 35 runs Parascopy once per row of the paralog catalog
# (assets/paralog_loci.tsv, or PARALOG_LOCI), with the hook below playing
# samtools and Parascopy:
#   - a catalog of three rows gives exactly three `parascopy cn-using` calls,
#     in catalog order, each with the row's model file of the chosen
#     population and its own output directory, and one table per row;
#   - PARASCOPY_POPULATION=AFR takes the AFR model files;
#   - a row whose region has no read is not assessed: no Parascopy call for it,
#     a table that says why, and the other rows still run;
#   - a row whose model file is not installed, and a row whose Parascopy run
#     fails, are not assessed the same way, and the step still exits 0;
#   - the catalog in the repository runs its SMN1 row into the table the
#     reports read.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
use_output_hook   # writes host-path.sh, which the hook below reads
D="${GENOME_DIR}/reference/parascopy-${PARASCOPY_DATA_VERSION}"
mkdir -p "${D}/homology_table"
printf 'x\n' > "${D}/homology_table/GRCh38.bed.gz"
for pop in EUR AFR; do
  mkdir -p "${D}/models_GRCh38_1KGP/${pop}"
  for m in SMN1 GBA PMS2; do printf 'x\n' > "${D}/models_GRCh38_1KGP/${pop}/${m}.gz"; done
done

cat > "${CASE_WORK}/hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
image=$1; shift
arg() {  # arg FLAG: the value after FLAG in the call
  local prev=""
  for a in "$@"; do [ "$prev" = "$FLAG" ] && { printf '%s' "$a"; return; }; prev=$a; done
}
case "$*" in
  "samtools view -c"*)
    # FAKE_NO_READS: a region with no read
    region=${*: -1}
    if [ "$region" = "${FAKE_NO_READS:-}" ]; then echo 0; else echo 25; fi ;;
  "parascopy depth"*)
    FLAG=-o; mkdir -p "$(host_path "$(arg "$@")")" ;;
  "parascopy cn-using"*)
    model=$2
    case "$model" in *"/${FAKE_CN_FAIL:-none}") echo "fake parascopy: failed on ${model}" >&2; exit 1 ;; esac
    FLAG=-o; o=$(host_path "$(arg "$@")"); mkdir -p "$o"
    name=$(basename "$model" .gz)
    printf '#chrom\tstart\tend\tlocus\tsample\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\tinfo\thomologous_regions\nchr1\t100\t200\t%s\tsample1\tPASS\t4\t40\tPASS\t2,2\t30,30\t.\tchr1:300-400:+\n' \
      "$name" | gzip -c > "${o}/res.samples.bed.gz" ;;
  *) exit 1 ;;   # the HBA depth check: not assessed, which this case does not read
esac
HOOK
chmod +x "${CASE_WORK}/hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hook"

L="${CASE_WORK}/loci.tsv"
{
  printf '# a test catalog\n'
  printf 'model\tgenes\tregion\treading\tpscn_reliable\ttag\tcross_check\ttable\tsource\n'
  printf 'SMN1.gz\tSMN1,SMN2\tchr5:70951946-70951946\tx\tyes\tcarrier\tparaphase\tsmn\ts\n'
  printf 'GBA.gz\tGBA1,GBAP1\tchr1:155234452-155244699\tx\tno\tcarrier\tparaphase\tgba\ts\n'
  printf 'PMS2.gz\tPMS2,PMS2CL\tchr7:5973239-6009125\tx\tno\tSF\tparaphase\tpms2\ts\n'
} > "$L"
OUT="${GENOME_DIR}/sample1/paralogs"
M="/genome/reference/parascopy-${PARASCOPY_DATA_VERSION}/models_GRCh38_1KGP"
calls() { grep -c 'parascopy cn-using' "$FAKE_DOCKER_LOG" || true; }
call_models() { grep -o "parascopy cn-using [^ ]*" "$FAKE_DOCKER_LOG" | awk '{print $3}' | tr '\n' ' '; }
na() { awk -F'\t' 'NR > 1 && $5 ~ /^not assessed/ {print $5}' "$1" 2>/dev/null; }

# --- three rows, three calls ------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
PARALOG_LOCI="$L" run_expect 0 three "${SCRIPTS}/35-paralogs.sh" sample1
[ "$(calls)" = 3 ] || fail "three catalog rows gave $(calls) parascopy cn-using calls, want 3"
want="${M}/EUR/SMN1.gz ${M}/EUR/GBA.gz ${M}/EUR/PMS2.gz "
[ "$(call_models)" = "$want" ] || fail "cn-using models: '$(call_models)', want '${want}'"
for t in smn gba pms2; do
  grep -q "parascopy cn-using .* -o /genome/sample1/paralogs/sample1_parascopy\.part/${t} " "$FAKE_DOCKER_LOG" \
    || fail "no cn-using call wrote to sample1_parascopy.part/${t}"
  [ -s "${OUT}/sample1_${t}_copy_number.tsv" ] || fail "no table sample1_${t}_copy_number.tsv"
  [ -s "${OUT}/sample1_parascopy/${t}/res.samples.bed.gz" ] || fail "no Parascopy output under sample1_parascopy/${t}"
done
[ "$(awk -F'\t' 'NR == 2 {print $4 "/" $6}' "${OUT}/sample1_gba_copy_number.tsv")" = GBA/4 ] \
  || fail "the GBA table does not carry Parascopy's row: $(cat "${OUT}/sample1_gba_copy_number.tsv")"
[ "$(head -n 1 "${OUT}/sample1_smn_copy_number.tsv")" = "$(printf 'chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions')" ] \
  || fail "the SMN table header changed: $(head -n 1 "${OUT}/sample1_smn_copy_number.tsv")"

# --- another population -----------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
PARALOG_LOCI="$L" PARASCOPY_POPULATION=AFR run_expect 0 afr "${SCRIPTS}/35-paralogs.sh" sample1
want="${M}/AFR/SMN1.gz ${M}/AFR/GBA.gz ${M}/AFR/PMS2.gz "
[ "$(call_models)" = "$want" ] || fail "AFR cn-using models: '$(call_models)', want '${want}'"

# --- a region without reads: not assessed, no call --------------------------------
: > "$FAKE_DOCKER_LOG"
PARALOG_LOCI="$L" FAKE_NO_READS=chr1:155234452-155244699 run_expect 0 no-reads "${SCRIPTS}/35-paralogs.sh" sample1
[ "$(calls)" = 2 ] || fail "with no read over GBA: $(calls) cn-using calls, want 2"
grep -q 'cn-using [^ ]*/GBA\.gz' "$FAKE_DOCKER_LOG" && fail "Parascopy ran on GBA without a read there"
[ "$(na "${OUT}/sample1_gba_copy_number.tsv")" = "not assessed: no read over chr1:155234452-155244699" ] \
  || fail "the GBA table does not say not assessed: $(cat "${OUT}/sample1_gba_copy_number.tsv")"
[ -z "$(na "${OUT}/sample1_pms2_copy_number.tsv")" ] || fail "PMS2 was not assessed too"
[ ! -e "${OUT}/sample1_parascopy/gba" ] || fail "a Parascopy directory from the earlier run is left for GBA"

# --- a model that is not installed, and a Parascopy run that fails -----------------
rm "${D}/models_GRCh38_1KGP/EUR/PMS2.gz"
: > "$FAKE_DOCKER_LOG"
PARALOG_LOCI="$L" FAKE_CN_FAIL=GBA.gz run_expect 0 missing-and-failed "${SCRIPTS}/35-paralogs.sh" sample1
[ "$(calls)" = 2 ] || fail "with PMS2's model missing: $(calls) cn-using calls, want 2 (SMN1, GBA)"
case "$(na "${OUT}/sample1_pms2_copy_number.tsv")" in
  "not assessed: no model PMS2.gz in the installed Parascopy data"*) ;;
  *) fail "the PMS2 table does not name the missing model: $(cat "${OUT}/sample1_pms2_copy_number.tsv")" ;;
esac
case "$(na "${OUT}/sample1_gba_copy_number.tsv")" in
  "not assessed: parascopy cn-using failed"*) ;;
  *) fail "the GBA table does not say Parascopy failed: $(cat "${OUT}/sample1_gba_copy_number.tsv")" ;;
esac
output_has missing-and-failed 'WARNING: GBA\.gz: not assessed'
[ "$(awk -F'\t' 'NR == 2 {print $6}' "${OUT}/sample1_smn_copy_number.tsv")" = 4 ] || fail "the SMN row did not run"
printf 'x\n' > "${D}/models_GRCh38_1KGP/EUR/PMS2.gz"

# --- the catalog in the repository --------------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 repo-catalog "${SCRIPTS}/35-paralogs.sh" sample1
grep -q "parascopy cn-using ${M}/EUR/SMN1\.gz " "$FAKE_DOCKER_LOG" || fail "the repository's catalog did not run SMN1"
[ -s "${OUT}/sample1_smn_copy_number.tsv" ] || fail "the repository's catalog wrote no SMN table"

# --- a catalog with a bad row stops the step before any container ---------------------
printf 'GBA.gz\tGBA1\tchr1:1-2 x\tx\tno\tcarrier\tnone\tgba2\ts\n' >> "$L"
: > "$FAKE_DOCKER_LOG"
PARALOG_LOCI="$L" run_expect 1 bad-row "${SCRIPTS}/35-paralogs.sh" sample1
[ "$(calls)" = 0 ] || fail "a bad catalog row still ran Parascopy"
echo "paralog-catalog-loop: all checks passed"
