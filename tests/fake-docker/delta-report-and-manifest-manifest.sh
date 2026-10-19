#!/usr/bin/env bash
# bin/write_manifest.sh's ancestry_panel and hla_database rows.
#   ancestry_panel: "not installed (setup.sh --ancestry-panel)" without the
#     panel, "<name> (installed)" with the panel and its sites list,
#     "not used" with ANCESTRY_PANEL=none;
#   hla_database: the release in the hla.dat given as HLA_DAT, even outside
#     GENOME_DIR and with a T1K index folder of another release present;
#     without HLA_DAT the hla.dat setup.sh installs; then the release in the
#     index folder name, said to be from there; "unknown" when the sample has
#     HLA calls and nothing else is found; "not installed" otherwise.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
mkdir -p "${GENOME_DIR}/sample1"
M="${GENOME_DIR}/sample1/run_manifest.tsv"
PANEL_NAME=$(awk -F= '$1 == "PGSC_PANEL" {gsub(/"/, "", $2); print $2}' "${REPO_ROOT}/versions.env")
# shellcheck disable=SC2016  # the dollar signs are sed text, nothing expands
HLA_REL=$(sed -n 's/^HLA_DB_RELEASE=\${HLA_DB_RELEASE:-\([0-9.]*\)}.*/\1/p' "${REPO_ROOT}/scripts/lib/common.sh")
[ -n "$PANEL_NAME" ] && [ -n "$HLA_REL" ] || fail "could not read PGSC_PANEL (${PANEL_NAME}) or HLA_DB_RELEASE (${HLA_REL})"
row() { awk -F'\t' -v k="$1" '$1 == "data" && $2 == k {print $3}' "$M"; }
manifest() { run_expect 0 "$1" env "${@:2}" bash "${REPO_ROOT}/bin/write_manifest.sh" sample1 test; }
want() {  # want KEY TEXT WHY
  local got
  got=$(row "$1")
  [ "$got" = "$2" ] || fail "$3: ${1} reads '${got}', expected '${2}'"
}

# ancestry_panel
manifest no-panel
want ancestry_panel "not installed (setup.sh --ancestry-panel)" "no panel file"
mkdir -p "${GENOME_DIR}/reference/pgsc_calc"
printf 'x' > "${GENOME_DIR}/reference/pgsc_calc/${PANEL_NAME}.tar.zst"
manifest panel-no-sites
want ancestry_panel "not installed (setup.sh --ancestry-panel)" "the panel without its sites list"
printf 'x' > "${GENOME_DIR}/reference/pgsc_calc/${PANEL_NAME}_GRCh38_sites.tsv"
manifest panel
want ancestry_panel "${PANEL_NAME} (installed)" "the panel and its sites list"
manifest panel-none ANCESTRY_PANEL=none
want ancestry_panel "not used (ANCESTRY_PANEL=none)" "ANCESTRY_PANEL=none"

# hla_database: nothing, then calls only
manifest hla-none
want hla_database "not installed (setup.sh installs it)" "no hla.dat, no index, no calls"
mkdir -p "${GENOME_DIR}/sample1/hla_t1k"
printf 'HLA-A\t2\tA*01:01\t10\t60\t.\t0\t-1\n' > "${GENOME_DIR}/sample1/hla_t1k/sample1_hla_genotype.tsv"
manifest hla-calls-only
want hla_database "unknown (HLA calls exist, but no hla.dat or T1K index was found)" "HLA calls with no database found"

# a T1K index folder of release 3.60.0, and an hla.dat of 3.66.0 outside GENOME_DIR
mkdir -p "${GENOME_DIR}/t1k_idx/t1k-1.0.9_imgt-3.60.0_gencode-46"
OUTSIDE="${CASE_WORK}/elsewhere/hla.dat"
mkdir -p "$(dirname "$OUTSIDE")"
printf 'ID   HLA00001; standard; DNA; HUM; 3503 BP.\nDT   01/01/2026 (Rel. 3.65.0, Created)\nDT   01/04/2026 (Rel. 3.66.0, Last updated)\n' > "$OUTSIDE"
manifest hla-index-only
want hla_database "IPD-IMGT/HLA 3.60.0 (from the T1K index folder name; hla.dat not found)" "only the T1K index folder"
manifest hla-given HLA_DAT="$OUTSIDE"
want hla_database "IPD-IMGT/HLA 3.66.0" "HLA_DAT outside GENOME_DIR"

# the hla.dat setup.sh installs, with no DT line: its folder names the release
mkdir -p "${GENOME_DIR}/hla/IPD-IMGT-HLA_${HLA_REL}"
printf 'ID   HLA00001; standard; DNA; HUM; 3503 BP.\n' > "${GENOME_DIR}/hla/IPD-IMGT-HLA_${HLA_REL}/hla.dat"
manifest hla-installed
want hla_database "IPD-IMGT/HLA ${HLA_REL}" "the hla.dat setup.sh installs"
echo "ancestry_panel and hla_database say what is installed and which release the run used."
