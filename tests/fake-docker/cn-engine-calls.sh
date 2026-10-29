#!/usr/bin/env bash
# Steps 35 (the HBA1/HBA2 block), 41 (RhD) and 46 (LPA KIV-2) count copies
# through one engine, bin/copy_number.py, with a hook that plays Parascopy and
# samtools bedcov (made-up depths: 30x everywhere, the KIV-2 block at 20
# repeats) and runs the engine for real on the host:
#   - each step calls copy_number.py bed once and step 35 and 46 call
#     copy_number.py call once (step 41's report imports the engine), and no
#     step calls hba_copy_number.py or lpa_kiv2.py any more;
#   - the tables say what the depths give: four alpha-globin genes, RhD D+
#     with two copies, 20.0 KIV-2 repeats;
#   - a control segment with no read (FAKE_ZERO=flank_e19, the others at 30x)
#     fails step 46 with BROKEN_CONTROL, never a count.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n' \
  | gzip -c > "${GENOME_DIR}/sample1/vcf/sample1.vcf.gz"
P="${GENOME_DIR}/reference/parascopy-${PARASCOPY_DATA_VERSION}"
mkdir -p "${P}/homology_table" "${P}/models_GRCh38_1KGP/EUR"
printf 'placeholder\n' > "${P}/homology_table/GRCh38.bed.gz"
printf 'placeholder\n' > "${P}/models_GRCh38_1KGP/EUR/SMN1.gz"
use_output_hook   # writes host-path.sh, which the hook below reads

# fakecov.py BED [q1]: samtools bedcov over BED (any number of columns; the
# name is column 6 of a marker row, else column 4), one summed depth appended.
cat > "${CASE_WORK}/fakecov.py" <<'PY'
import os, sys
q1 = len(sys.argv) > 2
kiv = 30.0 * (float(os.environ.get("FAKE_KIV2", "20")) + 2) / 14
for line in open(sys.argv[1]):
    c = line.rstrip("\n").split("\t")
    if not line.strip() or line.startswith("#") or c[0] == "chrom":
        continue
    name = c[5] if len(c) >= 6 else c[3]
    mean = kiv if name.startswith("kiv2_") else 30.0
    if q1 and name.startswith("kiv2_"):
        mean *= 0.1
    elif q1 and name in ("HBA1", "HBA2"):
        mean *= 0.75
    if name == os.environ.get("FAKE_ZERO"):
        mean = 0.0
    print("\t".join(c + [str(int(mean * (int(c[2]) - int(c[1]))))]))
PY

cat > "${CASE_WORK}/hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
image=$1; shift
arg() {  # arg FLAG ARGS...: the value after FLAG
  local f=$1 prev=""; shift
  for a in "$@"; do [ "$prev" = "$f" ] && { printf '%s' "$a"; return; }; prev=$a; done
}
map() {  # a container path on the host; anything else as it is
  case "$1" in
    /pgp-bin/*) printf '%s' "${REPO_ROOT}/bin/${1#/pgp-bin/}" ;;
    /pgp-assets/*) printf '%s' "${REPO_ROOT}/assets/${1#/pgp-assets/}" ;;
    /genome*) host_path "$1" ;;
    *) printf '%s' "$1" ;;
  esac
}
cov() { python3 "${CASE_WORK}/fakecov.py" "$@"; }
case "$*" in
  "parascopy depth"*) : > "$(host_path "$(arg -o "$@")")" ;;
  "parascopy cn-using"*)
    o=$(host_path "$(arg -o "$@")"); mkdir -p "$o"
    printf 'chr5\t70925087\t70953015\tSMN1\t.\tPASS\t4\t40\tPASS\t2,2\t30\t.\tchr5:70049669-70077595\n' | gzip -c > "${o}/res.samples.bed.gz" ;;
  "python3 /pgp-bin/"*)
    args=(); for a in "$@"; do args+=("$(map "$a")"); done
    exec "${args[@]}" ;;
  "sh -c "*)
    script=${3//\/pgp-bin\//${REPO_ROOT}/bin/}
    args=(); for a in "${@:4}"; do args+=("$(map "$a")"); done
    exec sh -c "$script" "${args[@]}" ;;
  "samtools bedcov -Q 1 "*) cov "$(map "$5")" q1 ;;
  "samtools bedcov "*) cov "$(map "$3")" ;;
  "bash -euo pipefail -c "*)
    case "$5" in
      *targets.bed*)   # step 41: the marker sites, then the RHD regions
        o=$(map "$8")
        { grep -v '^#' "$(map "$7")" | tail -n +2; cat "${o}/rhd.bed"; } > "${o}/targets.bed"
        cov "${o}/targets.bed" > "$(map "${10}")" ;;
      *)               # step 46: all MAPQ, then MAPQ >= 1
        cov "$(map "$7")" > "$(map "$9")"
        cov "$(map "$7")" q1 > "$(map "${10}")" ;;
    esac ;;
esac
HOOK
chmod +x "${CASE_WORK}/hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hook"
S="${GENOME_DIR}/sample1"
# calls TEXT: docker calls whose command holds TEXT (the log quotes with %q: backslashes dropped)
calls() { awk -v p="$1" '{gsub(/\\/, "")} index($0, p) {n++} END {print n + 0}' "$FAKE_DOCKER_LOG"; }
val() { awk -F'\t' -v k="$2" '$1 == k {print $2; exit}' "$1"; }

# --- step 35 ---------------------------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 paralogs "${SCRIPTS}/35-paralogs.sh" sample1
[ "$(calls 'copy_number.py bed --locus HBA')" = 1 ] || fail "step 35: copy_number.py bed --locus HBA was not called once"
[ "$(calls 'copy_number.py call --locus HBA')" = 1 ] || fail "step 35: copy_number.py call --locus HBA was not called once"
[ "$(calls 'hba_copy_number.py')" = 0 ] || fail "step 35 still calls hba_copy_number.py"
H="${S}/paralogs/sample1_hba_copy_number.tsv"
[ "$(val "$H" alpha_globin_genes)" = 4 ] || fail "step 35: four alpha-globin genes expected: $(cat "$H")"
[ "$(val "$H" status)" = called ] || fail "step 35: HBA not called: $(cat "$H")"

# --- step 41 ---------------------------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 bloodtype "${SCRIPTS}/41-bloodtype.sh" sample1
[ "$(calls 'copy_number.py bed --locus RHD')" = 1 ] || fail "step 41: copy_number.py bed --locus RHD was not called once"
[ "$(awk -F'\t' '$4 == "RHD_GENE" {n++} END {print n + 0}' "${S}/bloodtype/sample1_bloodtype_depth.tsv")" = 1 ] \
  || fail "step 41: the depth table has no RHD_GENE row"
R=$(awk -F'\t' '$1 == "RhD"' "${S}/bloodtype/sample1_bloodtype.tsv")
case "$R" in *"D+"*"RHD/RHD (RHD copies by read depth 2.00)"*) ;; *) fail "step 41: RhD with two copies expected: ${R}" ;; esac

# --- step 46 ---------------------------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 lpa "${SCRIPTS}/46-lpa.sh" sample1
[ "$(calls 'copy_number.py bed --locus LPA_KIV2')" = 1 ] || fail "step 46: copy_number.py bed --locus LPA_KIV2 was not called once"
[ "$(calls 'copy_number.py call --locus LPA_KIV2')" = 1 ] || fail "step 46: copy_number.py call --locus LPA_KIV2 was not called once"
[ "$(calls 'lpa_kiv2.py')" = 0 ] || fail "step 46 still calls lpa_kiv2.py"
L="${S}/lpa/sample1_lpa_kiv2.tsv"
[ "$(val "$L" kiv2_copies)" = 20.0 ] || fail "step 46: 20.0 KIV-2 repeats expected: $(cat "$L")"

# A control exon with no read while the others have 30x: a named error, no table
FAKE_ZERO=flank_e19 run_expect 3 lpa-broken "${SCRIPTS}/46-lpa.sh" sample1
output_has lpa-broken 'BROKEN_CONTROL: control segment flank_e19 of LPA_KIV2 has no read'
[ ! -s "$L" ] || fail "step 46 wrote a table from a broken control: $(cat "$L")"
