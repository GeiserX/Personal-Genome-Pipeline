#!/usr/bin/env bash
# Step 32 (pypgx) takes its genes from assets/pgx_genes.tsv. This test checks:
#  1. the table against pypgx's own gene table at the pinned tag: every bam
#     row is a gene pypgx marks SV=TRUE, and every SV=TRUE gene is a bam row
#     (a planted SLC28A3 bam row and UGT2B17 set to vcf must both be reported);
#  2. scripts/32-pypgx.sh under a fake docker, with a fake pypgx that logs its
#     calls: the genes the step ran before the table (four bam, 19 vcf), given
#     as a table, give the same calls in the same order; each bam row added to
#     that table adds exactly one BAM-mode call; the shipped table gives one
#     BAM-mode call per bam row;
#  3. the not-callable line: from the callable_without_alt column, printed for
#     a BAM whose header lacks the gene's ALT contig and not for one that has
#     it, and not printed when the column says yes;
#  4. the PYPGX module reads the same table and spells no gene list;
#  5. the reports: every bam gene has exactly one row in the PGx copy number
#     block (its call, or why there is none), pypgx's, even beside a CNVpytor
#     section, and a not-callable gene gets its not-assessable line (two
#     planted summaries must be reported);
#  6. one owner per gene: no row of the depth engine's catalog
#     (bin/cn_loci.tsv) is a gene pypgx calls from the BAM (a planted UGT2B17
#     row must be reported).
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
CAT="${REPO}/assets/pgx_genes.tsv"
FAILS=0
bad() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
ok() { echo "ok: $*"; }

# The genes and modes of scripts/32-pypgx.sh before the table, as it spelled
# them (BAM_GENES, VCF_GENES).
OLD_BAM="CYP2D6 CYP2A6 GSTM1 GSTT1"
OLD_VCF="CYP1A2 CYP2B6 CYP2C9 CYP2C19 CYP3A4 CYP3A5 CYP4F2 DPYD TPMT NUDT15 UGT1A1 SLCO1B1 VKORC1 NAT2 COMT MTHFR ABCB1 G6PD IFNL3"

rows() { awk -F'\t' '!/^#/ && $1 != "gene" && NF' "$1"; }
row_of() { awk -F'\t' -v g="$1" '$1 == g' "$CAT"; }
# set_col GENE N VALUE: the shipped table with column N of GENE's row set to VALUE
set_col() { awk -F'\t' -v OFS='\t' -v g="$1" -v n="$2" -v v="$3" '$1 == g {$n = v} {print}' "$CAT"; }

# --- 1. the table against pypgx's gene table ---------------------------------
# shellcheck source=../versions.env
. "${REPO}/versions.env"
PYPGX_VERSION=${PYPGX_IMAGE##*:}
PYPGX_VERSION=${PYPGX_VERSION%%--*}
[ "$PYPGX_VERSION" = "$PYPGX_BUNDLE_VERSION" ] || bad "PYPGX_IMAGE is ${PYPGX_VERSION}, PYPGX_BUNDLE_VERSION ${PYPGX_BUNDLE_VERSION}"
GT="${WORK}/gene-table.csv"
curl -sS --fail --retry 4 --retry-delay 5 --max-time 60 -o "$GT" \
  "https://raw.githubusercontent.com/sbslee/pypgx/${PYPGX_VERSION}/pypgx/api/data/gene-table.csv"

# check_table TABLE: 0 when the bam rows are exactly pypgx's SV=TRUE genes and
# every row is a pypgx target gene
check_table() {
  python3 - "$GT" "$1" <<'PY'
import csv, sys
gt = {r["Gene"]: r for r in csv.DictReader(open(sys.argv[1]))}
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[2]) if l.strip() and not l.startswith("#")][1:]
sv = {g for g, r in gt.items() if r["SV"] == "TRUE"}
bam = {r[0] for r in rows if r[1] == "bam"}
errs = [f"{g} is not a pypgx target gene" for g in (r[0] for r in rows) if gt.get(g, {}).get("Target") != "TRUE"]
errs += [f"{g} is a bam row but pypgx has no copy-number model for it (SV=FALSE)" for g in sorted(bam - sv)]
errs += [f"{g} has a copy-number model (SV=TRUE) but is not a bam row" for g in sorted(sv - bam)]
print("\n".join(errs) if errs else f"{len(bam)} bam rows = pypgx's {len(sv)} SV genes")
sys.exit(1 if errs else 0)
PY
}
if out=$(check_table "$CAT"); then ok "$out"; else bad "the shipped table: $out"; fi
{ cat "$CAT"; printf 'SLC28A3\tbam\tyes\t.\tnone\tno\n'; } > "${WORK}/planted-slc28a3.tsv"
if out=$(check_table "${WORK}/planted-slc28a3.tsv"); then bad "a planted SLC28A3 bam row passed"; else ok "control, planted SLC28A3 bam row reported: $out"; fi
set_col UGT2B17 2 vcf > "${WORK}/planted-ugt2b17-vcf.tsv"
if out=$(check_table "${WORK}/planted-ugt2b17-vcf.tsv"); then bad "UGT2B17 as a vcf row passed"; else ok "control, UGT2B17 as vcf reported: $out"; fi

# --- 2 and 3. scripts/32-pypgx.sh under a fake docker -------------------------
G="${WORK}/genome" S=sample1
mkdir -p "${G}/${S}/aligned" "${G}/${S}/vcf" "${WORK}/bin" "${WORK}/pypgx-bin"
for f in "${S}/aligned/${S}_sorted.bam.bai" "${S}/vcf/${S}.vcf.gz" "${S}/vcf/${S}.vcf.gz.tbi"; do echo placeholder > "${G}/${f}"; done
B="${G}/reference/pypgx-bundle"
mkdir -p "$B"
git -C "$B" init -q
git -C "$B" -c core.hooksPath=/dev/null -c user.name=test -c user.email=test@example.invalid commit -q --allow-empty -m bundle
git -C "$B" tag "$PYPGX_BUNDLE_VERSION"
ln -s "${REPO}/scripts/ci/fake-docker/docker" "${WORK}/bin/docker"
# The bash running this test runs the step too (it needs bash 4.4 or later).
BASH_DIR=$(dirname "$BASH")

# A fake pypgx: logs each call.
cat > "${WORK}/pypgx-bin/pypgx" <<'SH'
#!/usr/bin/env bash
echo "$*" >> "$PYPGX_LOG"
SH
chmod +x "${WORK}/pypgx-bin/pypgx"
# The docker run hook: the contig check reads the @SQ lines of the BAM
# placeholder (the header fixture); the pypgx container's bash -c runs with
# the fake pypgx; the summary call logs its arguments and writes a header.
cat > "${WORK}/hook" <<'SH'
#!/usr/bin/env bash
img=$1; shift
host() { printf '%s' "${GENOME_DIR}${1#/genome}"; }
if [ "$1" = python3 ] && [ "$2" = -c ] && [[ "$3" == *pysam.AlignmentFile* ]]; then
  bam=$(host "$4"); shift 4; miss=()
  for c in "$@"; do awk -F'\t' -v c="SN:$c" '$1 == "@SQ" && $2 == c {f = 1} END {exit !f}' "$bam" || miss+=("$c"); done
  echo "${miss[*]}"
elif [ "$1" = bash ] && [ "$2" = -c ]; then
  PATH="${PYPGX_BIN}:$PATH" bash -c "$3"
elif [ "$1" = python3 ] && [ "$2" = /pgp-bin/pgx_outside_calls.py ]; then
  printf '%s\n' "$*" >> "$SUMMARY_LOG"
  while [ $# -gt 0 ]; do [ "$1" = --out ] && printf 'Gene\n' > "$(host "$2")"; shift; done
fi
exit 0
SH
chmod +x "${WORK}/hook"

# run32 TABLE HEADER NAME: run the step with TABLE and a BAM whose header is
# HEADER (alt or noalt); leaves NAME.out, NAME.pypgx and NAME.summary
run32() {
  local hdr=$'@HD\tVN:1.6\n@SQ\tSN:chr22\tLN:50818468'
  [ "$2" = noalt ] || hdr+=$'\n@SQ\tSN:chr22_KI270879v1_alt\tLN:304548'
  printf '%s\n' "$hdr" > "${G}/${S}/aligned/${S}_sorted.bam"
  : > "${WORK}/$3.pypgx"; : > "${WORK}/$3.summary"
  local rc=0
  env -i HOME="$WORK" PATH="${WORK}/bin:${BASH_DIR}:/usr/bin:/bin" GENOME_DIR="$G" PGX_GENES="$1" \
    FAKE_DOCKER_LOG="${WORK}/docker.log" FAKE_DOCKER_RUN_HOOK="${WORK}/hook" PYPGX_BIN="${WORK}/pypgx-bin" \
    PYPGX_LOG="${WORK}/$3.pypgx" SUMMARY_LOG="${WORK}/$3.summary" \
    "$BASH" "${REPO}/scripts/32-pypgx.sh" "$S" > "${WORK}/$3.out" 2>&1 || rc=$?
  RC=$rc
}
# bam_calls NAME, vcf_calls NAME: the genes of the BAM-mode and VCF-only
# run-ngs-pipeline calls, in order; prep NAME: the genes depth was prepared for
bam_calls() { awk '$1 == "run-ngs-pipeline" && /--depth-of-coverage/ {print $2}' "${WORK}/$1.pypgx" | paste -sd ' ' -; }
vcf_calls() { awk '$1 == "run-ngs-pipeline" && !/--depth-of-coverage/ {print $2}' "${WORK}/$1.pypgx" | paste -sd ' ' -; }
prep() { sed -n 's/^prepare-depth-of-coverage .* --genes //p' "${WORK}/$1.pypgx"; }
sum_arg() { sed -n "s/.* --$2 \([^-]*[^ -]\) --.*/\1/p" "${WORK}/$1.summary" | head -n 1; }
expect() { if [ "$2" = "$3" ]; then ok "$1"; else bad "$1: got '$2', expected '$3'"; fi; }

# The table of the genes the step ran before: their rows of the shipped table
# with the mode they had, in the order they had.
old_table() {
  grep '^#' "$CAT"; grep '^gene	' "$CAT"
  local g
  for g in $OLD_BAM; do row_of "$g" | awk -F'\t' -v OFS='\t' '{$2 = "bam"; print}'; done
  for g in $OLD_VCF; do row_of "$g" | awk -F'\t' -v OFS='\t' '{$2 = "vcf"; print}'; done
}
old_table > "${WORK}/old.tsv"
[ "$(rows "${WORK}/old.tsv" | wc -l)" -eq 23 ] || bad "the old genes are not all in the shipped table"

run32 "${WORK}/old.tsv" alt old
expect "old table: exit 0" "$RC" 0
expect "old table: BAM-mode calls as before" "$(bam_calls old)" "$OLD_BAM"
expect "old table: VCF-only calls as before" "$(vcf_calls old)" "$OLD_VCF"
expect "old table: depth prepared for the old bam genes" "$(prep old)" "$OLD_BAM"
expect "old table: summary --bam-genes as before" "$(sum_arg old bam-genes)" "$OLD_BAM"
expect "old table: summary --vcf-genes as before" "$(sum_arg old vcf-genes)" "$OLD_VCF"

# Each bam row the shipped table adds, alone on top of the old table.
for g in $(rows "$CAT" | awk -F'\t' '$2 == "bam" {print $1}'); do
  case " $OLD_BAM " in *" $g "*) continue ;; esac
  { awk -F'\t' -v g="$g" '$1 != g' "${WORK}/old.tsv"; row_of "$g"; } > "${WORK}/plus-${g}.tsv"
  run32 "${WORK}/plus-${g}.tsv" alt "plus-${g}"
  n_old=$(wc -w <<<"$OLD_BAM")
  expect "+${g}: exactly one more BAM-mode call, for ${g}" "$(bam_calls "plus-${g}")" "${OLD_BAM} ${g}"
  if vcf_calls "plus-${g}" | tr ' ' '\n' | grep -qxF "$g"; then bad "+${g}: also called from the VCF only"; fi
  [ "$(wc -w <<<"$(bam_calls "plus-${g}")")" -eq $((n_old + 1)) ] || bad "+${g}: not one more BAM-mode call"
done

# The shipped table: one BAM-mode call per bam row, one VCF-only call per vcf row.
WANT_BAM=$(rows "$CAT" | awk -F'\t' '$2 == "bam" {print $1}' | paste -sd ' ' -)
WANT_VCF=$(rows "$CAT" | awk -F'\t' '$2 == "vcf" {print $1}' | paste -sd ' ' -)
run32 "$CAT" alt shipped
expect "shipped table, ALT header: BAM-mode calls" "$(bam_calls shipped)" "$WANT_BAM"
expect "shipped table, ALT header: VCF-only calls" "$(vcf_calls shipped)" "$WANT_VCF"
expect "shipped table: 13 BAM-mode calls" "$(wc -w <<<"$(bam_calls shipped)" | tr -d ' ')" 13
if grep -q 'not callable on this reference' "${WORK}/shipped.out"; then bad "ALT header: a not-callable line"; else ok "ALT header: no not-callable line"; fi
if grep -q -- '--not-callable' "${WORK}/shipped.summary"; then bad "ALT header: --not-callable passed"; fi

run32 "$CAT" noalt noalt
NOALT_BAM=$(tr ' ' '\n' <<<"$WANT_BAM" | grep -vxF GSTT1 | paste -sd ' ' -)
expect "no-ALT header: exit 0" "$RC" 0
expect "no-ALT header: GSTT1 left out of the BAM-mode calls" "$(bam_calls noalt)" "$NOALT_BAM"
expect "no-ALT header: GSTT1 left out of depth preparation" "$(prep noalt)" "$NOALT_BAM"
expect "no-ALT header: the not-callable line" \
  "$(grep -c '^NOTICE: the BAM has no chr22_KI270879v1_alt contig (reference without ALT contigs); GSTT1 is not callable on this reference and is skipped$' "${WORK}/noalt.out")" 1
if grep -q -- '--not-callable GSTT1' "${WORK}/noalt.summary"; then ok "no-ALT header: --not-callable GSTT1 to the summary"; else bad "no-ALT header: no --not-callable GSTT1"; fi
expect "no-ALT header: the summary still lists GSTT1" "$(sum_arg noalt bam-genes)" "$WANT_BAM"

set_col GSTT1 3 yes > "${WORK}/gstt1-yes.tsv"
run32 "${WORK}/gstt1-yes.tsv" noalt gstt1-yes
if grep -q 'not callable on this reference' "${WORK}/gstt1-yes.out"; then bad "callable_without_alt=yes: a not-callable line"; else ok "control: with callable_without_alt=yes no not-callable line, so it comes from the column"; fi

{ cat "$CAT"; printf 'CYP2C9\tbam\tyes\t.\tnone\tno\n'; } > "${WORK}/dup.tsv"
run32 "${WORK}/dup.tsv" alt dup
expect "a gene twice stops the step" "$RC" 1
{ cat "$CAT"; printf 'ABCG2\tBAM\tyes\t.\tnone\tno\n'; } > "${WORK}/mode.tsv"
run32 "${WORK}/mode.tsv" alt mode
expect "an unknown mode stops the step" "$RC" 1

# --- 4. the module reads the same table ---------------------------------------
NF="${REPO}/modules/local/pypgx/main.nf"
# shellcheck disable=SC2016 # the ${projectDir} is matched, not expanded
grep -qF 'file("${projectDir}/assets/pgx_genes.tsv"' "$NF" || bad "the PYPGX module does not read assets/pgx_genes.tsv"
if grep -qE '"(CYP[0-9A-Z]+ ){2,}' "$NF"; then bad "the PYPGX module spells a gene list"; else ok "the PYPGX module spells no gene list"; fi

# --- 5. the reports ------------------------------------------------------------
# A summary with every gene of the table: GSTT1 not callable, CYP2E1 failed in
# pypgx, every other bam gene with a copy-number call. Every bam gene must have
# a row in the PGx copy number block (its call, or "not called: pypgx failed"),
# and GSTT1 its not-assessable line. Two planted summaries must be reported: one
# where GSTT1 reads as called (no not-assessable line), and one where a bam
# gene's row says VCF (it drops out of the block).
python3 - "$REPO" "$CAT" <<'PY' || FAILS=$((FAILS + 1))
import os, sys, tempfile
repo, cat = sys.argv[1], sys.argv[2]
sys.path.insert(0, os.path.join(repo, "bin"))
import collect_summary, render_report
rows = [l.rstrip("\n").split("\t") for l in open(cat) if l.strip() and not l.startswith("#")][1:]
bam = [g for g, mode, *_ in rows if mode == "bam"]
hdr = ["Gene", "Diplotype", "Phenotype", "CNV_call", "Source", "Status", "Defining_variants", "Min_GQ",
       "Weakest_variant", "Defining_filter", "Quality"]


def summary(gstt1_status="not callable on this reference", vcf_source_for=None):
    lines = []
    for g, mode, *_ in rows:
        src = "VCF" if g == vcf_source_for else mode.upper()
        if g == "GSTT1":
            lines.append([g, "N/A", "N/A", "N/A", src, gstt1_status] + ["N/A"] * 5)
        elif g == "CYP2E1":
            lines.append([g, "FAILED", "N/A", "N/A", src, "failed"] + ["N/A"] * 5)
        else:
            lines.append([g, "*1/*1", "Normal Metabolizer", f"CNV-{g}" if mode == "bam" else "N/A", src, "called"]
                         + ["none", "N/A", "N/A", "N/A", "N/A"])
    return lines


def problems(lines):
    with tempfile.TemporaryDirectory() as d:
        os.makedirs(os.path.join(d, "pypgx"))
        with open(os.path.join(d, "pypgx", "s1_pypgx_summary.tsv"), "w") as f:
            f.write("\n".join("\t".join(r) for r in [hdr] + lines) + "\n")
        _, data = collect_summary.sec_pypgx(d, "s1")
    # The rows of the block (gene, text), with CNVpytor's section beside pypgx's:
    # pypgx owns every bam gene, so each has one row, and it is pypgx's.
    block = render_report.pgx_copy_number({"pypgx": {"state": "ok", "data": data},
                                           "cnvpytor": {"state": "ok", "data": {}}})
    cn = dict(block)
    text = "\n".join(render_report.pypgx_lines(data))
    bad = [f"{g}: {n} rows in the copy-number block" for g in sorted({g for g, _ in block})
           if (n := sum(x == g for x, _ in block)) > 1]
    for g in bam:
        want = ("not called: pypgx failed" if g == "CYP2E1" else f"CNV-{g}") + " (pypgx, read depth)"
        if g == "GSTT1":
            if "Not assessable on this reference: GSTT1" not in text:
                bad.append("GSTT1 has no not-assessable line")
        elif cn.get(g) != want:
            bad.append(f"{g}: copy-number block says {cn.get(g)!r}, not {want!r}")
    return bad


bad = problems(summary())
for b in bad:
    print("FAIL:", b)
if not bad:
    print(f"ok: the reports show the copy-number call of {len(bam) - 2} bam genes, CYP2E1's failed call "
          "and GSTT1's not-assessable line")
for name, planted in (("GSTT1 read as called", summary(gstt1_status="no call")),
                      ("SLC22A2's row says VCF", summary(vcf_source_for="SLC22A2"))):
    got = problems(planted)
    if got:
        print(f"ok: control, {name} reported: {'; '.join(got)}")
    else:
        print(f"FAIL: control, {name} passed")
        bad.append(name)
sys.exit(1 if bad else 0)
PY

# --- 6. one owner per gene ------------------------------------------------------
# A gene pypgx calls from the BAM (a bam row) is pypgx's: the depth engine's
# catalog has no row for it. A copy of the catalog with a UGT2B17 row planted
# back must be reported.
CN="${REPO}/bin/cn_loci.tsv"
cn_owner_clash() {   # prints the catalog loci that are bam rows of the table
  awk -F'\t' 'NR == FNR { if ($0 !~ /^#/ && $2 == "bam") bam[$1] = 1; next }
              $0 !~ /^#/ && $1 != "locus" && ($1 in bam) { print $1 }' "$CAT" "$1"
}
if [ ! -s "$CN" ]; then
  bad "bin/cn_loci.tsv is missing or empty"
else
  clash=$(cn_owner_clash "$CN")
  if [ -n "$clash" ]; then bad "bin/cn_loci.tsv has a row for a gene pypgx calls from the BAM: $clash"
  else ok "bin/cn_loci.tsv has no row for any of the $(awk -F'\t' '$2 == "bam"' "$CAT" | wc -l | tr -d ' ') bam genes"; fi
  { cat "$CN"; printf 'UGT2B17\tgene\tchr4\tUGT2B17:68516300-68627039\n'; } > "${WORK}/planted-cn-ugt2b17.tsv"
  clash=$(cn_owner_clash "${WORK}/planted-cn-ugt2b17.tsv")
  if [ "$clash" = "UGT2B17" ]; then ok "control, a planted UGT2B17 catalog row reported"
  else bad "control, a planted UGT2B17 catalog row passed (got '${clash}')"; fi
fi

if [ "$FAILS" -gt 0 ]; then
  echo "${FAILS} check(s) failed"
  exit 1
fi
echo "all checks passed"
