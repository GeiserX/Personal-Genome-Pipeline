#!/usr/bin/env bash
# Three tables reach both reports, with the real images of step 24 (the
# report's Python) and of the Nextflow HTML_REPORT process:
#   - step 42's genotype rules (traits/<id>_genotype_rules.tsv): the TAS2R38
#     line, PAV/PAV, under the traits rows;
#   - step 05's rules: a class 4 call over half inside an IG locus is moved
#     to class 3, listed apart with the rule and the reason, and left out of
#     the class 4-5 count;
#   - the repeat catalog step 09 ran on (expansion_hunter/<id>_eh_catalog.json)
#     gives MUC1 no thresholds, so it is not graded, never normal; and step
#     9b's motif table of the other flagged motif loci
#     (expansion_hunter/<id>_str_motifs.tsv) gives SAMD12 its verdict.
# The tables are made with the repository's own writers (bin/traits_report.py,
# bin/annotsv_postfilter.py) on made-up genotypes and calls; the Stranger
# VCF, the catalog and the motif table are made up in their real layouts.
#   1. Step 24 and generate-report.sh on the sample: both reports carry the
#      three. Control: the same sample without the three tables (and step 05
#      run without the rule table) carries none of them.
#   2. HTML_REPORT, run for real through a two-line pipeline on the same
#      files: it places each table (a name it has no place for stops the
#      task) and the report carries the three. Control: the bare set.
#   3. main.nf passes the three tables on: a stub run of every tool
#      (-profile test_all) stages each one in every HTML_REPORT task that has
#      the step's main output. Control: a copy without the three joins
#      stages none.
# What this cannot show: a real carrier of a repeat expansion or a real
# class 4 or 5 SV; the HG002 slice has neither.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

G="${CASE_TMP}/genome"
T=RFU
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
check "the fixture VCF is there" test -s "$VCF"
rm -rf "$G"
mkdir -p "$G"

# --- the tables: DIR/full has the three, DIR/bare has none of them ------------------------
python3 - "$REPO" "$G" "$T" <<'PY'
import json, os, shutil, sys
repo, g, s = sys.argv[1:4]
sys.path.insert(0, os.path.join(repo, "bin"))
import annotsv_postfilter, traits_report  # noqa: E401


def put(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)


for kind in ("full", "bare"):
    d = os.path.join(g, kind, s)
    put(f"{d}/run_manifest.tsv", "run\tscript\te2e\n")
    # Step 42: TAS2R38 homozygous for PAV (G, G, C on the plus strand).
    geno = os.path.join(g, kind, "genotypes.tsv")
    put(geno, "chr7\t141972804\tT\tC\tC/C\t40\nchr7\t141972905\tG\t.\tG/G\t40\nchr7\t141973545\tC\tG\tG/G\t40\n")
    rows = traits_report.read_table(os.path.join(repo, "assets", "traits.tsv"))
    rules = traits_report.read_rules(os.path.join(repo, "assets", "genotype_rules.tsv"))
    recs = traits_report.read_genotypes(geno)
    lines = traits_report.rule_lines(rules, recs)
    os.makedirs(f"{d}/traits", exist_ok=True)
    traits_report.write_tsv(f"{d}/traits/{s}_traits.tsv", traits_report.report(rows, recs, rules, None, lines))
    if kind == "full":
        traits_report.write_tsv(f"{d}/traits/{s}_genotype_rules.tsv", lines, traits_report.RULE_OUT_COLUMNS)
    # Step 05: two class 4 deletions; the rule table puts the second over half inside an IG locus.
    head = ("AnnotSV_ID\tSV_chrom\tSV_start\tSV_end\tSV_length\tSV_type\tID\tINFO\tAnnotation_mode\tGene_name\t"
            "B_loss_source\tB_loss_AFmax\tACMG_class\tAnnotSV_ranking_score\n")
    calls = ("20_3000001_3450000_DEL_1\t20\t3000000\t3450000\t-450000\tDEL\tKEEP\tEND=3450000;SVTYPE=DEL;"
             "SVCALLER=manta\tfull\tFIXK1\t\t\t4\t0.90\n"
             "20_20000001_20450000_DEL_1\t20\t20000000\t20450000\t-450000\tDEL\tMOVED\tEND=20450000;SVTYPE=DEL;"
             "SVCALLER=manta\tfull\tFIXM1\t\t\t4\t0.90\n")
    put(f"{d}/annotsv/{s}_sv_class4_5.tsv", head + calls)
    args = ["--class45", f"{d}/annotsv/{s}_sv_class4_5.tsv", "--out", f"{d}/annotsv/{s}_sv_class4_5_checked.tsv"]
    if kind == "full":
        put(os.path.join(g, "rules.tsv"), "#chrom\tstart\tend\trule\taction\tname\tsource\n"
            "chr20\t19990000\t20500000\tig_tr_locus\tclass_3\tIGFIX\tfixture\n")
        args += ["--rules", os.path.join(g, "rules.tsv")]
    annotsv_postfilter.main(args)
    os.remove(f"{d}/annotsv/{s}_sv_class4_5.tsv")
    # Steps 09 and 9b: MUC1 without a status, SAMD12 flagged.
    vhead = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    eh = ("chr1\t155188000\t.\tC\t<STR25>\t.\tPASS\tEND=155188001;REPID=MUC1;VARID=MUC1\tGT:REPCN\t0/1:20/25\n"
          "chr8\t118366913\t.\tA\t<STR40>\t.\tPASS\tEND=118366918;REPID=SAMD12;VARID=SAMD12{st}\tGT:REPCN\t0/1:1/40\n")
    put(f"{d}/expansion_hunter/{s}_eh.vcf", vhead + eh.format(st=""))
    t = os.path.getmtime(f"{d}/expansion_hunter/{s}_eh.vcf") - 60
    os.utime(f"{d}/expansion_hunter/{s}_eh.vcf", (t, t))
    put(f"{d}/expansion_hunter/{s}_eh_stranger.vcf", vhead + eh.format(st=";STR_STATUS=pre_mutation"))
    if kind == "full":
        put(f"{d}/expansion_hunter/{s}_eh_catalog.json", json.dumps([
            {"LocusId": "MUC1", "NoThresholds": "the reference repeat is 1200 bases, longer than 300: "
                                                "short reads cannot size it"},
            {"LocusId": "SAMD12", "NormalMax": 1, "PathologicMin": 105}], indent=1))
        put(f"{d}/expansion_hunter/{s}_str_motifs.tsv",
            "locus\tchrom\tstart\tend\trepeat_reads\tdominant_motif\tdominant_class\tdominant_share\tmotif_shares"
            "\tverdict\tnote\nSAMD12\tchr8\t118366812\t118366918\t5\tTGAAA\tpathogenic\t1.0\tTGAAA=1.0"
            "\tpathogenic_motif\tmotif TGAAA (pathogenic) in 5 of 5 repeat reads: confirm with repeat-primed PCR\n")
PY
check "the full set has the three tables" test -s "${G}/full/${T}/traits/${T}_genotype_rules.tsv" \
  -a -s "${G}/full/${T}/expansion_hunter/${T}_eh_catalog.json" -a -s "${G}/full/${T}/expansion_hunter/${T}_str_motifs.tsv"
check "the bare set has none of them" test ! -e "${G}/bare/${T}/traits/${T}_genotype_rules.tsv" \
  -a ! -e "${G}/bare/${T}/expansion_hunter/${T}_eh_catalog.json" -a ! -e "${G}/bare/${T}/expansion_hunter/${T}_str_motifs.tsv"
check_eq "step 05's rule moved MOVED to class 3 in the full set" \
  "$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next } $h["ID"] == "MOVED" { print $h["AnnotSV_class"] ">" $h["ACMG_class"] " " $h["Rule"] }' \
     "${G}/full/${T}/annotsv/${T}_sv_class4_5_checked.tsv")" "4>3 ig_tr_locus"

# The three sections, as patterns on a report (grep -E).
RULES_HTML='<h3>Genotype rules \(several sites read together\)</h3>'
RULES_ROW='<td>G/G, G/G, C/C</td><td>PAV/PAV \(no phase needed\)</td>'
MOVED_HTML="<h3>Moved to class 3 by step 05's rules</h3>"
MOVED_ROW='<td>DEL 20:20000000-20450000</td><td>4</td><td>3</td><td>FIXM1</td><td>ig_tr_locus</td>'
NG_ROW='<td>MUC1</td><td>20/25</td><td>not graded</td>'
MOTIF='SAMD12 motif count: pathogenic_motif'
# html_checks LABEL FILE WANT: WANT 1 means each section is there, 0 that none is
html_checks() {
  local label=$1 h n
  h=$(cat "$2" 2>/dev/null)
  for p in "$RULES_HTML" "$RULES_ROW" "$MOVED_HTML" "$MOVED_ROW" "$NG_ROW" "$MOTIF"; do
    n=$(grep -Ec -- "$p" <<< "$h" || true)
    if [ "$3" = 1 ]; then check_ge "${label}: ${p}" "$n" 1; else check_eq "${label}: no ${p}" "$n" 0; fi
  done
  if [ "$3" = 1 ]; then
    check "${label}: the class 4-5 table no longer lists MOVED" lacks '<tr><td>DEL 20:20000000-20450000</td><td>3</td>' "$h"
  fi
}

# --- 1. the bash reports ----------------------------------------------------------------
for kind in full bare; do
  GENOME_DIR="${G}/${kind}" run_step 24-html-report.sh "$T"
  check_step_exit 24-html-report.sh
  GENOME_DIR="${G}/${kind}" run_step generate-report.sh "$T"
  check_step_exit generate-report.sh
done
html_checks "bash HTML" "${G}/full/${T}/${T}_report.html" 1
html_checks "control: bash HTML without the tables" "${G}/bare/${T}/${T}_report.html" 0
TXT=$(cat "${G}/full/${T}/${T}_report.txt" 2>/dev/null)
check "text: the TAS2R38 rule line" has '\(TAS2R38 rs713598, rs1726866, rs10246939\): PAV/PAV \(no phase needed\)\. ' "$TXT"
check "text: the class 4-5 count leaves MOVED out" has 'Class 4-5 SVs: 1 ' "$TXT"
check "text: MOVED listed apart with the rule" has 'DEL 20:20000000-20450000  FIXM1  AnnotSV class 4 -> 3  \(ig_tr_locus: ' "$TXT"
check "text: MUC1 not graded with the catalog's reason" \
  has 'Not graded \(no thresholds in the catalog\): 1 \(MUC1: the reference repeat is 1200 bases' "$TXT"
check "text: the SAMD12 motif count" has "$MOTIF" "$TXT"
TXT0=$(cat "${G}/bare/${T}/${T}_report.txt" 2>/dev/null)
check "control: text without the tables has no rule line, no moved call, no not-graded locus" \
  lacks 'PAV/PAV \(no phase needed\)|Moved to class 3|Not graded \(no thresholds' "$TXT0"
check "control: text without the rule table counts both calls" has 'Class 4-5 SVs: 2 ' "$TXT0"

# --- 2. HTML_REPORT on the same files ------------------------------------------------------
NP="${CASE_TMP}/project"
rm -rf "$NP"
mkdir -p "$NP"
cp -R "${REPO}/bin" "${NP}/bin"
report_run() {
  local kind=$1 out="${CASE_TMP}/nf-$1" files="" f
  rm -rf "$out"
  while IFS= read -r f; do files="${files}${files:+, }file('${f}')"; done < <(
    find "${G}/${kind}/${T}" -type f \( -path '*/traits/*' -o -path '*/annotsv/*' -o -path '*/expansion_hunter/*' \) | sort)
  cat > "${NP}/main.nf" <<NF
include { HTML_REPORT } from '${REPO}/modules/local/html_report/main'

workflow {
    HTML_REPORT(channel.of([[id: '${T}', sex: 'male'], file('${VCF}'), [${files}]]))
}
NF
  cat > "${NP}/nextflow.config" <<NFCONF
params.outdir = '${out}/results'
params.publish_dir_mode = 'copy'
docker.enabled = true
docker.runOptions = '--user $(id -u):$(id -g) --network none'
NFCONF
  echo "+ nextflow run (HTML_REPORT, ${kind})"
  (
    cd "$CASE_TMP" && nextflow run "${NP}/main.nf" -ansi-log false \
      -c "${REPO}/conf/containers.config" -work-dir "${CASE_TMP}/work-${kind}"
  ) > "${CASE_TMP}/nf-${kind}.log" 2>&1
  local rc=$?
  check_eq "HTML_REPORT (${kind}) exits 0" "$rc" 0
  [ "$rc" = 0 ] || grep -E 'ERROR|Error executing process|Command error' -A 8 "${CASE_TMP}/nf-${kind}.log" | head -40
}
report_run full
report_run bare
html_checks "Nextflow HTML" "${CASE_TMP}/nf-full/results/${T}/${T}_report.html" 1
html_checks "control: Nextflow HTML without the tables" "${CASE_TMP}/nf-bare/results/${T}/${T}_report.html" 0

# --- 3. main.nf stages the three tables in HTML_REPORT (stub run) ---------------------------
# PAIRS: each new table and the step output that says the step ran for the sample
PAIRS="genotype_rules.tsv:traits.tsv str_motifs.tsv:eh_stranger.vcf eh_catalog.json:eh.vcf"
stub_staged() {
  local d id p new old hits got
  for p in $PAIRS; do
    new=${p%%:*} old=${p##*:} hits=0 got=0
    while IFS= read -r d; do
      d=$(dirname "$d")
      id=$(sed -n "s/^### name: [\"'].*HTML_REPORT (\(.*\))[\"']\$/\1/p" "${d}/.command.run")
      [ -e "${d}/input/${id}_${old}" ] || continue
      hits=$((hits + 1))
      [ ! -e "${d}/input/${id}_${new}" ] || got=$((got + 1))
    done < <(grep -l "^### name: [\"'].*HTML_REPORT (" -r "$1/work" --include=.command.run 2>/dev/null)
    echo "$new $hits $got"
  done
}
stub_run() {
  local d="${CASE_TMP}/$1"
  echo "+ nextflow run main.nf -profile test_all -stub ($1)"
  (cd "$d" && nextflow run main.nf -profile test_all -stub -ansi-log false --outdir "${d}/out" -work-dir "${d}/work") \
    > "${d}.log" 2>&1
  check_eq "stub run ($1) exits 0" "$?" 0
}
for n in stub stub-nojoin; do
  rm -rf "${CASE_TMP:?}/${n}"
  mkdir -p "${CASE_TMP}/${n}"
  cp -R "${REPO}/main.nf" "${REPO}/nextflow.config" "${REPO}/nextflow_schema.json" "${REPO}/versions.env" \
    "${REPO}/conf" "${REPO}/modules" "${REPO}/workflows" "${REPO}/bin" "${REPO}/assets" "${REPO}/scripts" "${CASE_TMP}/${n}/"
  mkdir -p "${CASE_TMP}/${n}/tests/fixtures"  # test_all's consumer_sites
  cp -R "${REPO}/tests/fixtures/consumer_export" "${CASE_TMP}/${n}/tests/fixtures/"
done
JOINS='BAM_ANALYSIS.out.str_motifs.map|BAM_ANALYSIS.out.eh_catalog.map|ch_genotype_rules.map'
check_eq "main.nf joins the three tables into the report inputs" "$(grep -cE "$JOINS" "${CASE_TMP}/stub/main.nf")" 3
sed -i -E "/${JOINS}/d" "${CASE_TMP}/stub-nojoin/main.nf"
check_eq "control copy: the joins are gone" "$(grep -cE "$JOINS" "${CASE_TMP}/stub-nojoin/main.nf")" 0
stub_run stub
stub_run stub-nojoin
while read -r new hits got; do
  check_ge "stub run: HTML_REPORT tasks with the step output before ${new}" "$hits" 1
  check_eq "stub run: each of them gets ${new}" "$got" "$hits"
done < <(stub_staged "${CASE_TMP}/stub")
while read -r new hits got; do
  check_ge "control: HTML_REPORT tasks with the step output before ${new}" "$hits" 1
  check_eq "control: without the join none gets ${new}" "$got" 0
done < <(stub_staged "${CASE_TMP}/stub-nojoin")

rm -rf "$G" "$NP" "${CASE_TMP:?}/stub" "${CASE_TMP:?}/stub-nojoin" "${CASE_TMP:?}"/work-* "${CASE_TMP:?}"/nf-*
finish
