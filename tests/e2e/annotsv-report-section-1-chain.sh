#!/usr/bin/env bash
# The class 4 and 5 SVs reach both reports, with the real images of steps
# 05 (AnnotSV), 22 (SURVIVOR, bcftools) and 24 (the report's Python), on the
# made-up AnnotSV data and calls of tests/fixtures/annotsv/class45.sh, in a
# GENOME_DIR of their own:
#   - step 05's checked table: EXON (exon 1 end to the base before exon 3)
#     class 4 on exon edges with the benign-set AF 0.05; NEAR, 10 bp inside
#     the intron, not on exon edges;
#   - step 22's consensus holds DEL2, which Manta and Delly both call;
#   - the HTML report (step 24) and the text report (generate-report.sh):
#     DEL2 with 2 callers (manta, delly), DEL3 a single caller, EXON's flags;
#   - red control: with the consensus emptied, DEL2 is a single caller.
# The HG002 slice has no class 4 or 5 SV, and the full AnnotSV data (5 GB)
# is not on the runner: a real class 4 or 5 call is a real-data item.
. "$(dirname "$0")/lib.sh"
# shellcheck source=../fixtures/annotsv/class45.sh
. "${REPO}/tests/fixtures/annotsv/class45.sh"

G="${CASE_TMP}/genome"
T="${SAMPLE}c45"
rm -rf "$G"
class45_bundle "${G}/annotsv_annotations"
class45_calls "${G}/${T}" "$T"
# The reports write a run manifest when there is none, from the reference
# this GENOME_DIR does not have.
printf 'run\tscript\te2e\n' > "${G}/${T}/run_manifest.tsv"

GENOME_DIR="$G" run_step 05-annotsv.sh "$T"
check_step_exit 05-annotsv.sh
C="${G}/${T}/annotsv/${T}_sv_class4_5_checked.tsv"
cat "$C" 2>/dev/null
col() {
  awk -F'\t' -v id="$1" -v c="$2" 'NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
    $h["ID"] == id { print (c in h) ? $h[c] : "no column " c; exit }' "$C" 2>/dev/null
}
check_eq "EXON: class" "$(col EXON ACMG_class)" 4
check_eq "EXON: on exon edges" "$(col EXON Exon_edge_detail)" "NM_900001 (FIXG1): exon 1 and exon 3 edges"
check_eq "EXON: benign-set AF" "$(col EXON Population_AF) $(col EXON Population_source) $(col EXON Common)" \
  "0.05 fixture_common_loss yes"
check_eq "NEAR: not on exon edges" "$(col NEAR Exon_edge)" no

GENOME_DIR="$G" run_step 22-survivor-merge.sh "$T"
check_step_exit 22-survivor-merge.sh
CONS="${G}/${T}/sv_merged/${T}_sv_consensus.vcf.gz"
check_eq "consensus records" "$(gzip -cd "$CONS" 2>/dev/null | grep -vc '^#' || true)" 1

GENOME_DIR="$G" run_step 24-html-report.sh "$T"
check_step_exit 24-html-report.sh
HTML=$(cat "${G}/${T}/${T}_report.html" 2>/dev/null)
check "HTML: the AnnotSV card" has 'Structural Variants, Class 4 and 5 \(AnnotSV\)' "$HTML"
check "HTML: DEL2 with its two callers" has '<td>DEL 20:200000[0-9]-2450000</td>.*<td>2 callers \(manta, delly\)</td>' "$HTML"
check "HTML: DEL3 a single caller" has '<td>DEL 20:300000[0-9]-3450000</td>.*<td>single caller \(manta\)</td>' "$HTML"

GENOME_DIR="$G" run_step generate-report.sh "$T"
check_step_exit generate-report.sh
TXT=$(cat "${G}/${T}/${T}_report.txt" 2>/dev/null)
sed -n '/Class 4 and 5 (AnnotSV)/,/^$/p' <<< "$TXT"
check "text: DEL2 with its two callers" has 'DEL  20:200000[0-9]-2450000 .*support: 2 callers \(manta, delly\)$' "$TXT"
check "text: DEL3 a single caller" has 'DEL  20:300000[0-9]-3450000 .*support: single caller \(manta\)$' "$TXT"
check "text: EXON's exon edges and frequency" \
  has 'on exon edges: NM_900001 \(FIXG1\): exon 1 and exon 3 edges' "$TXT"
check "text: the counts" has 'Class 4-5 SVs: 4 \(single caller: 3; deletions on exon edges: 1; common: 2\)' "$TXT"

# Red control: the consensus emptied to its header.
gzip -cd "$CONS" | grep '^#' | gzip -c > "${CASE_TMP}/empty.vcf.gz" && mv "${CASE_TMP}/empty.vcf.gz" "$CONS"
GENOME_DIR="$G" run_step generate-report.sh "$T"
check "red control: with an empty consensus DEL2 is a single caller" \
  has 'DEL  20:200000[0-9]-2450000 .*support: single caller \(manta\)$' "$(cat "${G}/${T}/${T}_report.txt" 2>/dev/null)"
rm -rf "$G"

finish
