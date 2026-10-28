#!/usr/bin/env bash
# Step 39 and the Nextflow COUPLE_CARRIER wiring with AnnotSV's benign loss
# set, on synthetic Delly calls of two partners, CSV_A (female) and CSV_B
# (male), in a GENOME_DIR of their own. Coordinates are the real spans of
# assets/carrier_genes.tsv; the exons are 100 bp at each gene's middle.
#   - CSV_A: a 30 bp deletion over an ASPA exon, a 5 kb deletion over HEXA
#     that matches a common deletion of the benign set (breakpoints 100 to 200
#     bp off), and a 5 kb deletion over SMPD1 with no common match. CSV_B: the
#     same HEXA deletion, homozygous. Every long deletion has RDCN below the
#     normal count, so read depth confirms each loss;
#   - with the benign set under GENOME_DIR/annotsv_annotations (where setup.sh
#     puts it), step 39 reports SMPD1 as CSV_A's only, and neither ASPA (too
#     short) nor HEXA (common); the report counts both and says the filter ran;
#   - red control: without the set, HEXA is shared (1 in 2 per pregnancy, the
#     false risk the filter removes), ASPA still is not reported, and the
#     report says common deletions were not filtered;
#   - a stub run of the pipeline with --pedigree, clinvar and delly stages the
#     benign set into COUPLE_CARRIER when --annotsv_annotations holds it, and
#     stages none when the directory lacks it (red control).
# Cost: about 3 minutes, most of it the two stub runs.
. "$(dirname "$0")/lib.sh"

C="${E2E_WORK}/genome-couple-common-sv"
GENES="${REPO}/assets/carrier_genes.tsv"
BEN_REL="annotsv_annotations/Annotations_Human/SVincludedInFt/BenignSV/GRCh38"
in_e2e() { docker run --rm -v "${E2E_WORK}:/w" "$BCFTOOLS_IMAGE" "$@"; }
in_e2e rm -rf /w/genome-couple-common-sv
mkdir -p "${C}/reference" "${C}/${BEN_REL}"
export GENOME_DIR="$C"

# mid GENE: chrom and the middle of the gene's span
mid() { awk -F'\t' -v g="$1" '$1 == g {print $2, int(($3 + $4) / 2)}' "$GENES"; }
read -r A_C A_M <<< "$(mid ASPA)"
read -r H_C H_M <<< "$(mid HEXA)"
read -r S_C S_M <<< "$(mid SMPD1)"
check "ASPA, HEXA and SMPD1 are in the carrier list" test -n "${A_M:-}" -a -n "${H_M:-}" -a -n "${S_M:-}"

printf '%s\t%s\t%s\n' "$A_C" $((A_M - 50)) $((A_M + 50)) "$H_C" $((H_M - 50)) $((H_M + 50)) \
  "$S_C" $((S_M - 50)) $((S_M + 50)) | sort -k1,1 -k2,2n > "${C}/reference/exons.bed"

# del CHROM POS END GT RDCN: a PASS Delly 2.x deletion record
del() { printf '%s\t%s\t.\tN\t<DEL>\t60\tPASS\tIMPRECISE;SVTYPE=DEL;END=%s\tGT:FT:RDCN\t%s:PASS:%s\n' "$@"; }
vcf_head() { printf '##fileformat=VCFv4.2\n##source=Synthetic_e2e_couple_common_sv\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$1"; }
for s in CSV_A CSV_B; do
  mkdir -p "${C}/${s}/clinvar" "${C}/${s}/delly"
  vcf_head "$s" > "${C}/${s}/clinvar/${s}_clinvar_hits.vcf"
done
{
  vcf_head CSV_A
  del "$A_C" $((A_M - 15)) $((A_M + 15)) 0/1 2
  del "$H_C" $((H_M - 2500)) $((H_M + 2500)) 0/1 1
  del "$S_C" $((S_M - 2500)) $((S_M + 2500)) 0/1 1
} | gzip -n > "${C}/CSV_A/delly/CSV_A_sv.vcf.gz"
{ vcf_head CSV_B; del "$H_C" $((H_M - 2500)) $((H_M + 2500)) 1/1 0; } | gzip -n > "${C}/CSV_B/delly/CSV_B_sv.vcf.gz"

# The benign loss set as AnnotSV writes it: chromosome without "chr", 0-based
# start, source, coordinates, frequency. The HEXA record is common (31%); a
# rare one (0.4%) sits exactly on the SMPD1 deletion.
BEN="${C}/${BEN_REL}/benign_Loss_SV_GRCh38.sorted.bed"
printf '%s\t%s\t%s\tgnomAD-SV_e2e_%s\t%s\t%s\n' \
  "${H_C#chr}" $((H_M - 2400)) $((H_M + 2700)) DEL1 "${H_C#chr}:$((H_M - 2400))-$((H_M + 2700))" 0.3100 \
  "${S_C#chr}" $((S_M - 2500)) $((S_M + 2500)) DEL2 "${S_C#chr}:$((S_M - 2500))-$((S_M + 2500))" 0.0040 > "$BEN"

printf 'sample\tsex\tfather\tmother\tpartner\nCSV_A\tfemale\t.\t.\tCSV_B\nCSV_B\tmale\t.\t.\tCSV_A\n' > "${C}/pedigree.tsv"

col() { awk -F'\t' -v g="$2" -v k="$3" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $1 == g {print $c[k]; exit}' "$1" 2>/dev/null; }
genes_of() { awk -F'\t' 'NR > 1 {print $1 ":" $5}' "$1" 2>/dev/null | sort | paste -sd, -; }

# --- step 39 with the benign set ------------------------------------------------
export CARRIER_EXONS="${C}/reference/exons.bed"
run_step 39-couple.sh couple "${C}/pedigree.tsv"
check_step_exit 39-couple.sh
T="${C}/couple/couple_CSV_A_CSV_B.tsv"
MD="${C}/couple/couple_CSV_A_CSV_B.md"
cat "$T" 2>/dev/null | cut -f1,5,6,8,10
check "step 39 found the benign set where setup.sh puts it" grep -q "Common deletions: ${BEN}" "$STEP_LOG"
check_eq "with the benign set: SMPD1 is CSV_A's only, nothing else" "$(genes_of "$T")" SMPD1:one_partner
check "SMPD1: the 5 kb deletion with read depth and no common match" \
  has '^deletion chr11:[0-9]+-[0-9]+ \(5\.0 kb\) het, 1 exon of the gene \(Delly\)$' "$(col "$T" SMPD1 findings_a)"
check "the report counts the short and the common deletion of CSV_A, and the common one of CSV_B" \
  sh -c "grep -qF '| CSV_A | 0 / 0 | 1 over exons of listed genes; 1 shorter than 50 bp set aside; 1 set aside as common deletions |' '$MD' &&
         grep -qF '| CSV_B | 0 / 0 | 0 over exons of listed genes; 1 set aside as common deletions |' '$MD'"
check "the report says the filter ran on one common deletion of the set" \
  grep -qF "a deletion matching one of 1% or more in AnnotSV's benign set (1 record, 50% reciprocal overlap) is set aside" "$MD"
cp "$T" "${CASE_TMP}/with_benign.tsv"

# --- red control: no benign set -------------------------------------------------
CARRIER_BENIGN_SV="${C}/reference/none.bed" run_step 39-couple.sh couple "${C}/pedigree.tsv"
check_step_exit 39-couple.sh
cat "$T" 2>/dev/null | cut -f1,5,6,8,10
check_eq "red: without the set, HEXA (the common deletion in both) is shared" "$(col "$T" HEXA result)" shared
check_eq "red: and reported as a risk, 1 in 2 per pregnancy" "$(col "$T" HEXA offspring_risk)" "1 in 2 per pregnancy"
check "the 30 bp ASPA deletion is still not reported" test -z "$(col "$T" ASPA result)"
check "the report says common deletions were not filtered, for each partner" \
  sh -c "grep -q 'Common deletions: not filtered, no benign SV set was given' '$MD' &&
         grep -q 'Common deletions of CSV_A: not filtered' '$MD' && grep -q 'Common deletions of CSV_B: not filtered' '$MD'"
run_step 39-couple.sh couple "${C}/pedigree.tsv"
check "the step back with the set writes the first table again" cmp -s "$T" "${CASE_TMP}/with_benign.tsv"
echo "couple-common-sv: with the benign set $(genes_of "${CASE_TMP}/with_benign.tsv")" >> "$E2E_NOTES"

# --- the pipeline: --annotsv_annotations reaches COUPLE_CARRIER -----------------
# A stub run (no containers) of the stub samplesheet, whose pedigree marks two
# samples as partners. The task's staged inputs show whether the benign set
# came through.
STUB="${CASE_TMP}/stub"
rm -rf "$STUB"; mkdir -p "$STUB"
git -C "$REPO" ls-files -z | (cd "$REPO" && tar --null -T - -cf -) | tar -xf - -C "$STUB"
echo '{"kir": false}' > "${CASE_TMP}/no-kir.json"
# couple_runs WORK: the .command.run of each COUPLE_CARRIER task; staged WORK:
# how many of them stage the benign set (a staging line names the file)
couple_runs() { find "${STUB}/$1" -name .command.run -exec grep -l '^### name: .*COUPLE_CARRIER' {} + 2>/dev/null; }
ran() { couple_runs "$1" | wc -l | tr -d ' '; }
staged() { couple_runs "$1" | while read -r f; do grep -l 'benign_Loss_SV_GRCh38\.sorted\.bed' "$f"; done | wc -l | tr -d ' '; }
stub_run() {  # stub_run WORK ANNOTATIONS
  ( cd "$STUB" && nextflow run main.nf -profile test_all -stub -ansi-log false --outdir "${STUB}/out-$1" -work-dir "${STUB}/$1" \
      --pedigree assets/stub/pedigree.tsv --tools clinvar,delly,html_report --annotsv_annotations "$2" \
      -params-file "${CASE_TMP}/no-kir.json" ) > "${STUB}/$1.log" 2>&1
}
cp -R "${C}/annotsv_annotations" "${STUB}/ann-with"
mkdir -p "${STUB}/ann-without"
stub_run work-with "${STUB}/ann-with"
RC=$?
tail -n 8 "${STUB}/work-with.log"
check_eq "stub run with the benign set: exit" "$RC" 0
check_ge "it ran COUPLE_CARRIER" "$(ran work-with)" 1
check_ge "COUPLE_CARRIER got the benign set" "$(staged work-with)" 1
stub_run work-without "${STUB}/ann-without"
RC=$?
tail -n 8 "${STUB}/work-without.log"
check_eq "stub run with an annotations directory that lacks the set: exit" "$RC" 0
check_ge "it ran COUPLE_CARRIER" "$(ran work-without)" 1
check_eq "red: without the file COUPLE_CARRIER gets no benign set" "$(staged work-without)" 0
rm -rf "$STUB"
in_e2e rm -rf /w/genome-couple-common-sv
finish
