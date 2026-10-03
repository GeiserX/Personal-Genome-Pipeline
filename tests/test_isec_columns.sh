#!/usr/bin/env bash
# Step 06 + both reports on a tiny synthetic genome.
#
# Builds a 200-base reference, a three-record sample VCF and a three-record
# ClinVar VCF, runs the real scripts/06-clinvar-screen.sh, then
# scripts/24-html-report.sh and scripts/generate-report.sh, and checks:
#   1. every ClinVar hit row shows ClinVar's gene and significance. The old
#      step read bcftools isec's 0002.vcf, which is the sample's side of the
#      intersection, so every row showed ".|.";
#   2. a pathogenic allele inside a multiallelic record (GT 1/2) is found;
#   3. a VCF whose FILTER is '.' everywhere gives the same hits, with a notice;
#   4. a ClinVar file named 1,2,... against chr1,... stops with an error, and
#      a sample VCF named 1,2,... stops with a hint about the sample;
#   5. zero hits print a single 0 in both reports, with no
#      "integer expression expected" (the `grep -c ... || echo 0` bug);
#   6. a ClinVar allele the person does not carry is no hit: GT 0/0, GT ./.,
#      the 0/0 half of a multiallelic record (ALT T,G with GT 0/1, ClinVar G)
#      and a hom-ref ALT '.' row at a ClinVar ALT '.' position were all listed;
#   7. clinvar/S1_clinvar_hits.tsv has one row per hit with genotype, gene,
#      significance and review status;
#   8. a sample record on a contig the reference lacks, in a header without
#      contig lengths, and a ClinVar record on NT_113889.1 (which the full
#      chr-renamed ClinVar file has) made bcftools norm exit 255; now both are
#      left out, with a notice, and the hits are unchanged.
# Needs docker: step 06 runs bcftools in staphb/bcftools:1.21.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
BCFTOOLS_IMAGE="staphb/bcftools:1.21"
WORK=$(mktemp -d)
# Files written by containers are root-owned; remove them from a container.
cleanup() {
  docker run --rm -v "${WORK}:/w" "$BCFTOOLS_IMAGE" sh -c 'rm -rf /w/*' >/dev/null 2>&1 || true
  rm -rf "$WORK" 2>/dev/null || true
}
trap cleanup EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

# chr1 = "ACGTTGCA" x 25. Base at position p is the (p mod 8)-th letter (0 -> 8th).
SEQ=$(printf 'ACGTTGCA%.0s' $(seq 1 25))

# make_genome <dir> <sample FILTER value> <clinvar contig> <clinvar mode: hit|nohit> [offref]
# hit mode adds four ClinVar alleles the person does not carry (34, 42, 58, 66);
# offref adds a record on a contig the reference lacks to each file.
make_genome() {
  local gd=$1 filter=$2 cv_contig=$3 mode=$4 offref=${5:-}
  # expansion_hunter/ exists so that the pre-fix step 24, which stopped when step 9
  # had never run, reaches the ClinVar table and the check sees its real output.
  mkdir -p "${gd}/reference" "${gd}/clinvar" "${gd}/S1/vcf" "${gd}/S1/expansion_hunter"
  printf '>chr1\n%s\n' "$SEQ" > "${gd}/reference/Homo_sapiens_assembly38.fasta"
  printf 'chr1\t200\t6\t200\t201\n' > "${gd}/reference/Homo_sapiens_assembly38.fasta.fai"

  {
    if [ "$offref" = offref ]; then
      # No contig lengths, as some vendor headers have
      printf '##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##contig=<ID=chrUn_test>\n'
    else
      printf '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=200>\n'
    fi
    printf '##FILTER=<ID=PASS,Description="All filters passed">\n'
    printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
    printf 'chr1\t10\t.\tC\tT\t50\t%s\t.\tGT\t0/1\n' "$filter"
    printf 'chr1\t20\t.\tT\tA,G\t50\t%s\t.\tGT\t1/2\n' "$filter"
    printf 'chr1\t30\t.\tG\tA\t50\t%s\t.\tGT\t1/1\n' "$filter"
    # Positions 34, 42, 58 and 66 are C. None of these carries the ClinVar allele.
    printf 'chr1\t34\t.\tC\tT\t50\t%s\t.\tGT\t0/0\n' "$filter"
    printf 'chr1\t42\t.\tC\tA\t50\t%s\t.\tGT\t./.\n' "$filter"
    printf 'chr1\t58\t.\tC\tT,G\t50\t%s\t.\tGT\t0/1\n' "$filter"
    printf 'chr1\t66\t.\tC\t.\t50\t%s\t.\tGT\t0/0\n' "$filter"
    if [ "$offref" = offref ]; then
      printf 'chrUn_test\t5\t.\tA\tG\t50\t%s\t.\tGT\t0/1\n' "$filter"
    fi
  } > "${gd}/S1/vcf/S1.vcf"

  local p1=10 p2=20
  if [ "$mode" = nohit ]; then p1=18; p2=26; fi   # 18 -> C, 26 -> C
  {
    printf '##fileformat=VCFv4.1\n##contig=<ID=%s,length=200>\n' "$cv_contig"
    if [ "$offref" = offref ]; then printf '##contig=<ID=NT_113889.1,length=1000>\n'; fi
    printf '##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene(s)">\n'
    printf '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Significance">\n'
    printf '##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="Review status">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    printf '%s\t%s\t1001\tC\tT\t.\t.\tGENEINFO=GENEA:11;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter\n' "$cv_contig" "$p1"
    if [ "$mode" = hit ]; then
      printf '%s\t20\t1002\tT\tG\t.\t.\tGENEINFO=GENEB:22;CLNSIG=Likely_pathogenic;CLNREVSTAT=reviewed_by_expert_panel\n' "$cv_contig"
    else
      printf '%s\t%s\t1002\tC\tG\t.\t.\tGENEINFO=GENEB:22;CLNSIG=Likely_pathogenic;CLNREVSTAT=reviewed_by_expert_panel\n' "$cv_contig" "$p2"
    fi
    if [ "$mode" = hit ]; then
      printf '%s\t34\t1004\tC\tT\t.\t.\tGENEINFO=GENED:44;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter\n' "$cv_contig"
      printf '%s\t42\t1005\tC\tA\t.\t.\tGENEINFO=GENEE:55;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter\n' "$cv_contig"
    fi
    printf '%s\t50\t1003\tC\tA\t.\t.\tGENEINFO=GENEC:33;CLNSIG=Pathogenic;CLNREVSTAT=no_assertion_criteria_provided\n' "$cv_contig"
    if [ "$mode" = hit ]; then
      printf '%s\t58\t1006\tC\tG\t.\t.\tGENEINFO=GENEF:66;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter\n' "$cv_contig"
      printf '%s\t66\t1007\tC\t.\t.\t.\tGENEINFO=GENEG:77;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter\n' "$cv_contig"
    fi
    if [ "$offref" = offref ]; then
      printf 'NT_113889.1\t5\t1008\tA\tG\t.\t.\tGENEINFO=GENEH:88;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter\n'
    fi
  } > "${gd}/clinvar/clinvar_pathogenic_chr.vcf"

  docker run --rm -v "${gd}:/g" "$BCFTOOLS_IMAGE" sh -c '
    set -e
    bcftools view -Oz -o /g/S1/vcf/S1.vcf.gz /g/S1/vcf/S1.vcf
    bcftools index -t /g/S1/vcf/S1.vcf.gz
    bcftools view -Oz -o /g/clinvar/clinvar_pathogenic_chr.vcf.gz /g/clinvar/clinvar_pathogenic_chr.vcf
    bcftools index -t /g/clinvar/clinvar_pathogenic_chr.vcf.gz'
}

run_step06() {
  GENOME_DIR=$1 bash "${REPO}/scripts/06-clinvar-screen.sh" S1
}

run_reports() {
  local gd=$1
  GENOME_DIR=$gd bash "${REPO}/scripts/24-html-report.sh" S1 > "${gd}/html.log" 2>&1 || true
  GENOME_DIR=$gd bash "${REPO}/scripts/generate-report.sh" S1 > "${gd}/text.log" 2>&1 || true
}

# Hit rows of the HTML ClinVar table (one <tr> per hit)
html_rows() { grep -o '<tr><td>chr1</td>.*</tr>' "$1/S1/S1_report.html" 2>/dev/null || true; }

# ---- Case 1-2: PASS VCF, two ClinVar alleles present (one inside a multiallelic record)
GD="${WORK}/pass"
make_genome "$GD" PASS chr1 hit
run_step06 "$GD" > "${GD}/step06.log" 2>&1 || fail "step 06 exited non-zero on a valid input: $(tail -3 "${GD}/step06.log")"
run_reports "$GD"
ROWS=$(html_rows "$GD")
echo "HTML ClinVar rows:"; printf '%s\n' "${ROWS:-<none>}" | sed 's/^/    /'
if printf '%s' "$ROWS" | grep -q '\.|\.'; then
  fail "HTML ClinVar table shows '.|.' instead of gene and significance"
fi
if printf '%s' "$ROWS" | grep -q '<td>GENEA</td><td>Pathogenic</td>'; then
  pass "HTML row for chr1:10 shows GENEA / Pathogenic"
else
  fail "HTML row for chr1:10 does not show GENEA / Pathogenic"
fi
if printf '%s' "$ROWS" | grep -q '<td>GENEB</td><td>Likely pathogenic</td>'; then
  pass "pathogenic allele inside a multiallelic record (GT 1/2) is reported"
else
  fail "pathogenic allele inside a multiallelic record (GT 1/2) is missing"
fi
if [ "$(printf '%s\n' "$ROWS" | grep -c '<td>het</td>' || true)" -eq 2 ]; then
  pass "both hits show genotype het"
else
  fail "hits do not show genotype het"
fi
if grep -q 'Pathogenic/Likely Pathogenic hits: 2$' "${GD}/S1/S1_report.txt" 2>/dev/null \
   && grep -q 'GENEA' "${GD}/S1/S1_report.txt" && grep -q 'GENEB' "${GD}/S1/S1_report.txt"; then
  pass "text report lists 2 hits with GENEA and GENEB"
else
  fail "text report does not list 2 hits with genes: $(grep -A4 'ClinVar' "${GD}/S1/S1_report.txt" 2>/dev/null | tr '\n' '|')"
fi

# ---- Case 6: alleles the person does not carry are no hits
HITS_VCF="${GD}/S1/clinvar/S1_clinvar_hits.vcf"
NONCARRIER=$(grep -v '^#' "$HITS_VCF" 2>/dev/null \
  | awk -F'\t' '$2 == 34 || $2 == 42 || $2 == 66 || ($2 == 58 && $5 == "G") {print $2 ":" $4 ">" $5 " GT " $10}' \
  | paste -sd' ' -)
if [ -z "$NONCARRIER" ] && [ "$(grep -c -v '^#' "$HITS_VCF" 2>/dev/null || true)" = 2 ]; then
  pass "GT 0/0, GT ./., the 0/0 half of a multiallelic record and a hom-ref ALT '.' row are not hits"
else
  fail "hits the person does not carry are listed: ${NONCARRIER:-none}; $(grep -c -v '^#' "$HITS_VCF" 2>/dev/null || true) hits, want 2"
fi
if grep -qE '^##bcftools_(isec|annotate)Command' "$HITS_VCF" 2>/dev/null; then
  fail "the hits VCF carries the screen's own isec/annotate command lines"
else
  pass "the hits VCF carries no isec/annotate command line"
fi

# ---- Case 7: the hits as a table
TSV="${GD}/S1/clinvar/S1_clinvar_hits.tsv"
if [ "$(head -n 1 "$TSV" 2>/dev/null)" = "$(printf 'chrom\tpos\tref\talt\tgenotype\tclinvar_id\tgeneinfo\tclnsig\tclnrevstat')" ] \
   && [ "$(awk 'END {print NR}' "$TSV")" = 3 ] \
   && grep -qxF "$(printf 'chr1\t10\tC\tT\t0/1\t1001\tGENEA:11\tPathogenic\tcriteria_provided,_single_submitter')" "$TSV" \
   && grep -qxF "$(printf 'chr1\t20\tT\tG\t0/1\t1002\tGENEB:22\tLikely_pathogenic\treviewed_by_expert_panel')" "$TSV"; then
  pass "hits TSV: a header and one row per hit with genotype, gene, significance and review status"
else
  fail "hits TSV is missing or wrong: $(tr '\t\n' ' |' < "$TSV" 2>/dev/null || echo 'no file')"
fi

# ---- Case 8: a contig the reference lacks, on either side, is left out with a notice
GD="${WORK}/offref"
make_genome "$GD" PASS chr1 hit offref
if run_step06 "$GD" > "${GD}/step06.log" 2>&1; then
  pass "a record on a contig the reference lacks (sample chrUn_test, ClinVar NT_113889.1) no longer stops step 06"
else
  fail "step 06 exited non-zero with records on contigs the reference lacks: $(tail -3 "${GD}/step06.log" | tr '\n' '|')"
fi
if grep -q "left out.*chrUn_test" "${GD}/step06.log"; then
  pass "the left-out contig is named in a notice"
else
  fail "no notice names the left-out contig chrUn_test"
fi
N=$(grep -c -v '^#' "${GD}/S1/clinvar/S1_clinvar_hits.vcf" 2>/dev/null || true)
if [ "${N:-0}" = 2 ] && grep -q 'GENEA' "${GD}/S1/clinvar/S1_clinvar_hits.vcf" && grep -q 'GENEB' "${GD}/S1/clinvar/S1_clinvar_hits.vcf"; then
  pass "the same 2 hits with off-reference contigs present"
else
  fail "off-reference contigs changed the hits: ${N:-no} hits, want 2"
fi

# ---- Case 3: FILTER '.' everywhere gives the same hits, with a notice
GD="${WORK}/dotfilter"
make_genome "$GD" . chr1 hit
run_step06 "$GD" > "${GD}/step06.log" 2>&1 || fail "step 06 exited non-zero on an all-'.' FILTER VCF"
if grep -q "no PASS record" "${GD}/step06.log"; then
  pass "all-'.' FILTER VCF prints the filter-mode notice"
else
  fail "all-'.' FILTER VCF prints no notice"
fi
N=$(grep -c -v '^#' "${GD}/S1/clinvar/S1_clinvar_hits.vcf" 2>/dev/null || true)
if [ "${N:-0}" = 2 ]; then
  pass "all-'.' FILTER VCF gives the same 2 hits"
else
  fail "all-'.' FILTER VCF gives ${N:-no} hits, expected 2"
fi

# ---- Case 4: ClinVar named 1,2,... against a chr-prefixed VCF must fail, not report 0
GD="${WORK}/nochr"
make_genome "$GD" PASS 1 hit
if run_step06 "$GD" > "${GD}/step06.log" 2>&1; then
  fail "ClinVar without chr prefix exited 0 (zero hits would look like a clean result)"
elif grep -q 'no contig name in common' "${GD}/step06.log"; then
  pass "ClinVar without chr prefix exits non-zero with the contig message"
else
  fail "ClinVar without chr prefix failed without the contig message: $(tail -2 "${GD}/step06.log")"
fi

# ---- Case 4b: a sample VCF named 1,2,... is the one to blame, not ClinVar
GD="${WORK}/samplenochr"
make_genome "$GD" PASS chr1 hit
docker run --rm -v "${GD}:/g" "$BCFTOOLS_IMAGE" sh -c '
  set -e
  sed -e "s/^chr1\t/1\t/" -e "s/ID=chr1,/ID=1,/" /g/S1/vcf/S1.vcf | bcftools view -Oz -o /g/S1/vcf/S1.vcf.gz
  bcftools index -f -t /g/S1/vcf/S1.vcf.gz'
if run_step06 "$GD" > "${GD}/step06.log" 2>&1; then
  fail "a sample VCF named 1,2,... exited 0"
elif grep -q 'sample VCF is not chr-named' "${GD}/step06.log" && ! grep -q 'Rebuild the ClinVar file' "${GD}/step06.log"; then
  pass "a sample VCF named 1,2,... stops with a message about the sample, not ClinVar"
else
  fail "a sample VCF named 1,2,... got the wrong hint: $(grep -A4 'no contig name' "${GD}/step06.log" | tr '\n' '|')"
fi

# ---- Case 5: zero hits print a single 0 in both reports
GD="${WORK}/zero"
make_genome "$GD" PASS chr1 nohit
run_step06 "$GD" > "${GD}/step06.log" 2>&1 || fail "step 06 exited non-zero on a zero-hit input"
run_reports "$GD"
if grep -q 'Count: 0 pathogenic hits' "${GD}/step06.log" && ! grep -qx '0' "${GD}/step06.log"; then
  pass "step 06 prints a single 0"
else
  fail "step 06 count line is not a single 0: $(grep -A1 'Count:' "${GD}/step06.log" | tr '\n' '|')"
fi
if grep -q 'class="badge badge-green">0</span>' "${GD}/S1/S1_report.html" 2>/dev/null; then
  pass "HTML badge shows a single 0"
else
  fail "HTML badge is not a single 0: $(grep -A1 'ClinVar matches' "${GD}/S1/S1_report.html" 2>/dev/null | tr '\n' ' ')"
fi
if grep -q 'Pathogenic/Likely Pathogenic hits: 0$' "${GD}/S1/S1_report.txt" 2>/dev/null \
   && ! grep -A1 'Pathogenic/Likely Pathogenic hits' "${GD}/S1/S1_report.txt" | grep -qx '0'; then
  pass "text report shows a single 0"
else
  fail "text report hit count is not a single 0: $(grep -A1 'Pathogenic/Likely Pathogenic hits' "${GD}/S1/S1_report.txt" 2>/dev/null | tr '\n' '|')"
fi
if grep -q 'integer expression expected' "${GD}/text.log" "${GD}/html.log"; then
  fail "a report printed 'integer expression expected'"
else
  pass "no 'integer expression expected'"
fi

echo ""
if [ "$FAILS" -gt 0 ]; then
  echo "${FAILS} check(s) failed"
  exit 1
fi
echo "All checks passed"
