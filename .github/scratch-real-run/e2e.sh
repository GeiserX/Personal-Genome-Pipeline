#!/usr/bin/env bash
# Scratch helper (deleted before merge): real run of the pipeline on the
# HG002 chr20 slice, then the input-check controls and two stub checks.
set -euo pipefail

D=${1:?usage: e2e.sh <data_dir>}
cd "${GITHUB_WORKSPACE:?}"

common=(-profile docker -c "$D/ci.config" --reference "$D/ref.fa")

echo "::group::real run"
test ! -e "$D/ref.dict" && echo "no ref.dict next to the reference, and mito_variants is off"
nextflow run main.nf "${common[@]}" \
  --input "$D/samplesheet.csv" \
  --outdir results_e2e \
  --tools vcfanno,delly,clinvar,mosdepth,manta,duphold,survivor_merge,expansion_hunter \
  --clinvar "$D/clinvar_pathogenic_chr.vcf.gz" --clinvar_index "$D/clinvar_pathogenic_chr.vcf.gz.tbi" \
  --cadd_snv "$D/cadd_tiny.tsv.gz" --cadd_snv_index "$D/cadd_tiny.tsv.gz.tbi" \
  --spliceai_snv "$D/spliceai_scores.masked.snv.hg38.vcf.gz" --spliceai_snv_index "$D/spliceai_scores.masked.snv.hg38.vcf.gz.tbi" \
  --expansion_catalog "$D/eh_catalog_chr20.json" \
  --delly_exclude "$D/delly_exclude.tsv"
echo "::endgroup::"

R=results_e2e/HG002
fail=0
check() { if "$@"; then echo "PASS: $*"; else echo "FAIL: $*"; fail=1; fi; }

echo "== published files"
find results_e2e -type f | sort

echo "== VCF_PRECHECK"
cat "$(grep -l 'FILTER counts' work/*/*/.command.out | head -1)"

echo "== VCFANNO: indexed annotated VCF that carries the score tags"
check test -s "$R/vep/HG002_annotated.vcf.gz"
check test -s "$R/vep/HG002_annotated.vcf.gz.tbi"
bcftools view -h "$R/vep/HG002_annotated.vcf.gz" | grep -E '^##INFO=<ID=(CADD_PHRED|SpliceAI),' || true
n_cadd_in=$(zcat "$D/cadd_tiny.tsv.gz" | grep -vc '^#')
n_cadd_out=$(bcftools query -f '%INFO/CADD_PHRED\n' "$R/vep/HG002_annotated.vcf.gz" | grep -c '^25\.5$' || true)
n_spl_in=$(zcat "$D/spliceai_scores.masked.snv.hg38.vcf.gz" | grep -vc '^#')
n_spl_out=$(bcftools query -f '%INFO/SpliceAI\n' "$R/vep/HG002_annotated.vcf.gz" | grep -c 'TESTGENE' || true)
echo "CADD rows in the score file: ${n_cadd_in}; records annotated with CADD_PHRED: ${n_cadd_out}"
echo "SpliceAI rows in the score file: ${n_spl_in}; records annotated with SpliceAI: ${n_spl_out}"
bcftools view -H "$R/vep/HG002_annotated.vcf.gz" | grep -m2 'CADD_PHRED=' | cut -f1-8 || true
bcftools view -H "$R/vep/HG002_annotated.vcf.gz" | grep -m2 'SpliceAI=' | cut -f1-8 || true
check test "$n_cadd_out" -eq "$n_cadd_in"
check test "$n_spl_out" -eq "$n_spl_in"
check test "$(bcftools view -H "$R/vep/HG002_annotated.vcf.gz" | wc -l)" -eq "$(bcftools view -H "$D/HG002.vcf.gz" | wc -l)"

echo "== DELLY: indexed VCF, exclude map passed"
check test -s "$R/delly/HG002_sv.vcf.gz"
check test -s "$R/delly/HG002_sv.vcf.gz.tbi"
echo "delly records: $(bcftools view -H "$R/delly/HG002_sv.vcf.gz" | wc -l)"
grep -h -A4 'delly call' work/*/*/.command.sh
check grep -q -- '-x delly_exclude.tsv' "$(grep -l 'delly call' work/*/*/.command.sh | head -1)"

echo "== DUPHOLD: annotated and filtered outputs"
check test -s "$R/sv_duphold/HG002_sv_duphold.vcf"
check test -s "$R/sv_filtered/HG002_sv_filtered.vcf.gz"
check test -s "$R/sv_filtered/HG002_sv_filtered.vcf.gz.tbi"
cat "$R/sv_filtered/HG002_sv_filtered.log"
n_ann=$(grep -vc '^#' "$R/sv_duphold/HG002_sv_duphold.vcf" || true)
n_fil=$(bcftools view -H "$R/sv_filtered/HG002_sv_filtered.vcf.gz" | wc -l)
echo "duphold annotated: ${n_ann}; filtered: ${n_fil}"
if [ "$n_fil" -lt "$n_ann" ]; then
  echo "PASS: filtered file has fewer records"
elif [ "$n_fil" -eq "$n_ann" ] && grep -q 'nothing removed' "$R/sv_filtered/HG002_sv_filtered.log"; then
  echo "PASS: equal counts, reason logged"
else
  echo "FAIL: duphold filter counts"; fail=1
fi
bcftools query -f '%CHROM:%POS %INFO/SVTYPE [DHFFC=%DHFFC DHBFC=%DHBFC]\n' "$R/sv_duphold/HG002_sv_duphold.vcf" | head -20

echo "== SURVIVOR_MERGE, CLINVAR, MOSDEPTH"
check test -s "$R/sv_merged/HG002_sv_consensus.vcf.gz"
check test -d "$R/clinvar"
check test -d "$R/coverage"

echo "== EXPANSION_HUNTER: --sex male in the command"
grep -h -A8 'ExpansionHunter' work/*/*/.command.sh | head -12
check grep -q -- '--sex male' "$(grep -l 'ExpansionHunter' work/*/*/.command.sh | head -1)"

echo "== input-check controls (each must stop with the quoted message)"
expect_fail() {
  local msg=$1; shift
  if nextflow run main.nf "$@" > control.log 2>&1; then
    echo "FAIL (run succeeded): $msg"; tail -20 control.log; fail=1; return 0
  fi
  if grep -F -- "$msg" control.log; then echo "PASS: stopped with: $msg"; else echo "FAIL (message missing): $msg"; tail -30 control.log; fail=1; fi
}
expect_fail "unknown tool clinvar_screen" "${common[@]}" --input "$D/samplesheet.csv" --outdir ctl1 --tools vep,clinvar_screen
expect_fail "Sample 'HG002' appears more than once" "${common[@]}" --input "$D/samplesheet_dup.csv" --outdir ctl2 --tools roh
expect_fail "Tool 'cnvpytor' is enabled in --tools but --cnvpytor_resources is not set" "${common[@]}" --input "$D/samplesheet.csv" --outdir ctl3 --tools cnvpytor
expect_fail "Tool 'survivor_merge' needs at least two SV callers" "${common[@]}" --input "$D/samplesheet.csv" --outdir ctl4 --tools manta,survivor_merge
expect_fail "Tool 'annotsv' is enabled in --tools but --annotsv_annotations is not set" "${common[@]}" --input "$D/samplesheet.csv" --outdir ctl5 --tools manta,duphold,annotsv
expect_fail "expansion_hunter needs the sample's sex" "${common[@]}" --input "$D/samplesheet_nosex.csv" --outdir ctl6 --tools expansion_hunter --expansion_catalog "$D/eh_catalog_chr20.json"
bgzip -c "$D/ref.fa" > "$D/ref.fa.gz"
expect_fail "is compressed" -profile docker -c "$D/ci.config" --reference "$D/ref.fa.gz" --input "$D/samplesheet.csv" --outdir ctl7 --tools roh
expect_fail "Sample 'HG002_nofilter': no record in HG002_nofilter.vcf.gz has FILTER=PASS" "${common[@]}" --input "$D/samplesheet_nofilter.csv" --outdir ctl8 --tools roh

echo "== --allow_unfiltered positive control"
if nextflow run main.nf "${common[@]}" --input "$D/samplesheet_nofilter.csv" --outdir ctl9 --tools roh --allow_unfiltered > control.log 2>&1; then
  if grep -F "allow_unfiltered is set" control.log; then echo "PASS: run went on with the warning"; else echo "FAIL: warning missing"; fail=1; fi
else
  echo "FAIL: --allow_unfiltered run failed"; tail -30 control.log; fail=1
fi

echo "== default stub run: two skip warnings, no vep/ directory"
nextflow run main.nf -profile test,docker -stub -c "$D/ci.config" --outdir stub_default > stub.log 2>&1 || { tail -40 stub.log; fail=1; }
check grep -q -F "vcfanno skipped: no score file is set" stub.log
check grep -q -F "prs skipped: --pgs_scoring is not set" stub.log
check test ! -e stub_default/stub_sample/vep
check test ! -e stub_default/stub_sample/prs

echo "== stub run with manta,duphold: both duphold directories"
nextflow run main.nf -profile test,docker -stub -c "$D/ci.config" --outdir stub_sv --tools manta,duphold > stub_sv.log 2>&1 || { tail -40 stub_sv.log; fail=1; }
find stub_sv/stub_sample -maxdepth 2 | sort
check test -d stub_sv/stub_sample/sv_duphold
check test -d stub_sv/stub_sample/sv_filtered

echo "== memory closures"
nextflow config main.nf -flat -profile docker | grep -E 'process\.(memory|time)|withLabel:process_(single|low|medium|high)\.(memory|time)' || true
check bash -c "nextflow config main.nf -flat -profile docker | grep -q 'task.attempt'"

exit "$fail"
