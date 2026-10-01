#!/usr/bin/env bash
# Scratch helper (deleted before merge): real runs of CLINICAL_FILTER and
# SLIVAR in their containers on a small hand-made VEP-style VCF.
# usage: modules.sh <repo_root_with_modules> <work_dir> <expect: new|main>
set -euo pipefail

SRC=${1:?}; T=${2:?}; EXPECT=${3:?}
mkdir -p "$T"; cd "$T"

# --- inputs ---------------------------------------------------------------
if [ ! -s constraint.tsv ]; then
  curl -fsSL -o constraint.tsv https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv
  curl -fsSL -o slivar https://github.com/brentp/slivar/releases/download/v0.3.1/slivar
fi
mkdir -p fake
printf '#!/bin/sh\necho "fake slivar: crashing on purpose" >&2\nexit 1\n' > fake/slivar

# VEP-style CSQ with the fields slivar needs (Gene, Feature). test.vcf adds a
# SpliceAI INFO record for CLINICAL_FILTER; slivar.vcf has none, so SLIVAR's
# predictor branch does not matter for the slivar checks.
header() {
  echo '##fileformat=VCFv4.2'
  echo '##contig=<ID=chr2,length=242193529>'
  echo '##contig=<ID=chr17,length=83257441>'
  echo '##FILTER=<ID=PASS,Description="All filters passed">'
  echo '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|Existing_variation|gnomADe_AF|CLIN_SIG">'
  if [ "${1:-}" = spliceai ]; then
    echo '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3 variant annotation. Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">'
  fi
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
}
records() {
  printf 'chr2\t165310000\t.\tC\tT\t50\tPASS\tCSQ=T|stop_gained|HIGH|SCN2A|ENSG00000136531|Transcript|ENST00000375437|protein_coding|.|0.00001|.\tGT\t0/1\n'
  printf 'chr2\t178527000\t.\tA\tG\t50\tPASS\tCSQ=G|missense_variant|MODERATE|TTN|ENSG00000155657|Transcript|ENST00000589042|protein_coding|.|0.001|.\tGT\t0/1\n'
  printf 'chr2\t178528000\t.\tT\tC\t50\tPASS\tCSQ=C|missense_variant|MODERATE|TTN|ENSG00000155657|Transcript|ENST00000589042|protein_coding|.|0.002|.\tGT\t0/1\n'
  printf 'chr17\t43045712\t.\tG\tA\t50\tPASS\tCSQ=A|stop_gained|HIGH|BRCA1|ENSG00000012048|Transcript|ENST00000357654|protein_coding|.|0.0001|pathogenic\tGT\t0/1\n'
}
{ header spliceai; records
  printf 'chr17\t43045800\t.\tC\tT\t50\tPASS\tCSQ=T|intron_variant|MODIFIER|BRCA1|ENSG00000012048|Transcript|ENST00000357654|protein_coding|.|0.0002|.;SpliceAI=T|BRCA1|0.85|0.00|0.00|0.00|1|2|3|4\tGT\t0/1\n'
} > test.vcf
{ header; records; } > slivar.vcf
for v in test slivar; do
  bcftools view "$v.vcf" -Oz -o "$v.vcf.gz"
  bcftools index -f -t "$v.vcf.gz"
done

cat > ci.config <<'CICONF'
docker.enabled = true
docker.runOptions = '-u 0:0'
trace.enabled = false
process.errorStrategy = 'finish'
CICONF

cat > test.nf <<NF
nextflow.enable.dsl = 2
params.outdir = 'out'
params.publish_dir_mode = 'copy'
params.slivar_bin = null
include { CLINICAL_FILTER } from '${SRC}/modules/local/clinical_filter/main'
include { SLIVAR }          from '${SRC}/modules/local/slivar/main'
workflow {
    CLINICAL_FILTER(Channel.of([[id: 'S1'], file("${T}/test.vcf.gz"), file("${T}/test.vcf.gz.tbi")]))
    SLIVAR(Channel.of([[id: 'S1'], file("${T}/slivar.vcf.gz"), file("${T}/slivar.vcf.gz.tbi")]), file("${T}/constraint.tsv"), file(params.slivar_bin))
}
NF

cat > test_slivar.nf <<NF
nextflow.enable.dsl = 2
params.outdir = 'out_fake'
params.publish_dir_mode = 'copy'
params.slivar_bin = null
include { SLIVAR } from '${SRC}/modules/local/slivar/main'
workflow {
    SLIVAR(Channel.of([[id: 'S1'], file("${T}/slivar.vcf.gz"), file("${T}/slivar.vcf.gz.tbi")]), file("${T}/constraint.tsv"), file(params.slivar_bin))
}
NF

fail=0

echo "== CLINICAL_FILTER + SLIVAR with the real slivar 0.3.1 and the gnomAD v4.1 constraint table"
if nextflow run test.nf -c ci.config --slivar_bin "${T}/slivar" > run.log 2>&1; then
  echo "run: succeeded"
  ok=1
else
  echo "run: FAILED"; ok=0
  grep -v -E 'Pulling|Waiting|Verifying|Download complete|Pull complete|Digest:|Status:|Unable to find image' run.log | tail -40
  for d in work/*/*/; do
    echo "--- $d exit=$(cat "$d/.exitcode" 2>/dev/null || echo none) $(grep -m1 -oE 'CLINICAL_FILTER|SLIVAR' "$d/.command.run" || true)"
    grep -v -E 'Pulling|Waiting|Verifying|Download complete|Pull complete|Digest:|Status:|Unable to find image' "$d/.command.err" | tail -15
  done
fi

if [ "$EXPECT" = "new" ]; then
  [ "$ok" -eq 1 ] || fail=1
  echo "== clinical VCF holds the SpliceAI-high record (written via the _spliceai_high tier)"
  bcftools view -H out/S1/clinical/S1_clinical.vcf.gz | cut -f1-5
  if bcftools view -H out/S1/clinical/S1_clinical.vcf.gz | grep -q $'chr17\t43045800'; then echo "PASS: SpliceAI tier present"; else echo "FAIL: SpliceAI tier missing"; fail=1; fi
  grep -h -B2 -A3 'spliceai_high' work/*/*/.command.sh | head -12
  echo "== slivar summary"
  column -t -s $'\t' out/S1/slivar/S1_slivar_summary.tsv
  misz=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next } $6 == "SCN2A" { print $c["mis_z"] "\t" $c["CONSTRAINED"]; exit }' out/S1/slivar/S1_slivar_summary.tsv)
  echo "SCN2A mis_z / CONSTRAINED: ${misz}"
  if [ -n "$misz" ] && [ "${misz%%$'\t'*}" != "." ]; then echo "PASS: SCN2A mis_z is a number"; else echo "FAIL: SCN2A mis_z"; fail=1; fi
  grep -h 'gnomAD constraint' work/*/*/.command.err || true
  echo "== compound hets VCF"
  bcftools view -H out/S1/slivar/S1_compound_hets.vcf.gz | cut -f1-5 || true
fi

echo "== SLIVAR with a slivar binary that exits 1"
if nextflow run test_slivar.nf -c ci.config --slivar_bin "${T}/fake/slivar" > run_fake.log 2>&1; then
  echo "fake slivar run: SUCCEEDED"
  fake_ok=1
else
  echo "fake slivar run: FAILED"
  fake_ok=0
  grep -E "terminated with an error exit status|fake slivar" run_fake.log | head -5
fi
if [ "$EXPECT" = "new" ] && [ "$fake_ok" -eq 1 ]; then echo "FAIL: a crashing slivar did not fail the task"; fail=1; fi
if [ "$EXPECT" = "new" ] && [ "$fake_ok" -eq 0 ]; then echo "PASS: a crashing slivar fails the task"; fi
if [ "$EXPECT" = "main" ]; then
  echo "main: real run ok=${ok}; crashing slivar run ok=${fake_ok} (1 means the crash was swallowed)"
  echo "main: slivar summary header (no LOEUF/pLI/mis_z means the constraint join was skipped):"
  head -1 out/S1/slivar/S1_slivar_summary.tsv 2>/dev/null || echo "(no summary)"
  echo "main: SLIVAR log with the crashing binary:"
  grep -h -E 'Compound het detection returned no results|fake slivar' work/*/*/.command.out work/*/*/.command.err 2>/dev/null | sort -u || true
fi

exit "$fail"
