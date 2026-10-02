#!/usr/bin/env bash
# `setup.sh --refresh clinvar` replaces an installed ClinVar for real:
#   - the raw file is NCBI's current one, checked against NCBI's .md5;
#   - clinvar_chr and clinvar_pathogenic_chr (and their indexes) are rebuilt
#     from it, and step 06's normalised copy is removed so step 06 rebuilds it;
#   - clinvar/RELEASE holds the new ##fileDate;
#   - a failed download leaves every installed file as it was.
# validate-setup.sh prints the release and warns when it is over 35 days old.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
C="${GENOME_DIR}/clinvar"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
for f in clinvar_pathogenic_chr.norm.vcf.gz clinvar_pathogenic_chr.norm.vcf.gz.tbi; do
  printf 'old\n' > "${C}/${f}"
done
echo 2025-01-01 > "${C}/RELEASE"
use_output_hook

# What NCBI serves: a ClinVar VCF of 2026-09-28.
SERVED="${CASE_WORK}/served"
mkdir -p "$SERVED"
printf '##fileformat=VCFv4.1\n##fileDate=2026-09-28\n##source=ClinVar\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t1\tA\tG\t.\t.\tCLNSIG=Pathogenic\n' \
  | gzip -c > "${SERVED}/clinvar.vcf.gz"
printf 'new index\n' > "${SERVED}/clinvar.vcf.gz.tbi"

run_expect 1 refresh-bad-name "${SCRIPTS}/setup.sh" --refresh dbsnp "$GENOME_DIR"
output_has refresh-bad-name 'Usage: .* --refresh clinvar <genome_dir>'

: > "$FAKE_DOCKER_LOG"
FAKE_DOWNLOAD_DIR="$SERVED" run_expect 0 refresh "${SCRIPTS}/setup.sh" --refresh clinvar "$GENOME_DIR"
docker_log_has '^curl [^ ]*/clinvar\.vcf\.gz\.md5 ' "the refresh did not read NCBI's md5"
cmp -s "${SERVED}/clinvar.vcf.gz" "${C}/clinvar.vcf.gz" || fail "clinvar.vcf.gz is not the served release"
cmp -s "${SERVED}/clinvar.vcf.gz.tbi" "${C}/clinvar.vcf.gz.tbi" || fail "clinvar.vcf.gz.tbi is not the served index"
for f in clinvar_chr.vcf.gz clinvar_pathogenic_chr.vcf.gz; do
  gzip -cd "${C}/${f}" 2>/dev/null | grep -q '^##fileformat=VCF' || fail "${f} was not rebuilt"
  [ -f "${C}/${f}.tbi" ] && ! grep -q placeholder "${C}/${f}.tbi" || fail "${f}.tbi was not rebuilt"
done
docker_log_has '^run .*--rename-chrs .*/clinvar/\.refresh/chr_rename\.txt /genome/clinvar/\.refresh/clinvar\.vcf\.gz ' \
  "clinvar_chr.vcf.gz was not built from the new download"
[ ! -e "${C}/clinvar_pathogenic_chr.norm.vcf.gz" ] || fail "step 06's normalised copy of the old release was kept"
[ "$(cat "${C}/RELEASE")" = 2026-09-28 ] || fail "clinvar/RELEASE is '$(cat "${C}/RELEASE")', want 2026-09-28"
[ ! -e "${C}/.refresh" ] || fail "the refresh left clinvar/.refresh behind"
output_has refresh 'ClinVar release 2026-09-28'

# A refresh whose download fails changes nothing.
BEFORE=$(cd "$C" && cksum clinvar.vcf.gz clinvar_chr.vcf.gz clinvar_pathogenic_chr.vcf.gz RELEASE)
FAKE_DOWNLOAD_DIR="$SERVED" FAKE_DOWNLOAD_FAIL='clinvar\.vcf\.gz$' run_rc refresh-offline "${SCRIPTS}/setup.sh" --refresh clinvar "$GENOME_DIR"
[ "$RC" -ne 0 ] || fail "a refresh whose download failed exited 0"
[ "$(cd "$C" && cksum clinvar.vcf.gz clinvar_chr.vcf.gz clinvar_pathogenic_chr.vcf.gz RELEASE)" = "$BEFORE" ] \
  || fail "a failed refresh changed the installed ClinVar files"

# validate-setup.sh: the release and its age.
echo 2025-01-01 > "${C}/RELEASE"
run_rc validate-old "${SCRIPTS}/validate-setup.sh" sample1
output_has validate-old 'ClinVar release 2025-01-01 is [0-9]+ days old \(over 35\)'
output_has validate-old 'setup\.sh --refresh clinvar'
date -u +%Y-%m-%d > "${C}/RELEASE"
run_rc validate-new "${SCRIPTS}/validate-setup.sh" sample1
output_has validate-new "ClinVar release $(date -u +%Y-%m-%d) \(0 days old\)"
output_lacks validate-new 'over 35'
rm "${C}/RELEASE"
run_rc validate-unknown "${SCRIPTS}/validate-setup.sh" sample1
output_has validate-unknown 'ClinVar release date unknown'
