# shellcheck shell=sh
# PCGR_IMAGE row. CPSR (step 17) needs the PCGR data bundle (about 7 GB) and
# the VEP cache of PCGR_VEP_CACHE_RELEASE (about 24 GB), which do not fit a
# runner. This runs what the image itself carries: cpsr, the VEP it calls
# (with --database instead of the cache) and the pcgrr R package that writes
# the report. It also checks that cpsr still takes every flag step 17 passes:
# 2.3 dropped --classify_all, and an unknown flag stops cpsr before it reads
# the VCF.
set -e
cpsr --version > cpsr_version.txt 2>&1 || true
cat cpsr_version.txt
cpsr --help > cpsr_help.txt 2>&1 || true
: > cpsr_missing_flags.txt
# The flags of scripts/17-cpsr.sh and modules/local/cpsr/main.nf.
for f in --input_vcf --vep_dir --refdata_dir --output_dir --genome_assembly --sample_id \
         --panel_id --secondary_findings --force_overwrite; do
  grep -q -e "$f" cpsr_help.txt || echo "$f" >> cpsr_missing_flags.txt
done
echo "flags step 17 passes that cpsr --help does not list: $(tr '\n' ' ' < cpsr_missing_flags.txt)"
# vep is in the pcgr conda environment, on PATH; Rscript is in the pcgrr one.
find_bin() {
  command -v "$1" 2>/dev/null && return 0
  for b in /opt/*/envs/*/bin/"$1"; do
    [ -x "$b" ] && { echo "$b"; return 0; }
  done
  echo "ERROR: $1 is not on PATH nor in /opt/*/envs/*/bin" >&2
  return 1
}
VEP=$(find_bin vep)
RSCRIPT=$(find_bin Rscript)
echo "vep: ${VEP}; Rscript: ${RSCRIPT}"
"$VEP" --input_file /in/sample50.vcf --output_file vep.vcf --vcf --database --assembly GRCh38 \
  --symbol --force_overwrite --no_stats
"$RSCRIPT" -e 'suppressMessages(library(pcgrr)); cat(as.character(packageVersion("pcgrr")), "\n")' > pcgrr_version.txt
cat pcgrr_version.txt
