# shellcheck shell=sh
# PCGR_IMAGE row. CPSR (step 17) needs the PCGR data bundle (about 8 GB) and a
# VEP 113 cache (about 26 GB), which do not fit a runner. This runs what the
# image itself carries: cpsr, the VEP it calls (with --database instead of the
# cache) and the pcgrr R package that writes the report.
set -e
cpsr --version > cpsr_version.txt 2>&1 || true
cat cpsr_version.txt
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
