#!/usr/bin/env bash
# Step 31 phases each compound-het pair from the reads, and the Nextflow
# COMPHET_PHASE process writes the same pairs table.
#
# The fixture is tests/fixtures/comphet/make_fixture.py: read pairs on a
# made-up contig (chrP) and a VEP VCF with three pairs, all rare and HIGH:
#   GENEC  SNV at 1000 and a deletion at 1003 on one copy   -> cis, 22 reads
#   GENET  SNVs at 5000 and 5400, one on each copy; only a read and its mate
#          cover both                                       -> trans, 20 read pairs
#   GENEF  SNVs at 20000 and 70000, 50 kb apart             -> unknown
# CHPS is the same fixture with GENEC's deletion on the other copy: the same
# positions now read trans, so the cis check rests on the alleles in the reads.
# The fixture sample (case 42) and MANEQ (case clinical-filter-mane-and-quality-1,
# no BAM) are checked too: every pair gets a verdict, and without a BAM it is
# unknown with the reason.
. "$(dirname "$0")/lib.sh"

COLS=$'GENE\tCHROM\tPOS_1\tREF_1\tALT_1\tPOS_2\tREF_2\tALT_2\tDISTANCE\tPHASE\tCIS_READS\tTRANS_READS\tPHASE_NOTE'
# pair FILE GENE: "PHASE CIS_READS TRANS_READS" of the gene's pair.
pair() { awk -F'\t' -v g="$2" '$1 == g {print $10, $11, $12; exit}' "$1" 2>/dev/null; }
# cell FILE POS COLUMN: the value of COLUMN in the row of POS.
cell() {
  awk -F'\t' -v p="$2" -v k="$3" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    $c["POS"] == p && (k in c) {print $c[k]; exit}' "$1" 2>/dev/null
}

# make_sample S [--swap]: the fixture as sample S, its reads as a sorted, indexed BAM.
make_sample() {
  local s=$1 f="${CASE_TMP}/fixture-$1"
  rm -rf "$f" "${GENOME_DIR:?}/${s}"
  python3 "${REPO}/tests/fixtures/comphet/make_fixture.py" "$f" "${@:2}"
  mkdir -p "${GENOME_DIR}/${s}/vep" "${GENOME_DIR}/${s}/aligned"
  sed "s/CHPH/${s}/g" "${f}/vep.vcf" > "${GENOME_DIR}/${s}/vep/${s}_vep.vcf"
  sed "s/CHPH/${s}/g" "${f}/reads.sam" > "${GENOME_DIR}/${s}/aligned/${s}.sam"
  sam sort -o "${s}/aligned/${s}_sorted.bam" "${s}/aligned/${s}.sam"
  sam index "${s}/aligned/${s}_sorted.bam"
  check "${s}: the fixture BAM is indexed" nonempty "${s}/aligned/${s}_sorted.bam.bai"
}

# --- step 31 on CHPH -----------------------------------------------------------
make_sample CHPH
run_step 31-slivar.sh CHPH
check_step_exit 31-slivar.sh
LOG=$(cat "$STEP_LOG")
P="${GENOME_DIR}/CHPH/slivar/CHPH_compound_het_pairs.tsv"
CH="${GENOME_DIR}/CHPH/slivar/CHPH_compound_hets.tsv"
cat "$P" "$CH" 2>/dev/null
check_eq "step 31: the pairs table's header" "$(head -1 "$P" 2>/dev/null)" "$COLS"
check_eq "step 31: three pairs" "$(tail -n +2 "$P" 2>/dev/null | wc -l | tr -d ' ')" 3
check_eq "step 31: GENEC (SNV and deletion 3 bp apart, one copy) is cis on 22 reads" "$(pair "$P" GENEC)" "cis 22 0"
check_eq "step 31: GENET (400 bp apart, seen only across mates) is trans on 20 read pairs" "$(pair "$P" GENET)" "trans 0 20"
check_eq "step 31: GENEF (50 kb apart) is unknown" "$(pair "$P" GENEF)" "unknown 0 0"
check "step 31: the log counts them" has 'phase from the reads: 1 trans, 1 cis, 1 unknown$' "$LOG"
check "step 31: the per-variant table has PHASE after GT" has $'\tGT\tPHASE\tMAX_AF\t' "$(head -1 "$CH" 2>/dev/null)"
check_eq "step 31: each variant's PHASE" \
  "$(for p in 1000 1003 5000 5400 20000 70000; do printf '%s ' "$(cell "$CH" "$p" PHASE)"; done)" \
  "cis cis trans trans unknown unknown "

# --- control: the deletion on the other copy -----------------------------------
make_sample CHPS --swap
run_step 31-slivar.sh CHPS
check_step_exit 31-slivar.sh
check_eq "control CHPS: GENEC with the deletion on the other copy reads trans" \
  "$(pair "${GENOME_DIR}/CHPS/slivar/CHPS_compound_het_pairs.tsv" GENEC)" "trans 0 22"

# --- the fixture sample of case 42 and MANEQ without a BAM --------------------
H="${GENOME_DIR}/${SAMPLE}/slivar/${SAMPLE}_compound_het_pairs.tsv"
check_eq "${SAMPLE} (case 42): step 31 wrote the pairs table" "$(head -1 "$H" 2>/dev/null)" "$COLS"
check_eq "${SAMPLE}: every pair is trans, cis or unknown" \
  "$(awk -F'\t' 'NR > 1 && $10 !~ /^(trans|cis|unknown)$/' "$H" 2>/dev/null | wc -l | tr -d ' ')" 0
M="${GENOME_DIR}/MANEQ/slivar/MANEQ_compound_het_pairs.tsv"
if [ -s "$M" ]; then
  check_eq "MANEQ (no BAM): the GENEH pair is unknown, with the reason" \
    "$(awk -F'\t' '$1 == "GENEH" {print $10 "|" $13}' "$M")" "unknown|no BAM (aligned/MANEQ_sorted.bam with its .bai)"
else
  fail "${M} is missing: case clinical-filter-mane-and-quality-1-tiers runs step 31 on MANEQ"
fi

# --- Nextflow COMPHET_PHASE on step 31's compound-het VCF and the BAM ---------
command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
O="${GENOME_DIR}/nf-comphet"
rm -rf "$O"
NP="${CASE_TMP}/project"
rm -rf "$NP"
mkdir -p "$NP" "$O"
cp -R "${REPO}/bin" "${NP}/bin"
V="${GENOME_DIR}/CHPH/slivar/CHPH_compound_hets.vcf.gz"
B="${GENOME_DIR}/CHPH/aligned/CHPH_sorted.bam"
cat > "${NP}/main.nf" <<NF
include { COMPHET_PHASE } from '${REPO}/modules/local/slivar/main'

workflow {
    COMPHET_PHASE(channel.of(
        [[id: 'CHPH'], file('${V}'), file('${B}'), file('${B}.bai')],
        [[id: 'NOBAM'], file('${V}'), [], []]))
}
NF
cat > "${NP}/nextflow.config" <<NFCONF
params.outdir = '${O}/results'
params.publish_dir_mode = 'copy'
docker.enabled = true
docker.runOptions = '--user $(id -u):$(id -g) --network none'
NFCONF
echo "+ nextflow run (COMPHET_PHASE on CHPH, with and without the BAM)"
(
  cd "$CASE_TMP" && nextflow run "${NP}/main.nf" -ansi-log false \
    -c "${REPO}/conf/containers.config" -work-dir "${CASE_TMP}/work"
) 2>&1 | tee "$STEP_LOG"
NF_RC=${PIPESTATUS[0]}
check_eq "nextflow run exits 0" "$NF_RC" 0
[ "$NF_RC" = 0 ] || grep -E 'Error executing process|Command exit status|Command error' -A 8 "${CASE_TMP}/.nextflow.log" | head -60
NFP="${O}/results/CHPH/slivar/CHPH_compound_het_pairs.tsv"
cat "$NFP" 2>/dev/null
check "COMPHET_PHASE writes the pairs table step 31 writes, byte for byte" cmp -s "$P" "$NFP"
check_eq "COMPHET_PHASE: GENET is trans" "$(pair "$NFP" GENET)" "trans 0 20"
check_eq "COMPHET_PHASE without a BAM: every pair unknown, with the reason" \
  "$(awk -F'\t' 'NR > 1 {print $10 "|" $13}' "${O}/results/NOBAM/slivar/NOBAM_compound_het_pairs.tsv" 2>/dev/null | sort -u)" \
  "unknown|no BAM for this sample in the samplesheet"

finish
