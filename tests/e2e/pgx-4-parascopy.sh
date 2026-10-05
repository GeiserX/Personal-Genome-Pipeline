#!/usr/bin/env bash
# Step 35 (Parascopy, opt-in) gives HG002's SMN1/SMN2 locus a copy-number
# line with its quality, on the fixture's chr5 slice. Its data comes from
# `setup.sh --parascopy-data`. Three stand-ins for what the fixture lacks, for
# this test only:
#   - the SMN1 model also lists a copy of one of its region groups on chr11,
#     which the fixture reference does not have, and Parascopy refuses a model
#     naming a contig the reference lacks: the case gives step 35 a copy of
#     the reference with an N-filled chr11 of the real length (REF_FASTA);
#   - so the BAM must name chr11 too: the case uses the fixture's GIAB slice
#     BAM, whose header has every GRCh38 contig, as sample HG002par;
#   - the BAM holds reads only in the slices, so the background depth comes
#     from 100 bp windows over the chr20 slice (PARASCOPY_DEPTH_BED, which
#     makes step 35 add --no-gc) instead of Parascopy's genome-wide windows.
# The copy numbers on such a slice are not a measurement of HG002; the case
# checks that the step runs end to end and reports a quality.
. "$(dirname "$0")/lib.sh"

echo "+ scripts/setup.sh --parascopy-data"
"${REPO}/scripts/setup.sh" --parascopy-data "$GENOME_DIR"
check "setup.sh --parascopy-data installed the homology table" \
  test -s "${GENOME_DIR}/reference/parascopy-${PARASCOPY_DATA_VERSION}/homology_table/GRCh38.bed.gz"

P=HG002par
mkdir -p "${GENOME_DIR}/${P}/aligned"
cp "${FIXTURE_DIR}/${SAMPLE}_slice.bam" "${GENOME_DIR}/${P}/aligned/${P}_sorted.bam"
cp "${FIXTURE_DIR}/${SAMPLE}_slice.bam.bai" "${GENOME_DIR}/${P}/aligned/${P}_sorted.bam.bai"
REF="${GENOME_DIR}/reference/parascopy_test_ref.fasta"
python3 - "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta" "$REF" <<'PY'
import shutil, sys
with open(sys.argv[2], "wb") as out, open(sys.argv[1], "rb") as src:
    shutil.copyfileobj(src, out)
    out.write(b">chr11\n")
    left = 135086622            # GRCh38 chr11, as the homology table lists it
    while left:
        n = min(left, 60)
        out.write(b"N" * n + b"\n")
        left -= n
PY
sam faidx "reference/$(basename "$REF")"
check "the stand-in reference has chr11" grep -q $'^chr11\t135086622\t' "${REF}.fai"

BED="${GENOME_DIR}/${P}/parascopy_windows.bed"
awk 'BEGIN { for (s = 10050000; s < 10450000; s += 100) printf "chr20\t%d\t%d\n", s, s + 100 }' > "$BED"
REF_FASTA="$REF" PARASCOPY_DEPTH_BED="$BED" run_step 35-paralogs.sh "$P"
check_step_exit 35-paralogs.sh

TSV="${GENOME_DIR}/${P}/paralogs/${P}_smn_copy_number.tsv"
cat "$TSV" 2>/dev/null
check_ge "SMN1/SMN2 rows on chr5" "$(awk -F'\t' 'NR > 1 && $1 == "chr5"' "$TSV" 2>/dev/null | wc -l | tr -d ' ')" 1
# shellcheck disable=SC2016  # awk fields
check "a row with an aggregate copy number and its quality" \
  awk -F'\t' 'NR > 1 && $1 == "chr5" && $6 != "" && $6 != "*" && $7 ~ /^[0-9.]+$/ {found = 1} END {exit !found}' "$TSV"
rm -f "$REF" "${REF}.fai"

finish
