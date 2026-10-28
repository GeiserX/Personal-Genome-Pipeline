#!/usr/bin/env bash
# Stranger — Annotate ExpansionHunter STR VCF with clinical pathogenicity status
# Adds STR_STATUS (normal/pre_mutation/full_mutation), disease name, OMIM number,
# inheritance mode, and normal/pathogenic repeat ranges to each locus.
# Input:  ExpansionHunter VCF produced by step 09 (and the sample's BAM, for RFC1)
# Output: Annotated VCF in $GENOME_DIR/<sample>/expansion_hunter/, and
#         <sample>_rfc1_motif.tsv when Stranger flags RFC1 (bin/rfc1_motif.py),
#         <sample>_str_motifs.tsv when it flags another locus whose motif the
#         catalog names (the FAME loci, SCA31, SCA37)
# Catalog: STRANGER_REPEATS=FILE, STRANGER_REPEATS=bundled for Stranger's
# bundled GRCh38 one; by default the catalog step 09 ran on
# (<sample>_eh_catalog.json, the STRchive one) when it carries Stranger's
# thresholds, else the bundled GRCh38 one.
#
# This step exits cleanly (exit 0) when the ExpansionHunter VCF does not exist
# so the pipeline can treat it as optional without failing run-all.sh.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
OUTPUT_DIR="${SAMPLE_DIR}/expansion_hunter"
EH_VCF="${OUTPUT_DIR}/${SAMPLE}_eh.vcf"
OUT_VCF="${OUTPUT_DIR}/${SAMPLE}_eh_stranger.vcf"
MOTIF_TSV="${OUTPUT_DIR}/${SAMPLE}_rfc1_motif.tsv"
OTHER_TSV="${OUTPUT_DIR}/${SAMPLE}_str_motifs.tsv"
EH_CAT="${OUTPUT_DIR}/${SAMPLE}_eh_catalog.json"
ALIGN_DIR=${ALIGN_DIR:-aligned}
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"

echo "=== Stranger STR Annotation: ${SAMPLE} ==="

# Exit cleanly if the ExpansionHunter VCF does not exist.
# Step 09 may have been skipped or run separately — this is not an error.
if [ ! -f "$EH_VCF" ]; then
  echo "INFO: ExpansionHunter VCF not found: ${EH_VCF}"
  echo "INFO: Run scripts/09-expansion-hunter.sh first, or skip this annotation step."
  exit 0
fi

# The catalog Stranger grades with: STRANGER_REPEATS, else the one step 09
# ran on when it carries Stranger's thresholds (NormalMax), else none (the
# bundled GRCh38 one below).
REPEATS=""
case "${STRANGER_REPEATS:-}" in
  bundled) ;;
  "") if [ -s "$EH_CAT" ] && grep -q '"NormalMax"' "$EH_CAT"; then REPEATS=$EH_CAT; fi ;;
  *) REPEATS=$STRANGER_REPEATS ;;
esac
# The motif count reads its motifs and loci from that catalog when it names them.
CAT_ARGS=()
if [ -n "$REPEATS" ] && [ -f "$REPEATS" ] && grep -q '"PathogenicMotifs"' "$REPEATS"; then
  CAT_ARGS=(--catalog "$(cpath "$REPEATS")")
fi

stranger_run() {
  echo "Input:  ${EH_VCF}"
  # Stranger annotates each STR locus with clinical pathogenicity thresholds from a
  # repeat catalog. Its default is the GRCh37 catalog it bundles, while step 09
  # calls the GRCh38 loci, so the bundled GRCh38 catalog is passed explicitly
  # (found next to the installed package). A custom catalog (TSV or JSON) can be
  # supplied via STRANGER_REPEATS.
  if [ -n "$REPEATS" ]; then
    if [ ! -f "$REPEATS" ]; then
      echo "ERROR: STRANGER_REPEATS catalog not found: ${REPEATS}" >&2
      exit 1
    fi
    REPEATS_C=$(cpath "$REPEATS") || exit 2
    if [ "$REPEATS" = "$EH_CAT" ]; then
      echo "Repeat catalog: the one step 09 ran on, ${EH_CAT} ($(grep -o '"SourceDisplay": "[^"]*"' "$EH_CAT" | head -n 1 | cut -d'"' -f4))"
    else
      echo "Repeat catalog: ${REPEATS} (custom)"
    fi
    atomic_out "$OUT_VCF" run_in \
      --cpus 1 --memory 1g \
      "${STRANGER_IMAGE}" \
      stranger \
        --repeats-file "$REPEATS_C" \
        "/genome/${SAMPLE}/expansion_hunter/${SAMPLE}_eh.vcf"
  else
    echo "Repeat catalog: bundled GRCh38 clinical catalog (variant_catalog_grch38.json)"
    # shellcheck disable=SC2016  # the dollar signs are for the container's sh
    atomic_out "$OUT_VCF" run_in \
      --cpus 1 --memory 1g \
      "${STRANGER_IMAGE}" \
      sh -c 'catalog=$(python -c "import os, stranger; print(os.path.join(os.path.dirname(stranger.__file__), \"resources\", \"variant_catalog_grch38.json\"))") &&
        [ -f "$catalog" ] || { echo "ERROR: no GRCh38 catalog in the Stranger image (${catalog:-not found})" >&2; exit 1; }
        exec stranger --repeats-file "$catalog" "$1"' \
      stranger "/genome/${SAMPLE}/expansion_hunter/${SAMPLE}_eh.vcf"
  fi
}

# Skip only a finished output: written through a temporary name, so an empty
# or cut-short file from a failed run is never taken for one.
STRANGER_RAN=""
if have_output "$OUT_VCF"; then
  echo "Stranger output already exists: ${OUT_VCF}"
  echo "Delete to re-run: rm ${OUT_VCF}"
else
  stranger_run
  STRANGER_RAN=1
fi

# RFC1: Stranger's catalog sizes the degenerate AARRG unit and flags a low
# count as full_mutation, so the flag says nothing until the motif is known.
# When RFC1 is flagged, bin/rfc1_motif.py counts the pentamer motif of the
# reads at the locus (soft-clipped bases and mates mapped elsewhere
# included) and writes a benign, pathogenic or uncertain verdict. It runs
# again whenever Stranger did, so its table is never older than the VCF.
RFC1_FLAG=$(awk -F'\t' '!/^#/ && $8 ~ /(^|;)REPID=RFC1(;|$)/ && $8 ~ /(^|;)STR_STATUS=[^;]*(pre|full)_mutation/ {print "yes"; exit}' "$OUT_VCF")
if [ -z "$RFC1_FLAG" ]; then
  rm -f "$MOTIF_TSV"
elif [ -z "$STRANGER_RAN" ] && [ -s "$MOTIF_TSV" ] && [ "$MOTIF_TSV" -nt "$OUT_VCF" ] \
     && ! grep -q $'\tnot_checked\t' "$MOTIF_TSV"; then
  # A not_checked table (no BAM last time) is redone, so a BAM that appeared
  # since gets counted.
  echo "RFC1 motif count already exists: ${MOTIF_TSV}"
  echo "Delete to re-run: rm ${MOTIF_TSV}"
else
  if [ -f "$BAM" ] && [ -f "${BAM}.bai" ]; then
    MOTIF_ARGS=(--bam "$(cpath "$BAM")")
  else
    MOTIF_ARGS=(--skip "no BAM (${ALIGN_DIR}/${SAMPLE}_sorted.bam with its .bai)")
  fi
  echo "RFC1 is flagged: counting the repeat motif of its reads (bin/rfc1_motif.py)"
  # pysam is in PYPGX_IMAGE.
  run_in --cpus 1 --memory 2g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYPGX_IMAGE" \
    python3 /pgp-bin/rfc1_motif.py --vcf "$(cpath "$OUT_VCF")" --out "$(cpath "${MOTIF_TSV}.part")" "${MOTIF_ARGS[@]}" ${CAT_ARGS[@]+"${CAT_ARGS[@]}"}
  mv "${MOTIF_TSV}.part" "$MOTIF_TSV"
  echo "RFC1 motif: ${MOTIF_TSV}"
fi

# The other loci whose motif decides (catalog fields PathogenicMotifs): one
# row per locus Stranger flags, with the catalog's motifs. Only with a
# catalog that names motifs, and again whenever Stranger ran; a not_checked
# table (no BAM last time) is redone, as RFC1's is.
if [ "${#CAT_ARGS[@]}" -eq 0 ]; then
  rm -f "$OTHER_TSV"
elif [ -n "$STRANGER_RAN" ] || [ ! -e "$OTHER_TSV" ] || [ "$OUT_VCF" -nt "$OTHER_TSV" ] \
     || grep -q $'\tnot_checked\t' "$OTHER_TSV"; then
  if [ -f "$BAM" ] && [ -f "${BAM}.bai" ]; then
    MOTIF_ARGS=(--bam "$(cpath "$BAM")")
  else
    MOTIF_ARGS=(--skip "no BAM (${ALIGN_DIR}/${SAMPLE}_sorted.bam with its .bai)")
  fi
  run_in --cpus 1 --memory 2g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYPGX_IMAGE" \
    python3 /pgp-bin/rfc1_motif.py --vcf "$(cpath "$OUT_VCF")" "${CAT_ARGS[@]}" --flagged \
    --out "$(cpath "${OTHER_TSV}.part")" "${MOTIF_ARGS[@]}"
  if [ -e "${OTHER_TSV}.part" ]; then
    mv "${OTHER_TSV}.part" "$OTHER_TSV"
    echo "Motif count of the other flagged motif loci: ${OTHER_TSV}"
  else
    rm -f "$OTHER_TSV"
  fi
fi

echo "=== Stranger complete ==="
echo "Annotated VCF: ${OUT_VCF}"
echo ""
echo "Each locus now has:"
echo "  STR_STATUS: normal / pre_mutation / full_mutation"
echo "  Disease name, OMIM number, inheritance mode, and repeat size ranges"
echo "  See docs/09b-stranger.md for interpretation guidance."
