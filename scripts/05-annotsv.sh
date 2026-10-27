#!/usr/bin/env bash
# AnnotSV — Annotate structural variants with ACMG pathogenicity classification
# Input: step 15's duphold-filtered calls (duphold/<sample>_sv_filtered.vcf.gz),
#   the input of the Nextflow ANNOTSV; without them, Manta's diploidSV.vcf.gz,
#   with a note. SV_VCF overrides both. Before AnnotSV, bin/sv_svlen_fix.py
#   gives sequence-resolved records their net length change, and, unless
#   SV_VCF is set, bin/cnvpytor_to_vcf.py adds the Delly (step 19) and
#   CNVpytor (step 18) deletions and duplications no Manta call covers. The
#   sample's sex (male or female, second argument, when known) sets the
#   normal copies on chrX and chrY for the read-depth check of long Delly calls.
# Output: *_sv_annotated.tsv (ACMG class 1-5 for each SV),
#   *_sv_class4_5.tsv (its class 4 and 5 rows), *_sv_class4_5_checked.tsv
#   (each of those with its caller, an exon-edge flag for a deletion, its
#   benign-set frequency, its read-depth support and the class after the
#   rule table assets/sv_rule_regions.tsv, bin/annotsv_postfilter.py) and
#   *_annotsv_input.vcf
#   (what AnnotSV read, each record with its caller in INFO SVCALLER)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> [male|female]}
SEX=${2:-}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
case "$SEX" in
  ''|male|female) ;;
  *) echo "ERROR: sex must be 'male' or 'female' (or left out), got '${SEX}'." >&2; exit 2 ;;
esac
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
# Allow SV_VCF env var override (e.g., for Sniffles2 long-read SVs)
if [ -n "${SV_VCF:-}" ]; then
  MANTA_VCF="$SV_VCF"
else
  MANTA_VCF="${SAMPLE_DIR}/manta/results/variants/diploidSV.vcf.gz"
  # Fall back to manta2/ if a second Manta run was used
  [ ! -f "$MANTA_VCF" ] && MANTA_VCF="${SAMPLE_DIR}/manta2/results/variants/diploidSV.vcf.gz"
  # Step 15's depth-filtered calls, unless Manta's calls are newer than them
  FILTERED="${SAMPLE_DIR}/duphold/${SAMPLE}_sv_filtered.vcf.gz"
  if [ -f "$FILTERED" ] && [ -f "${FILTERED}.tbi" ] && ! [ "$MANTA_VCF" -nt "$FILTERED" ]; then
    MANTA_VCF="$FILTERED"
  elif [ -f "$FILTERED" ]; then
    echo "NOTE: ${FILTERED} is unfinished or older than Manta's calls; annotating Manta's calls."
    echo "  Run scripts/15-duphold.sh ${SAMPLE} first for the depth-filtered list the pipeline annotates."
  else
    echo "NOTE: no duphold-filtered calls (scripts/15-duphold.sh); annotating Manta's calls without the depth filter."
  fi
fi
OUTPUT_DIR="${SAMPLE_DIR}/annotsv"
# The AnnotSV image holds code only; its annotation data is a separate
# download that setup.sh puts here.
ANNOTATIONS_DIR="${GENOME_DIR}/annotsv_annotations"

echo "=== AnnotSV: ${SAMPLE} ==="
echo "Input: ${MANTA_VCF}"
echo "Annotations: ${ANNOTATIONS_DIR}/"
echo "Output: ${OUTPUT_DIR}/"

if [ ! -d "${ANNOTATIONS_DIR}/Annotations_Human/Genes/GRCh38" ]; then
  echo "ERROR: AnnotSV annotation data not found in ${ANNOTATIONS_DIR}/Annotations_Human/" >&2
  echo "  Run ./scripts/setup.sh ${GENOME_DIR} to download it (~5.3 GB), or see docs/05-annotsv.md." >&2
  exit 1
fi

# setup.sh adds two files to the bundle: the MANE Select transcripts, which
# AnnotSV reports in the Tx, Frameshift, Location and exon columns instead of
# the transcript with the longest overlap (often a predicted XM_ one), and
# gnomAD-SV's common SVs in its benign set. Without them the ACMG classes are
# wrong for the reasons in docs/05-annotsv.md, so the step stops.
MANE_TX="${ANNOTATIONS_DIR}/mane_select_refseq.txt"
GNOMAD_SV_STAMP="${ANNOTATIONS_DIR}/Annotations_Human/SVincludedInFt/BenignSV/GRCh38/gnomad-sv.installed"
for f in "$MANE_TX" "$GNOMAD_SV_STAMP"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: ${f} not found: AnnotSV's data lacks the MANE Select list or gnomAD-SV's common SVs." >&2
    echo "  Run ./scripts/setup.sh --annotsv-data ${GENOME_DIR} to add them, or see docs/05-annotsv.md." >&2
    exit 1
  fi
done

if [ ! -f "$MANTA_VCF" ]; then
  echo "ERROR: SV VCF not found: ${MANTA_VCF}" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# AnnotSV classifies only the records of the VCF it is given. The other
# callers' deletions and duplications that no Manta call covers (reciprocal
# overlap of half the longer call, same type) join Manta's: Delly's PASS
# calls, and CNVpytor's with e-val below 0.01 and q0 below 0.5. A Delly call
# longer than 1,000 bp also needs its read-depth copy number (RDCN) to show
# the change, judged against 1 copy on a male's chrX and chrY outside the
# PARs and 2 elsewhere. Each file is used when it exists. With SV_VCF set,
# that file alone is annotated.
EXTRA_ARGS=()
CALLER_NAME=sv_vcf
if [ -z "${SV_VCF:-}" ]; then
  CALLER_NAME=manta
  DELLY_VCF="${SAMPLE_DIR}/delly/${SAMPLE}_sv.vcf.gz"
  CNVPYTOR_VCF="${SAMPLE_DIR}/cnvpytor/${SAMPLE}_cnvs.vcf.gz"
  if [ -f "$DELLY_VCF" ]; then
    EXTRA_ARGS+=(--delly "$(cpath "$DELLY_VCF")" --sex "$SEX" --par-bed /pgp/par_grch38.bed)
    echo "Delly: long calls need read depth (RDCN) to agree; sex ${SEX:-not given (1 or 2 copies of chrX and chrY outside the PARs both count as normal)}."
  else
    echo "Delly: no ${DELLY_VCF#"$GENOME_DIR"/} (scripts/19-delly.sh), so no Delly calls are added."
  fi
  if [ -f "$CNVPYTOR_VCF" ]; then
    EXTRA_ARGS+=(--cnvpytor "$(cpath "$CNVPYTOR_VCF")")
  else
    echo "CNVpytor: no ${CNVPYTOR_VCF#"$GENOME_DIR"/} (scripts/18-cnvpytor.sh), so no CNVpytor calls are added."
  fi
fi

# The input AnnotSV reads. bin/sv_svlen_fix.py first: AnnotSV measures a
# sequence-resolved SV over POS-1..END, so a deletion that also inserts bases
# (REF 137 bp, ALT 2 bp) reads as 137 deleted coding bases, a frameshift
# (criterion 2E-1), although the net change, -135, keeps the frame.
INPUT_VCF="${OUTPUT_DIR}/${SAMPLE}_annotsv_input.vcf"
SV_VCF_C=$(cpath "$MANTA_VCF")
# shellcheck disable=SC2016  # the $ expressions expand inside the container
run_in --cpus 1 --memory 2g -v "${PGP_ROOT}/bin:/pgp-bin:ro" \
  -v "${PGP_ROOT}/assets/par_grch38.bed:/pgp/par_grch38.bed:ro" "${ANNOTSV_IMAGE}" bash -c '
  set -euo pipefail
  python3 /pgp-bin/sv_svlen_fix.py "$1" | python3 /pgp-bin/cnvpytor_to_vcf.py --base - --base-name "$2" "${@:3}"
' _ "$SV_VCF_C" "$CALLER_NAME" ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"} > "${INPUT_VCF}.part"
mv "${INPUT_VCF}.part" "$INPUT_VCF"

# AnnotSV builds sorted copies of its annotation files inside the annotations
# directory the first time it runs, so that directory is writable here.
# -bracketedAltMode DEL,DUP,INS: AnnotSV would otherwise read a breakend pair
# on one chromosome up to 10 Mb apart as a whole-span inversion, which is what
# step 04 left as breakends on purpose (a lone inversion junction of 1 Mb or
# more); it annotates them as breakends.
TSV="${OUTPUT_DIR}/${SAMPLE}_sv_annotated.tsv"
run_in --rw "$ANNOTATIONS_DIR" --cpus 4 --memory 8g \
  "${ANNOTSV_IMAGE}" \
  AnnotSV \
    -SVinputFile "$(cpath "$INPUT_VCF")" \
    -outputFile "$(cpath "$TSV")" \
    -genomeBuild GRCh38 \
    -annotationMode both \
    -annotationsDir /genome/annotsv_annotations \
    -txFile /genome/annotsv_annotations/mane_select_refseq.txt \
    -bracketedAltMode DEL,DUP,INS

# The class 4 and 5 SVs, one "full" row each, found by column name: AnnotSV
# drops the columns of annotation files a bundle lacks, so a column number
# points elsewhere on another bundle.
CLASS45="${OUTPUT_DIR}/${SAMPLE}_sv_class4_5.tsv"
if [ -s "$TSV" ]; then
  awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i
                         if (!("ACMG_class" in c) || !("Annotation_mode" in c)) { print "no ACMG_class or Annotation_mode column" > "/dev/stderr"; exit 2 }
                         print; next }
               $c["Annotation_mode"] == "full" && ($c["ACMG_class"] == "4" || $c["ACMG_class"] == "5")' \
    "$TSV" > "${CLASS45}.part"
  mv "${CLASS45}.part" "$CLASS45"
  N45=$(($(wc -l < "$CLASS45") - 1))
else
  echo "WARNING: AnnotSV wrote no ${TSV} (no SV it could annotate)."
  : > "$CLASS45"
  N45=0
fi

# The checks a class 4 or 5 call needs (bin/annotsv_postfilter.py): its
# caller, whether a deletion sits on the exon edges of one transcript of
# AnnotSV's genes file (a processed transcript aligned across an intron), the
# benign-set frequency of its type, and whether its read depth shows the
# change. The rule table moves to class 3 a call inside a coding tandem
# repeat, the IG and TR loci or the MHC, one the depth contradicts, and one
# whose gene symbol is a pseudogene's alias. How many callers saw it comes
# from step 22's consensus, which the reports join when they are written.
CHECKED="${OUTPUT_DIR}/${SAMPLE}_sv_class4_5_checked.tsv"
GENES_BED="${ANNOTATIONS_DIR}/Annotations_Human/Genes/GRCh38/genes.RefSeq.sorted.bed"
GENES_ARGS=()
if [ -s "$GENES_BED" ]; then
  GENES_ARGS=(--genes "$(cpath "$GENES_BED")")
else
  echo "Exon edges: no ${GENES_BED#"$GENOME_DIR"/} in the annotation data, so no deletion is checked against them."
fi
run_in --cpus 1 --memory 2g -v "${PGP_ROOT}/bin:/pgp-bin:ro" \
  -v "${PGP_ROOT}/assets/sv_rule_regions.tsv:/pgp/sv_rule_regions.tsv:ro" "${ANNOTSV_IMAGE}" \
  python3 /pgp-bin/annotsv_postfilter.py --class45 "$(cpath "$CLASS45")" ${GENES_ARGS[@]+"${GENES_ARGS[@]}"} \
    --rules /pgp/sv_rule_regions.tsv --out "$(cpath "$CHECKED")"

echo "=== AnnotSV complete ==="
echo "Input: ${INPUT_VCF}"
echo "Results: ${TSV}"
echo "Class 4-5 (pathogenic or likely pathogenic): ${N45} SVs in ${CLASS45}"
echo "Their caller, exon-edge, population, depth and rule-table checks (ACMG_class after them): ${CHECKED}"
