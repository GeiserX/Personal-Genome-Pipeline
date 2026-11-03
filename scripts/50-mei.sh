#!/usr/bin/env bash
# 50-mei.sh — Mobile element insertions (Alu, L1, SVA) with MEGAnE
# Usage: ./scripts/50-mei.sh <sample_name> [male|female]
#
# Opt-in. MEGAnE (Kojima et al. 2023, Nat Genet, PMID 37169872) finds the
# copies of Alu, L1, SVA and HERV elements a genome carries and the reference
# does not, from read pairs with one end in the element and from reads that
# cross its edge, and genotypes each insertion. bin/mei_report.py then lists
# the insertions it calls and reads a short catalog of insertions known to
# cause disease (assets/mei_sites.tsv: the FKTN SVA of Fukuyama muscular
# dystrophy, the TAF1 SVA of X-linked dystonia-parkinsonism, the BBS1 SVA of
# Bardet-Biedl syndrome) as seen, not seen, or not callable (too few reads
# there to tell). With AnnotSV's data in GENOME_DIR (annotsv_annotations/, as
# step 05 reads it), AnnotSV annotates the insertions the table lists: the
# genes each falls in, where in them, the frequency of a benign insertion at
# that place, and AnnotSV's class (NA for every insertion: AnnotSV 3.5 scores
# deletions and duplications only). The same commands as the MEI_DEPTH, MEI
# and MEI_REPORT Nextflow processes (MEI_REPORT runs AnnotSV too).
#
# Input:  aligned/<sample>_sorted.bam with its .bai (ALIGN_DIR as for the other
#         BAM steps), paired-end short reads; the reference; MEGAnE's k-mer
#         set of the reference (setup.sh --megane-kmers:
#         reference/megane/<reference name>.mk and .mi; MEGANE_KMERS overrides);
#         optional: AnnotSV's data, annotsv_annotations/ (setup.sh, as for step 05)
# Output: mei/<sample>_megane.vcf        MEGAnE's genotyped insertions, as it wrote them
#         mei/<sample>_mei.tsv           one row per Alu, L1 or SVA insertion MEGAnE identified,
#                                        with at least one copy, and AnnotSV's genes, location,
#                                        benign frequency and class ("." without AnnotSV's data)
#         mei/<sample>_mei_annotsv.tsv   AnnotSV's table of those insertions (with its data only)
#         mei/<sample>_mei_sites.tsv     the catalog with seen / not seen / not callable
#         mei/<sample>_mei_depth.tsv     samtools bedcov over the catalog windows
#         mei/<sample>_megane.log        MEGAnE's log
# MEGAnE reads the whole BAM; about an hour for a 30x genome on 4 threads, in
# its authors' measure. It runs with 16 GB of memory, the figure the MEI
# process asks for too (conf/base.config).
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> [male|female]}
SEX=${2:-}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image MEGANE_IMAGE SAMTOOLS_IMAGE PYTHON_IMAGE
case "$SEX" in
  male|female) ;;
  "") SEX=unknown ;;
  *) echo "ERROR: sex must be 'male' or 'female' (or left out), got '${SEX}'." >&2; exit 2 ;;
esac

S="${GENOME_DIR}/${SAMPLE}"
BAM="${S}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam"
n=$(basename "$REF_FASTA"); n=${n%.gz}; n=${n%.*}
MK=${MEGANE_KMERS:-${GENOME_DIR}/reference/megane/${n}.mk}
case "$MK" in /*) ;; *) MK="${GENOME_DIR}/${MK}" ;; esac
SITES="${PGP_ROOT}/assets/mei_sites.tsv"
ANNOTATIONS_DIR="${GENOME_DIR}/annotsv_annotations"
OUTDIR="${S}/mei"
P="${OUTDIR}/${SAMPLE}"
WORK="${OUTDIR}/megane_work"

echo "=== Mobile element insertions (MEGAnE): ${SAMPLE} ==="
for f in "$BAM" "${BAM}.bai" "$REF_FASTA" "${REF_FASTA}.fai" "$SITES"; do
  [ -f "$f" ] || { echo "ERROR: File not found: ${f}" >&2; exit 1; }
done
for f in "$MK" "${MK%.mk}.mi"; do
  [ -s "$f" ] || { echo "ERROR: MEGAnE's k-mer set not found: ${f} (run: ./scripts/setup.sh --megane-kmers ${GENOME_DIR})" >&2; exit 1; }
done
# AnnotSV runs when its data is there, and then with the two files step 05
# requires (setup.sh --annotsv-data): without them its genes and benign
# frequencies are not the ones docs/05-annotsv.md describes.
ANNOTSV=false
if [ -d "${ANNOTATIONS_DIR}/Annotations_Human/Genes/GRCh38" ]; then
  for f in "${ANNOTATIONS_DIR}/mane_select_refseq.txt" \
           "${ANNOTATIONS_DIR}/Annotations_Human/SVincludedInFt/BenignSV/GRCh38/gnomad-sv.installed"; do
    if [ ! -s "$f" ]; then
      echo "ERROR: ${f} not found: AnnotSV's data lacks the MANE Select list or gnomAD-SV's common SVs." >&2
      echo "  Run ./scripts/setup.sh --annotsv-data ${GENOME_DIR} to add them, or see docs/05-annotsv.md." >&2
      exit 1
    fi
  done
  require_image ANNOTSV_IMAGE
  ANNOTSV=true
fi
echo "BAM: ${BAM}"
echo "k-mer set: ${MK}"
echo "Sex: ${SEX}"
if $ANNOTSV; then
  echo "AnnotSV: ${ANNOTATIONS_DIR}/"
else
  echo "AnnotSV: no ${ANNOTATIONS_DIR#"$GENOME_DIR"/}/ (./scripts/setup.sh), so the insertions get no gene or frequency (\".\")."
fi
mkdir -p "$OUTDIR"
rm -rf "$WORK"
rm -f "${P}_megane.vcf" "${P}_mei.tsv" "${P}_mei_sites.tsv" "${P}_mei_depth.tsv" "${P}_megane.log" \
  "${P}_mei_mainchr.txt" "${P}_mei_sites.bed" "${P}_mei_annotsv.tsv" "${P}_mei_annotsv_input.vcf"
mkdir -p "$WORK"
B=$(cpath "$BAM") O=$(cpath "$P") WC=$(cpath "$WORK")

# The catalog windows as a BED, on the contigs the BAM has (samtools bedcov
# stops on one it does not have), and the chromosomes MEGAnE reads: chr1 to
# chr22, chrX and chrY as the BAM names them. MEGAnE's own list holds every
# GRCh38 contig and stops when one is missing from the BAM.
echo "[1/4] Depth over the catalog windows, and the chromosomes to read (samtools)..."
# shellcheck disable=SC2016  # $1 and $2 belong to the inner bash
run_in -v "${SITES}:/pgp-sites.tsv:ro" --cpus 1 --memory 2g "$SAMTOOLS_IMAGE" bash -euo pipefail -c '
    samtools view -H "$1" | awk -F "\t" "\$1 == \"@SQ\" {sub(/^SN:/, \"\", \$2); print \$2}" > "$2_bam_contigs.txt"
    grep -E "^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$" "$2_bam_contigs.txt" > "$2_mei_mainchr.txt" || true
    awk -F "\t" "NR == FNR {c[\$1]; next} !/^#/ && \$1 != \"gene\" && (\$2 in c) {print \$2 \"\t\" \$3 - 1 \"\t\" \$4 \"\t\" \$1}" \
      "$2_bam_contigs.txt" /pgp-sites.tsv > "$2_mei_sites.bed"
    if [ -s "$2_mei_sites.bed" ]; then samtools bedcov -c "$2_mei_sites.bed" "$1" > "$2_mei_depth.tsv"; else : > "$2_mei_depth.tsv"; fi
    rm -f "$2_bam_contigs.txt" "$2_mei_sites.bed"' _ "$B" "$O"
[ -s "${P}_mei_mainchr.txt" ] || { echo "ERROR: the BAM has none of chr1-chr22, chrX, chrY: nothing for MEGAnE to read." >&2; exit 1; }

# MEGAnE refuses a .fai that is not newer than its FASTA, which a copied or
# unpacked reference can have: it reads the reference through a link with a
# fresh copy of the .fai beside it. -only_ins: insertions only, not the
# reference copies the sample lacks. Its -mainchr comes after the wrapper's
# own, so it wins. MEGAnE's single-sample VCF is MEI_final_percentile_genotyped.vcf;
# when it finds too few insertions to set its filter thresholds (a BAM of one
# chromosome), it writes MEI_final_failed_genotyped.vcf instead, with fixed ones.
echo "[2/4] MEGAnE call_genotype_38 (${THREADS} threads)..."
# MEGAnE stops when chrX is missing from the BAM unless told the sample has no
# sex chromosomes (-no_sex_chr: every chromosome read as an autosome); a BAM cut
# to some autosomes needs it.
NOSEX=()
grep -qE '^(chr)?X$' "${P}_mei_mainchr.txt" || NOSEX=(-no_sex_chr)
ln -s "$REF_FASTA_C" "${WORK}/ref.fa"
cp "${REF_FASTA}.fai" "${WORK}/ref.fa.fai"
run_in -w "$WC" --cpus "$THREADS" --memory 16g "$MEGANE_IMAGE" call_genotype_38 -i "$B" -fa "${WC}/ref.fa" -mk "$(cpath "$MK")" \
  -mainchr "${O}_mei_mainchr.txt" -sex "$SEX" ${NOSEX[@]+"${NOSEX[@]}"} -sample_name "$SAMPLE" -p "$THREADS" -only_ins -outdir "${WC}/out"
VCF="${WORK}/out/MEI_final_percentile_genotyped.vcf"
if [ ! -s "$VCF" ] && [ -s "${WORK}/out/MEI_final_failed_genotyped.vcf" ]; then
  VCF="${WORK}/out/MEI_final_failed_genotyped.vcf"
  echo "WARNING: MEGAnE found too few insertions to set its filter thresholds and used fixed ones ($(basename "$VCF"))."
fi
cat "${WORK}/out"/*.log > "${P}_megane.log" 2>/dev/null || true
[ -s "$VCF" ] || { echo "ERROR: MEGAnE wrote no MEI_final_percentile_genotyped.vcf (log: ${P}_megane.log)" >&2; exit 1; }
cp "$VCF" "${P}_megane.vcf"

# The insertions the calls table lists, one <INS:ME> record each with an ID
# the table is joined on. AnnotSV gives an insertion its genes and their
# location, and the benign set's frequency of an insertion there; it scores
# no insertion, so their class is NA. It writes no table when it had nothing
# to annotate. It builds sorted copies of its data inside the annotations
# directory the first time it runs, as in step 05, so that is writable.
ANNOTSV_ARGS=()
if $ANNOTSV; then
  echo "[3/4] AnnotSV on the insertions..."
  run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 2g "$PYTHON_IMAGE" \
    python3 /pgp-bin/mei_report.py --sample "$SAMPLE" --vcf "${O}_megane.vcf" --annotsv-input "${O}_mei_annotsv_input.vcf"
  AC=$(cpath "$ANNOTATIONS_DIR")
  run_in --rw "$ANNOTATIONS_DIR" --cpus 4 --memory 8g "$ANNOTSV_IMAGE" \
    AnnotSV -SVinputFile "${O}_mei_annotsv_input.vcf" -outputFile "${O}_mei_annotsv.tsv" -genomeBuild GRCh38 \
      -annotationMode both -annotationsDir "$AC" -txFile "${AC}/mane_select_refseq.txt"
  [ -s "${P}_mei_annotsv.tsv" ] || : > "${P}_mei_annotsv.tsv"
  rm -f "${P}_mei_annotsv_input.vcf"
  ANNOTSV_ARGS=(--annotsv "${O}_mei_annotsv.tsv")
else
  echo "[3/4] AnnotSV: skipped (no annotation data)."
fi

echo "[4/4] Calls and the catalog (bin/mei_report.py)..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" -v "${SITES}:/pgp-sites.tsv:ro" --cpus 1 --memory 2g "$PYTHON_IMAGE" \
  python3 /pgp-bin/mei_report.py --sample "$SAMPLE" --vcf "${O}_megane.vcf" --sites /pgp-sites.tsv \
    --depth "${O}_mei_depth.tsv" ${ANNOTSV_ARGS[@]+"${ANNOTSV_ARGS[@]}"} --out-calls "${O}_mei.tsv" --out-sites "${O}_mei_sites.tsv"
rm -rf "$WORK" "${P}_mei_mainchr.txt"

for f in "${P}_mei.tsv" "${P}_mei_sites.tsv"; do
  [ -s "$f" ] || { echo "ERROR: not written: ${f}" >&2; exit 1; }
done
echo ""
echo "Insertions: ${P}_mei.tsv ($(($(wc -l < "${P}_mei.tsv") - 1)) called)"
echo "Known disease insertions: ${P}_mei_sites.tsv"
if $ANNOTSV; then echo "AnnotSV's table of the insertions: ${P}_mei_annotsv.tsv"; fi
