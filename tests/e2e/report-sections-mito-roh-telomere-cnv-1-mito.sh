#!/usr/bin/env bash
# Step 20 annotates its chrM calls with ClinVar and gnomAD v3.1's chrM sites,
# validate-setup.sh warns without the pinned gnomAD file and fails with another
# file under its name, and the text report shows the annotation, the poly-C tract
# count, the autosomal ROH lines and the telomere read count.
#
# Red first: the earlier step-20 cases ran without the gnomAD file, so the
# table must say it is not annotated, and validate-setup.sh must warn; a
# file that is not the pinned one must fail. Then the pinned file goes in and
# both turn.
. "$(dirname "$0")/lib.sh"

G="${GENOME_DIR}/gnomad_mito/${GNOMAD_MITO_URL##*/}"
T="${GENOME_DIR}/${SAMPLE}/mito/${SAMPLE}_chrM_annotated.tsv"
# Alleles with a value in column N of the annotation table (11 gnomad_af_hom, 7 clinvar_id).
filled() { awk -F'\t' -v c="$1" '!/^#/ && $1 != "chrom" && $c != "." {n++} END {print n + 0}' "$T" 2>/dev/null; }
validate_line() { grep -E 'gnomAD v3\.1 mitochondrial sites' "$STEP_LOG" | head -n 1; }

# --- 1. without the gnomAD file -------------------------------------------------
rm -f "$G"
run_step 20-mtoolbox.sh "$SAMPLE"
check_step_exit 20-mtoolbox.sh
check "the step wrote the annotation table" test -s "$T"
head -n 3 "$T"
check "without the gnomAD file the table says it is not annotated" has $'^# gnomad\tnot annotated' "$(head -n 3 "$T")"
check_eq "alleles with a gnomAD value without the file" "$(filled 11)" 0
run_step validate-setup.sh "$SAMPLE"
check "validate-setup.sh warns without the gnomAD file ($(validate_line))" has '\[WARN\].*not found at' "$(validate_line)"

# --- 2. another file under its name ---------------------------------------------
mkdir -p "$(dirname "$G")"
printf 'not the pinned file\n' > "$G"
run_step validate-setup.sh "$SAMPLE"
check "validate-setup.sh fails on a file that is not the pinned one ($(validate_line))" has '\[FAIL\].*sha256 differs' "$(validate_line)"

# --- 3. the pinned file ---------------------------------------------------------
rm -f "$G"
curl -fsSL --retry 3 -o "${G}.part" "$GNOMAD_MITO_URL"
check_eq "the download has the pinned sha256" "$(sha256sum "${G}.part" 2>/dev/null | awk '{print $1}')" "$GNOMAD_MITO_SHA256"
mv -f "${G}.part" "$G"
run_step validate-setup.sh "$SAMPLE"
check "validate-setup.sh passes the pinned file ($(validate_line))" has '\[OK\].*present' "$(validate_line)"

run_step 20-mtoolbox.sh "$SAMPLE"
check_step_exit 20-mtoolbox.sh
head -n 3 "$T"
check "the table names the gnomAD file" has $'^# gnomad\tgnomad\\.genomes\\.v3\\.1\\.sites\\.chrM\\.vcf\\.bgz' "$(head -n 3 "$T")"
check "the table names the ClinVar file and its date" has $'^# clinvar\tclinvar_chr\\.vcf\\.gz\t[0-9]{4}-?[0-9]{2}-?[0-9]{2}' "$(head -n 3 "$T")"
check_ge "alleles with a gnomAD value" "$(filled 11)" 1
check_ge "alleles with a ClinVar record" "$(filled 7)" 1
# One ClinVar match, read back from ClinVar's own file: same id, position and alleles.
read -r POS REF ALT ID < <(awk -F'\t' '!/^#/ && $1 != "chrom" && $7 != "." {print $2, $3, $4, $7; exit}' "$T")
BACK=$(bcf query -r "chrM:${POS:-0}" -f '%ID %REF %ALT\n' clinvar/clinvar_chr.vcf.gz 2>/dev/null | awk -v i="${ID:-x}" '$1 == i')
check "ClinVar record ${ID:-none} at m.${POS:-?} is in ClinVar's file with the same alleles (${BACK:-none})" \
  has "^${ID:-x} ${REF:-x} ${ALT:-x}\$" "${BACK:-}"
awk -F'\t' '!/^#/ && ($7 != "." || $1 == "chrom")' "$T" | head -n 6

# --- 4. the text report ---------------------------------------------------------
run_step generate-report.sh "$SAMPLE"
check_step_exit generate-report.sh
TXT=$(cat "${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.txt" 2>/dev/null)
grep -A 12 '^## Mitochondrial Variants' <<< "$TXT" | head -n 20
check "the report counts heteroplasmy outside the poly-C tracts" \
  has '^  Heteroplasmic \(AF 0\.05 to 0\.95, outside the D-loop poly-C tracts\): [0-9]+$' "$TXT"
check "the report has the poly-C tract count" \
  has '^  Length heteroplasmy in the D-loop poly-C tracts \(m\.303-315, m\.16184-16193; common, usually technical\): [0-9]+$' "$TXT"
check "the report names the annotation and both sources" \
  has "^  Annotation \\(mito/${SAMPLE}_chrM_annotated\\.tsv\\): ClinVar clinvar_chr\\.vcf\\.gz \\([0-9]{4}-?[0-9]{2}-?[0-9]{2}\\); gnomAD gnomad\\.genomes\\.v3\\.1\\.sites\\.chrM\\.vcf\\.bgz\$" "$TXT"
check "the report counts ClinVar P/LP chrM calls" has '^  ClinVar pathogenic or likely pathogenic: [0-9]+$' "$TXT"
check "the ROH section says autosomes only and no population AFs" has 'No population allele frequencies were used' "$TXT"
S11="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}_roh_summary.txt"
check "step 11's summary keeps its columns" has $'^chrom\tstart\tend\tlength_bp\tlength_mb$' "$(cat "$S11" 2>/dev/null)"

# --- 5. the telomere read floor on TelomereHunter's real summary -------------------
# Step 10 on a copy of the sample (the BAM hard-linked), as bash-step-10 runs
# it. The slice has no telomeric read, so the report must warn.
TS="${SAMPLE}telfloor"
mkdir -p "${GENOME_DIR}/${TS}/aligned"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${GENOME_DIR}/${TS}/aligned/${TS}_sorted.bam"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" "${GENOME_DIR}/${TS}/aligned/${TS}_sorted.bam.bai"
run_step 10-telomere-hunter.sh "$TS"
check_step_exit 10-telomere-hunter.sh
head -n 2 "${GENOME_DIR}/${TS}/telomere/${TS}/${TS}/${TS}_summary.tsv" 2>/dev/null | cut -f1-11
TEL=$(python3 - "${REPO}/bin" "${GENOME_DIR}/${TS}" "$TS" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
import collect_summary, render_report
s = collect_summary.collect(sys.argv[3], sys.argv[2])
d = s["sections"]["telomere"]["data"]
print(f"reads={d.get('intratelomeric_reads')} depth={d.get('depth')} floor={d.get('read_floor')} low={d.get('low_reads')}")
print(render_report.text_report(s))
PY
)
grep -E '^reads=|Intratelomeric|WARNING' <<< "$TEL"
check "the summary reads TelomereHunter's intratel_reads, total_reads and read_length" \
  has '^reads=[0-9]+ depth=[0-9.]+ floor=[0-9]+ low=(True|False)$' "$TEL"
check "the slice's telomeric reads are under the floor, and the report warns" \
  has '^  WARNING: too few telomeric reads for a stable estimate' "$TEL"
in_genome "$BCFTOOLS_IMAGE" rm -rf "$TS"

finish
