#!/usr/bin/env bash
# Step 25 and the PRS processes on HG002's gVCF (case 21), declared male, with
# the real pgsc_calc:
#   - strand-ambiguous rows: a synthetic score on the chr20 slice has five
#     hom-ref rows (effect allele = reference base, the other allele neither
#     it nor its complement), three hom-ref rows whose other allele is the
#     complement (A/T, C/G: ambiguous) and three chrX hom-ref rows. By
#     default the summary says "3 dropped" and pgsc_calc matches 8 of 11;
#     with PRS_KEEP_AMBIGUOUS=true it says "3 kept", matches all 11 and the
#     sum grows by exactly 2 x (32 + 64 + 128). One score (PGS000055) also
#     has a row that agrees with the reference only after a strand flip: it
#     refuses the flag ("3 dropped, strand check failed", its 3 ambiguous
#     rows out of the file);
#   - sex: PGS000004 (breast cancer, female) is "not applicable to the
#     declared sex" with its raw sum and no percentile, PGS000662 (prostate
#     cancer, male) "applies";
#   - the rerun keeps the first run: its summary under its time, with
#     "3 dropped", and its pgsc_calc folder as pgsc_calc.<time>;
#   - the four lipid and blood-pressure scores of assets/pgs_scores.tsv
#     (PGS002257, PGS000065, PGS000064, PGS000066) are the real PGS Catalog
#     files: each gets a row with its total, and the strand check's verdict
#     on them and the step's time and pgsc_calc's task peaks go to the
#     notes. The fixture covers a few Mb, so their match rates, sums and
#     percentiles say nothing about a genome;
#   - Nextflow with --prs_keep_ambiguous on the same scores and genotypes:
#     the same sums, "3 kept" and the refused score, --keep_ambiguous on
#     pgsc_calc's command, the dated summary published beside the plain one;
#     and a row with a VCF alone stops the PRS task (WGS only).
# Removes what it made.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
D="${GENOME_DIR}/${SAMPLE}"
REF="${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta"
check "case 21 wrote a gVCF" vcf_ok "${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz"
NEW_IDS="PGS002257 PGS000065 PGS000064 PGS000066"
for id in $NEW_IDS; do
  check "${id} is in assets/pgs_scores.tsv" grep -q "^${id}	" "${REPO}/assets/pgs_scores.tsv"
done

# --- sites -----------------------------------------------------------------------
# chrom, pos, effect allele, other allele, weight, kind
SITES="${CASE_TMP}/sites.tsv"
python3 - "${FIXTURE_DIR}/HG002_truth_chr20.bed" "${FIXTURE_DIR}/HG002_truth_chr20.vcf.gz" "$REF" > "$SITES" <<'PY'
import gzip, sys
bed, truth, fasta = sys.argv[1:4]
iv = [tuple(int(x) for x in l.split()[1:3]) for l in open(bed) if l.startswith("chr20\t")]
var = [int(l.split("\t")[1]) for l in gzip.open(truth, "rt") if l.startswith("chr20\t")]
fai = {l.split()[0]: [int(x) for x in l.split()[1:5]] for l in open(fasta + ".fai")}
length, offset, bases, width = fai["chr20"]
fh = open(fasta, "rb")
def base(pos):
    i = pos - 1
    fh.seek(offset + (i // bases) * width + i % bases)
    return fh.read(1).decode().upper()
comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
homref, p = [], 10_050_000
while len(homref) < 9 and p < 10_450_000:
    inside = any(s + 200 < p <= e - 200 for s, e in iv)
    clear = all(abs(v - p) > 200 for v in var)
    b = base(p)
    if inside and clear and b in comp:
        homref.append((p, b))
    p += 7_919
for i, (p, b) in enumerate(homref[:5]):
    other = next(x for x in "ACGT" if x not in (b, comp[b]))
    print(f"20\t{p}\t{b}\t{other}\t{2 ** i}\thomref")
for (p, b), w in zip(homref[5:8], (32, 64, 128)):
    print(f"20\t{p}\t{b}\t{comp[b]}\t{w}\tambiguous")
# A row given on the other strand: its alleles are the complements of the
# reference base and of a third base.
p, b = homref[8]
third = next(x for x in "ACGT" if x not in (b, comp[b]))
print(f"20\t{p}\t{comp[b]}\t{comp[third]}\t0.5\tflipped")
PY
cat "$SITES"
check_eq "score sites (5 hom-ref, 3 ambiguous hom-ref, 1 flipped)" "$(wc -l < "$SITES" | tr -d ' ')" 9
XR=chrX:73700001-74000000
bcf query -r "$XR" -f '%POS\t%INFO/END\t[%GT]\n' "${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz" > "${CASE_TMP}/x_blocks.tsv"
bcf query -r "$XR" -i 'GT="alt"' -f '%POS\n' "${SAMPLE}/vcf/${SAMPLE}.vcf.gz" > "${CASE_TMP}/x_vars.tsv"
X_SITES="${CASE_TMP}/x_sites.tsv"
python3 - "${CASE_TMP}/x_blocks.tsv" "${CASE_TMP}/x_vars.tsv" "$REF" > "$X_SITES" <<'PY2'
import sys
blocks, variants, fasta = sys.argv[1:4]
var = [int(l) for l in open(variants) if l.strip()]
fai = {l.split()[0]: [int(x) for x in l.split()[1:5]] for l in open(fasta + ".fai")}
length, offset, bases, width = fai["chrX"]
fh = open(fasta, "rb")
def base(pos):
    i = pos - 1
    fh.seek(offset + (i // bases) * width + i % bases)
    return fh.read(1).decode().upper()
comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
rows = [l.rstrip("\n").split("\t") for l in open(blocks) if l.strip()]
refs = [(int(r[0]), int(r[1])) for r in rows if len(r) == 3 and r[1] not in ("", ".") and r[2] in ("0/0", "0|0")]
out, last = [], 0
for start, end in refs:
    p = (start + end) // 2
    if end - start >= 50 and p - last > 2000 and all(abs(v - p) > 100 for v in var) and base(p) in comp:
        out.append((p, base(p)))
        last = p
for (p, b), w in zip(out[:3], (0.25, 0.25, 0.25)):
    other = next(x for x in "ACGT" if x not in (b, comp[b]))
    print(f"X\t{p}\t{b}\t{other}\t{w}\tchrx")
PY2
cat "$X_SITES"
check_eq "chrX hom-ref rows" "$(wc -l < "$X_SITES" | tr -d ' ')" 3

# --- the scores: synthetic for the list's first nine, the real files for the four new ones ---
SCORES="${GENOME_DIR}/prs_scores"
rm -rf "$SCORES" "${D}/prs"
mkdir -p "$SCORES"
PGS_BASE=https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores
for id in $(awk -F'\t' '$1 ~ /^PGS[0-9]+$/ {print $1}' "${REPO}/assets/pgs_scores.tsv"); do
  case " ${NEW_IDS} " in
    *" ${id} "*)
      url="${PGS_BASE}/${id}/ScoringFiles/Harmonized/${id}_hmPOS_GRCh38.txt.gz"
      curl -fsSL --retry 3 -o "${SCORES}/${id}.txt.gz" "$url"
      want=$(curl -fsSL --retry 3 "${url}.md5" | awk '{print $1}')
      check_eq "${id}: the PGS Catalog's GRCh38-harmonised file, md5 checked" \
        "$( (md5sum "${SCORES}/${id}.txt.gz" 2>/dev/null || echo x) | awk '{print $1}')" "${want:-none}"
      continue ;;
  esac
  { printf '#pgs_id=%s\n#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n' "$id"
    awk -F'\t' -v OFS='\t' -v id="$id" '$6 != "flipped" || id == "PGS000055" {print $1, $2, $3, $4, $5}' "$SITES"
    cut -f1-5 "$X_SITES"; } | gzip -c > "${SCORES}/${id}.txt.gz"
done
# Rows of each real file pgsc_calc gets: GRCh38-placed, on chr1-22 or chrX (the sample is male).
for id in $NEW_IDS; do
  gzip -dc "${SCORES}/${id}.txt.gz" | awk -F'\t' '/^#/ {next} !h {for (i = 1; i <= NF; i++) c[$i] = i; h = 1; next}
    $c["hm_chr"] ~ /^([1-9]|1[0-9]|2[0-2]|X)$/ && $c["hm_pos"] != "" && $c["effect_weight"] != "" {n++} END {print n + 0}' \
    > "${CASE_TMP}/${id}.rows"
done

SUMMARY="${D}/prs/${SAMPLE}_prs_summary.tsv"
col() {  # col FILE PGS_ID COLUMN
  awk -F'\t' -v id="$2" -v k="$3" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $c["PGS_ID"] == id {print $c[k]}' "$1" 2>/dev/null
}

# --- step 25, defaults -----------------------------------------------------------------
ANCESTRY_PANEL=none run_step 25-prs.sh "$SAMPLE" male
check_step_exit "25-prs.sh male"
cat "$SUMMARY" 2>/dev/null
SUM1=$(col "$SUMMARY" PGS000018 Score_SUM)
check_eq "default: PGS000018's 3 ambiguous rows are dropped" "$(col "$SUMMARY" PGS000018 Ambiguous)" "3 dropped"
check_eq "default: 11 rows (8 autosomal, 3 chrX)" "$(col "$SUMMARY" PGS000018 Variants_Total)" 11
check_eq "default: pgsc_calc matched 8 (the ambiguous 3 excluded)" "$(col "$SUMMARY" PGS000018 Variants_Matched)" 8
check_eq "male: PGS000004 (breast cancer) is female-specific" "$(col "$SUMMARY" PGS000004 Applies_To)" female
check_eq "male: PGS000004 is not applicable to the declared sex" "$(col "$SUMMARY" PGS000004 Sex_Check)" "not applicable to the declared sex"
check "male: PGS000004 keeps its raw sum" has '^-?[0-9.e+-]+$' "$(col "$SUMMARY" PGS000004 Score_SUM)"
check_eq "male: PGS000004 has no percentile" "$(col "$SUMMARY" PGS000004 Percentile)" NA
check_eq "male: PGS000662 (prostate cancer) applies" "$(col "$SUMMARY" PGS000662 Sex_Check)" applies
check_eq "male: PGS000018 is not sex-specific" "$(col "$SUMMARY" PGS000018 Sex_Check)" ""
check "the step log marks PGS000004 not applicable" has 'PGS000004\): sum .*female-specific: not applicable to the declared sex' "$(cat "$STEP_LOG")"
for id in $NEW_IDS; do
  check_eq "${id}: a row with every placed row of the real file in its total" "$(col "$SUMMARY" "$id" Variants_Total)" "$(cat "${CASE_TMP}/${id}.rows")"
done
DATED1=$(find "${D}/prs" -maxdepth 1 -name "${SAMPLE}_prs_summary.*.tsv" | head -n 1)
check "the first run's summary is also kept under its time" test -s "${DATED1:-/nonexistent}"

# --- step 25, PRS_KEEP_AMBIGUOUS=true (a minute later, so the first run keeps its own name) -----
touch -d '-2 minutes' "$SUMMARY" "$DATED1" 2>/dev/null
mv "$DATED1" "${CASE_TMP}/first_dated.tsv"
STAMP1=$(date -r "$SUMMARY" +%Y%m%d-%H%M)
t0=$SECONDS
PRS_KEEP_AMBIGUOUS=true ANCESTRY_PANEL=none run_step 25-prs.sh "$SAMPLE" male
T_STEP=$((SECONDS - t0))
check_step_exit "25-prs.sh male, PRS_KEEP_AMBIGUOUS=true"
cat "$SUMMARY" 2>/dev/null
SUM2=$(col "$SUMMARY" PGS000018 Score_SUM)
check_eq "keep: PGS000018's 3 ambiguous rows are kept" "$(col "$SUMMARY" PGS000018 Ambiguous)" "3 kept"
check_eq "keep: pgsc_calc matched all 11" "$(col "$SUMMARY" PGS000018 Variants_Matched)" 11
check "keep: the sum grew by 2 x (32 + 64 + 128) = 448 (${SUM1:-none} to ${SUM2:-none})" \
  awk -v a="${SUM1:-x}" -v b="${SUM2:-x}" 'BEGIN { if (a !~ /^-?[0-9.e+-]+$/ || b !~ /^-?[0-9.e+-]+$/) exit 1; d = b - a - 448; exit !(d < 1e-6 && d > -1e-6) }'
check_eq "keep: PGS000055, with a flipped row, refuses the flag" "$(col "$SUMMARY" PGS000055 Ambiguous)" "3 dropped, strand check failed"
check_eq "keep: PGS000055's ambiguous rows are out of its file (12 - 3)" "$(col "$SUMMARY" PGS000055 Variants_Total)" 9
check "keep: pgsc_calc ran with --keep_ambiguous" grep -q -- '--keep_ambiguous' "${D}/prs/pgsc_calc/nextflow.log"
check "control: the default run's pgsc_calc had no --keep_ambiguous" bash -c '! grep -q -- --keep_ambiguous "$1"' _ "${D}/prs/pgsc_calc.${STAMP1}/nextflow.log"
check "rerun: the first run's pgsc_calc folder is kept as pgsc_calc.${STAMP1}" test -s "${D}/prs/pgsc_calc.${STAMP1}/results/sample/match/sample_summary.csv"
check_eq "rerun: the first run's summary is kept under its time, with its own Ambiguous" \
  "$(col "${D}/prs/${SAMPLE}_prs_summary.${STAMP1}.tsv" PGS000018 Ambiguous)" "3 dropped"
STRAND=$(grep -E '^  PGS0[0-9]+: .*strand check' "$STEP_LOG")
echo "$STRAND"
for id in $NEW_IDS; do
  check "keep: prs-format checked ${id}'s strand against the reference" has "^  ${id}: .*strand check: [0-9]+ rows agree with the reference, [0-9]+ flipped" "$STRAND"
done
TRACE="${D}/prs/pgsc_calc/pipeline_info/execution_trace.txt"
TOP=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
  { v = $c["peak_rss"]; n = v + 0; if (v ~ /GB/) n *= 1024; else if (v ~ /KB/) n /= 1024; else if (v !~ /MB/) n /= 1048576
    printf "%.0f\t%s\t%s\n", n, $c["name"], $c["realtime"] }' "$TRACE" 2>/dev/null | sort -rn | head -4)
check "pgsc_calc's trace has its task peaks" test -n "$TOP"
{
  printf '#### PRS: the four lipid and blood-pressure scores on the fixture (step 25, male, PRS_KEEP_AMBIGUOUS=true)\n\n```\n'
  printf 'step 25 wall time with the 13 scores: %s s\n' "$T_STEP"
  for id in $NEW_IDS; do
    printf '%s: total %s, matched %s, sum %s, ambiguous %s\n' "$id" "$(col "$SUMMARY" "$id" Variants_Total)" \
      "$(col "$SUMMARY" "$id" Variants_Matched)" "$(col "$SUMMARY" "$id" Score_SUM)" "$(col "$SUMMARY" "$id" Ambiguous)"
  done
  echo "$STRAND"
  echo "pgsc_calc tasks with the highest peak RSS (MB, name, time):"; echo "$TOP"
  printf '```\n\n'
} >> "$E2E_NOTES"
cp "$SUMMARY" "${CASE_TMP}/bash_summary.tsv"

# --- Nextflow, --prs_keep_ambiguous ---------------------------------------------------------
NF_OUT="${GENOME_DIR}/nf-prs-strand"
WORK="${E2E_WORK}/nf-work-prs-strand"
SHEET="${CASE_TMP}/samplesheet.csv"
CALC="${GENOME_DIR}/tools/pgsc_calc-${PGSC_CALC_VERSION}"
rm -rf "$NF_OUT" "$WORK"
printf 'sample,vcf,vcf_index,gvcf,gvcf_index,sex\n%s,%s,%s,%s,%s,male\n' "$SAMPLE" \
  "${D}/vcf/${SAMPLE}.vcf.gz" "${D}/vcf/${SAMPLE}.vcf.gz.tbi" "${D}/vcf/${SAMPLE}.g.vcf.gz" "${D}/vcf/${SAMPLE}.g.vcf.gz.tbi" > "$SHEET"
nf_run() {  # nf_run SHEET OUT: --tools prs --prs_keep_ambiguous; exit code in NF_RC
  ( cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -work-dir "$WORK" \
      --input "$1" --reference "$REF" --sex_check warn --tools prs --prs_keep_ambiguous \
      --pgs_scoring "$SCORES" --pgsc_calc "$CALC" --outdir "$2" --max_cpus 4 --max_memory 14.GB ) 2>&1 | tee "$STEP_LOG"
  NF_RC=${PIPESTATUS[0]}
}
echo "+ nextflow run main.nf --tools prs --prs_keep_ambiguous (VCF + gVCF row, male)"
nf_run "$SHEET" "$NF_OUT"
check_eq "nextflow run exits 0" "$NF_RC" 0
NFS="${NF_OUT}/${SAMPLE}/prs/${SAMPLE}_prs_summary.tsv"
cat "$NFS" 2>/dev/null
check_eq "nextflow: PGS000018's 3 ambiguous rows are kept" "$(col "$NFS" PGS000018 Ambiguous)" "3 kept"
check_eq "nextflow: PGS000018's sum equals step 25's" "$(col "$NFS" PGS000018 Score_SUM)" "${SUM2:-none}"
check_eq "nextflow: PGS000055 refuses the flag" "$(col "$NFS" PGS000055 Ambiguous)" "3 dropped, strand check failed"
check_eq "nextflow: PGS000004 is not applicable to the declared sex" "$(col "$NFS" PGS000004 Sex_Check)" "not applicable to the declared sex"
check "nextflow: every score's sum, matched count and Ambiguous equal step 25's" \
  diff <(cut -f2,3,4,5,10,11,12 "${CASE_TMP}/bash_summary.tsv") <(cut -f2,3,4,5,10,11,12 "$NFS")
check_ge "nextflow: the summary is also published under its time" \
  "$(find "${NF_OUT}/${SAMPLE}/prs" -maxdepth 1 -name "${SAMPLE}_prs_summary.*.tsv" | wc -l | tr -d ' ')" 1
CMD=$(grep -l 'nextflow -log pgsc_calc.nextflow.log run' "$WORK"/*/*/.command.sh 2>/dev/null | head -n 1)
check "nextflow: the PRS task gave pgsc_calc --keep_ambiguous" grep -q -- '--keep_ambiguous' "${CMD:-/dev/null}"

# The reverse: a row with a VCF alone is refused, before pgsc_calc runs.
SHEET2="${CASE_TMP}/samplesheet_vcf.csv"
printf 'sample,vcf,vcf_index,sex\n%s,%s,%s,male\n' "$SAMPLE" "${D}/vcf/${SAMPLE}.vcf.gz" "${D}/vcf/${SAMPLE}.vcf.gz.tbi" > "$SHEET2"
echo "+ nextflow run main.nf --tools prs --prs_keep_ambiguous (VCF row, no gVCF)"
nf_run "$SHEET2" "${GENOME_DIR}/nf-prs-strand-vcf"
check "nextflow: a VCF-only row with --prs_keep_ambiguous fails the run" test "$NF_RC" -ne 0
check "nextflow: and says it needs the gVCF" has 'prs_keep_ambiguous needs the sample.s gVCF' "$(cat "$STEP_LOG" "${CASE_TMP}/.nextflow.log" 2>/dev/null)"

printf '#### PRS: strand-ambiguous rows and sex-specific scores\n\nstep 25 default sum %s, keep_ambiguous %s (448 more); Nextflow %s\n\n' \
  "${SUM1:-none}" "${SUM2:-none}" "$(col "$NFS" PGS000018 Score_SUM)" >> "$E2E_NOTES"
rm -rf "$SCORES" "${D:?}/prs" "$NF_OUT" "${GENOME_DIR}/nf-prs-strand-vcf" "$WORK"
finish
