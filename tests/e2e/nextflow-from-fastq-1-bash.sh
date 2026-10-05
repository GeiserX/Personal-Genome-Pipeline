#!/usr/bin/env bash
# The bash leg of the parity check (scripts/ci/parity-diff.sh): the single
# step scripts on the fixture reads, as sample HG002P: 01b (fastp), 02
# (minimap2, markdup), 03 (DeepVariant, male), 06, 08 and 36 (T1K's HLA types
# as PharmCAT's outside calls), 07, 11 and 25. The next case runs the
# Nextflow pipeline on the same reads under the same sample name; the E2E
# workflow then compares the two.
#
# PRS scores a synthetic file (the nine PGS ids of step 25 so nothing is
# downloaded): five chr20 sites where HG002 matches the reference (GIAB
# v4.2.1 truth, no variant within 200 bp; effect allele = reference base) and
# three truth SNVs (effect allele = ALT). The next case scores the same file.
. "$(dirname "$0")/lib.sh"

P=HG002P
D="${GENOME_DIR}/${P}"
rm -rf "$D"
mkdir -p "${D}/fastq"
for r in R1 R2; do
  ln -f "${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_${r}.fastq.gz" "${D}/fastq/${P}_${r}.fastq.gz"
done

run_step 01b-fastp-qc.sh "$P"
check_step_exit 01b-fastp-qc.sh
run_step 02-alignment.sh "$P"
check_step_exit 02-alignment.sh
check "step 02 aligned the trimmed reads" has 'Using trimmed FASTQs from fastp' "$(cat "$STEP_LOG")"
# The fixture's slices only, as case 21 (the rest of the reference has no reads)
INTERVALS=$(awk '{printf "%s%s:%d-%d", (NR > 1 ? " " : ""), $1, $2 + 1, $3}' "${FIXTURE_DIR}/regions.bed")
INTERVALS="$INTERVALS" run_step 03-deepvariant.sh "$P" male
check_step_exit 03-deepvariant.sh
run_step 06-clinvar-screen.sh "$P"
check_step_exit 06-clinvar-screen.sh
# HLA types first: PharmCAT reads them as outside calls (step 36), as the
# pipeline's PHARMCAT reads PGX_CONSENSUS's.
run_step 08-hla-typing.sh "$P"
check_step_exit 08-hla-typing.sh
run_step 36-pgx-consensus.sh "$P"
check_step_exit 36-pgx-consensus.sh
run_step 07-pharmacogenomics.sh "$P"
check_step_exit 07-pharmacogenomics.sh
run_step 11-roh-analysis.sh "$P"
check_step_exit 11-roh-analysis.sh

# --- synthetic PRS scores ---------------------------------------------------------
SITES="${CASE_TMP}/score_sites.tsv"
python3 - "${FIXTURE_DIR}/HG002_truth_chr20.bed" "${FIXTURE_DIR}/HG002_truth_chr20.vcf.gz" \
  "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta" > "$SITES" <<'PY'
import gzip, sys
bed, truth, fasta = sys.argv[1:4]
iv = [tuple(int(x) for x in l.split()[1:3]) for l in open(bed) if l.startswith("chr20\t")]
var, snv = [], []
for l in gzip.open(truth, "rt"):
    if not l.startswith("chr20\t"):
        continue
    f = l.rstrip("\n").split("\t")
    var.append(int(f[1]))
    gt = f[9].split(":")[0].replace("|", "/")
    if len(f[3]) == 1 and len(f[4]) == 1 and 10_050_000 < int(f[1]) < 10_450_000 and gt in ("0/1", "1/1"):
        snv.append((int(f[1]), f[4]))
fai = {l.split()[0]: [int(x) for x in l.split()[1:5]] for l in open(fasta + ".fai")}
length, offset, bases, width = fai["chr20"]
def base(pos):
    i = pos - 1
    with open(fasta, "rb") as fh:
        fh.seek(offset + (i // bases) * width + i % bases)
        return fh.read(1).decode().upper()
out, p = [], 10_050_000
while len(out) < 5 and p < 10_450_000:
    inside = any(s + 200 < p <= e - 200 for s, e in iv)
    clear = all(abs(v - p) > 200 for v in var)
    b = base(p)
    if inside and clear and b in "ACGT":
        out.append((p, b))
    p += 7_919
out += snv[:: max(1, len(snv) // 3)][:3]
for p, a in out:
    print(f"chr20\t{p}\t{a}")
PY
cat "$SITES"
check_eq "score sites (5 hom-ref, 3 truth SNVs)" "$(wc -l < "$SITES" | tr -d ' ')" 8
PGS="${E2E_WORK}/parity-pgs"
rm -rf "$PGS" "${GENOME_DIR}/prs_scores"
mkdir -p "$PGS" "${GENOME_DIR}/prs_scores"
for id in PGS000018 PGS000014 PGS000004 PGS000662 PGS000016 PGS000334 PGS000027 PGS000017 PGS000055; do
  { printf '#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\teffect_weight\n'
    awk -F'\t' -v OFS='\t' '{print "20", $2, $3, 2 ^ (NR - 1)}' "$SITES"; } | gzip -c > "${PGS}/${id}.txt.gz"
done
cp "${PGS}"/*.txt.gz "${GENOME_DIR}/prs_scores/"
run_step 25-prs.sh "$P"
check_step_exit 25-prs.sh
rm -rf "${GENOME_DIR}/prs_scores"

# --- what the parity check reads ------------------------------------------------
for f in "aligned/${P}_sorted.bam" "vcf/${P}.vcf.gz" "vcf/${P}.g.vcf.gz" "clinvar/${P}_clinvar_hits.tsv" \
         "vcf/${P}_roh.txt" "vcf/${P}.report.json" "prs/${P}_prs_summary.tsv"; do
  check "bash leg wrote ${f}" nonempty "${P}/${f}"
done
check_eq "PRS input of the bash leg" "$(awk -F'\t' 'NR == 2 {print $NF}' "${D}/prs/${P}_prs_summary.tsv" 2>/dev/null)" gvcf

finish
