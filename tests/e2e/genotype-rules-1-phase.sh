#!/usr/bin/env bash
# Step 42's genotype rules (assets/genotype_rules.tsv) on synthetic data, with
# the real bcftools, WhatsHap and Python images, and the Nextflow mirror:
#   - scripts/ci/check-rule-alleles.sh --online: every rule site matches
#     dbSNP's GRCh38 record, and the dbSNP table pinned in the script is what
#     dbSNP's VCF of the pinned build holds today; controls: a rules copy with
#     rs1726866's ref and alt swapped fails, naming it, and a copy of the
#     script whose pinned table differs from dbSNP fails the online compare;
#   - TAS2R38 on a made-up chr7 (N everywhere except a 2.6 kb window around the
#     gene, its three sites carrying the GRCh38 bases), one gVCF heterozygous at
#     the three sites, four samples:
#       GRSPAN    read pairs that span the three sites, from PAV and AVI copies
#                 -> WhatsHap phases them in one block: PAV/AVI, read_backed
#       GRSWAP    the same gVCF, read pairs from G,A,C and C,G,T copies
#                 -> other_haplotype, so the result rests on the reads
#       GRNOSPAN  the same gVCF, single reads that each cover one site
#                 -> unphased: no read or read pair joins the sites
#       GRREF     a gVCF that is one reference block (GQ 42), no BAM
#                 -> every site homozygous reference with the block's GQ;
#                    the reference haplotype AAI is not one the rule reads
#     and the traits table keeps one row per table row plus the APOE line;
#   - nextflow run main.nf --tools traits on GRSPAN (VCF, gVCF and BAM) writes
#     the same rules table, traits table and genotypes as the bash step, byte
#     for byte; control: the comparison tells GRSPAN's rules table from
#     GRNOSPAN's.
# Writes four sample folders and a reference folder under GENOME_DIR, all
# removed at the end; no other case reads them.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
CHECK="${REPO}/scripts/ci/check-rule-alleles.sh"
RULES="${REPO}/assets/genotype_rules.tsv"
TABLE="${REPO}/assets/traits.tsv"

# --- the dbSNP orientation check, online ---------------------------------------------
bash "$CHECK" --online > "${CASE_TMP}/dbsnp.log" 2>&1
RC=$?
cat "${CASE_TMP}/dbsnp.log"
check_eq "check-rule-alleles.sh --online exits 0" "$RC" 0
check "the pinned dbSNP table is what dbSNP holds today" has '^ok: the pinned table is what dbSNP' "$(cat "${CASE_TMP}/dbsnp.log")"
awk -F'\t' -v OFS='\t' '$1 == "site" && $3 == "rs1726866" { t = $6; $6 = $7; $7 = t } { print }' "$RULES" > "${CASE_TMP}/swapped.tsv"
check_eq "control: the swapped copy differs from the rules in one line" \
  "$(diff "$RULES" "${CASE_TMP}/swapped.tsv" | grep -c '^>' || true)" 1
bash "$CHECK" --online "${CASE_TMP}/swapped.tsv" > "${CASE_TMP}/swapped.log" 2>&1
check "control: the swapped copy fails (exit $?), naming rs1726866" has 'FAIL: rs1726866: ref and alt are swapped' \
  "$(cat "${CASE_TMP}/swapped.log")"
mkdir -p "${CASE_TMP}/root/scripts/ci"
cp "${REPO}/versions.env" "${CASE_TMP}/root/"
sed 's/^#: rs7412\tchr19\t44908822\tC\tT$/#: rs7412\tchr19\t44908822\tC\tG/' "$CHECK" > "${CASE_TMP}/root/scripts/ci/check-rule-alleles.sh"
check "control: the script copy's pinned table differs in one line" \
  test "$(diff "$CHECK" "${CASE_TMP}/root/scripts/ci/check-rule-alleles.sh" | grep -c '^>')" -eq 1
bash "${CASE_TMP}/root/scripts/ci/check-rule-alleles.sh" --online "$RULES" > "${CASE_TMP}/table.log" 2>&1
check "control: a pinned table that differs from dbSNP fails the online compare" \
  has '^FAIL: the pinned table differs from dbSNP' "$(cat "${CASE_TMP}/table.log")"

# --- synthetic chr7, reads and gVCFs --------------------------------------------------
RD="${G}/reference-genotype-rules"
rm -rf "$RD"
mkdir -p "$RD"
for s in GRSPAN GRSWAP GRNOSPAN GRREF; do rm -rf "${G:?}/${s}"; mkdir -p "${G}/${s}/vcf" "${G}/${s}/aligned"; done
python3 - "$RD" "$G" <<'PY'
import random, sys
rd, g = sys.argv[1], sys.argv[2]
random.seed(42)
LEN, W0, W1 = 141975000, 141972000, 141974600   # chr7 length; the window [W0, W1), 1-based
S3, S2, S1 = 141972804, 141972905, 141973545     # rs10246939 (T), rs1726866 (G), rs713598 (C)
REF = {S3: "T", S2: "G", S1: "C"}
win = [random.choice("ACGT") for _ in range(W1 - W0)]
for p, b in REF.items():
    win[p - W0] = b
seq = lambda p: win[p - W0]                       # noqa: E731
with open(f"{rd}/chr7.fa", "w") as f:
    f.write(">chr7\n")
    line = 60
    buf = []
    for start in range(1, LEN + 1, line):
        end = min(start + line, LEN + 1)
        if end <= W0 or start >= W1:
            buf.append("N" * (end - start))
        else:
            buf.append("".join(win[p - W0] if W0 <= p < W1 else "N" for p in range(start, end)))
        if len(buf) >= 100000:
            f.write("\n".join(buf) + "\n"); buf = []
    f.write("\n".join(buf) + "\n")

def gvcf(path, het):
    hdr = ("##fileformat=VCFv4.2\n##contig=<ID=chr7,length=%d>\n" % LEN
           + '##ALT=<ID=*,Description="Any other allele">\n'
           + '##INFO=<ID=END,Number=1,Type=Integer,Description="End of the block">\n'
           + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
           + '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
           + '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
           + '##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum depth">\n'
           + '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">\n'
           + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
    rec, start = [], W0
    for p in (sorted(het) if het else []):
        rec.append(f"chr7\t{start}\t.\t{seq(start)}\t<*>\t0\t.\tEND={p - 1}\tGT:GQ:MIN_DP:PL\t0/0:42:30:0,90,900")
        rec.append(f"chr7\t{p}\t.\t{seq(p)}\t{het[p]},<*>\t40\tPASS\t.\tGT:GQ:DP:PL\t0/1:35:40:400,0,400,990,990,990")
        start = p + 1
    rec.append(f"chr7\t{start}\t.\t{seq(start)}\t<*>\t0\t.\tEND={W1 - 1}\tGT:GQ:MIN_DP:PL\t0/0:42:30:0,90,900")
    with open(path, "w") as f:
        f.write(hdr + "\n".join(rec) + "\n")

HET = {S3: "C", S2: "A", S1: "G"}
for s in ("GRSPAN", "GRSWAP", "GRNOSPAN"):
    gvcf(f"{g}/{s}/vcf/{s}.g.vcf", HET)
    with open(f"{g}/{s}/vcf/{s}.vcf", "w") as f:   # the variant-only VCF the samplesheet needs
        f.write("##fileformat=VCFv4.2\n##contig=<ID=chr7,length=%d>\n" % LEN
                + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
                + "".join(f"chr7\t{p}\t.\t{seq(p)}\t{HET[p]}\t40\tPASS\t.\tGT\t0/1\n" for p in sorted(HET)))
gvcf(f"{g}/GRREF/vcf/GRREF.g.vcf", {})

def hap(alleles):
    """The window with these alleles at S1, S2, S3."""
    h = list(win)
    for p, b in zip((S1, S2, S3), alleles):
        h[p - W0] = b
    return h

def comp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))

def sam(path, frags):
    with open(path, "w") as f:
        f.write(f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr7\tLN:{LEN}\n")
        for line in frags:
            f.write(line + "\n")

def pairs(h, tag, n=20):
    """Read pairs: read 1 covers S3 and S2, its mate S1."""
    out = []
    for i in range(n):
        a, b = S3 - 40 + i, S1 - 60 + i
        s1, s2 = "".join(h[a - W0:a - W0 + 150]), "".join(h[b - W0:b - W0 + 150])
        tlen = b + 150 - a
        name = f"{tag}{i}"
        out.append(f"{name}\t99\tchr7\t{a}\t60\t150M\t=\t{b}\t{tlen}\t{s1}\t{'I' * 150}")
        out.append(f"{name}\t147\tchr7\t{b}\t60\t150M\t=\t{a}\t{-tlen}\t{s2}\t{'I' * 150}")
    return out

def singles(h, tag, n=20):
    """Single 60-base reads centred on one site each: none covers two."""
    out = []
    for p in (S3, S2, S1):
        for i in range(n):
            a = p - 30 + (i % 10) - 5
            s = "".join(h[a - W0:a - W0 + 60])
            out.append(f"{tag}{p}_{i}\t{16 if i % 2 else 0}\tchr7\t{a}\t60\t60M\t*\t0\t0\t{s}\t{'I' * 60}")
    return out

PAV, AVI = hap("GGC"), hap("CAT")
sam(f"{g}/GRSPAN/aligned/GRSPAN.sam", pairs(PAV, "pav") + pairs(AVI, "avi"))
sam(f"{g}/GRSWAP/aligned/GRSWAP.sam", pairs(hap("GAC"), "gac") + pairs(hap("CGT"), "cgt"))
sam(f"{g}/GRNOSPAN/aligned/GRNOSPAN.sam", singles(PAV, "pav") + singles(AVI, "avi"))
PY
check "the synthetic chr7 was written" test -s "${RD}/chr7.fa"
sam faidx reference-genotype-rules/chr7.fa
check_eq "chr7.fa.fai: one contig of 141975000 bases" "$(cut -f1,2 "${RD}/chr7.fa.fai" 2>/dev/null | tr '\t' ' ')" "chr7 141975000"
for s in GRSPAN GRSWAP GRNOSPAN GRREF; do
  in_genome "$BCFTOOLS_IMAGE" bash -c "bcftools view -Oz -o '${s}/vcf/${s}.g.vcf.gz' '${s}/vcf/${s}.g.vcf' && bcftools index -t '${s}/vcf/${s}.g.vcf.gz'"
  check "${s}: the gVCF is indexed" test -s "${G}/${s}/vcf/${s}.g.vcf.gz.tbi"
done
for s in GRSPAN GRSWAP GRNOSPAN; do
  in_genome "$BCFTOOLS_IMAGE" bash -c "bcftools view -Oz -o '${s}/vcf/${s}.vcf.gz' '${s}/vcf/${s}.vcf' && bcftools index -t '${s}/vcf/${s}.vcf.gz'"
  in_genome "$SAMTOOLS_IMAGE" bash -c "samtools sort -o '${s}/aligned/${s}_sorted.bam' '${s}/aligned/${s}.sam' && samtools index '${s}/aligned/${s}_sorted.bam'"
  check "${s}: the BAM is indexed" test -s "${G}/${s}/aligned/${s}_sorted.bam.bai"
done
check_eq "GRSPAN: 80 reads in the BAM" "$(sam view -c GRSPAN/aligned/GRSPAN_sorted.bam 2>/dev/null)" 80
check_eq "GRSPAN: the gVCF genotypes at the three sites" \
  "$(bcf query -r chr7:141972804,chr7:141972905,chr7:141973545 -f '[%GT]\n' GRSPAN/vcf/GRSPAN.g.vcf.gz 2>/dev/null | tr '\n' ' ')" "0/1 0/1 0/1 "

# --- the bash step ---------------------------------------------------------------------
export REF_FASTA="${RD}/chr7.fa"
rule() { awk -F'\t' -v r="$2" -v k="$3" 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
  $1 == r { print $c[k]; exit }' "${G}/$1/traits/$1_genotype_rules.tsv" 2>/dev/null; }
ROWS=$(grep -v '^#' "$TABLE" | awk 'NR > 1' | wc -l | tr -d ' ')
for s in GRSPAN GRSWAP GRNOSPAN GRREF; do
  run_step 42-traits.sh "$s"
  check_step_exit 42-traits.sh
  cp "$STEP_LOG" "${CASE_TMP}/${s}.log"
  check_eq "${s}: one traits row per table row plus the APOE line" \
    "$(awk 'NR > 1' "${G}/${s}/traits/${s}_traits.tsv" 2>/dev/null | wc -l | tr -d ' ')" "$((ROWS + 1))"
  check_eq "${s}: one rules line per rule" \
    "$(awk 'NR > 1' "${G}/${s}/traits/${s}_genotype_rules.tsv" 2>/dev/null | wc -l | tr -d ' ')" \
    "$(awk -F'\t' '$1 == "rule"' "$RULES" | wc -l | tr -d ' ')"
  check_eq "${s}: APOE, whose sites are not on the synthetic reference, is not covered" "$(rule "$s" APOE status)" not_covered
done
unset REF_FASTA
for s in GRSPAN GRSWAP GRNOSPAN; do
  check_eq "${s}: the three TAS2R38 genotypes, reference allele first" "$(rule "$s" TAS2R38 genotypes)" "C/G,G/A,T/C"
  check_eq "${s}: the TAS2R38 traits rows read het" \
    "$(awk -F'\t' '$2 == "TAS2R38" { printf "%s:%s ", $6, $7 }' "${G}/${s}/traits/${s}_traits.tsv" 2>/dev/null)" \
    "C/G:called G/A:called T/C:called "
done
check_eq "GRSPAN: WhatsHap phased the pairs: PAV/AVI" "$(rule GRSPAN TAS2R38 result) $(rule GRSPAN TAS2R38 status)" "PAV/AVI called"
check_eq "GRSPAN: phase read_backed" "$(rule GRSPAN TAS2R38 phase)" read_backed
check_eq "GRSPAN: the GQ of the three records" "$(rule GRSPAN TAS2R38 gq)" "35,35,35"
check_eq "GRSWAP: the same genotypes phased into G,A,C and C,G,T: other_haplotype" \
  "$(rule GRSWAP TAS2R38 status) $(rule GRSWAP TAS2R38 phase)" "other_haplotype read_backed"
check "GRSWAP: the reading names both copies" has 'G,A,C and C,G,T' "$(rule GRSWAP TAS2R38 reading)"
check_eq "GRNOSPAN: no read joins the sites: unphased" "$(rule GRNOSPAN TAS2R38 status) $(rule GRNOSPAN TAS2R38 result)" "unphased ."
check "GRNOSPAN: the reading says no read or read pair covered them" has 'no read or read pair covered them together' \
  "$(rule GRNOSPAN TAS2R38 reading)"
check_eq "GRREF: reference blocks only: homozygous reference at the three sites" "$(rule GRREF TAS2R38 genotypes)" "C/C,G/G,T/T"
check_eq "GRREF: each site with the block's GQ" "$(rule GRREF TAS2R38 gq)" "42,42,42"
check_eq "GRREF: the reference haplotype AAI is not one the rule reads" \
  "$(rule GRREF TAS2R38 status) $(rule GRREF TAS2R38 phase)" "other_haplotype not_needed"
check_eq "GRREF: the TAS2R38 traits rows read homozygous reference" \
  "$(awk -F'\t' '$2 == "TAS2R38" { printf "%s:%s ", $6, $7 }' "${G}/GRREF/traits/GRREF_traits.tsv" 2>/dev/null)" \
  "C/C:called G/G:called T/T:called "
check "GRREF: the step says there is no BAM to phase from" has 'No BAM' "$(cat "${CASE_TMP}/GRREF.log")"
echo "TAS2R38 lines:"
for s in GRSPAN GRSWAP GRNOSPAN GRREF; do
  echo "  ${s}: $(rule "$s" TAS2R38 result) $(rule "$s" TAS2R38 status) $(rule "$s" TAS2R38 phase)"
done

# --- the Nextflow mirror ---------------------------------------------------------------
if command -v nextflow >/dev/null; then
  NF_OUT="${G}/nf-genotype-rules"
  SHEET="${CASE_TMP}/samplesheet.csv"
  S="${G}/GRSPAN"
  printf 'sample,vcf,vcf_index,gvcf,gvcf_index,bam,bam_index\nGRSPAN,%s,%s,%s,%s,%s,%s\n' \
    "${S}/vcf/GRSPAN.vcf.gz" "${S}/vcf/GRSPAN.vcf.gz.tbi" "${S}/vcf/GRSPAN.g.vcf.gz" "${S}/vcf/GRSPAN.g.vcf.gz.tbi" \
    "${S}/aligned/GRSPAN_sorted.bam" "${S}/aligned/GRSPAN_sorted.bam.bai" > "$SHEET"
  (
    cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
      -work-dir "${E2E_WORK}/nf-work-genotype-rules" \
      --input "$SHEET" --reference "${RD}/chr7.fa" --tools traits --sex_check warn \
      --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB
  ) 2>&1 | tee "$STEP_LOG"
  NF_RC=${PIPESTATUS[0]}
  check_eq "nextflow run --tools traits on GRSPAN exits 0" "$NF_RC" 0
  check "the run had a TRAITS_PHASE task" has 'TRAITS_PHASE' "$(cat "$STEP_LOG")"
  for f in GRSPAN_genotype_rules.tsv GRSPAN_traits.tsv GRSPAN_traits_genotypes.tsv; do
    check "Nextflow published traits/${f}" test -s "${NF_OUT}/GRSPAN/traits/${f}"
    check "traits/${f}: Nextflow and the bash step agree byte for byte" \
      cmp -s "${NF_OUT}/GRSPAN/traits/${f}" "${S}/traits/${f}"
  done
  check "control: the comparison tells GRSPAN's rules table from GRNOSPAN's" \
    test "$(cmp -s "${NF_OUT}/GRSPAN/traits/GRSPAN_genotype_rules.tsv" "${G}/GRNOSPAN/traits/GRNOSPAN_genotype_rules.tsv"; echo $?)" -ne 0
  in_genome "$BCFTOOLS_IMAGE" rm -rf nf-genotype-rules
else
  fail "nextflow is not on PATH"
fi

in_genome "$BCFTOOLS_IMAGE" rm -rf GRSPAN GRSWAP GRNOSPAN GRREF reference-genotype-rules
finish
