#!/usr/bin/env bash
# Step 35's HBA1/HBA2 count from read depth (bin/hba_copy_number.py) on
# synthetic reads, with a planted one-gene deletion. The fixture has no reads
# at the alpha-globin cluster, so the case writes its own: 150 bp reads cut
# from the fixture reference every 5 bp over chr16:95,000-235,000 (30x),
# with the read group of the GIAB slice BAM. A read whose 150 bases occur
# more than once in that stretch (either strand) gets MAPQ 0, as an aligner
# gives a read that fits HBA1 and HBA2 equally; the others MAPQ 60. Two
# samples:
#   HG002hba4  the reads as cut: four alpha-globin genes
#   HG002hba3  every second read over HBA2 left out (one HBA2 copy lost, as
#              -a4.2): three
# Each is merged with the chr5 and chr20 reads of the GIAB slice BAM, so the
# Parascopy part of step 35 runs as in pgx-4-parascopy.sh (its stand-in
# reference with an N-filled chr11, background windows on chr20).
# Then the HBA block of modules/local/parascopy/main.nf, cut from the file,
# runs in the same image on the same BAM: its table must equal step 35's.
# Synthetic reads sit at their true place and carry no GC bias, so this
# shows the counting and the wiring, not the accuracy on a real genome.
# Last, the same block runs on HG002's real Illumina reads at the cluster.
# Paraphase on HG002's HiFi reads calls two HBA1 and two HBA2
# (paralogs-paraphase-1), so the truth is four.
#   HG002novo  GIAB's 60x novoalign BAM (the fixture's source), sampled to
#              about 30x. novoalign left out the reads that fit both genes:
#              no read over them has MAPQ 0, HBA2 sits at about half depth,
#              and before the MAPQ 0 check this counted three (estimate
#              2.76). It must now say not assessed, naming the cause.
#   HG002bwa   Google's 35x NovaSeq PCR-free BAM of HG002, aligned with
#              BWA-MEM to GRCh38 without ALT contigs: four.
#   HG002mm2   the reads of HG002bwa over the cluster (duplicates left out)
#              realigned with step 02's aligner and preset (minimap2 -a -x
#              sr) to chr16:1-300,000 of the fixture reference: four.
. "$(dirname "$0")/lib.sh"

echo "+ scripts/setup.sh --parascopy-data"
"${REPO}/scripts/setup.sh" --parascopy-data "$GENOME_DIR"

REF="${GENOME_DIR}/reference/parascopy_hba_ref.fasta"
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

W="hba_e2e"
mkdir -p "${GENOME_DIR}/${W}"
cp "${FIXTURE_DIR}/${SAMPLE}_slice.bam" "${GENOME_DIR}/${W}/slice.bam"
cp "${FIXTURE_DIR}/${SAMPLE}_slice.bam.bai" "${GENOME_DIR}/${W}/slice.bam.bai"
sam view -b -o "${W}/slice_chr5_chr20.bam" "${W}/slice.bam" chr5 chr20
sam view -H "${W}/slice.bam" > "${GENOME_DIR}/${W}/header.sam"
sam faidx reference/GRCh38_no_alt_analysis_set.fasta chr16:95001-235000 > "${GENOME_DIR}/${W}/chr16.fa"

for P in HG002hba4 HG002hba3; do
  python3 - "${GENOME_DIR}/${W}/chr16.fa" "${GENOME_DIR}/${W}/header.sam" "$P" "${GENOME_DIR}/${W}/${P}_syn.sam" <<'PY'
import sys
fa, header, sample, out = sys.argv[1:]
seq = "".join(l.strip() for l in open(fa) if not l.startswith(">")).upper()
start0 = 95000                      # 0-based start of chr16:95001-235000
rg = next((l.split("\t")[1][3:] for l in open(header) if l.startswith("@RG")), None)
hba2 = (172875, 173710)             # bin/hba_copy_number.py's HBA2, 0-based
L, step, n = 150, 5, 0
comp = str.maketrans("ACGTN", "TGCAN")
seen = {}
for i in range(len(seq) - L + 1):
    seen[seq[i:i + L]] = seen.get(seq[i:i + L], 0) + 1
def mapq(r):                        # 0 for a read that fits more than one place
    return 0 if seen.get(r, 0) + seen.get(r.translate(comp)[::-1], 0) > 1 else 60
with open(out, "w") as f:
    f.write(open(header).read())
    for i in range(0, len(seq) - L + 1, step):
        s = start0 + i
        if sample.endswith("3") and s < hba2[1] and s + L > hba2[0]:
            n += 1
            if n % 2:
                continue
        f.write(f"syn{s}\t0\tchr16\t{s + 1}\t{mapq(seq[i:i + L])}\t{L}M\t*\t0\t0\t{seq[i:i + L]}\t{'I' * L}"
                + (f"\tRG:Z:{rg}" if rg else "") + "\n")
PY
  mkdir -p "${GENOME_DIR}/${P}/aligned"
  sam sort -o "${W}/${P}_syn.bam" "${W}/${P}_syn.sam"
  sam merge -f "${P}/aligned/${P}_sorted.bam" "${W}/slice_chr5_chr20.bam" "${W}/${P}_syn.bam"
  sam index "${P}/aligned/${P}_sorted.bam"
done

BED="${GENOME_DIR}/${W}/parascopy_windows.bed"
awk 'BEGIN { for (s = 10050000; s < 10450000; s += 100) printf "chr20\t%d\t%d\n", s, s + 100 }' > "$BED"

# The module's HBA block, as written in the .nf file, made a plain shell
# script: Groovy's \\ is the shell's \, ${bam} and ${meta.id}
# become the sample's BAM and name.
MOD="${REPO}/modules/local/parascopy/main.nf"
BLOCK="${CASE_TMP}/module_hba.sh"
sed -n '/# HBA1\/HBA2 from depth (as scripts\/35-paralogs.sh)/,/rm -f hba_regions.bed hba_bedcov.tsv/p' "$MOD" \
  | sed 's/\\\\/\\/g' > "$BLOCK"
check_ge "the module's HBA block was cut from main.nf" "$(wc -l < "$BLOCK" | tr -d ' ')" 8

hba_value() { awk -F'\t' -v k="$2" '$1 == k {print $2}' "$1" 2>/dev/null; }

for P in HG002hba4 HG002hba3; do
  want=${P#HG002hba}
  REF_FASTA="$REF" PARASCOPY_DEPTH_BED="$BED" run_step 35-paralogs.sh "$P"
  check_step_exit 35-paralogs.sh
  T="${GENOME_DIR}/${P}/paralogs/${P}_hba_copy_number.tsv"
  cat "$T" 2>/dev/null
  check_eq "${P}: step 35's HBA table says called" "$(hba_value "$T" status)" called
  check_eq "${P}: alpha-globin genes from depth" "$(hba_value "$T" alpha_globin_genes)" "$want"
  check_ge "${P}: reads of MAPQ 0 over the genes, per thousand of the depth" \
    "$(awk -v x="$(hba_value "$T" mapq0_share)" 'BEGIN { printf "%d", x * 1000 }')" 100
  check "${P}: the SMN table is still written" test -s "${GENOME_DIR}/${P}/paralogs/${P}_smn_copy_number.tsv"

  # The module's block in the same image, on the same BAM
  M="${GENOME_DIR}/${W}/${P}_module"
  rm -rf "$M"; mkdir -p "$M"
  sed -e "s|\${bam}|/genome/${P}/aligned/${P}_sorted.bam|g" -e "s|\${meta.id}|${P}|g" "$BLOCK" > "${M}/block.sh"
  docker run --rm --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -v "${REPO}/bin:/pgp-bin:ro" -w "/genome/${W}/${P}_module" \
    "$PARASCOPY_IMAGE" sh -c 'PATH="/pgp-bin:$PATH"; sh -e block.sh' 2>&1 | tail -n 3
  check "${P}: the module's HBA block writes the same table as step 35" cmp -s "$T" "${M}/${P}_hba_copy_number.tsv"
  echo "hba: ${P} $(hba_value "$T" estimate) ($(hba_value "$T" alpha_globin_genes) genes)" >> "$E2E_NOTES"
done

# --- HG002's real Illumina reads at the alpha-globin cluster ------------------
GIAB_BAM="https://giab.s3.amazonaws.com/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam"
BWA_BAM="https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"
sam_https() {  # the image's libcurl has no CA bundle, so lend it the host's
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${GENOME_DIR}:/genome" -w "/genome/${W}" \
    "$SAMTOOLS_IMAGE" samtools "$@"
}
retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && sleep $((i * 15)); done; return 1; }
# real_block NAME: the module's block on hba_e2e/NAME.bam; its table lands in hba_e2e/NAME_module/
real_block() {
  local m="${GENOME_DIR}/${W}/$1_module"
  sam index "${W}/$1.bam"
  check_ge "$1: HG002 Illumina reads at the cluster" "$(sam view -c "${W}/$1.bam" 2>/dev/null)" 10000
  mkdir -p "$m"
  sed -e "s|\${bam}|/genome/${W}/$1.bam|g" -e "s|\${meta.id}|$1|g" "$BLOCK" > "${m}/block.sh"
  docker run --rm --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -v "${REPO}/bin:/pgp-bin:ro" -w "/genome/${W}/$1_module" \
    "$PARASCOPY_IMAGE" sh -c 'PATH="/pgp-bin:$PATH"; sh -e block.sh' 2>&1 | tail -n 3
  cat "${m}/$1_hba_copy_number.tsv" 2>/dev/null
  local t="${m}/$1_hba_copy_number.tsv"
  echo "hba: $1 real reads, $(hba_value "$t" status) $(hba_value "$t" alpha_globin_genes), estimate $(hba_value "$t" estimate);" \
    "HBA2 $(hba_value "$t" hba2_depth)x, HBA1 $(hba_value "$t" hba1_depth)x," \
    "flanks $(hba_value "$t" flank_left_depth)x and $(hba_value "$t" flank_right_depth)x, MAPQ 0 share $(hba_value "$t" mapq0_share)" >> "$E2E_NOTES"
}
real_table() { echo "${GENOME_DIR}/${W}/$1_module/$1_hba_copy_number.tsv"; }

R=HG002novo
retry sam_https view -b -s 42.5 -o "${R}.bam" "$GIAB_BAM" chr16:95001-235000
rm -f "${GENOME_DIR}/${W}/$(basename "$GIAB_BAM").bai"
real_block "$R"
check_eq "${R} (novoalign, no read of MAPQ 0 over the genes): not assessed, never three" \
  "$(hba_value "$(real_table "$R")" status)" "not assessed"
check_eq "${R}: MAPQ 0 share over the genes" "$(hba_value "$(real_table "$R")" mapq0_share)" "0.000"
check "${R}: the message names the cause" grep -q "mapping quality 0.*left out the reads that fit both genes" "$(real_table "$R")"

R=HG002bwa
retry sam_https view -b -o "${R}.bam" "$BWA_BAM" chr16:95001-235000
rm -f "${GENOME_DIR}/${W}/$(basename "$BWA_BAM").bai"
real_block "$R"
check_eq "${R} (BWA-MEM, 35x): called" "$(hba_value "$(real_table "$R")" status)" called
check_eq "${R}: alpha-globin genes from depth (truth: 4)" "$(hba_value "$(real_table "$R")" alpha_globin_genes)" 4

# The same reads through step 02's aligner. A reference of chr16:1-300,000
# keeps GRCh38's coordinates; HBA1 and HBA2 both lie in it, so a read that
# fits both still gets MAPQ 0.
R=HG002mm2
sam faidx reference/GRCh38_no_alt_analysis_set.fasta chr16:1-300000 | sed '1s/.*/>chr16/' > "${GENOME_DIR}/${W}/chr16_300k.fa"
sam collate -u -o "${W}/bwa_collated.bam" "${W}/HG002bwa.bam" "${W}/collate_tmp"
sam fastq -F 0xD00 -n -1 "${W}/mm2_R1.fq" -2 "${W}/mm2_R2.fq" -0 /dev/null -s /dev/null "${W}/bwa_collated.bam"
in_genome "$MINIMAP2_IMAGE" minimap2 -t 2 -a -x sr "${W}/chr16_300k.fa" "${W}/mm2_R1.fq" "${W}/mm2_R2.fq" \
  > "${GENOME_DIR}/${W}/${R}.sam" 2>/dev/null
sam sort -o "${W}/${R}.bam" "${W}/${R}.sam"
real_block "$R"
check_eq "${R} (minimap2 -a -x sr, as step 02): called" "$(hba_value "$(real_table "$R")" status)" called
check_eq "${R}: alpha-globin genes from depth (truth: 4)" "$(hba_value "$(real_table "$R")" alpha_globin_genes)" 4

rm -f "$REF" "${REF}.fai"
rm -rf "${GENOME_DIR:?}/${W}"
finish
