#!/usr/bin/env bash
# Read groups, sample name, unaligned BAM input (pgp-9ms.2, pgp-9ms.10). Temporary.
# shellcheck source=lib.sh
. "$(dirname "$0")/lib.sh"
prelude

G=/mnt/scratch/gnew
GO=/mnt/scratch/gold
sudo mkdir -p /mnt/scratch && sudo chown "$(id -u):$(id -g)" /mnt/scratch
mkdir -p "$G/reference" "$G/s1/fastq"
python3 "$SCRATCH/gen_data.py" ref "$G/reference/Homo_sapiens_assembly38.fasta" chr20:120000 chrM:16569
samtools faidx "$G/reference/Homo_sapiens_assembly38.fasta"
python3 "$SCRATCH/gen_data.py" pairs "$G/reference/Homo_sapiens_assembly38.fasta" chr20 4000 /tmp/a1.fq.gz /tmp/a2.fq.gz
python3 "$SCRATCH/gen_data.py" pairs "$G/reference/Homo_sapiens_assembly38.fasta" chrM 2000 /tmp/m1.fq.gz /tmp/m2.fq.gz
cat /tmp/a1.fq.gz /tmp/m1.fq.gz > "$G/s1/fastq/s1_R1.fastq.gz"
cat /tmp/a2.fq.gz /tmp/m2.fq.gz > "$G/s1/fastq/s1_R2.fastq.gz"
# Ten long reads as an unaligned BAM (the usual PacBio HiFi delivery)
python3 "$SCRATCH/gen_data.py" long "$G/reference/Homo_sapiens_assembly38.fasta" chr20 10 /tmp/lr.fq
samtools import -0 /tmp/lr.fq -o "$G/s1/fastq/s1.bam"
echo "unaligned BAM records: $(samtools view -c "$G/s1/fastq/s1.bam")"
cp -r "$G" "$GO"

hdr_rg() { samtools view -H "$1" | grep '^@RG' || echo "(no @RG line)"; }

# --- Step 02 (minimap2) ---
expect_ok "02 new: align FASTQ" env GENOME_DIR="$G" THREADS=4 bash "$NEW/scripts/02-alignment.sh" s1
RG=$(hdr_rg "$G/s1/aligned/s1_sorted.bam"); echo "$RG"
check "02 new: BAM has @RG with SM:s1" "$RG" "@RG ... SM:s1" "$(grep -q $'\tSM:s1' <<< "$RG" && echo 1 || echo 0)"
expect_ok "02 old (origin-main): align FASTQ" env GENOME_DIR="$GO" THREADS=4 bash "$OLD/scripts/02-alignment.sh" s1
RG=$(hdr_rg "$GO/s1/aligned/s1_sorted.bam"); echo "$RG"
check "02 old (red-first): BAM has no @RG" "$RG" "no @RG line" "$(grep -q '^@RG' <<< "$RG" && echo 0 || echo 1)"

# --- Step 02a (BWA-MEM2) ---
expect_ok "02a new: align FASTQ with BWA-MEM2" env GENOME_DIR="$G" THREADS=4 bash "$NEW/scripts/02a-alignment-bwamem2.sh" s1
RG=$(hdr_rg "$G/s1/aligned_bwamem2/s1_sorted.bam"); echo "$RG"
check "02a new: BAM has @RG with SM:s1" "$RG" "@RG ... SM:s1" "$(grep -q $'\tSM:s1' <<< "$RG" && echo 1 || echo 0)"

# --- Step 02b (long reads from an unaligned BAM) ---
expect_ok "02b new: align 10-read unaligned BAM" env GENOME_DIR="$G" PLATFORM=hifi THREADS=4 bash "$NEW/scripts/02b-alignment-longread.sh" s1
N=$(samtools view -c "$G/s1/aligned_longread/s1_sorted.bam" 2>/dev/null || echo 0)
check "02b new: output BAM records" "$N" "10" "$([ "$N" = 10 ] && echo 1 || echo 0)"
RG=$(hdr_rg "$G/s1/aligned_longread/s1_sorted.bam"); echo "$RG"
check "02b new: @RG PL:PACBIO" "$RG" "PL:PACBIO" "$(grep -q 'PL:PACBIO' <<< "$RG" && echo 1 || echo 0)"
observe "02b old (origin-main): align 10-read unaligned BAM" env GENOME_DIR="$GO" PLATFORM=hifi THREADS=4 bash "$OLD/scripts/02b-alignment-longread.sh" s1
N=$(samtools view -c "$GO/s1/aligned_longread/s1_sorted.bam" 2>/dev/null || echo "no BAM")
check "02b old (red-first): output BAM records" "$N" "not 10" "$([ "$N" != 10 ] && echo 1 || echo 0)"

# --- Step 03 (DeepVariant sample name) ---
expect_ok "03 new: DeepVariant" env GENOME_DIR="$G" bash "$NEW/scripts/03-deepvariant.sh" s1
SN=$(bcftools query -l "$G/s1/vcf/s1.vcf.gz" 2>&1)
check "03 new: VCF sample column" "$SN" "s1" "$([ "$SN" = s1 ] && echo 1 || echo 0)"
echo "records: $(bcftools view -H "$G/s1/vcf/s1.vcf.gz" | wc -l)"
expect_ok "03 old (origin-main): DeepVariant" env GENOME_DIR="$GO" bash "$OLD/scripts/03-deepvariant.sh" s1
SN=$(bcftools query -l "$GO/s1/vcf/s1.vcf.gz" 2>&1)
check "03 old (red-first): VCF sample column" "$SN" "not s1" "$([ "$SN" != s1 ] && echo 1 || echo 0)"

# --- Step 20 (GATK Mutect2 on chrM needs a read group) ---
expect_ok "20 new: Mutect2 chrM on BAM with @RG" env GENOME_DIR="$G" bash "$NEW/scripts/20-mtoolbox.sh" s1
observe "20 old (origin-main): Mutect2 chrM on BAM without @RG" env GENOME_DIR="$GO" bash "$OLD/scripts/20-mtoolbox.sh" s1
grep -m3 -iE 'read group|SM tag|A USER ERROR' "$LOGS/20_old_(origin-main):_Mutect2_chrM_on_BAM_without_@RG.log" || true

finish
