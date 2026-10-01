#!/usr/bin/env bash
# Step 21 (Cyrius) on a public GIAB HG002 slice: the CYP2D6 locus plus the
# normalisation regions Cyrius reads (pgp-9ms.5). Temporary.
# shellcheck source=lib.sh
. "$(dirname "$0")/lib.sh"
prelude

G=/mnt/scratch/cy
GO=/mnt/scratch/cyold
sudo mkdir -p /mnt/scratch && sudo chown "$(id -u):$(id -g)" /mnt/scratch
mkdir -p "$G/s1/aligned" /tmp/cyrwhl
BAM_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam"

# The regions Cyrius 1.1.1 reads, from its own wheel
curl -sfL -o /tmp/cyrwhl/c.whl https://files.pythonhosted.org/packages/55/b8/83b8fc9ad78718b417905a34380269886e1cb4a5ce7b725a08b45c9990c7/cyrius-1.1.1-py3-none-any.whl
(cd /tmp/cyrwhl && unzip -q -o c.whl cyrius/data/CYP2D6_region_38.bed)
BED=/tmp/cyrwhl/cyrius/data/CYP2D6_region_38.bed
echo "regions: $(wc -l < "$BED")"
cd /tmp && samtools view -H "$BAM_URL" | grep -m3 '^@SQ'
START=$(date +%s)
timeout 2400 samtools view -M -L "$BED" -b -o "$G/s1/aligned/s1_sorted.bam" "$BAM_URL"
samtools index "$G/s1/aligned/s1_sorted.bam"
echo "slice: $(samtools view -c "$G/s1/aligned/s1_sorted.bam") reads, $(du -h "$G/s1/aligned/s1_sorted.bam" | cut -f1), $(( $(date +%s) - START )) s"
mkdir -p "$GO/s1"
cp -r "$G/s1/aligned" "$GO/s1/"

expect_ok "21 new: Cyrius on the HG002 slice" env GENOME_DIR="$G" bash "$NEW/scripts/21-cyrius.sh" s1
TSV="$G/s1/cyrius/s1_cyp2d6.tsv"
check "21 new: result TSV" "$(wc -l < "$TSV" 2>/dev/null || echo missing) lines" "header + 1 row" "$([ "$(wc -l < "$TSV" 2>/dev/null || echo 0)" -ge 2 ] && echo 1 || echo 0)"
expect_fail "21 new: no network" env DOCKER_RUN_EXTRA="--network none" GENOME_DIR="$G" bash "$NEW/scripts/21-cyrius.sh" s1
L="$LOGS/21_new:_no_network.log"
check "21 new: pip's own error is in the log" "$(grep -m1 -E '^ERROR|Could not|Failed to establish' "$L" || echo '<none>')" "a pip error line" "$(grep -qE '^ERROR|Could not|Failed to establish' "$L" && echo 1 || echo 0)"
expect_fail "21 old (origin-main, red-first): Cyrius" env GENOME_DIR="$GO" bash "$OLD/scripts/21-cyrius.sh" s1
grep -m2 -E 'star_caller|not found|127' "$LOGS/21_old_(origin-main,_red-first):_Cyrius.log" || true

finish
