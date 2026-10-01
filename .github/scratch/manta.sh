#!/usr/bin/env bash
# Manta and Strelka2 run twice on the same sample (pgp-9ms.7). Temporary.
# shellcheck source=lib.sh
. "$(dirname "$0")/lib.sh"
prelude

G=/mnt/scratch/mnew
GO=/mnt/scratch/mold
sudo mkdir -p /mnt/scratch && sudo chown "$(id -u):$(id -g)" /mnt/scratch
mkdir -p "$G/reference" "$G/s1/fastq"
python3 "$SCRATCH/gen_data.py" ref "$G/reference/Homo_sapiens_assembly38.fasta" chr20:200000
samtools faidx "$G/reference/Homo_sapiens_assembly38.fasta"
python3 "$SCRATCH/gen_data.py" pairs "$G/reference/Homo_sapiens_assembly38.fasta" chr20 40000 "$G/s1/fastq/s1_R1.fastq.gz" "$G/s1/fastq/s1_R2.fastq.gz"
expect_ok "02 new: align for Manta" env GENOME_DIR="$G" THREADS=4 bash "$NEW/scripts/02-alignment.sh" s1
mkdir -p "$GO/s1"
cp -r "$G/reference" "$GO/"
cp -r "$G/s1/aligned" "$GO/s1/"

expect_ok "04 new: Manta first run" env GENOME_DIR="$G" bash "$NEW/scripts/04-manta.sh" s1
check "04 new: diploidSV index after first run" "$(ls "$G/s1/manta/results/variants/" 2>&1 | tr '\n' ' ')" "diploidSV.vcf.gz.tbi present" "$([ -f "$G/s1/manta/results/variants/diploidSV.vcf.gz.tbi" ] && echo 1 || echo 0)"
expect_ok "04 new: Manta second run" env GENOME_DIR="$G" bash "$NEW/scripts/04-manta.sh" s1
L="$LOGS/04_new:_Manta_second_run.log"
check "04 new: second run logs already done" "$(grep -m1 'already done' "$L" || echo '<none>')" "Manta already done" "$(grep -q 'already done' "$L" && echo 1 || echo 0)"
# Interrupted run: results gone, workflow still there -> resume without reconfiguring
sudo rm -rf "$G/s1/manta/results"
observe "04 new: resume after results removed" env GENOME_DIR="$G" bash "$NEW/scripts/04-manta.sh" s1
L="$LOGS/04_new:_resume_after_results_removed.log"
check "04 new: resume skips configManta" "$(grep -m1 'resuming' "$L" || echo '<none>')" "resuming" "$(grep -q 'resuming' "$L" && echo 1 || echo 0)"
tail -3 "$L"
# Leftover runDir with neither workflow nor results -> cleared, then configured
sudo rm -rf "$G/s1/manta" && mkdir -p "$G/s1/manta/workspace" && touch "$G/s1/manta/workspace/junk"
expect_ok "04 new: leftover runDir is cleared" env GENOME_DIR="$G" bash "$NEW/scripts/04-manta.sh" s1

expect_ok "04 old (origin-main): Manta first run" env GENOME_DIR="$GO" bash "$OLD/scripts/04-manta.sh" s1
expect_fail "04 old (origin-main, red-first): Manta second run" env GENOME_DIR="$GO" bash "$OLD/scripts/04-manta.sh" s1
grep -m2 -iE 'error|already' "$LOGS/04_old_(origin-main,_red-first):_Manta_second_run.log" || true

expect_ok "03c new: Strelka2 first run" env GENOME_DIR="$G" THREADS=4 bash "$NEW/scripts/03c-strelka2-germline.sh" s1
expect_ok "03c new: Strelka2 second run" env GENOME_DIR="$G" THREADS=4 bash "$NEW/scripts/03c-strelka2-germline.sh" s1
L="$LOGS/03c_new:_Strelka2_second_run.log"
check "03c new: second run logs already done" "$(grep -m1 'already done' "$L" || echo '<none>')" "Strelka2 already done" "$(grep -q 'already done' "$L" && echo 1 || echo 0)"
expect_ok "03c old (origin-main): Strelka2 first run" env GENOME_DIR="$GO" THREADS=4 bash "$OLD/scripts/03c-strelka2-germline.sh" s1
expect_fail "03c old (origin-main, red-first): Strelka2 second run" env GENOME_DIR="$GO" THREADS=4 bash "$OLD/scripts/03c-strelka2-germline.sh" s1

finish
