#!/usr/bin/env bash
# setup.sh / validate-setup.sh under a fake docker, then AnnotSV for real
# (pgp-9ms.1, .9, .11). Temporary.
# shellcheck source=lib.sh
. "$(dirname "$0")/lib.sh"
prelude

sudo mkdir -p /mnt/scratch && sudo chown "$(id -u):$(id -g)" /mnt/scratch
G=/mnt/scratch/sg
mkdir -p "$G/reference" "$G/clinvar"
echo ">chr1" > "$G/reference/Homo_sapiens_assembly38.fasta"
printf 'chr1\t1\t6\t60\t61\n' > "$G/reference/Homo_sapiens_assembly38.fasta.fai"
echo "placeholder" > "$G/clinvar/clinvar.vcf.gz"   # raw file present, its .tbi and both subsets missing

# Fake docker: logs every argv (empty arguments show as '') and succeeds
mkdir -p /tmp/fakedocker
cat > /tmp/fakedocker/docker <<'EOF'
#!/usr/bin/env bash
{ printf 'docker'; for a in "$@"; do printf " %q" "$a"; done; printf '\n'; } >> /tmp/fakedocker.log
exit 0
EOF
chmod +x /tmp/fakedocker/docker
# A PATH without wget (stock macOS has curl only)
mkdir -p /tmp/nowget
for d in /usr/local/sbin /usr/local/bin /usr/sbin /usr/bin /sbin /bin; do
  for f in "$d"/*; do
    b=$(basename "$f"); [ "$b" = wget ] && continue
    [ -e "/tmp/nowget/$b" ] || ln -s "$f" "/tmp/nowget/$b"
  done
done
NOWGET_PATH="/tmp/fakedocker:/tmp/nowget"
check "PATH without wget" "$(PATH=$NOWGET_PATH command -v wget || echo 'wget: absent'), $(PATH=$NOWGET_PATH command -v curl)" "wget absent, curl present" "$(PATH=$NOWGET_PATH command -v wget >/dev/null && echo 0 || echo 1)"

# --- setup.sh, first run: no wget, raw ClinVar only; AnnotSV data downloaded for real ---
: > /tmp/fakedocker.log
expect_ok "setup.sh new: no wget, raw ClinVar only, real AnnotSV download" env PATH="$NOWGET_PATH" bash "$NEW/scripts/setup.sh" "$G"
cp /tmp/fakedocker.log "$LOGS/fakedocker-setup1.log"
check "setup.sh new: ClinVar .tbi fetched with curl" "$(ls -l "$G/clinvar/clinvar.vcf.gz.tbi" 2>&1 | awk '{print $5" bytes"}')" "non-empty" "$([ -s "$G/clinvar/clinvar.vcf.gz.tbi" ] && echo 1 || echo 0)"
N=$(grep -c 'bcftools annotate --rename-chrs\|bcftools view -i' /tmp/fakedocker.log)
check "setup.sh new: derived ClinVar builds launched" "$N" "2" "$([ "$N" = 2 ] && echo 1 || echo 0)"
E=$(grep -E "^docker run .* '' " /tmp/fakedocker.log | wc -l)
check "setup.sh new: docker argv with an empty argument" "$E" "0" "$([ "$E" = 0 ] && echo 1 || echo 0)"
grep -E 'Docker images:|AnnotSV' "$LOGS/setup.sh_new:_no_wget,_raw_ClinVar_only,_real_AnnotSV_download.log" || true
check "setup.sh new: AnnotSV annotations in place" "$(ls "$G/annotsv_annotations/Annotations_Human/Genes/" 2>&1 | tr '\n' ' ')" "GRCh38 present" "$([ -d "$G/annotsv_annotations/Annotations_Human/Genes/GRCh38" ] && echo 1 || echo 0)"
echo "AnnotSV annotations on disk: $(du -sh "$G/annotsv_annotations" | cut -f1)"
check "setup.sh new: no tarball or .part left" "$(ls "$G" | tr '\n' ' ')" "no .tar.gz or .part" "$(ls "$G" | grep -qE 'tar.gz|\.part' && echo 0 || echo 1)"

# --- setup.sh, rerun: raw ClinVar + .tbi present, subsets still missing -> both rebuilt ---
: > /tmp/fakedocker.log
expect_ok "setup.sh new: rerun with raw ClinVar and .tbi present" env PATH="$NOWGET_PATH" bash "$NEW/scripts/setup.sh" "$G"
N=$(grep -c 'bcftools annotate --rename-chrs\|bcftools view -i' /tmp/fakedocker.log)
check "setup.sh new rerun: derived ClinVar builds launched" "$N" "2" "$([ "$N" = 2 ] && echo 1 || echo 0)"

# --- validate-setup.sh under the fake docker ---
expect_fail "validate-setup.sh new: fake docker (placeholder FASTA is too small)" env PATH="/tmp/fakedocker:$PATH" GENOME_DIR="$G" bash "$NEW/scripts/validate-setup.sh" s1
L="$LOGS/validate-setup.sh_new:_fake_docker_(placeholder_FASTA_is_too_small).log"
check "validate-setup.sh new: reaches the image list and summary" "$(grep -cE '=== (Docker Images|Summary) ===' "$L") sections, $(grep -c 'unbound variable' "$L") unbound" "2 sections, 0 unbound" "$([ "$(grep -cE '=== (Docker Images|Summary) ===' "$L")" = 2 ] && ! grep -q 'unbound variable' "$L" && echo 1 || echo 0)"
grep -E 'AnnotSV|Docker image' "$L" || true

# --- red-first: origin/main setup.sh and validate-setup.sh ---
expect_fail "setup.sh old (origin-main, red-first): wget present" env PATH="/tmp/fakedocker:$PATH" bash "$OLD/scripts/setup.sh" "$G"
grep -m1 'unbound variable' "$LOGS/setup.sh_old_(origin-main,_red-first):_wget_present.log" || true
rm -f "$G/clinvar/clinvar.vcf.gz.tbi"
expect_fail "setup.sh old (origin-main, red-first): no wget" env PATH="$NOWGET_PATH" bash "$OLD/scripts/setup.sh" "$G"
grep -m2 -E 'wget|ERROR' "$LOGS/setup.sh_old_(origin-main,_red-first):_no_wget.log" || true
expect_fail "validate-setup.sh old (origin-main, red-first)" env PATH="/tmp/fakedocker:$PATH" GENOME_DIR="$G" bash "$OLD/scripts/validate-setup.sh" s1
grep -m1 'unbound variable' "$LOGS/validate-setup.sh_old_(origin-main,_red-first).log" || true

# --- AnnotSV for real on a ten-record Manta-style VCF ---
mkdir -p "$G/s1/manta/results/variants"
{
  echo '##fileformat=VCFv4.1'
  echo '##ALT=<ID=DEL,Description="Deletion">'
  echo '##ALT=<ID=DUP,Description="Duplication">'
  echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
  echo '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">'
  echo '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">'
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n'
  while read -r c s e t; do
    len=$(( e - s )); [ "$t" = DEL ] && len=-$len
    printf '%s\t%s\tMantaSV%s\tN\t<%s>\t100\tPASS\tSVTYPE=%s;END=%s;SVLEN=%s\tGT\t0/1\n' "$c" "$s" "$s" "$t" "$t" "$e" "$len"
  done <<'SV'
chr1 155234000 155240000 DEL
chr5 70925000 70955000 DEL
chr7 117480000 117490000 DEL
chr11 5225000 5227000 DEL
chr13 32320000 32330000 DUP
chr15 48400000 48410000 DEL
chr16 2050000 2060000 DUP
chr17 43050000 43060000 DEL
chr22 42126000 42131000 DEL
chrX 31100000 31200000 DEL
SV
} | bgzip -c > "$G/s1/manta/results/variants/diploidSV.vcf.gz"
tabix -p vcf "$G/s1/manta/results/variants/diploidSV.vcf.gz"
echo "SV records: $(bcftools view -H "$G/s1/manta/results/variants/diploidSV.vcf.gz" | wc -l)"

expect_ok "05 new: AnnotSV with -annotationsDir" env GENOME_DIR="$G" bash "$NEW/scripts/05-annotsv.sh" s1
TSV="$G/s1/annotsv/s1_sv_annotated.tsv"
ls -la "$G/s1/annotsv/" || true
A=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="ACMG_class") c=i; next} c && $c!="" {n++} END{print n+0}' "$TSV" 2>/dev/null || echo 0)
check "05 new: rows with a non-empty ACMG_class" "$A" ">0" "$([ "${A:-0}" -gt 0 ] && echo 1 || echo 0)"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="ACMG_class")c=i; if($i=="Gene_name")g=i}} NR>1 && c{print $2":"$3"-"$4, $g, "ACMG_class="$c}' "$TSV" 2>/dev/null | head -12
mv "$G/annotsv_annotations" "$G/annotsv_annotations.away"
expect_fail "05 new: without annotations" env GENOME_DIR="$G" bash "$NEW/scripts/05-annotsv.sh" s1
grep -m2 -E 'ERROR|setup.sh' "$LOGS/05_new:_without_annotations.log" || true
expect_fail "05 old (origin-main, red-first): AnnotSV without -annotationsDir" env GENOME_DIR="$G" bash "$OLD/scripts/05-annotsv.sh" s1
grep -m3 -iE 'error|annotation|exit' "$LOGS/05_old_(origin-main,_red-first):_AnnotSV_without_-annotationsDir.log" || true
mv "$G/annotsv_annotations.away" "$G/annotsv_annotations"

finish
