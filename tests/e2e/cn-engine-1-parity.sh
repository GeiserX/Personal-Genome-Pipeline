#!/usr/bin/env bash
# The copy-number engine (bin/copy_number.py, bin/cn_loci.tsv) against the
# per-gene code it replaced, on real reads: HG002's Illumina reads at the
# alpha-globin cluster (chr16:95,001-235,000), the RH locus
# (chr1:25.15-25.45 Mb) and LPA (chr6:160.53-160.67 Mb) from the 35x BWA-MEM
# GRCh38 no-ALT BAM of DeepVariant's case study, written again under the
# fixture's header, at full depth and at an eighth of it. For each BAM:
#   - step 46 of this tree and of main write the same lpa table, byte for byte;
#   - step 41 of this tree and of main write the same blood group table and
#     the same evidence file;
#   - step 35's HBA block (bedcov all MAPQ and MAPQ >= 1, then the verdict)
#     writes the same table through copy_number.py as main's
#     hba_copy_number.py does.
# main is fetched at run time. Once main itself has bin/copy_number.py there
# is no older code left to compare with: the case says so and checks only
# that this tree's tables were written.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR" T="$CASE_TMP"
SRC="https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"
FIX="${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
dk() { docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$1" "${@:2}"; }
retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && { echo "  try ${i} of 4 failed; again in $((i * 15)) s"; sleep $((i * 15)); }; done; return 1; }

# --- main's tree --------------------------------------------------------------------
BASE="${T}/main"
rm -rf "$BASE"; mkdir -p "$BASE"
check "main fetches" retry git -C "$REPO" fetch --quiet --depth 1 origin main
git -C "$REPO" archive FETCH_HEAD bin scripts assets versions.env | tar -x -C "$BASE"
check "main's tree has the step scripts" test -x "${BASE}/scripts/46-lpa.sh"
OLD=true
[ -e "${BASE}/bin/copy_number.py" ] && OLD=false
$OLD || echo "cn-engine: main already has bin/copy_number.py; no older code to compare with" >> "$E2E_NOTES"

# --- HG002's reads at the three loci, on the fixture's header --------------------------
fetch_reads() {
  rm -f "${T}/win.bam" "${T}/$(basename "$SRC").bai"
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$SAMTOOLS_IMAGE" \
    samtools view -b -o "${T}/win.bam" "$SRC" chr16:95001-235000 chr1:25150000-25450000 chr6:160530000-160670000
}
check "HG002's reads at the three loci download" retry fetch_reads
{ dk "$SAMTOOLS_IMAGE" samtools view -H "$FIX" | grep -v -e '^@RG' -e '^@PG'
  dk "$SAMTOOLS_IMAGE" samtools view -H "${T}/win.bam" | grep '^@RG'
  dk "$SAMTOOLS_IMAGE" samtools view "${T}/win.bam" \
    | awk -F'\t' '$7 == "=" || $7 == "*" || $7 == "chr1" || $7 == "chr6" || $7 == "chr16"'
} > "${T}/cn.sam"
dk "$SAMTOOLS_IMAGE" samtools sort -o "${T}/full.bam" "${T}/cn.sam"
dk "$SAMTOOLS_IMAGE" samtools view -b -s 42.125 -o "${T}/thin.bam" "${T}/full.bam"
rm -f "${T}/cn.sam" "${T}/win.bam"
check_ge "HG002 reads at the three loci" "$(dk "$SAMTOOLS_IMAGE" samtools view -c "${T}/full.bam" 2>/dev/null || echo 0)" 50000

# dkb BIN IMAGE CMD...: dk with a tree's bin/ at its own path (this checkout is outside the work folder).
dkb() { docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${E2E_WORK}:${E2E_WORK}" -v "$1:$1:ro" -w "$T" "$2" "${@:3}"; }
# hba_block TREE OUT BAM: step 35's HBA block with TREE's bin (main: hba_copy_number.py).
hba_block() {
  local bin="$1/bin" out=$2 bam=$3 d="${T}/hba.$$"
  rm -rf "$d"; mkdir -p "$d"
  if [ -e "${bin}/copy_number.py" ]; then
    dkb "$bin" "$PARASCOPY_IMAGE" python3 "${bin}/copy_number.py" bed --locus HBA > "${d}/r.bed"
  else
    dkb "$bin" "$PARASCOPY_IMAGE" python3 "${bin}/hba_copy_number.py" bed > "${d}/r.bed"
  fi
  dkb "$bin" "$PARASCOPY_IMAGE" samtools bedcov "${d}/r.bed" "$bam" > "${d}/c.tsv"
  dkb "$bin" "$PARASCOPY_IMAGE" samtools bedcov -Q 1 "${d}/r.bed" "$bam" > "${d}/q.tsv"
  if [ -e "${bin}/copy_number.py" ]; then
    dkb "$bin" "$PARASCOPY_IMAGE" python3 "${bin}/copy_number.py" call --locus HBA --bedcov "${d}/c.tsv" \
      --bedcov-mapq1 "${d}/q.tsv" --table "$out"
  else
    dkb "$bin" "$PARASCOPY_IMAGE" python3 "${bin}/hba_copy_number.py" check --bedcov "${d}/c.tsv" \
      --bedcov-mapq1 "${d}/q.tsv" --out "$out"
  fi
  rm -rf "$d"
}

for V in full thin; do
  N="${SAMPLE}cn${V}"
  rm -rf "${G:?}/${N}"; mkdir -p "${G}/${N}/aligned" "${G}/${N}/vcf"
  cp "${T}/${V}.bam" "${G}/${N}/aligned/${N}_sorted.bam"
  dk "$SAMTOOLS_IMAGE" samtools index "${G}/${N}/aligned/${N}_sorted.bam"
  cp "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${G}/${N}/vcf/${N}.vcf.gz"
  cp "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" "${G}/${N}/vcf/${N}.vcf.gz.tbi"
  NEW="${T}/new-${V}" OLDO="${T}/old-${V}"
  rm -rf "$NEW" "$OLDO"; mkdir -p "$NEW" "$OLDO"

  run_step 46-lpa.sh "$N"; check_step_exit 46-lpa.sh
  cp "${G}/${N}/lpa/${N}_lpa_kiv2.tsv" "${NEW}/lpa.tsv"
  run_step 41-bloodtype.sh "$N"; check_step_exit 41-bloodtype.sh
  cp "${G}/${N}/bloodtype/${N}_bloodtype.tsv" "${NEW}/bloodtype.tsv"
  cp "${G}/${N}/bloodtype/${N}_bloodtype_evidence.tsv" "${NEW}/evidence.tsv"
  hba_block "$REPO" "${NEW}/hba.tsv" "${G}/${N}/aligned/${N}_sorted.bam"
  for f in lpa bloodtype evidence hba; do check "${V}: this tree wrote its ${f} table" test -s "${NEW}/${f}.tsv"; done
  cat "${NEW}/lpa.tsv" "${NEW}/hba.tsv"; grep '^RhD' "${NEW}/bloodtype.tsv"
  echo "cn-engine: ${V} LPA $(awk -F'\t' '$1 == "kiv2_copies" {print $2}' "${NEW}/lpa.tsv")," \
    "HBA $(awk -F'\t' '$1 == "alpha_globin_genes" {print $2}' "${NEW}/hba.tsv")," \
    "RhD $(awk -F'\t' '$1 == "RhD" {print $3}' "${NEW}/bloodtype.tsv")" >> "$E2E_NOTES"
  $OLD || continue

  echo "+ main's scripts/46-lpa.sh ${N}"
  "${BASE}/scripts/46-lpa.sh" "$N" 2>&1 | tail -n 3
  cp "${G}/${N}/lpa/${N}_lpa_kiv2.tsv" "${OLDO}/lpa.tsv"
  echo "+ main's scripts/41-bloodtype.sh ${N}"
  "${BASE}/scripts/41-bloodtype.sh" "$N" 2>&1 | tail -n 3
  cp "${G}/${N}/bloodtype/${N}_bloodtype.tsv" "${OLDO}/bloodtype.tsv"
  cp "${G}/${N}/bloodtype/${N}_bloodtype_evidence.tsv" "${OLDO}/evidence.tsv"
  hba_block "$BASE" "${OLDO}/hba.tsv" "${G}/${N}/aligned/${N}_sorted.bam"
  for f in lpa bloodtype evidence hba; do
    check "${V}: main wrote its ${f} table" test -s "${OLDO}/${f}.tsv"
    check "${V}: the ${f} table is main's, byte for byte" cmp -s "${NEW}/${f}.tsv" "${OLDO}/${f}.tsv"
    cmp -s "${NEW}/${f}.tsv" "${OLDO}/${f}.tsv" || diff "${OLDO}/${f}.tsv" "${NEW}/${f}.tsv" | head -20
  done
done
check "control: the comparison sees a difference (full against thin)" bash -c '! cmp -s "$1" "$2"' _ \
  "${T}/new-full/lpa.tsv" "${T}/new-thin/lpa.tsv"

rm -rf "${G:?}/${SAMPLE}cnfull" "${G:?}/${SAMPLE}cnthin" "${T:?}"/*.bam "$BASE"
finish
