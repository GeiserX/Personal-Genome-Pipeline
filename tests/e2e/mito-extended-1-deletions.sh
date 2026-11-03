#!/usr/bin/env bash
# Step 20's large mtDNA deletions (eKLIPse), mtDNA copy number and confirmed
# variant list, on both paths.
#   1. The image: EKLIPSE_IMAGE was pulled by its digest, and its eKLIPse code
#      is the code of Bioconda's eklipse 1.8 build 6, whose recipe builds
#      upstream commit 3606cb2edac983d2623ddc667b49206c3d01373c.
#   2. assets/mito_confirmed.tsv is what scripts/ci/build-mito-confirmed.sh
#      writes from NCBI's ClinVar of today; a copy with one id changed fails
#      the check (red control).
#   3. Two synthetic chrM BAMs under the fixture's header: 2000x of reads from
#      the fixture's chrM, one with 10% of the molecules carrying the common
#      4977 bp deletion (m.8470 to m.13446, junction reads soft-clipped as an
#      aligner writes them), one without (red control). Step 20 must call one
#      deletion within 50 bp of both breakpoints in the first and none in the
#      second, the report must name it as the common deletion, and the copy
#      number of a BAM with no autosomal reads must be "not assessed".
#   4. Nextflow: MITO_DELETIONS, run for real in its image on the first BAM,
#      writes the same tables as step 20, and HTML_REPORT places them (a name
#      it has no place for stops the task) and shows the deletion.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR" T="$CASE_TMP"
dk() { docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$1" "${@:2}"; }
retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && { echo "  try ${i} of 4 failed; again in $((i * 15)) s"; sleep $((i * 15)); }; done; return 1; }
COMMIT=3606cb2edac983d2623ddc667b49206c3d01373c

# --- 1. the image ---------------------------------------------------------------------
DIGEST=${EKLIPSE_IMAGE##*@}
check "EKLIPSE_IMAGE is pinned by digest (${DIGEST})" has '^sha256:[0-9a-f]{64}$' "$DIGEST"
check "the image was pulled by that digest" \
  has "@${DIGEST}" "$(docker image inspect --format '{{join .RepoDigests " "}}' "$EKLIPSE_IMAGE" 2>/dev/null)"
PKG="${T}/eklipse-1.8-hdfd78af_6.conda"
TAG=${EKLIPSE_IMAGE%@*}; TAG=${TAG##*:}
check_eq "the image is Bioconda's build 1.8--hdfd78af_6" "$TAG" "1.8--hdfd78af_6"
check "Bioconda's package of that build downloads" \
  retry curl -fsSL -o "$PKG" "https://conda.anaconda.org/bioconda/noarch/eklipse-1.8-hdfd78af_6.conda"
rm -rf "${T}/pkg"; mkdir -p "${T}/pkg"
python3 -c 'import sys, zipfile; zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])' "$PKG" "${T}/pkg"
zstd -dc "${T}"/pkg/info-*.tar.zst | tar -x -C "${T}/pkg" info/recipe/meta.yaml
zstd -dc "${T}"/pkg/pkg-*.tar.zst | tar -x -C "${T}/pkg" bin/eKLIPse_sc.py bin/eKLIPse.py
URL=$(awk '$1 == "url:" {print $2; exit}' "${T}/pkg/info/recipe/meta.yaml" 2>/dev/null)
check_eq "the recipe builds upstream commit ${COMMIT:0:7}" "$URL" "https://github.com/dooguypapua/eKLIPse/archive/${COMMIT}.zip"
for f in eKLIPse.py eKLIPse_sc.py; do
  IMG=$(docker run --rm "$EKLIPSE_IMAGE" sh -c "sha256sum \"\$(dirname \"\$(command -v eKLIPse.py)\")/${f}\"" 2>/dev/null | cut -d' ' -f1)
  check_eq "the image's ${f} is the package's" "${IMG:-none}" "$(sha256sum "${T}/pkg/bin/${f}" 2>/dev/null | cut -d' ' -f1)"
done

# --- 2. the confirmed list ---------------------------------------------------------------
CV="${T}/clinvar"
mkdir -p "$CV"
fetch_clinvar() {
  curl -fsSL -o "${CV}/clinvar.vcf.gz" https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz &&
  curl -fsSL -o "${CV}/clinvar.vcf.gz.tbi" https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
}
check "NCBI's ClinVar downloads" retry fetch_clinvar
"${REPO}/scripts/ci/build-mito-confirmed.sh" --clinvar "${CV}/clinvar.vcf.gz" --check > "${T}/build.log" 2>&1
BUILD_RC=$?
cat "${T}/build.log"
check_eq "build-mito-confirmed.sh --check passes on today's ClinVar" "$BUILD_RC" 0
sed 's/\t9589\t/\t9590\t/' "${REPO}/assets/mito_confirmed.tsv" > "${T}/tampered.tsv"
"${REPO}/scripts/ci/build-mito-confirmed.sh" --clinvar "${CV}/clinvar.vcf.gz" --check --out "${T}/tampered.tsv" \
  > "${T}/build-red.log" 2>&1
check_eq "red control: a copy with one ClinVar id changed fails the check" "$?" 1
rm -rf "$CV"

# --- 3. step 20 on two synthetic chrM BAMs ---------------------------------------------
sam faidx reference/GRCh38_no_alt_analysis_set.fasta chrM > "${T}/chrM.fa"
check_ge "the fixture's chrM is whole (bases)" "$(grep -v '^>' "${T}/chrM.fa" | tr -d '\n' | wc -c | tr -d ' ')" 16569
sam view -H "${SAMPLE}/aligned/${SAMPLE}_sorted.bam" > "${T}/header.sam"
DEL="${SAMPLE}mtdel" WT="${SAMPLE}mtwt"
for N in "$DEL" "$WT"; do
  args=()
  [ "$N" = "$DEL" ] && args=(--deletion 8470-13446 --level 0.10)
  python3 "${REPO}/tests/test_mito_extended.py" synth-sam --fasta "${T}/chrM.fa" --header "${T}/header.sam" \
    --sample "$N" --depth 2000 "${args[@]}" --out "${T}/${N}.sam"
  dk "$SAMTOOLS_IMAGE" samtools sort -o "${T}/${N}.bam" "${T}/${N}.sam"
  rm -rf "${G:?}/${N}"; mkdir -p "${G}/${N}/aligned"
  cp "${T}/${N}.bam" "${G}/${N}/aligned/${N}_sorted.bam"
  dk "$SAMTOOLS_IMAGE" samtools index "${G}/${N}/aligned/${N}_sorted.bam"
  rm -f "${T}/${N}.sam"
done
check_ge "soft-clipped junction reads in the deletion BAM" \
  "$(dk "$SAMTOOLS_IMAGE" samtools view "${G}/${DEL}/aligned/${DEL}_sorted.bam" | awk '$6 ~ /S/' | wc -l | tr -d ' ')" 100
check_eq "no soft-clipped read in the control BAM" \
  "$(dk "$SAMTOOLS_IMAGE" samtools view "${G}/${WT}/aligned/${WT}_sorted.bam" | awk '$6 ~ /S/' | wc -l | tr -d ' ')" 0

# near N: deletions of N's table within 50 bp of the common deletion's breakpoints.
near() { awk -F'\t' 'NR > 1 && ($1 - 8469) ^ 2 <= 2500 && ($2 - 13447) ^ 2 <= 2500' "${G}/$1/mito/$1_chrM_deletions.tsv" 2>/dev/null | wc -l | tr -d ' '; }
summary_text() {
  python3 - "${REPO}/bin" "${G}/$1" "$1" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
import collect_summary, render_report
s = collect_summary.collect(sys.argv[3], sys.argv[2])
print(render_report.text_report(s))
PY
}
for N in "$DEL" "$WT"; do
  run_step 20-mtoolbox.sh "$N"
  check_step_exit 20-mtoolbox.sh
  check "${N}: step 20 wrote the deletion table" nonempty "${N}/mito/${N}_chrM_deletions.tsv"
  check "${N}: step 20 wrote the copy-number table" nonempty "${N}/mito/${N}_mtdna_copies.tsv"
  cat "${G}/${N}/mito/${N}_chrM_deletions.tsv" 2>/dev/null
done
check_eq "one deletion within 50 bp of both breakpoints, with 10% of the molecules deleted" "$(near "$DEL")" 1
LEVEL=$(awk -F'\t' 'NR > 1 && ($1 - 8469) ^ 2 <= 2500 && ($2 - 13447) ^ 2 <= 2500 {print $4}' \
  "${G}/${DEL}/mito/${DEL}_chrM_deletions.tsv" 2>/dev/null | head -n 1)
check "its level is above 0 and at most 20% (${LEVEL:-none})" \
  python3 -c "import sys; v = float(sys.argv[1]); sys.exit(not 0 < v <= 20)" "${LEVEL:-x}"
check_eq "red control: no deletion in the BAM without one" \
  "$(awk 'NR > 1' "${G}/${WT}/mito/${WT}_chrM_deletions.tsv" 2>/dev/null | wc -l | tr -d ' ')" 0
check_eq "the copy number of a BAM with no autosomal reads is not assessed" \
  "$(awk -F'\t' '$1 == "reason" {print $2}' "${G}/${DEL}/mito/${DEL}_mtdna_copies.tsv" 2>/dev/null)" low_autosomal_depth
TXT=$(summary_text "$DEL")
grep -E 'Large deletions|common 4977|mtDNA copy number' <<< "$TXT"
check "the report names the common deletion with its level" \
  has '^ +m\.[0-9]+_[0-9]+del \([0-9]+ bp\)  level [0-9.]+%  The common 4977 bp deletion' "$TXT"
check "the report says blood often lacks a deletion" has 'often found in muscle and not in blood' "$TXT"
check "the report says the copy number was not assessed on this BAM" \
  has 'mtDNA copy number: not assessed \(too few autosomal reads' "$TXT"
TXT_WT=$(summary_text "$WT")
check "red control: the report of the BAM without a deletion names none" \
  has 'Large deletions \(eKLIPse\): 0, of them 0 named' "$TXT_WT"

# --- 4. Nextflow: MITO_DELETIONS and HTML_REPORT -------------------------------------------
if command -v nextflow >/dev/null; then
  NP="${T}/project"
  rm -rf "$NP" "${T}/nf"; mkdir -p "$NP"
  cp -R "${REPO}/bin" "${NP}/bin"
  BAM="${G}/${DEL}/aligned/${DEL}_sorted.bam"
  VCF="${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
  cat > "${NP}/main.nf" <<NF
include { MITO_DELETIONS } from '${REPO}/modules/local/mito_variants/main'
include { HTML_REPORT } from '${REPO}/modules/local/html_report/main'

workflow {
    MITO_DELETIONS(channel.of([[id: '${DEL}', sex: 'female'], file('${BAM}'), file('${BAM}.bai')]))
    tables = MITO_DELETIONS.out.deletions.join(MITO_DELETIONS.out.copies)
    HTML_REPORT(tables.map { meta, d, c -> [meta, file('${VCF}'),
        [d, c, file('${G}/${DEL}/mito/${DEL}_chrM_filtered.vcf.gz'), file('${REPO}/assets/mito_confirmed.tsv')]] })
}
NF
  cat > "${NP}/nextflow.config" <<NFCONF
params.outdir = '${T}/nf/results'
params.publish_dir_mode = 'copy'
docker.enabled = true
docker.runOptions = '--user $(id -u):$(id -g) --network none'
NFCONF
  echo "+ nextflow run (MITO_DELETIONS, HTML_REPORT)"
  (cd "$T" && nextflow run "${NP}/main.nf" -ansi-log false -c "${REPO}/conf/containers.config" -work-dir "${T}/nf/work") \
    > "${T}/nf.log" 2>&1
  NF_RC=$?
  check_eq "the Nextflow run exits 0" "$NF_RC" 0
  [ "$NF_RC" = 0 ] || grep -E 'ERROR|Error executing process|Command error' -A 12 "${T}/nf.log" | head -60
  R="${T}/nf/results/${DEL}"
  for f in chrM_deletions.tsv mtdna_copies.tsv; do
    check "MITO_DELETIONS writes the same ${f} as step 20" cmp -s "${R}/mito/${DEL}_${f}" "${G}/${DEL}/mito/${DEL}_${f}"
  done
  check "HTML_REPORT names the common deletion" grep -q 'The common 4977 bp deletion' "${R}/${DEL}_report.html"
else
  fail "nextflow is not on PATH"
fi

rm -rf "${G:?}/${DEL}" "${G:?}/${WT}" "${T}/pkg" "$PKG"
finish
