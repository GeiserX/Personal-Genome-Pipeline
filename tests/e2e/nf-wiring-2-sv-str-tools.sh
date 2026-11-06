#!/usr/bin/env bash
# The Nextflow modules a default run-all.sh run plans for the BAM, with their
# real images: MANTA, DUPHOLD (and DUPHOLD_FILTER), EXPANSION_HUNTER and
# CNVPYTOR (and CNVPYTOR_VCF), through main.nf on the fixture's VCF and BAM.
# Each must publish a readable output with content: Manta's diploid SVs,
# duphold's depth tags on them, ExpansionHunter's repeat calls, and CNVpytor's
# VCF with its index. ExpansionHunter gets the bundled catalog cut to the
# fixture's contigs (as bash-step-09-str does); CNVpytor the pinned v1.3.2
# GC and mask files docs/00-reference-setup.md downloads.
# Control: a copy of the pipeline whose EXPANSION_HUNTER passes a wrong flag
# fails the run, so a broken module command line cannot ship green.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

G="$GENOME_DIR"
REF="${G}/reference/GRCh38_no_alt_analysis_set.fasta"
VCF="${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
BAM="${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
for f in "$REF" "${REF}.fai" "$VCF" "${VCF}.tbi" "$BAM" "${BAM}.bai"; do
  check "the base cases left ${f#"$G"/}" test -s "$f"
done

# ExpansionHunter's bundled GRCh38 catalog, cut to the loci on the fixture's
# contigs: it stops on a locus whose contig the reference lacks.
CAT="${CASE_TMP}/eh_catalog.json"
docker run --rm "$EXPANSIONHUNTER_IMAGE" cat /usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json \
  > "${CASE_TMP}/catalog_full.json"
python3 - "${CASE_TMP}/catalog_full.json" "${REF}.fai" "$CAT" <<'PY'
import json, sys
loci = json.load(open(sys.argv[1]))
contigs = {line.split("\t")[0] for line in open(sys.argv[2])}
def regions(locus):
    out = []
    for key in ("ReferenceRegion", "OfftargetRegions"):
        r = locus.get(key, [])
        out += r if isinstance(r, list) else [r]
    return out
def known(c):
    return c in contigs or "chr" + c in contigs
kept = [l for l in loci if all(known(r.split(":")[0]) for r in regions(l))]
print(f"catalog loci: {len(loci)}, on the fixture's contigs: {len(kept)}")
json.dump(kept, open(sys.argv[3], "w"), indent=1)
PY
check_ge "catalog loci on the fixture's contigs" \
  "$(python3 -c 'import json, sys; print(len(json.load(open(sys.argv[1]))))' "$CAT" 2>/dev/null || echo 0)" 5

# CNVpytor's pinned resources (docs/00-reference-setup.md): hg38's GC and mask
# files are read; the others only have to exist.
RES="${CASE_TMP}/cnvpytor"
mkdir -p "$RES"
BASE=https://github.com/abyzovlab/CNVpytor/raw/v1.3.2/cnvpytor/data
for f in gc_hg38.pytor mask_hg38.pytor gc_hg19.pytor mask_hg19.pytor gc_chm13v2.0.pytor gc_chm13v1.1.pytor gc_kn99.pytor; do
  curl -fsSL --retry 3 -o "${RES}/${f}" "${BASE}/${f}" || fail "download of ${f}"
done
for f in gc_hg38.pytor mask_hg38.pytor; do
  check "${f} is an HDF5 file" sh -c "head -c 8 '${RES}/${f}' | od -An -c | tr -d ' ' | grep -q '^211HDF'"
done

SHEET="${CASE_TMP}/sheet.csv"
printf 'sample,vcf,vcf_index,bam,bam_index,sex\n%s,%s,%s,%s,%s,male\n' "$SAMPLE" "$VCF" "${VCF}.tbi" "$BAM" "${BAM}.bai" > "$SHEET"
OUT_REL="nf-wiring-sv"
OUT="${G}/${OUT_REL}"
# nf_run DIR NAME TOOLS: main.nf of DIR with the real images on the fixture
nf_run() {
  rm -rf "$OUT"
  echo "+ nextflow run main.nf -profile docker --tools $3 ($2)"
  (
    cd "$CASE_TMP" && nextflow run "$1/main.nf" -profile docker -ansi-log false \
      -work-dir "${E2E_WORK}/nf-work-$2" --input "$SHEET" --reference "$REF" --tools "$3" \
      --expansion_catalog "$CAT" --cnvpytor_resources "$RES" \
      --outdir "$OUT" --max_cpus 4 --max_memory 14.GB
  ) > "${CASE_TMP}/$2.log" 2>&1
}
nf_run "$REPO" tools manta,duphold,expansion_hunter,cnvpytor
NF_RC=$?
grep -vE 'Pulling|Waiting|Verifying|Download complete|Pull complete|Already exists' "${CASE_TMP}/tools.log" | tail -n 30
check_eq "nextflow run exits 0" "$NF_RC" 0
[ "$NF_RC" = 0 ] || grep -E 'Error executing process|Command exit status|Command error' -A 10 "${CASE_TMP}/tools.log" | head -60
cp "${CASE_TMP}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
for p in MANTA DUPHOLD DUPHOLD_FILTER EXPANSION_HUNTER CNVPYTOR CNVPYTOR_VCF; do
  check "${p} ran" grep -Eq "Submitted process > [A-Z_:]*:${p} \(${SAMPLE}\)" "${CASE_TMP}/tools.log"
done

R="${OUT_REL}/${SAMPLE}"
MANTA_VCF="${R}/manta/results/variants/diploidSV.vcf.gz"
check "MANTA writes a readable diploid SV VCF" vcf_ok "$MANTA_VCF"
N_SV=$(vcf_count "$MANTA_VCF")
echo "Manta SVs: ${N_SV}"
check_ge "Manta calls SVs on the fixture" "$N_SV" 1
DH="${G}/${R}/sv_duphold/${SAMPLE}_sv_duphold.vcf"
check "DUPHOLD declares its depth tag DHFFC" grep -q '^##FORMAT=<ID=DHFFC' "$DH"
check_eq "DUPHOLD keeps every Manta record" "$(grep -vc '^#' "$DH" 2>/dev/null || echo 0)" "$N_SV"
check_ge "and tags them with DHFFC" "$(grep -v '^#' "$DH" 2>/dev/null | grep -c 'DHFFC' || true)" 1
check "DUPHOLD_FILTER writes a readable VCF" vcf_ok "${R}/sv_filtered/${SAMPLE}_sv_filtered.vcf.gz"
check "and its index" nonempty "${R}/sv_filtered/${SAMPLE}_sv_filtered.vcf.gz.tbi"
EH="${G}/${R}/expansion_hunter/${SAMPLE}_eh.vcf"
check_ge "EXPANSION_HUNTER calls repeats" "$(grep -vc '^#' "$EH" 2>/dev/null || true)" 1
check_ge "with a repeat count (REPCN) on them" "$(grep -v '^#' "$EH" 2>/dev/null | grep -c 'REPCN' || true)" 1
check "EXPANSION_HUNTER publishes the catalog it ran" cmp -s "$CAT" "${G}/${R}/expansion_hunter/${SAMPLE}_eh_catalog.json"
CNV="${R}/cnvpytor/${SAMPLE}_cnvs.vcf.gz"
check "CNVPYTOR_VCF writes a readable VCF" vcf_ok "$CNV"
check "and its index" nonempty "${CNV}.tbi"
check "CNVPYTOR writes its call table" test -e "${G}/${R}/cnvpytor/${SAMPLE}_cnvs.txt"
check "the CNV VCF names the fixture's contigs" grep -q 'contig=<ID=chr20' <<< "$(bcf view -h "$CNV" 2>/dev/null)"
echo "CNVpytor calls: $(vcf_count "$CNV"), table rows: $(wc -l < "${G}/${R}/cnvpytor/${SAMPLE}_cnvs.txt" 2>/dev/null | tr -d ' ')" >> "$E2E_NOTES"
echo "Manta SVs: ${N_SV}, EH loci: $(grep -vc '^#' "$EH" 2>/dev/null || true)" >> "$E2E_NOTES"

# --- control: a broken EXPANSION_HUNTER command line fails the run ---------------
BROKEN="${CASE_TMP}/broken"
rm -rf "$BROKEN"
mkdir -p "$BROKEN"
cp -R "${REPO}/main.nf" "${REPO}/nextflow.config" "${REPO}/nextflow_schema.json" "${REPO}/versions.env" \
  "${REPO}/conf" "${REPO}/modules" "${REPO}/workflows" "${REPO}/bin" "${REPO}/assets" "$BROKEN/"
sed -i 's/--reads \${bam}/--readz ${bam}/' "${BROKEN}/modules/local/expansion_hunter/main.nf"
check_eq "control copy: EXPANSION_HUNTER passes --readz" \
  "$(grep -c -- '--readz \${bam}' "${BROKEN}/modules/local/expansion_hunter/main.nf")" 1
nf_run "$BROKEN" broken expansion_hunter
BR=$?
check "control: the run with the broken flag fails" test "$BR" -ne 0
check "control: at EXPANSION_HUNTER" grep -q "Error executing process > '.*EXPANSION_HUNTER" "${CASE_TMP}/broken.log"

rm -rf "$OUT" "$BROKEN" "$RES" "${E2E_WORK}/nf-work-tools" "${E2E_WORK}/nf-work-broken"
finish
