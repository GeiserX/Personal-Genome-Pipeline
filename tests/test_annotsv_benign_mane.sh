#!/usr/bin/env bash
# Step 05 (AnnotSV) with the two files setup.sh adds to AnnotSV's data, on a
# made-up annotation set small enough for the unit-test job, with the real
# AnnotSV image (ANNOTSV_IMAGE).
#
# The annotation set holds what AnnotSV needs to rank losses: genes (RefSeq and
# Ensembl), one benign and one pathogenic record per SV type, and a gene table
# with ClinGen's HI score. On chr20:
#   FIXG1   a haploinsufficient gene (HI 3) with two transcripts: NM_900001,
#           its MANE Select, and XM_900002, which has one more coding exon of
#           100 bp inside NM_900001's first intron;
#   G01-G40 forty one-exon genes;
#   FARG1   far from both, under the benign records of the set.
# Two deletions:
#   DEL1  chr20:1003001-1004000, inside NM_900001's first intron and over
#         XM_900002's extra exon. On XM_900002 it is a frameshift (2E-1,
#         +0.90, class 4); on NM_900001 it is intronic;
#   DEL2  chr20:2000001-2450000 over the forty genes (3C, +0.90, class 4).
# A gnomAD-SV-style sites file lists a common deletion around each one, and a
# MANE-style summary names NM_900001 as FIXG1's MANE Select.
#
# Checks:
#   1. setup.sh --annotsv-data turns the two files into what AnnotSV reads:
#      it runs AnnotSV once, as the caller, so the three records AnnotSV can
#      use out of nine (PASS, DEL/DUP/INS/INV, AN >= 1000, AF >= 0.001,
#      N_HOMALT >= 5) are already merged into the benign files, owned by the
#      caller, before any sample run; and the one MANE Select RefSeq id; a
#      second run changes nothing;
#   2. step 05: DEL1 is reported on NM_900001 with no frameshift and gets the
#      gnomAD frequency, class 1; DEL2 gets the gnomAD frequency and leaves
#      class 4 for class 3;
#   3. red controls with the same image and data: without -txFile DEL1 is
#      back on XM_900002 as a frameshift in class 4, and without gnomAD-SV
#      DEL2 stays class 4 with no frequency;
#   4. validate-setup.sh passes both files and fails on each one removed, and
#      step 05 stops when either is missing.
# Needs docker.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
# shellcheck source=../versions.env
. "${REPO}/versions.env"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }
# expect DESC GOT WANT: GOT matches the extended regex WANT.
expect() {
  if [[ "$2" =~ $3 ]]; then pass "$1 ($2)"; else fail "$1: got '$2', want /$3/"; fi
}

command -v docker >/dev/null 2>&1 || { echo "FAIL: docker is needed to run this test"; exit 1; }

G="${WORK}/genome"
ANN="${G}/annotsv_annotations"
H="${ANN}/Annotations_Human"

# build_annotations DIR: the made-up annotation set in DIR/Annotations_Human.
build_annotations() {
  local h="$1/Annotations_Human" b t i s
  mkdir -p "${h}/Genes/GRCh38" "${h}/Gene-based/Fixture"
  # chrom txStart txEnd strand gene tx cdsStart cdsEnd exonStarts exonEnds
  {
    printf '20\t1000000\t1010000\t+\tFIXG1\tNM_900001\t1000100\t1009900\t1000000,1006000,1009800,\t1000200,1006100,1010000,\n'
    printf '20\t1000000\t1010000\t+\tFIXG1\tXM_900002\t1000100\t1009900\t1000000,1003500,1006000,1009800,\t1000200,1003600,1006100,1010000,\n'
    for i in $(seq 1 40); do
      s=$((2000000 + i * 10000))
      printf '20\t%d\t%d\t+\tG%02d\tNM_91%04d\t%d\t%d\t%d,\t%d,\n' "$s" $((s + 2000)) "$i" "$i" $((s + 100)) $((s + 1900)) "$s" $((s + 2000))
    done
    printf '20\t50000000\t50020000\t+\tFARG1\tNM_920001\t50000100\t50019900\t50000000,\t50020000,\n'
  } > "${h}/Genes/GRCh38/genes.RefSeq.sorted.bed"
  sed 's/\t\(NM\|XM\)_9/\tENST9/' "${h}/Genes/GRCh38/genes.RefSeq.sorted.bed" > "${h}/Genes/GRCh38/genes.ENSEMBL.sorted.bed"
  cut -f6 "${h}/Genes/GRCh38/genes.RefSeq.sorted.bed" | awk '{print $1 "\t1"}' > "${h}/Genes/GRCh38/transcript_version.RefSeq.tsv"
  cut -f6 "${h}/Genes/GRCh38/genes.ENSEMBL.sorted.bed" | awk '{print $1 "\t1"}' > "${h}/Genes/GRCh38/transcript_version.ENSEMBL.tsv"
  printf 'genes\tHI\tTS\tGnomAD_pLI\tLOEUF_bin\tDDD_HI_percent\tOMIM_morbid\nFIXG1\t3\t0\t0.10\t5\t50\tno\n' \
    > "${h}/Gene-based/Fixture/fixture_genes.tsv"
  # AnnotSV reads gene aliases from this table before it annotates anything.
  mkdir -p "${h}/Gene-based/NCBIandHGNCgeneID"
  printf 'genes\tNCBI_gene_ID\tHGNC_gene_ID\nFIXG1\t900001\tHGNC:900001\n' \
    > "${h}/Gene-based/NCBIandHGNCgeneID/geneSymbol_NCBIandHGNCgeneID.tsv"
  # One record per SV type, far from both deletions. AnnotSV removes the
  # benign and pathogenic columns (and so the ranking) unless the files of
  # every type exist, and it checks them in the CHM13 folder. Each benign
  # record overlaps a gene (FARG1): AnnotSV 3.5.10 stops on a benign file
  # that overlaps none.
  for b in GRCh38 CHM13; do
    mkdir -p "${h}/SVincludedInFt/BenignSV/${b}" "${h}/FtIncludedInSV/PathogenicSV/${b}"
    for t in Loss Gain Ins Inv; do
      printf '20\t50001000\t50002000\tfixture_benign_%s\t20:50001000-50002000\t0.2000\n' "$t" \
        > "${h}/SVincludedInFt/BenignSV/${b}/benign_${t}_SV_${b}.sorted.bed"
      printf '20\t50100000\t50100500\tfixture phenotype\t\tfixture_pathogenic_%s\t20:50100000-50100500\n' "$t" \
        > "${h}/FtIncludedInSV/PathogenicSV/${b}/pathogenic_${t}_SV_${b}.sorted.bed"
    done
  done
}

# The gnomAD-SV-style sites: three records AnnotSV can use, six it cannot.
gnomad_sites() {
  local i='AN=100000;AC=5000;AF=0.05;N_HOMALT=200'
  printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  printf 'chr20\t1002900\tgnomAD-SV_fixture_DEL1\tN\t<DEL>\t999\tPASS\tALGORITHMS=manta;END=1004100;SVLEN=1200;SVTYPE=DEL;%s;AN_afr=100;AF_afr=0.9\n' "$i"
  printf 'chr20\t1999000\tgnomAD-SV_fixture_DEL2\tN\t<DEL>\t999\tPASS\tEND=2460000;SVTYPE=DEL;AN=90000;AC=1800;AF=0.02;N_HOMALT=40\n'
  printf 'chr20\t3000000\tgnomAD-SV_fixture_INS\tN\t<INS>\t999\tPASS\tEND=3000001;SVTYPE=INS;%s\n' "$i"
  printf 'chr20\t3100000\tfiltered_out\tN\t<DEL>\t999\tHIGH_NCR\tEND=3101000;SVTYPE=DEL;%s\n' "$i"
  printf 'chr20\t3200000\tbreakend\tN\t<BND>\t999\tPASS\tEND=3200000;SVTYPE=BND;%s\n' "$i"
  printf 'chr20\t3300000\tmulti_cnv\tN\t<CNV>\t999\tPASS\tEND=3310000;SVTYPE=CNV;CN_NUMBER=1000\n'
  printf 'chr20\t3400000\trare\tN\t<DEL>\t999\tPASS\tEND=3401000;SVTYPE=DEL;AN=100000;AC=50;AF=0.0005;N_HOMALT=5\n'
  printf 'chr20\t3500000\tfew_homalt\tN\t<DEL>\t999\tPASS\tEND=3501000;SVTYPE=DEL;AN=100000;AC=500;AF=0.005;N_HOMALT=2\n'
  printf 'chr20\t3600000\tsmall_an\tN\t<DEL>\t999\tPASS\tEND=3601000;SVTYPE=DEL;AN=800;AC=80;AF=0.1;N_HOMALT=9\n'
}

mane_summary() {
  printf '#NCBI_GeneID\tEnsembl_Gene\tHGNC_ID\tsymbol\tname\tRefSeq_nuc\tRefSeq_prot\tEnsembl_nuc\tEnsembl_prot\tMANE_status\tGRCh38_chr\tchr_start\tchr_end\tchr_strand\n'
  printf 'GeneID:1\tENSG9.1\tHGNC:1\tFIXG1\tfixture gene\tNM_900001.3\tNP_900001.3\tENST900001.2\tENSP900001.2\tMANE Select\tNC_000020.11\t1000001\t1010000\t+\n'
  printf 'GeneID:1\tENSG9.1\tHGNC:1\tFIXG1\tfixture gene\tNM_900003.1\tNP_900003.1\tENST900003.1\tENSP900003.1\tMANE Plus Clinical\tNC_000020.11\t1000001\t1010000\t+\n'
}

sv_vcf() {
  printf '##fileformat=VCFv4.2\n##contig=<ID=chr20,length=64444167>\n'
  printf '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">\n'
  printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n'
  printf 'chr20\t1003001\tDEL1\tN\t<DEL>\t100\tPASS\tEND=1004000;SVTYPE=DEL;SVLEN=-1000\tGT\t0/1\n'
  printf 'chr20\t2000001\tDEL2\tN\t<DEL>\t100\tPASS\tEND=2450000;SVTYPE=DEL;SVLEN=-450000\tGT\t0/1\n'
}

# field TSV END GENE COLUMN: COLUMN of the row of the SV ending at END; GENE
# empty picks the "full" row, a gene name its "split" row.
field() {
  awk -F'\t' -v e="$2" -v g="$3" -v c="$4" '
    NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
    $h["SV_end"] == e && ((g == "" && $h["Annotation_mode"] == "full") || (g != "" && $h["Annotation_mode"] == "split" && $h["Gene_name"] == g)) {
      print (c in h) ? $h[c] : "no column " c; exit
    }' "$1"
}

# annotsv NAME ANN_DIR [OPTION...]: the image on the fixture SVs, as step 05
# runs it, into ${WORK}/runs/NAME. AnnotSV also writes beside its output
# directory (NAME.SV_RE_intersect.tmp), so the parent is mounted.
annotsv() {
  local name=$1 ann=$2
  shift 2
  mkdir -p "${WORK}/runs/${name}"
  docker run --rm --network none --user "$(id -u):$(id -g)" -e HOME=/tmp \
    -v "${G}/s1/sv:/in:ro" -v "${ann}:/ann" -v "${WORK}/runs:/runs" "$ANNOTSV_IMAGE" \
    AnnotSV -SVinputFile /in/fixture_sv.vcf -outputFile "/runs/${name}/annotated.tsv" -outputDir "/runs/${name}" \
      -genomeBuild GRCh38 -annotationMode both -annotationsDir /ann "$@" > "${WORK}/runs/${name}.log" 2>&1 \
    || { fail "AnnotSV exited non-zero for ${name}"; tail -n 30 "${WORK}/runs/${name}.log"; }
}

# --- 1. setup.sh --annotsv-data ---------------------------------------------------------
build_annotations "$ANN"
cp -R "$ANN" "${WORK}/ann-plain"   # the set without the two files, for the red controls
mkdir -p "${G}/s1/sv"
sv_vcf > "${G}/s1/sv/fixture_sv.vcf"
# The two downloads, already in place: setup.sh keeps a download only once its
# checksum matches, so it takes a file it finds as whole.
gnomad_sites | gzip -c > "${ANN}/gnomad.v${GNOMAD_SV_VERSION}.sv.sites.vcf.gz"
mane_summary | gzip -c > "${ANN}/MANE.GRCh38.v${MANE_RELEASE}.summary.txt.gz"

rc=0
GENOME_DIR="$G" "${REPO}/scripts/setup.sh" --annotsv-data "$G" > "${WORK}/setup1.log" 2>&1 || rc=$?
expect "setup.sh --annotsv-data exits 0" "$rc" '^0$'
BEN="${H}/SVincludedInFt/BenignSV/GRCh38"
ADAPTER="${BEN}/gnomad.v4.0.sv.chrAll.from-v${GNOMAD_SV_VERSION}-sites.vcf.gz"
expect "setup.sh had AnnotSV merge and remove the gnomAD-SV file" "$([ -e "$ADAPTER" ] && echo present || echo gone)" '^gone$'
expect "gnomAD-SV records in AnnotSV's benign files" \
  "$(cat "${BEN}"/benign_*_SV_GRCh38.sorted.bed 2>/dev/null | cut -f4 | grep -v '^fixture_benign_' | sort | tr '\n' ' ')" \
  '^gnomAD-SV_fixture_DEL1 gnomAD-SV_fixture_DEL2 gnomAD-SV_fixture_INS $'
expect "DEL1 in the benign loss file with its frequency" \
  "$(awk -F'\t' '$4 == "gnomAD-SV_fixture_DEL1" {print $1 ":" $2 "-" $3 " " $6}' "${BEN}/benign_Loss_SV_GRCh38.sorted.bed" 2>/dev/null)" \
  '^(chr)?20:1002899-1004100 0\.0500$'
expect "the benign files belong to the caller" \
  "$(find "$BEN" -name 'benign_*_SV_GRCh38.sorted.bed' ! -user "$(id -u)" | wc -l | tr -d ' ')" '^0$'
expect "setup.sh leaves no merge run behind" "$(find "$ANN" -maxdepth 1 -name '.setup-merge*' | wc -l | tr -d ' ')" '^0$'
expect "the MANE Select list" "$(tr '\n' ' ' < "${ANN}/mane_select_refseq.txt" 2>/dev/null)" '^NM_900001\.3 $'
expect "the gnomAD-SV stamp" "$(cat "${H}/SVincludedInFt/BenignSV/GRCh38/gnomad-sv.installed" 2>/dev/null)" \
  "^gnomAD-SV ${GNOMAD_SV_VERSION} sha256 ${GNOMAD_SV_SHA256}\$"
expect "the downloads are removed once used" \
  "$(find "$ANN" -maxdepth 1 \( -name 'gnomad.v*.sv.sites.vcf.gz' -o -name 'MANE.*.summary.txt.gz' \) | wc -l | tr -d ' ')" '^0$'

# --- 2. step 05 -------------------------------------------------------------------------
rc=0
GENOME_DIR="$G" SV_VCF="${G}/s1/sv/fixture_sv.vcf" "${REPO}/scripts/05-annotsv.sh" s1 > "${WORK}/step05.log" 2>&1 || rc=$?
expect "step 05 exits 0" "$rc" '^0$'
OUT="${G}/s1/annotsv/s1_sv_annotated.tsv"
if [ -s "$OUT" ]; then
  expect "DEL1: FIXG1 is reported on its MANE Select transcript" "$(field "$OUT" 1004000 FIXG1 Tx)" '^NM_900001$'
  expect "DEL1: no frameshift on NM_900001 (intronic)" "$(field "$OUT" 1004000 FIXG1 Frameshift)" '^(NA|no)$'
  expect "DEL1: gnomAD-SV frequency" "$(field "$OUT" 1004000 '' B_loss_AFmax)" '^0\.05'
  expect "DEL1: class 1 (2F, common in gnomAD-SV)" "$(field "$OUT" 1004000 '' ACMG_class)" '^1'
  expect "DEL2: gnomAD-SV frequency" "$(field "$OUT" 2450000 '' B_loss_AFmax)" '^0\.02'
  expect "DEL2: out of class 4" "$(field "$OUT" 2450000 '' ACMG_class)" '^3'
else
  fail "step 05 wrote no ${OUT}"
  tail -n 40 "${WORK}/step05.log"
fi

# A second setup run finds both in place and adds nothing.
rc=0
GENOME_DIR="$G" "${REPO}/scripts/setup.sh" --annotsv-data "$G" > "${WORK}/setup2.log" 2>&1 || rc=$?
expect "second setup.sh --annotsv-data exits 0" "$rc" '^0$'
expect "second run: both already present" "$(grep -c 'already present' "${WORK}/setup2.log")" '^3$'
expect "second run: no gnomAD-SV file written again" "$([ -e "$ADAPTER" ] && echo present || echo absent)" '^absent$'

# --- 3. red controls --------------------------------------------------------------------
# Without -txFile, on the same data (gnomAD-SV merged): XM_900002 is back.
cp -R "$ANN" "${WORK}/ann-merged"
annotsv no-tx "${WORK}/ann-merged"
NT="${WORK}/runs/no-tx/annotated.tsv"
expect "control without -txFile: DEL1 on XM_900002" "$(field "$NT" 1004000 FIXG1 Tx 2>/dev/null)" '^XM_900002$'
expect "control without -txFile: a frameshift" "$(field "$NT" 1004000 FIXG1 Frameshift 2>/dev/null)" '^yes$'
expect "control without -txFile: DEL1 class 4" "$(field "$NT" 1004000 '' ACMG_class 2>/dev/null)" '^4'
# Without gnomAD-SV, with the MANE list: DEL2 stays class 4.
cp "${ANN}/mane_select_refseq.txt" "${WORK}/ann-plain/"
annotsv no-gnomad "${WORK}/ann-plain" -txFile /ann/mane_select_refseq.txt
NG="${WORK}/runs/no-gnomad/annotated.tsv"
expect "control without gnomAD-SV: DEL2 has no frequency" "$(field "$NG" 2450000 '' B_loss_AFmax 2>/dev/null)" '^$'
expect "control without gnomAD-SV: DEL2 class 4" "$(field "$NG" 2450000 '' ACMG_class 2>/dev/null)" '^4'

# --- 4. validate-setup.sh and step 05 without the files ---------------------------------
GENOME_DIR="$G" "${REPO}/scripts/validate-setup.sh" > "${WORK}/validate1.log" 2>&1 || true
expect "validate-setup: gnomAD-SV passes" "$(grep -c '\[OK\].*AnnotSV benign set (step 5): gnomAD-SV' "${WORK}/validate1.log")" '^1$'
expect "validate-setup: MANE passes" "$(grep -c '\[OK\].*AnnotSV transcripts (step 5): MANE' "${WORK}/validate1.log")" '^1$'
mv "${H}/SVincludedInFt/BenignSV/GRCh38/gnomad-sv.installed" "${WORK}/stamp.kept"
mv "${ANN}/mane_select_refseq.txt" "${WORK}/mane.kept"
GENOME_DIR="$G" "${REPO}/scripts/validate-setup.sh" > "${WORK}/validate2.log" 2>&1 || true
expect "validate-setup: gnomAD-SV removed fails" "$(grep -c '\[FAIL\].*AnnotSV benign set (step 5)' "${WORK}/validate2.log")" '^1$'
expect "validate-setup: MANE removed fails" "$(grep -c '\[FAIL\].*AnnotSV transcripts (step 5)' "${WORK}/validate2.log")" '^1$'
rc=0
GENOME_DIR="$G" SV_VCF="${G}/s1/sv/fixture_sv.vcf" "${REPO}/scripts/05-annotsv.sh" s1 > "${WORK}/step05-missing.log" 2>&1 || rc=$?
expect "step 05 without the MANE list stops" "$rc" '^1$'
expect "step 05 says how to add it" "$(grep -c 'setup.sh --annotsv-data' "${WORK}/step05-missing.log")" '^[1-9]'
mv "${WORK}/mane.kept" "${ANN}/mane_select_refseq.txt"
rc=0
GENOME_DIR="$G" SV_VCF="${G}/s1/sv/fixture_sv.vcf" "${REPO}/scripts/05-annotsv.sh" s1 > "${WORK}/step05-nostamp.log" 2>&1 || rc=$?
expect "step 05 without gnomAD-SV stops" "$rc" '^1$'
expect "step 05 names the gnomAD-SV stamp" "$(grep -c 'gnomad-sv.installed not found' "${WORK}/step05-nostamp.log")" '^1$'

if [ "$FAILS" -gt 0 ]; then
  echo "${FAILS} check(s) failed. Logs:"
  for f in setup1 step05 validate1 validate2; do echo "--- ${f}.log"; tail -n 25 "${WORK}/${f}.log" 2>/dev/null; done
  exit 1
fi
echo "All checks passed."
