#!/usr/bin/env bash
# Step 05 (AnnotSV) on what it now reads: Manta's calls after step 04, with
# the net length change of sequence-resolved records, plus the Delly and
# CNVpytor calls no Manta call covers. The real AnnotSV image (ANNOTSV_IMAGE)
# on a made-up annotation set small enough for the unit-test job.
#
# The annotation set, on chr20:
#   FIXG1   a haploinsufficient gene (HI 3), NM_900001 (its MANE Select):
#           exons 1000000-1001000, 1005000-1006000 and 1009000-1010000,
#           coding 1000500-1009500 (2,000 coding bases);
#   G01-G40 forty one-exon genes from 2,010,000;
#   FARG1   far from all of them, under the benign records of the set.
# Manta's calls (diploidSV.vcf.gz, as step 04 writes it):
#   PLANT   chr20:1005400, REF 137 bp, ALT 2 bp, END=1005536, SVLEN=-136:
#           inside FIXG1's second exon, a net in-frame loss of 135 bp
#           (6.75% of the coding sequence);
#   DEL2    chr20:2000001-2450000 over the forty genes (3C, class 4);
#   BND     a lone inversion junction, its two breakends 5 Mb apart
#           (chr20:20,000,000 and 25,000,000), as step 04 leaves one of
#           1 Mb or more.
# Delly: a PASS deletion matching DEL2 (Manta has it: not added).
# CNVpytor: a deletion chr20:30,000,001-30,200,000 with e-val 1e-8, q0 0,
#   that no Manta call covers.
#
# Checks:
#   1. step 05: PLANT reads SV_length -135, Frameshift "no" on NM_900001, no
#      2E-1, class 3; the CNVpytor deletion is in AnnotSV's table with
#      SVCALLER=cnvpytor and a negative SV_length (CNVpytor's own SVLEN is
#      positive); Delly's DEL2 is not there twice; the breakend pair
#      is annotated as breakends, no INV row; the class 4-5 table holds DEL2
#      and nothing else, found by column name;
#   2. red controls with the same image and data, AnnotSV on Manta's calls
#      as step 05 gave them before (no length fix, no other callers, AnnotSV's
#      default -bracketedAltMode): PLANT is a frameshift with 2E-1 in class
#      4, the CNVpytor deletion is absent, the breakend pair is a 5 Mb INV;
#      and the old printed filter ($120, $16) does not find DEL2 on this
#      bundle, whose columns are not a full bundle's.
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

# --- The annotation set ------------------------------------------------------------------
mkdir -p "${H}/Genes/GRCh38" "${H}/Gene-based/Fixture" "${H}/Gene-based/NCBIandHGNCgeneID"
# chrom txStart txEnd strand gene tx cdsStart cdsEnd exonStarts exonEnds
{
  printf '20\t1000000\t1010000\t+\tFIXG1\tNM_900001\t1000500\t1009500\t1000000,1005000,1009000,\t1001000,1006000,1010000,\n'
  for i in $(seq 1 40); do
    s=$((2000000 + i * 10000))
    printf '20\t%d\t%d\t+\tG%02d\tNM_91%04d\t%d\t%d\t%d,\t%d,\n' "$s" $((s + 2000)) "$i" "$i" $((s + 100)) $((s + 1900)) "$s" $((s + 2000))
  done
  printf '20\t50000000\t50020000\t+\tFARG1\tNM_920001\t50000100\t50019900\t50000000,\t50020000,\n'
} > "${H}/Genes/GRCh38/genes.RefSeq.sorted.bed"
sed 's/\tNM_9/\tENST9/' "${H}/Genes/GRCh38/genes.RefSeq.sorted.bed" > "${H}/Genes/GRCh38/genes.ENSEMBL.sorted.bed"
cut -f6 "${H}/Genes/GRCh38/genes.RefSeq.sorted.bed" | awk '{print $1 "\t1"}' > "${H}/Genes/GRCh38/transcript_version.RefSeq.tsv"
cut -f6 "${H}/Genes/GRCh38/genes.ENSEMBL.sorted.bed" | awk '{print $1 "\t1"}' > "${H}/Genes/GRCh38/transcript_version.ENSEMBL.tsv"
printf 'genes\tHI\tTS\tGnomAD_pLI\tLOEUF_bin\tDDD_HI_percent\tOMIM_morbid\nFIXG1\t3\t0\t0.10\t5\t50\tno\n' \
  > "${H}/Gene-based/Fixture/fixture_genes.tsv"
printf 'genes\tNCBI_gene_ID\tHGNC_gene_ID\nFIXG1\t900001\tHGNC:900001\n' \
  > "${H}/Gene-based/NCBIandHGNCgeneID/geneSymbol_NCBIandHGNCgeneID.tsv"
# One benign and one pathogenic record per SV type, far from the calls:
# AnnotSV drops the ranking unless the files of every type exist, in the
# GRCh38 and the CHM13 folders, and stops on a benign file that overlaps no gene.
for b in GRCh38 CHM13; do
  mkdir -p "${H}/SVincludedInFt/BenignSV/${b}" "${H}/FtIncludedInSV/PathogenicSV/${b}"
  for t in Loss Gain Ins Inv; do
    printf '20\t50001000\t50002000\tfixture_benign_%s\t20:50001000-50002000\t0.2000\n' "$t" \
      > "${H}/SVincludedInFt/BenignSV/${b}/benign_${t}_SV_${b}.sorted.bed"
    printf '20\t50100000\t50100500\tfixture phenotype\t\tfixture_pathogenic_%s\t20:50100000-50100500\n' "$t" \
      > "${H}/FtIncludedInSV/PathogenicSV/${b}/pathogenic_${t}_SV_${b}.sorted.bed"
  done
done
# The two files setup.sh adds, which step 05 requires.
printf 'NM_900001\n' > "${ANN}/mane_select_refseq.txt"
printf 'gnomAD-SV fixture\n' > "${H}/SVincludedInFt/BenignSV/GRCh38/gnomad-sv.installed"

# --- The calls ---------------------------------------------------------------------------
S="${G}/s1"
mkdir -p "${S}/manta/results/variants" "${S}/delly" "${S}/cnvpytor"
REF137=$(printf 'A%.0s' $(seq 1 137))
{
  printf '##fileformat=VCFv4.1\n##contig=<ID=chr20,length=64444167>\n'
  printf '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
  printf '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length">\n##INFO=<ID=MATEID,Number=.,Type=String,Description="Mate">\n'
  printf '##INFO=<ID=INV_JUNCTIONS,Number=1,Type=Integer,Description="Junctions">\n'
  printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n'
  printf 'chr20\t1005400\tPLANT\t%s\tAC\t500\tPASS\tEND=1005536;SVTYPE=DEL;SVLEN=-136\tGT\t0/1\n' "$REF137"
  printf 'chr20\t2000001\tDEL2\tN\t<DEL>\t500\tPASS\tEND=2450000;SVTYPE=DEL;SVLEN=-450000\tGT\t0/1\n'
  printf 'chr20\t20000000\tMantaBND:7:0:1:0:0:0:0\tN\tN]chr20:25000000]\t500\tPASS\tSVTYPE=BND;MATEID=MantaBND:7:0:1:0:0:0:1;INV_JUNCTIONS=1\tGT\t0/1\n'
  printf 'chr20\t25000000\tMantaBND:7:0:1:0:0:0:1\tN\tN]chr20:20000000]\t500\tPASS\tSVTYPE=BND;MATEID=MantaBND:7:0:1:0:0:0:0;INV_JUNCTIONS=1\tGT\t0/1\n'
} | gzip -c > "${S}/manta/results/variants/diploidSV.vcf.gz"
{
  printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n'
  printf 'chr20\t2000500\tDEL00000001\tN\t<DEL>\t60\tPASS\tIMPRECISE;SVTYPE=DEL;END=2449000\tGT\t0/1\n'
} | gzip -c > "${S}/delly/s1_sv.vcf.gz"
{
  printf '##fileformat=VCFv4.2\n##source=CNVpytor\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n'
  printf 'chr20\t30000001\tCNVpytor_del0\tN\t<DEL>\t.\tPASS\tEND=30200000;IMPRECISE;SVLEN=200000;SVTYPE=DEL;pytorRD=0.5;pytorP1=1e-08;pytorQ0=0.0\tGT:CN\t0/1:1\n'
} | gzip -c > "${S}/cnvpytor/s1_cnvs.vcf.gz"

# field TSV ID COLUMN [GENE]: COLUMN of the "full" row of the record named ID,
# or of its "split" row on GENE.
field() {
  awk -F'\t' -v id="$2" -v c="$3" -v g="${4:-}" '
    NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
    $h["ID"] == id && ((g == "" && $h["Annotation_mode"] == "full") || (g != "" && $h["Annotation_mode"] == "split" && $h["Gene_name"] == g)) {
      print (c in h) ? $h[c] : "no column " c; exit
    }' "$1"
}
# rows TSV ID: the number of "full" rows of the record named ID.
rows() { awk -F'\t' -v id="$2" 'NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
  $h["ID"] == id && $h["Annotation_mode"] == "full" { n++ } END { print n + 0 }' "$1"; }

# --- 1. step 05 ------------------------------------------------------------------------------
rc=0
GENOME_DIR="$G" "${REPO}/scripts/05-annotsv.sh" s1 > "${WORK}/step05.log" 2>&1 || rc=$?
expect "step 05 exits 0" "$rc" '^0$'
OUT="${S}/annotsv/s1_sv_annotated.tsv"
if [ -s "$OUT" ]; then
  expect "PLANT: SV_length is the net change" "$(field "$OUT" PLANT SV_length)" '^-135$'
  expect "PLANT: reported on NM_900001" "$(field "$OUT" PLANT Tx FIXG1)" '^NM_900001$'
  expect "PLANT: no frameshift (135 bp keep the frame)" "$(field "$OUT" PLANT Frameshift FIXG1)" '^no$'
  expect "PLANT: no 2E-1 in the ranking" "$(field "$OUT" PLANT AnnotSV_ranking_criteria | grep -c '2E-1' || true)" '^0$'
  expect "PLANT: class 3" "$(field "$OUT" PLANT ACMG_class)" '^3$'
  expect "CNVpytor's deletion is in AnnotSV's table" "$(rows "$OUT" CNVpytor_del0)" '^1$'
  expect "CNVpytor's deletion names its caller" "$(field "$OUT" CNVpytor_del0 INFO | grep -o 'SVCALLER=[a-z_]*')" '^SVCALLER=cnvpytor$'
  expect "CNVpytor's deletion: a negative SV_length, as Manta's deletions" "$(field "$OUT" CNVpytor_del0 SV_length)" '^-[0-9]+$'
  expect "Manta's DEL2 names its caller" "$(field "$OUT" DEL2 INFO | grep -o 'SVCALLER=[a-z_]*')" '^SVCALLER=manta$'
  expect "Delly's copy of DEL2 is not added" "$(rows "$OUT" DEL00000001)" '^0$'
  expect "the breakend pair is annotated as breakends" "$(field "$OUT" MantaBND:7:0:1:0:0:0:0 SV_type)" '^BND$'
  expect "no INV row" "$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next } $h["SV_type"] == "INV" { n++ } END { print n + 0 }' "$OUT")" '^0$'
  C45="${S}/annotsv/s1_sv_class4_5.tsv"
  expect "the class 4-5 table: DEL2 only" "$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next } { print $h["ID"] ":" $h["ACMG_class"] }' "$C45" 2>/dev/null | tr '\n' ' ')" '^DEL2:4 $'
  if [ "$(head -n 1 "$C45" 2>/dev/null)" = "$(head -n 1 "$OUT")" ]; then pass "the class 4-5 table has AnnotSV's header"
  else fail "the class 4-5 table does not start with AnnotSV's header"; fi
  expect "the log names what was added" "$(grep -o 'cnvpytor: 1 DEL/DUP calls with e-val < 0.01 and q0 < 0.5, 1 added' "${WORK}/step05.log")" 'added$'
else
  fail "step 05 wrote no ${OUT}"
  tail -n 40 "${WORK}/step05.log"
fi

# --- 2. red controls ---------------------------------------------------------------------------
mkdir -p "${WORK}/runs"
docker run --rm --network none --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "${S}/manta/results/variants:/in:ro" -v "${ANN}:/ann" -v "${WORK}/runs:/runs" "$ANNOTSV_IMAGE" \
  AnnotSV -SVinputFile /in/diploidSV.vcf.gz -outputFile /runs/before/annotated.tsv -outputDir /runs/before \
    -genomeBuild GRCh38 -annotationMode both -annotationsDir /ann -txFile /ann/mane_select_refseq.txt \
  > "${WORK}/runs/before.log" 2>&1 || { fail "AnnotSV exited non-zero on Manta's calls"; tail -n 30 "${WORK}/runs/before.log"; }
B="${WORK}/runs/before/annotated.tsv"
if [ -s "$B" ]; then
  expect "control: PLANT read on the REF span" "$(field "$B" PLANT SV_length)" '^-136$'
  expect "control: PLANT a frameshift" "$(field "$B" PLANT Frameshift FIXG1)" '^yes$'
  expect "control: PLANT gets 2E-1, class 4" "$(field "$B" PLANT AnnotSV_ranking_criteria | grep -o '2E-1' | head -n 1) $(field "$B" PLANT ACMG_class)" '^2E-1 4$'
  expect "control: no CNVpytor deletion" "$(rows "$B" CNVpytor_del0)" '^0$'
  expect "control: the breakend pair becomes a 5 Mb INV" \
    "$(field "$B" MantaBND:7:0:1:0:0:0:0 SV_type) $(field "$B" MantaBND:7:0:1:0:0:0:0 SV_length)" '^INV 5000000$'
else
  fail "AnnotSV wrote no ${B}"
fi
if [ -s "${OUT}" ]; then
  n=$(head -n 1 "$OUT" | tr '\t' '\n' | grep -nx ACMG_class | cut -d: -f1)
  expect "control: ACMG_class is another column than 120 on this bundle" "$([ "$n" != 120 ] && echo "$n")" '^[0-9]+$'
  expect "control: the old filter (\$120, \$16) does not find DEL2" \
    "$(awk -F'\t' 'NR==1 || ($120==4 || $120==5) && $16=="full"' "$OUT" | tail -n +2 | grep -c 'DEL2' || true)" '^0$'
fi

echo
if [ "$FAILS" -eq 0 ]; then echo "RESULT: ALL PASS"; else echo "RESULT: ${FAILS} FAILED"; exit 1; fi
