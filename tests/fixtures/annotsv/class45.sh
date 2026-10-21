# shellcheck shell=bash
# Made-up AnnotSV data and SV calls on chr20 for the class 4 and 5 checks
# (bin/annotsv_postfilter.py), read by tests/test_annotsv_postfilter.sh and
# tests/e2e/annotsv-report-section-1-chain.sh. The records describe no person.
#
#   class45_bundle ANN   AnnotSV's annotation data in ANN, small enough for
#                        the real AnnotSV image in a unit test:
#     FIXG1   a haploinsufficient gene (HI 3), NM_900001 (its MANE Select):
#             exons 1000000-1001000, 1005000-1006000, 1009000-1010000
#             (0-based starts), coding 1000500-1009500;
#     G01-G40 forty one-exon genes from 2,010,000, and H01-H40 from 3,010,000;
#     FARG1   far from all of them, under the fixture benign records;
#     a common benign deletion, AF 0.05, over all of FIXG1.
#   class45_calls DIR SAMPLE   Manta's and Delly's calls of SAMPLE in DIR
#                        (the step 04 and step 19 file names), bgzipped with
#                        BCFTOOLS_IMAGE:
#     EXON    Manta, chr20:1001000-1009000: from the last base of exon 1 to
#             the base before exon 3, as the reads of a processed
#             transcript aligned across both introns look (exon 2 is lost,
#             1,000 coding bases: a frameshift in an HI gene, class 4);
#     NEAR    Manta, chr20:1000990-1009000: the same, 10 bp inside intron 1;
#     DEL2    Manta and Delly (2 bp apart), chr20:2000001-2450000 over G01-G40
#             (criterion 3C, class 4);
#     DEL3    Manta alone, chr20:3000001-3450000 over H01-H40 (class 4).

class45_bundle() {
  local ann=$1 h="$1/Annotations_Human" b t i s
  mkdir -p "${h}/Genes/GRCh38" "${h}/Gene-based/Fixture" "${h}/Gene-based/NCBIandHGNCgeneID"
  # chrom txStart txEnd strand gene tx cdsStart cdsEnd exonStarts exonEnds
  {
    printf '20\t1000000\t1010000\t+\tFIXG1\tNM_900001\t1000500\t1009500\t1000000,1005000,1009000,\t1001000,1006000,1010000,\n'
    for i in $(seq 1 40); do
      s=$((2000000 + i * 10000))
      printf '20\t%d\t%d\t+\tG%02d\tNM_91%04d\t%d\t%d\t%d,\t%d,\n' "$s" $((s + 2000)) "$i" "$i" $((s + 100)) $((s + 1900)) "$s" $((s + 2000))
    done
    for i in $(seq 1 40); do
      s=$((3000000 + i * 10000))
      printf '20\t%d\t%d\t+\tH%02d\tNM_93%04d\t%d\t%d\t%d,\t%d,\n' "$s" $((s + 2000)) "$i" "$i" $((s + 100)) $((s + 1900)) "$s" $((s + 2000))
    done
    printf '20\t50000000\t50020000\t+\tFARG1\tNM_920001\t50000100\t50019900\t50000000,\t50020000,\n'
  } > "${h}/Genes/GRCh38/genes.RefSeq.sorted.bed"
  sed 's/\tNM_9/\tENST9/' "${h}/Genes/GRCh38/genes.RefSeq.sorted.bed" > "${h}/Genes/GRCh38/genes.ENSEMBL.sorted.bed"
  cut -f6 "${h}/Genes/GRCh38/genes.RefSeq.sorted.bed" | awk '{print $1 "\t1"}' > "${h}/Genes/GRCh38/transcript_version.RefSeq.tsv"
  cut -f6 "${h}/Genes/GRCh38/genes.ENSEMBL.sorted.bed" | awk '{print $1 "\t1"}' > "${h}/Genes/GRCh38/transcript_version.ENSEMBL.tsv"
  printf 'genes\tHI\tTS\tGnomAD_pLI\tLOEUF_bin\tDDD_HI_percent\tOMIM_morbid\nFIXG1\t3\t0\t0.10\t5\t50\tno\n' \
    > "${h}/Gene-based/Fixture/fixture_genes.tsv"
  printf 'genes\tNCBI_gene_ID\tHGNC_gene_ID\nFIXG1\t900001\tHGNC:900001\n' \
    > "${h}/Gene-based/NCBIandHGNCgeneID/geneSymbol_NCBIandHGNCgeneID.tsv"
  # One benign and one pathogenic record per SV type, far from the calls:
  # AnnotSV drops the ranking unless the files of every type exist, in the
  # GRCh38 and the CHM13 folders, and stops on a benign file that overlaps no
  # gene. The GRCh38 loss file also holds the common deletion over FIXG1.
  for b in GRCh38 CHM13; do
    mkdir -p "${h}/SVincludedInFt/BenignSV/${b}" "${h}/FtIncludedInSV/PathogenicSV/${b}"
    for t in Loss Gain Ins Inv; do
      {
        if [ "$b$t" = GRCh38Loss ]; then
          printf '20\t999000\t1011000\tfixture_common_loss\t20:999000-1011000\t0.0500\n'
        fi
        printf '20\t50001000\t50002000\tfixture_benign_%s\t20:50001000-50002000\t0.2000\n' "$t"
      } > "${h}/SVincludedInFt/BenignSV/${b}/benign_${t}_SV_${b}.sorted.bed"
      printf '20\t50100000\t50100500\tfixture phenotype\t\tfixture_pathogenic_%s\t20:50100000-50100500\n' "$t" \
        > "${h}/FtIncludedInSV/PathogenicSV/${b}/pathogenic_${t}_SV_${b}.sorted.bed"
    done
  done
  # The two files setup.sh adds, which step 05 requires.
  printf 'NM_900001\n' > "${ann}/mane_select_refseq.txt"
  printf 'gnomAD-SV fixture\n' > "${h}/SVincludedInFt/BenignSV/GRCh38/gnomad-sv.installed"
}

# _class45_vcf SAMPLE RECORDS: a VCF with the header both callers' files need.
_class45_vcf() {
  printf '##fileformat=VCFv4.2\n##contig=<ID=chr20,length=64444167>\n##FILTER=<ID=PASS,Description="All filters passed">\n'
  printf '##ALT=<ID=DEL,Description="Deletion">\n'
  printf '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
  printf '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n'
  printf '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">\n'
  printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$1"
  printf '%s' "$2"
}

class45_calls() {
  local dir=$1 s=$2 r
  mkdir -p "${dir}/manta/results/variants" "${dir}/delly"
  r=$'chr20\t1000990\tNEAR\tN\t<DEL>\t500\tPASS\tSVTYPE=DEL;END=1009000;SVLEN=-8010\tGT\t0/1\n'
  r+=$'chr20\t1001000\tEXON\tN\t<DEL>\t500\tPASS\tSVTYPE=DEL;END=1009000;SVLEN=-8000\tGT\t0/1\n'
  r+=$'chr20\t2000001\tDEL2\tN\t<DEL>\t500\tPASS\tSVTYPE=DEL;END=2450000;SVLEN=-449999\tGT\t0/1\n'
  r+=$'chr20\t3000001\tDEL3\tN\t<DEL>\t500\tPASS\tSVTYPE=DEL;END=3450000;SVLEN=-449999\tGT\t0/1\n'
  _class45_vcf "$s" "$r" > "${dir}/manta/results/variants/diploidSV.vcf"
  r=$'chr20\t2000003\tDEL00000001\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=2450002;SVLEN=-449999\tGT\t0/1\n'
  _class45_vcf "$s" "$r" > "${dir}/delly/${s}_sv.vcf"
  # shellcheck disable=SC2016,SC2154  # $1 expands inside the container; BCFTOOLS_IMAGE comes from versions.env
  docker run --rm --network none --user "$(id -u):$(id -g)" -v "${dir}:/d" "$BCFTOOLS_IMAGE" bash -euo pipefail -c '
    for f in /d/manta/results/variants/diploidSV.vcf "/d/delly/$1_sv.vcf"; do
      bcftools view -Oz -o "$f.gz" "$f" && bcftools index -t "$f.gz" && rm "$f"
    done' _ "$s"
}
