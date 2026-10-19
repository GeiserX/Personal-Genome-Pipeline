#!/usr/bin/env bash
# setup.sh builds both GENCODE files, the gene lines for T1K (gencode_genes)
# and the exon BED of the cDNA read check (gencode_exons), from one download
# of GENCODE's basic annotation, checked against GENCODE's MD5SUMS. The fake
# download serves a small real GTF and its MD5SUMS line.
#   1. A first setup fetches the GTF once and writes both files: the gene
#      lines, and each exon once as chrom, 0-based start, end, sorted.
#   2. validate-setup.sh reports both present.
#   3. With the exon BED deleted, validate-setup.sh names it missing and what
#      steps 06 and 17 do without it; the next setup fetches once more and
#      writes only the exon BED, leaving the gene file as it was.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
use_output_hook
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hook-outputs"

REL=$(sed -n 's/^GENCODE_RELEASE=//p' "${REPO_ROOT}/scripts/lib/common.sh")
GTF="gencode.v${REL}.basic.annotation.gtf.gz"
GENES="${GENOME_DIR}/reference/gencode.v${REL}.basic.genes.gtf"
EXONS="${GENOME_DIR}/reference/gencode.v${REL}.basic.exons.bed"
SERVED="${CASE_WORK}/served"
mkdir -p "$SERVED"
# Two genes on two contigs; the first gene's exon 1 is listed by two
# transcripts, so it must come out once.
{
  printf 'chr2\tHAVANA\tgene\t1001\t5000\t.\t+\t.\tgene_name "GENEB";\n'
  printf 'chr2\tHAVANA\texon\t4001\t5000\t.\t+\t.\tgene_name "GENEB";\n'
  printf 'chr2\tHAVANA\texon\t1001\t1200\t.\t+\t.\tgene_name "GENEB";\n'
  printf 'chr1\tHAVANA\tgene\t101\t900\t.\t+\t.\tgene_name "GENEA";\n'
  printf 'chr1\tHAVANA\ttranscript\t101\t900\t.\t+\t.\tgene_name "GENEA";\n'
  printf 'chr1\tHAVANA\texon\t101\t300\t.\t+\t.\tgene_name "GENEA";\n'
  printf 'chr1\tHAVANA\texon\t701\t900\t.\t+\t.\tgene_name "GENEA";\n'
  printf 'chr1\tHAVANA\texon\t101\t300\t.\t+\t.\tgene_name "GENEA";\n'
} | gzip -nc > "${SERVED}/${GTF}"
printf '%s  %s\n' "$(fake_md5 < "${SERVED}/${GTF}")" "$GTF" > "${SERVED}/MD5SUMS"
fetches() { grep -cE "^(curl|wget) [^ ]*/${GTF//./\\.} " "$FAKE_DOCKER_LOG" || true; }

# --- 1. one fetch, both files ---------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
FAKE_DOWNLOAD_DIR="$SERVED" run_expect 0 setup "${SCRIPTS}/setup.sh" "$GENOME_DIR"
[ "$(fetches)" -eq 1 ] || fail "setup.sh downloaded ${GTF} $(fetches) times, expected once"
[ -s "$GENES" ] || fail "setup.sh wrote no ${GENES}"
[ -s "$EXONS" ] || fail "setup.sh wrote no ${EXONS}"
[ "$(cut -f3 "$GENES" | sort -u | tr '\n' ' ')" = "gene " ] || fail "the gene file holds other lines: $(cut -f3 "$GENES" | sort -u | tr '\n' ' ')"
[ "$(wc -l < "$GENES" | tr -d ' ')" -eq 2 ] || fail "the gene file has $(wc -l < "$GENES") lines, expected 2"
want=$(printf 'chr1\t100\t300\nchr1\t700\t900\nchr2\t1000\t1200\nchr2\t4000\t5000')
[ "$(cat "$EXONS")" = "$want" ] || fail "the exon BED is not each exon once, 0-based and sorted: $(tr '\n\t' '| ' < "$EXONS")"
[ ! -e "${GENES}.download.gz" ] || fail "setup.sh kept the GTF download"
output_has setup "\[OK\] GENCODE ${REL} gene coordinates \(step 08\): "
output_has setup "\[OK\] GENCODE ${REL} exons \(steps 06 and 17\) already present\."

# --- 2. validate-setup.sh reports both -------------------------------------------------
run_rc validate "${SCRIPTS}/validate-setup.sh" sample1
output_has validate "GENCODE ${REL} gene coordinates: present"
output_has validate "GENCODE ${REL} exons: present"

# --- 3. the exon BED missing: named by validate, rebuilt alone by setup -------------------
rm -f "$EXONS"
touch -t 202001010000 "$GENES"
run_rc validate-missing "${SCRIPTS}/validate-setup.sh" sample1
output_has validate-missing "GENCODE ${REL} exons \(steps 06 and 17 skip their cDNA read check without them\) not found at: ${EXONS//./\\.}"
: > "$FAKE_DOCKER_LOG"
FAKE_DOWNLOAD_DIR="$SERVED" run_expect 0 setup-again "${SCRIPTS}/setup.sh" "$GENOME_DIR"
[ "$(fetches)" -eq 1 ] || fail "the second setup.sh downloaded ${GTF} $(fetches) times, expected once"
[ "$(cat "$EXONS" 2>/dev/null)" = "$want" ] || fail "the second setup.sh did not rebuild the exon BED"
[ -z "$(find "$GENES" -newer "${SERVED}/MD5SUMS")" ] || fail "the second setup.sh rewrote the gene file it already had"
output_has setup-again "\[OK\] GENCODE ${REL} gene coordinates \(step 08\) already present\."
output_has setup-again "\[OK\] GENCODE ${REL} exons \(steps 06 and 17\): "
