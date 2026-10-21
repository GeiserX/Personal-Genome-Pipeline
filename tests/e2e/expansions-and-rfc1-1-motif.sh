#!/usr/bin/env bash
# Step 09b's RFC1 motif count (bin/rfc1_motif.py) with the real Stranger and
# PYPGX_IMAGE. The fixture has no reads at RFC1, so the case makes three
# samples of its own: an ExpansionHunter VCF with one RFC1 record (the header
# of tests/smoke/stranger_input.vcf.in) and a BAM of made-up reads at the
# locus, the reads of tests/test_rfc1_motif.py.
#   path    RFC1 at 100/150 in-repeat; AAGGG reads (soft-clipped, unmapped
#           and placed at the mate, mates on another contig, reverse strand)
#           and three reference reads: Stranger flags full_mutation, the
#           count says pathogenic_motif, and the text report prints it.
#   benign  the same record; AAAAG reads: benign_motif.
#   normal  RFC1 at 5/5 on spanning reads: Stranger does not flag it, so
#           no motif table is written (the control that the count runs only
#           on a flag).
# What this cannot show: reads from a real carrier. Real AAGGG and AAAAG
# expansions, and how many of their reads the aligner soft-clips or places
# with the mate, need a real genome.
. "$(dirname "$0")/lib.sh"

make_sample() {  # make_sample <name> <case of tests/test_rfc1_motif.py> <REPCN> <SO> <REPCI>
  local s=$1 kind=$2 cn=$3 so=$4 ci=$5
  rm -rf "${GENOME_DIR:?}/${s}"
  mkdir -p "${GENOME_DIR}/${s}/expansion_hunter" "${GENOME_DIR}/${s}/aligned"
  python3 - "$REPO" "$kind" "$cn" "$so" "$ci" "${GENOME_DIR}/${s}/expansion_hunter/${s}_eh.vcf" \
    "${GENOME_DIR}/${s}/aligned/reads.sam" <<'PY'
import os, sys
repo, kind, cn, so, ci, vcf, sam = sys.argv[1:8]
sys.path.insert(0, os.path.join(repo, "tests"))
import test_rfc1_motif as t
head = [l for l in open(os.path.join(repo, "tests", "smoke", "stranger_input.vcf.in")) if l.startswith("#")]
alts = sorted({int(x) for x in cn.split("/")})
head[-1:-1] = [f'##ALT=<ID=STR{a},Description="Allele comprised of {a} repeat units">\n' for a in alts]
gt = "1/2" if len(alts) == 2 else "1/1"
with open(vcf, "w") as fh:
    fh.writelines(head)
    fh.write(f"{t.CHROM}\t{t.START + 1}\tRFC1\tA\t" + ",".join(f"<STR{a}>" for a in alts)
             + f"\t.\tPASS\tSVTYPE=STR;END={t.END};REF=11;RL=55;RU=AARRG;VARID=RFC1;REPID=RFC1"
             f"\tGT:SO:REPCN:REPCI:ADSP:ADFL:ADIR:LC\t{gt}:{so}:{cn}:{ci}:0/0:2/2:20/20:30.0\n")
at, elsewhere = t.case(kind)
with open(sam, "w") as fh:
    fh.write("@HD\tVN:1.6\tSO:unsorted\n")
    fh.write(f"@SQ\tSN:{t.CHROM}\tLN:190214555\n@SQ\tSN:{t.OTHER}\tLN:248956422\n")
    for r in at + list(elsewhere.values()):
        flag = (1 | (64 if r.is_read1 else 128) | (4 if r.is_unmapped else 0) | (16 if r.is_reverse else 0)
                | (256 if r.is_secondary else 0) | (512 if r.is_qcfail else 0) | (1024 if r.is_duplicate else 0))
        mate = "=" if r.mate_contig == r.contig else r.mate_contig
        fh.write("\t".join(map(str, [r.query_name, flag, r.contig, r.reference_start + 1,
                                     0 if r.is_unmapped else 60, r.cigar_text or "*", mate,
                                     r.next_reference_start + 1, 0, r.query_sequence,
                                     "I" * len(r.query_sequence)])) + "\n")
PY
  docker run --rm -u "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -w "/genome/${s}/aligned" "$SAMTOOLS_IMAGE" \
    sh -c "samtools sort -o ${s}_sorted.bam reads.sam && samtools index ${s}_sorted.bam"
  rm -f "${GENOME_DIR}/${s}/aligned/reads.sam"
}

motif_col() {  # motif_col <sample> <column>
  awk -F'\t' -v k="$2" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i} NR == 2 {print $c[k]}' \
    "${GENOME_DIR}/$1/expansion_hunter/$1_rfc1_motif.tsv" 2>/dev/null
}

stranger_status() {  # stranger_status <sample>
  awk -F'\t' '!/^#/ {n = split($8, kv, ";"); for (i = 1; i <= n; i++) if (kv[i] ~ /^STR_STATUS=/) print substr(kv[i], 12)}' \
    "${GENOME_DIR}/$1/expansion_hunter/$1_eh_stranger.vcf" 2>/dev/null
}

# --- path: AAGGG reads -------------------------------------------------------------
P="${SAMPLE}rfc1path"
make_sample "$P" pathogenic 100/150 INREPEAT/INREPEAT 80-120/130-170
check "the BAM of made-up reads is indexed" test -s "${GENOME_DIR}/${P}/aligned/${P}_sorted.bam.bai"
run_step 09b-stranger.sh "$P"
check_step_exit 09b-stranger.sh
check_eq "path: Stranger flags RFC1 at 100/150 (the catalog's low pathogenic minimum)" "$(stranger_status "$P")" full_mutation
check "path: the log says the motif count ran" has 'RFC1 is flagged: counting the repeat motif' "$(cat "$STEP_LOG")"
check_eq "path: verdict" "$(motif_col "$P" verdict)" pathogenic_motif
check_eq "path: dominant motif" "$(motif_col "$P" dominant_motif)" AAGGG
check_eq "path: 16 repeat reads (6 soft-clipped, 2 placed unmapped, 3 rescued mates, 2 reverse, 3 reference)" \
  "$(motif_col "$P" repeat_reads)" 16
check_eq "path: the shares of both motifs" "$(motif_col "$P" motif_shares)" "AAGGG=0.81,AAAAG=0.19"
check "path: no temporary table is left" test ! -e "${GENOME_DIR}/${P}/expansion_hunter/${P}_rfc1_motif.tsv.part"
run_step generate-report.sh "$P"
check_step_exit generate-report.sh
TXT="${GENOME_DIR}/${P}/${P}_report.txt"
check "path: the text report prints the verdict under the RFC1 flag" \
  grep -q 'RFC1 motif count: pathogenic_motif: AAGGG in 13 of 16 repeat reads' "$TXT"
check "path: the text report prints the in-repeat count as a lower bound" grep -q 'RFC1 .*>=80/>=130' "$TXT"
run_step 09b-stranger.sh "$P"
check "path: a second run keeps the table" has 'RFC1 motif count already exists' "$(cat "$STEP_LOG")"

# --- benign: AAAAG reads -------------------------------------------------------------
B="${SAMPLE}rfc1benign"
make_sample "$B" benign 100/150 INREPEAT/INREPEAT 80-120/130-170
run_step 09b-stranger.sh "$B"
check_step_exit 09b-stranger.sh
check_eq "benign: verdict" "$(motif_col "$B" verdict)" benign_motif
check_eq "benign: 10 AAAAG reads" "$(motif_col "$B" motif_shares)/$(motif_col "$B" repeat_reads)" "AAAAG=1.00/10"

# --- normal: no flag, no count -------------------------------------------------------
N="${SAMPLE}rfc1normal"
make_sample "$N" benign 5/5 SPANNING/SPANNING 5-5/5-5
run_step 09b-stranger.sh "$N"
check_step_exit 09b-stranger.sh
check_eq "normal: Stranger does not flag RFC1 at 5/5" "$(stranger_status "$N")" normal
check "normal: no motif table" test ! -e "${GENOME_DIR}/${N}/expansion_hunter/${N}_rfc1_motif.tsv"
check "normal: the log does not mention the motif count" lacks 'RFC1 is flagged' "$(cat "$STEP_LOG")"

rm -rf "${GENOME_DIR:?}/${P}" "${GENOME_DIR:?}/${B}" "${GENOME_DIR:?}/${N}"
finish
