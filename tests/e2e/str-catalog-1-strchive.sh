#!/usr/bin/env bash
# The STRchive catalog of steps 09 and 9b with the real tools.
#   1. setup.sh --strchive-data downloads the two STRchive files of
#      versions.env, checks their sha256 and writes catalog.json (82 loci).
#   2. ExpansionHunter 5.0.0 (step 09) runs on the catalog's loci on the
#      fixture's contigs: one record per repeat of each locus. Should it
#      stop, every locus is run alone to name the ones it refuses.
#   3. Stranger (step 9b) grades with the catalog step 09 ran on: records in
#      = records out, the STRchive source on every graded locus, and no
#      status on the loci the catalog leaves ungraded.
#   4. Made-up reads, in a sample of their own: RFC1 AAGGG reads and SAMD12
#      (FAME1) TTTCA reads, with an ExpansionHunter VCF that sizes both
#      above the catalog's normal maximum. Step 9b flags both with the real
#      Stranger, the RFC1 count reads the catalog's RFC1 (pathogenic_motif,
#      the shares of tests/test_rfc1_motif.py) and the SAMD12 count
#      pathogenic_motif; the same with TTTTA reads at SAMD12 gives
#      benign_motif.
#   5. STRANGER and STRANGER_RFC1 in Nextflow, given the same catalog, write
#      step 9b's two motif tables byte for byte.
#   Control: step 9b on the same made-up sample with STRANGER_REPEATS=bundled
#      writes no STRchive source and no table of the other motif loci.
# What this cannot show: reads of a real carrier at any of the 82 loci; the
# fixture's slice of HG002 has none.
. "$(dirname "$0")/lib.sh"

SETUP="${CASE_TMP}/setup"
rm -rf "$SETUP"
mkdir -p "$SETUP"
CAT="${SETUP}/reference/strchive-${STRCHIVE_VERSION}/catalog.json"

# --- 1. setup ---------------------------------------------------------------------------
"${REPO}/scripts/setup.sh" --strchive-data "$SETUP" > "${CASE_TMP}/setup.log" 2>&1
check_eq "setup.sh --strchive-data exits 0" "$?" 0
cat "${CASE_TMP}/setup.log"
check "setup.sh names the catalog it wrote" has "\[OK\] STRchive ${STRCHIVE_VERSION//./\\.} catalog" "$(cat "${CASE_TMP}/setup.log")"
check_eq "the catalog has 82 loci" "$(python3 -c 'import json, sys; print(len(json.load(open(sys.argv[1]))))' "$CAT" 2>/dev/null)" 82

# --- 2. ExpansionHunter on the loci of the fixture's contigs ------------------------------------
S="${SAMPLE}strc"
SD="${GENOME_DIR}/${S}"
rm -rf "$SD"
mkdir -p "${SD}/aligned"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${SD}/aligned/${S}_sorted.bam"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" "${SD}/aligned/${S}_sorted.bam.bai"
SUB="${GENOME_DIR}/reference/strchive_fixture.json"
python3 - "$CAT" "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta.fai" "$SUB" > "${CASE_TMP}/subset.txt" <<'PY'
import json, sys
loci = json.load(open(sys.argv[1]))
contigs = {line.split("\t")[0] for line in open(sys.argv[2])}
regions = lambda l: l["ReferenceRegion"] if isinstance(l["ReferenceRegion"], list) else [l["ReferenceRegion"]]
kept = [l for l in loci if all(r.split(":")[0] in contigs for r in regions(l))]
json.dump(kept, open(sys.argv[3], "w"), indent=1)
print(len(kept), sum(len(regions(l)) for l in kept), sum("NormalMax" in l for l in kept),
      ",".join(l["LocusId"] for l in kept if "NoThresholds" in l))
PY
read -r N_LOCI N_VAR N_GRADED UNGRADED < "${CASE_TMP}/subset.txt"
echo "loci on the fixture's contigs: ${N_LOCI} (${N_VAR} repeats, ${N_GRADED} graded; not graded: ${UNGRADED})"
check_ge "catalog loci on the fixture's contigs" "$N_LOCI" 40
EH_CATALOG="$SUB" run_step 09-expansion-hunter.sh "$S" male
check_step_exit 09-expansion-hunter.sh
if [ "$STEP_RC" != 0 ]; then
  echo "ExpansionHunter stopped: each locus alone"
  python3 -c 'import json, sys; [print(l["LocusId"]) for l in json.load(open(sys.argv[1]))]' "$SUB" | while read -r id; do
    python3 -c 'import json, sys; json.dump([l for l in json.load(open(sys.argv[1])) if l["LocusId"] == sys.argv[2]], open(sys.argv[3], "w"))' \
      "$SUB" "$id" "${GENOME_DIR}/reference/one_locus.json"
    docker run --rm -u "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" "$EXPANSIONHUNTER_IMAGE" ExpansionHunter \
      --reads "/genome/${S}/aligned/${S}_sorted.bam" --reference /genome/reference/GRCh38_no_alt_analysis_set.fasta \
      --variant-catalog /genome/reference/one_locus.json --output-prefix "/genome/${S}/one" --sex male \
      > "${CASE_TMP}/one.log" 2>&1 </dev/null || echo "  ${id} stops ExpansionHunter: $(tail -n 2 "${CASE_TMP}/one.log" | tr '\n' ' ')"
  done
fi
check "the log names the catalog" has 'Variant catalog: .*strchive_fixture\.json \(EH_CATALOG\)' "$(cat "$STEP_LOG")"
EH="${SD}/expansion_hunter/${S}_eh.vcf"
check_eq "one ExpansionHunter record per repeat of every locus" "$(grep -vc '^#' "$EH" 2>/dev/null || echo 0)" "$N_VAR"
check "step 09 copies the catalog beside the VCF" cmp -s "$SUB" "${SD}/expansion_hunter/${S}_eh_catalog.json"

# --- 3. Stranger on the catalog step 09 ran on --------------------------------------------------
run_step 09b-stranger.sh "$S"
check_step_exit 09b-stranger.sh
check "the log names the catalog of step 09" has "Repeat catalog: the one step 09 ran on, .*\(STRchive-${STRCHIVE_VERSION//./\\.}\)" "$(cat "$STEP_LOG")"
OUT="${SD}/expansion_hunter/${S}_eh_stranger.vcf"
check_eq "records in = records out" "$(grep -vc '^#' "$OUT" 2>/dev/null || echo 0)" "$N_VAR"
check_eq "every graded locus carries STRchive as its source" \
  "$(grep -v '^#' "$OUT" 2>/dev/null | grep -c "SourceDisplay=STRchive-${STRCHIVE_VERSION}" || true)" "$N_GRADED"
check_eq "and a status" "$(grep -v '^#' "$OUT" 2>/dev/null | grep -c 'STR_STATUS=' || true)" "$N_GRADED"
for id in ${UNGRADED//,/ }; do
  check "${id} (not graded in the catalog) has no status" \
    lacks "REPID=${id}(;|\$).*STR_STATUS|STR_STATUS.*REPID=${id}(;|\$)" "$(grep -v '^#' "$OUT")"
done

# --- 4. made-up reads at RFC1 and SAMD12 ------------------------------------------------------
make_sample() {  # make_sample <name> <SAMD12 unit>
  local s=$1 unit=$2
  rm -rf "${GENOME_DIR:?}/${s}"
  mkdir -p "${GENOME_DIR}/${s}/expansion_hunter" "${GENOME_DIR}/${s}/aligned"
  cp "$CAT" "${GENOME_DIR}/${s}/expansion_hunter/${s}_eh_catalog.json"
  python3 - "$REPO" "$unit" "${GENOME_DIR}/${s}/expansion_hunter/${s}_eh.vcf" "${GENOME_DIR}/${s}/aligned/reads.sam" <<'PY'
import os, random, sys
repo, unit, vcf, sam = sys.argv[1:5]
sys.path.insert(0, os.path.join(repo, "tests"))
import test_rfc1_motif as t
head = [l for l in open(os.path.join(repo, "tests", "smoke", "stranger_input.vcf.in")) if l.startswith("##")
        and not l.startswith("##contig")]
head += [f'##ALT=<ID=STR{a},Description="Allele comprised of {a} repeat units">\n' for a in (40, 100, 150)]
head += ["##contig=<ID=chr1,length=248956422>\n", "##contig=<ID=chr4,length=190214555>\n",
         "##contig=<ID=chr8,length=145138636>\n"]
fmt = "GT:SO:REPCN:REPCI:ADSP:ADFL:ADIR:LC"
expanded = "1/2:INREPEAT/INREPEAT:100/150:80-120/130-170:0/0:2/2:20/20:30.0"  # two long alleles
with open(vcf, "w") as fh:
    fh.writelines(head)
    fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
    # The catalog's RFC1 record is its AAGGG repeat, and SAMD12's its TGAAA one.
    fh.write(f"chr4\t39348480\tRFC1\tA\t<STR100>,<STR150>\t.\tPASS\tSVTYPE=STR;END=39348483;REF=0;RL=4;RU=AAGGG;"
             f"VARID=RFC1;REPID=RFC1\t{fmt}\t{expanded}\n")
    fh.write(f"chr8\t118366913\tSAMD12\tA\t<STR40>\t.\tPASS\tSVTYPE=STR;END=118366918;REF=1;RL=6;RU=TGAAA;"
             f"VARID=SAMD12;REPID=SAMD12\t{fmt}\t0/1:SPANNING/INREPEAT:1/40:1-1/30-50:10/0:2/2:0/10:30.0\n")
at, elsewhere = t.case("pathogenic")
rng = random.Random(3)
start = 118366812
fame = [(f"f{i}", "chr8", start, "150M", unit * 30) for i in range(5)]
left = "".join(rng.choice("ACGT") for _ in range(60))
ref_allele = left + "TAAAA" * 20 + "TGAAA" + left[:30]
fame += [(f"fr{i}", "chr8", start - 60, f"{len(ref_allele)}M", ref_allele) for i in range(6)]
with open(sam, "w") as fh:
    fh.write("@HD\tVN:1.6\tSO:unsorted\n")
    fh.write(f"@SQ\tSN:{t.CHROM}\tLN:190214555\n@SQ\tSN:{t.OTHER}\tLN:248956422\n@SQ\tSN:chr8\tLN:145138636\n")
    for r in at + list(elsewhere.values()):
        flag = (1 | (64 if r.is_read1 else 128) | (4 if r.is_unmapped else 0) | (16 if r.is_reverse else 0)
                | (256 if r.is_secondary else 0) | (512 if r.is_qcfail else 0) | (1024 if r.is_duplicate else 0))
        mate = "=" if r.mate_contig == r.contig else r.mate_contig
        fh.write("\t".join(map(str, [r.query_name, flag, r.contig, r.reference_start + 1,
                                     0 if r.is_unmapped else 60, r.cigar_text or "*", mate,
                                     r.next_reference_start + 1, 0, r.query_sequence,
                                     "I" * len(r.query_sequence)])) + "\n")
    for name, chrom, pos, cig, seq in fame:
        fh.write("\t".join(map(str, [name, 0, chrom, pos + 1, 60, cig, "*", 0, 0, seq, "I" * len(seq)])) + "\n")
PY
  docker run --rm -u "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -w "/genome/${s}/aligned" "$SAMTOOLS_IMAGE" \
    sh -c "samtools sort -o ${s}_sorted.bam reads.sam && samtools index ${s}_sorted.bam"
  rm -f "${GENOME_DIR}/${s}/aligned/reads.sam"
}
col() {  # col <table> <locus> <column>
  awk -F'\t' -v l="$2" -v k="$3" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $1 == l {print $c[k]}' "$1" 2>/dev/null
}
status() {  # status <sample> <REPID>
  awk -F'\t' -v id="$2" '!/^#/ && $8 ~ ("(^|;)REPID=" id "(;|$)") {n = split($8, kv, ";"); for (i = 1; i <= n; i++) if (kv[i] ~ /^STR_STATUS=/) print substr(kv[i], 12)}' \
    "${GENOME_DIR}/$1/expansion_hunter/$1_eh_stranger.vcf" 2>/dev/null
}

P="${SAMPLE}strcpath"
make_sample "$P" TTTCA
run_step 09b-stranger.sh "$P"
check_step_exit 09b-stranger.sh
check_eq "RFC1 at 100 and 150 AAGGG units: Stranger flags it with the catalog (normal to 11, pathogenic from 400)" \
  "$(status "$P" RFC1)" pre_mutation
check_eq "SAMD12 at 40 TGAAA units: flagged (normal to the reference's 1, pathogenic from 105)" "$(status "$P" SAMD12)" pre_mutation
R_TAB="${GENOME_DIR}/${P}/expansion_hunter/${P}_rfc1_motif.tsv"
O_TAB="${GENOME_DIR}/${P}/expansion_hunter/${P}_str_motifs.tsv"
check_eq "RFC1 through the catalog: pathogenic_motif" "$(col "$R_TAB" RFC1 verdict)" pathogenic_motif
check_eq "RFC1 through the catalog: the shares of the built-in count" "$(col "$R_TAB" RFC1 motif_shares)" "AAGGG=0.81,AAAAG=0.19"
check_eq "RFC1 through the catalog: the locus is both its repeats" "$(col "$R_TAB" RFC1 start)-$(col "$R_TAB" RFC1 end)" "39348425-39348483"
check_eq "SAMD12: TTTCA reads are pathogenic_motif" "$(col "$O_TAB" SAMD12 verdict)" pathogenic_motif
check_eq "SAMD12: five repeat reads, the six reference reads not counted" "$(col "$O_TAB" SAMD12 repeat_reads)" 5
check_eq "SAMD12: the motif is named as STRchive writes it (TGAAA, the TTTCA of the other strand)" \
  "$(col "$O_TAB" SAMD12 dominant_motif)" TGAAA
check_eq "the other table holds the flagged motif loci but RFC1" "$(awk 'NR > 1 {print $1}' "$O_TAB" 2>/dev/null | paste -sd, -)" SAMD12

B="${SAMPLE}strcbenign"
make_sample "$B" TTTTA
run_step 09b-stranger.sh "$B"
check_step_exit 09b-stranger.sh
check_eq "SAMD12: TTTTA reads are benign_motif" "$(col "${GENOME_DIR}/${B}/expansion_hunter/${B}_str_motifs.tsv" SAMD12 verdict)" benign_motif

# --- 5. the same in Nextflow -------------------------------------------------------------------
if command -v nextflow >/dev/null; then
  NP="${CASE_TMP}/project"
  rm -rf "$NP" "${CASE_TMP}/nf"
  mkdir -p "$NP"
  cp -R "${REPO}/bin" "${NP}/bin"
  cat > "${NP}/main.nf" <<NF
include { STRANGER; STRANGER_RFC1 } from '${REPO}/modules/local/stranger/main'

workflow {
    cat = channel.value(file('${CAT}'))
    STRANGER(channel.of([[id: '${P}'], file('${GENOME_DIR}/${P}/expansion_hunter/${P}_eh.vcf')]), cat)
    STRANGER_RFC1(STRANGER.out.vcf.map { meta, vcf ->
        [meta, vcf, file('${GENOME_DIR}/${P}/aligned/${P}_sorted.bam'), file('${GENOME_DIR}/${P}/aligned/${P}_sorted.bam.bai')] }, cat)
}
NF
  {
    echo "params.outdir = '${CASE_TMP}/nf/results'"
    echo "params.publish_dir_mode = 'copy'"
    echo "docker.enabled = true"
    echo "docker.runOptions = '--user $(id -u):$(id -g) --network none'"
  } > "${NP}/nextflow.config"
  ( cd "$CASE_TMP" && nextflow run "${NP}/main.nf" -ansi-log false -c "${REPO}/conf/containers.config" \
      -work-dir "${CASE_TMP}/nf/work" ) > "${CASE_TMP}/nf.log" 2>&1
  rc=$?
  check_eq "nextflow run of STRANGER and STRANGER_RFC1 exits 0" "$rc" 0
  [ "$rc" = 0 ] || grep -E 'ERROR|Error executing process|Command error' -A 8 "${CASE_TMP}/nf.log" | head -40
  NF_DIR="${CASE_TMP}/nf/results/${P}/expansion_hunter"
  check "STRANGER grades as step 9b, byte for byte" \
    cmp -s <(grep -v '^#' "${GENOME_DIR}/${P}/expansion_hunter/${P}_eh_stranger.vcf") <(grep -v '^#' "${NF_DIR}/${P}_eh_stranger.vcf")
  check "STRANGER_RFC1 writes step 9b's RFC1 table, byte for byte" cmp -s "$R_TAB" "${NF_DIR}/${P}_rfc1_motif.tsv"
  check "and its table of the other motif loci" cmp -s "$O_TAB" "${NF_DIR}/${P}_str_motifs.tsv"
else
  fail "nextflow is not on PATH"
fi

# --- control: the same sample with Stranger's bundled catalog -----------------------------------
# The checks above can fail: without the STRchive catalog no locus carries
# its source, SAMD12 gets the bundled thresholds and no motif table is written.
rm -f "${GENOME_DIR}/${P}/expansion_hunter/${P}_eh_stranger.vcf"
STRANGER_REPEATS=bundled run_step 09b-stranger.sh "$P"
check_step_exit 09b-stranger.sh
check "control: STRANGER_REPEATS=bundled names Stranger's GRCh38 catalog" has 'variant_catalog_grch38' "$(cat "$STEP_LOG")"
check_eq "control: no record carries STRchive as its source" \
  "$(grep -v '^#' "${GENOME_DIR}/${P}/expansion_hunter/${P}_eh_stranger.vcf" 2>/dev/null | grep -c 'SourceDisplay=STRchive' || true)" 0
check "control: no table of the other motif loci" test ! -e "$O_TAB"

rm -rf "${GENOME_DIR:?}/${S}" "${GENOME_DIR:?}/${P}" "${GENOME_DIR:?}/${B}" "$SUB" "${GENOME_DIR}/reference/one_locus.json"
finish
