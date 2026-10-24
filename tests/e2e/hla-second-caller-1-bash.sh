#!/usr/bin/env bash
# Step 08 with HLA_QC=true runs OptiType beside T1K, and step 36 passes HLA-A
# and HLA-B to PharmCAT only when the two callers agree at two fields.
#
# HG002ot: HG002's BAM (case 20), linked, so HG002's own results stay as the
# later cases read them.
#   - OptiType types the fixture's HLA reads; its HLA-A is HG002's published
#     pair (A*01:01 and A*26:01, as the T1K_IMAGE smoke row checks); the
#     concordance table compares the callers at HLA-A, HLA-B and HLA-C;
#   - where they agree step 36 passes the gene, naming both callers;
#   - a planted disagreement (OptiType's HLA-B1 rewritten) keeps HLA-B out
#     of the outside-call file and labels the row; HLA-A still passes;
#   - red control: the same inputs through a copy of bin/pgx_outside_calls.py
#     with the agreement test disabled put the disagreeing HLA-B in the
#     outside-call file, and the check above rejects it;
#   - HLA_QC=true at step 36 with no OptiType result holds both genes back;
#   - step 08 again without HLA_QC removes the old OptiType result.
#   - step 24 and the text report (render_report.py on the same summary)
#     print one line per gene of the concordance table, both calls named;
#     red control: a copy of bin/render_report.py with the OptiType block
#     disabled prints none, and the check rejects it.
# HG002sp, the positive control: HG002's reads over HLA-B replaced by reads
# simulated from two known alleles, B*57:01 and B*07:02 (genomic sequences of
# the pinned hla.dat), aligned with minimap2 as step 02 aligns. Both callers
# must call that pair, and step 36 passes it on.
. "$(dirname "$0")/lib.sh"

G=$GENOME_DIR
SRC=$SAMPLE
BAM_SRC="${G}/${SRC}/aligned/${SRC}_sorted.bam"
HLA_DAT=$(bash -c '. "$1/scripts/lib/common.sh" && install_data_file hla_dat >/dev/null && data_file hla_dat' _ "$REPO")
bash -c '. "$1/scripts/lib/common.sh" && install_data_file gencode_genes >/dev/null' _ "$REPO"
check "case 20's BAM is here" test -s "${BAM_SRC}.bai"
check "the pinned hla.dat is installed" test -s "$HLA_DAT"

# sample_from NAME BAM: a sample dir whose aligned BAM is BAM (linked when it can be).
sample_from() {
  rm -rf "${G:?}/$1"
  mkdir -p "${G}/$1/aligned"
  ln "$2" "${G}/$1/aligned/$1_sorted.bam" 2>/dev/null || cp "$2" "${G}/$1/aligned/$1_sorted.bam"
  ln "$2.bai" "${G}/$1/aligned/$1_sorted.bam.bai" 2>/dev/null || cp "$2.bai" "${G}/$1/aligned/$1_sorted.bam.bai"
}
# pair FILE GENE: the two-field pair of GENE in a T1K genotype (sorted, '/').
t1k_pair() {
  awk -F'\t' -v g="$2" '$1 == g {for (i = 3; i <= 6; i += 3) if ($i ~ /^HLA-/) {a = $i; sub(/^HLA-[A-Z0-9]+/, "", a)
    split(a, f, ":"); print f[1] ":" f[2]} if ($6 !~ /^HLA-/) {a = $3; sub(/^HLA-[A-Z0-9]+/, "", a); split(a, f, ":"); print f[1] ":" f[2]}}' \
    "$1" 2>/dev/null | sort | paste -sd/ -
}
# ot_pair FILE LETTER: OptiType's pair of one gene ('A', 'B'), as '*01:01/*26:01', sorted.
ot_pair() {
  awk -F'\t' -v l="$2" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i} NR == 2 {print $c[l "1"]; print $c[l "2"]}' "$1" 2>/dev/null \
    | sed "s/^$2//" | sort | paste -sd/ -
}
# held_back CALLS CONS GENE REASON: GENE is not in the outside-call file and its row is indeterminate for REASON.
held_back() {
  ! grep -q "^$3"$'\t' "$1" 2>/dev/null \
    && grep -Eq "^$3"$'\t''indeterminate'$'\t''no'$'\t'"$4" "$2" 2>/dev/null
}
not_held_back() { ! held_back "$@"; }
sorted_dip() { tr '/' '\n' <<< "$1" | sort | paste -sd/ -; }

# --- HG002ot: the fixture's own reads ---------------------------------------------
S=${SRC}ot
sample_from "$S" "$BAM_SRC"
OT="${G}/${S}/hla_optitype/${S}_optitype_result.tsv"
CONC="${G}/${S}/hla_optitype/${S}_hla_concordance.tsv"
GT="${G}/${S}/hla_t1k/${S}_hla_genotype.tsv"
CALLS="${G}/${S}/pgx_consensus/${S}_outside_calls.tsv"
CONS="${G}/${S}/pgx_consensus/${S}_pgx_consensus.tsv"

HLA_QC=true run_step 08-hla-typing.sh "$S"
check_step_exit 08-hla-typing.sh
check_ge "read pairs step 08 gave OptiType" "$(sed -n 's/^Read pairs at HLA-A, HLA-B and HLA-C: //p' "$STEP_LOG")" 100
cat "$OT" "$CONC" 2>/dev/null
check_eq "OptiType's HLA-A is HG002's published pair" "$(ot_pair "$OT" A)" '*01:01/*26:01'
check "OptiType called two HLA-B alleles" has '^\*[0-9]+:[0-9]+/\*[0-9]+:[0-9]+$' "$(ot_pair "$OT" B)"
check "the concordance table has HLA-A, HLA-B and HLA-C" \
  test "$(awk -F'\t' 'NR > 1 {print $1}' "$CONC" 2>/dev/null | paste -sd, -)" = HLA-A,HLA-B,HLA-C
for g in HLA-A HLA-B; do
  check "the concordance table: T1K and OptiType agree at ${g}" grep -q "^${g}"$'\t[^\t]*\t[^\t]*\tyes\t' "$CONC"
done
check_eq "the table's HLA-B is T1K's pair" "$(awk -F'\t' '$1 == "HLA-B" {print $2}' "$CONC" 2>/dev/null | tr '/' '\n' | sort | paste -sd/ -)" \
  "$(t1k_pair "$GT" HLA-B)"
cp "$CONC" "${E2E_WORK}/hla-qc-bash-concordance.tsv" 2>/dev/null   # case hla-second-caller-2 compares the Nextflow table
echo "- HLA second caller on the fixture: $(awk -F'\t' 'NR > 1 {printf "%s T1K %s OptiType %s agree %s; ", $1, $2, $3, $4}' "$CONC" 2>/dev/null)" >> "$E2E_NOTES"

# Both reports: one line per gene of the concordance table.
run_step 24-html-report.sh "$S"
check_step_exit 24-html-report.sh
HTML="${G}/${S}/${S}_report.html"
TXT="${CASE_TMP}/${S}_report.txt"
python3 -I "${REPO}/bin/render_report.py" --summary "${G}/${S}/summary.json" -o "$TXT" >/dev/null
cp "${G}/${S}/summary.json" "${E2E_WORK}/hla-qc-bash-summary.json" 2>/dev/null   # case hla-second-caller-2 compares the Nextflow one
sed -n '/OptiType check/,/^$/p' "$TXT"
# report_has TXT HTML: every row of the concordance table is a line of both reports.
report_has() {
  local g t o a note w x n=0 ok=0
  while IFS=$'\t' read -r g t o a note; do
    case $a in yes) w=agree ;; no) w=disagree ;; *) w=$a ;; esac
    x="${w}: T1K ${t}, OptiType ${o}"
    case $note in ''|-) ;; *) x+=" (${note})" ;; esac
    n=$((n + 1))
    grep -qF "$(printf '    %-8s ' "$g")${x}" "$1" \
      && grep -qF "${g}: T1K vs OptiType" "$2" && grep -qF "$x" "$2" && ok=$((ok + 1))
  done < <(tail -n +2 "$CONC")
  [ "$n" -eq 3 ] && [ "$ok" -eq 3 ]
}
not_report_has() { ! report_has "$@"; }
check "the text report heads the OptiType check" grep -qF 'OptiType check (HLA_QC=true, two fields):' "$TXT"
check "both reports print each gene's agreement and both calls, as the table says" report_has "$TXT" "$HTML"
# Red control: the report without its OptiType block.
BR="${CASE_TMP}/bin_no_optitype"
rm -rf "$BR"
cp -R "${REPO}/bin" "$BR"
sed -i 's/ot = S\["hla"\]\["data"\].get("optitype")/ot = None/; s/if hla.get("optitype") is not None:/if False:/' "${BR}/render_report.py"
check "the control's mutations applied" test "$(grep -c '^        ot = None$\|^    if False:$' "${BR}/render_report.py")" -eq 2
python3 -I "${BR}/render_report.py" --summary "${G}/${S}/summary.json" -o "${CASE_TMP}/ctl.txt" -o "${CASE_TMP}/ctl.html" >/dev/null
check "control: without the OptiType block the check above rejects both reports" \
  not_report_has "${CASE_TMP}/ctl.txt" "${CASE_TMP}/ctl.html"

run_step 36-pgx-consensus.sh "$S"
check_step_exit 36-pgx-consensus.sh
check "step 36 read OptiType's result" has '^  OptiType HLA types \(step 08, HLA_QC=true\): ' "$(cat "$STEP_LOG")"
cat "$CALLS" "$CONS" 2>/dev/null
for g in HLA-A HLA-B; do
  check "${g} reaches PharmCAT, both callers named" \
    grep -Eq "^${g}"$'\t''\*[0-9]+:[0-9]+/\*[0-9]+:[0-9]+'$'\t''yes'$'\t''T1K and OptiType \(step 08\), which agree at two fields'$'\t' "$CONS"
  check "${g} is in the outside-call file" grep -q "^${g}"$'\t' "$CALLS"
done

# The planted disagreement: OptiType's B1 becomes an allele T1K did not call.
cp "$OT" "${CASE_TMP}/optitype_real.tsv"
T1K_B=$(t1k_pair "$GT" HLA-B)
PLANT='B*58:01'; [[ "/${T1K_B}/" != *"/*58:01/"* ]] || PLANT='B*15:01'
awk -F'\t' -v OFS='\t' -v p="$PLANT" 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "B1") c = i} NR == 2 {$c = p} {print}' \
  "${CASE_TMP}/optitype_real.tsv" > "$OT"
check "the planted OptiType HLA-B (${PLANT}) differs from T1K's (${T1K_B})" test "$(ot_pair "$OT" B)" != "$T1K_B"
run_step 36-pgx-consensus.sh "$S"
check_step_exit 36-pgx-consensus.sh
cat "$CALLS" "$CONS" 2>/dev/null
check "the disagreeing HLA-B is held back and labelled" held_back "$CALLS" "$CONS" HLA-B 'T1K and OptiType disagree at two fields'
check "its row prints both calls" grep -Eq "^HLA-B"$'\t'".*T1K: .*; OptiType: .*\\${PLANT#B}" "$CONS"
check "HLA-A, where they agree, still reaches PharmCAT" grep -q $'^HLA-A\t' "$CALLS"

# Red control: the agreement test disabled, the same inputs.
BROKEN="${CASE_TMP}/pgx_outside_calls_broken.py"
sed 's/if not same_pair(names, ot):/if False:/' "${REPO}/bin/pgx_outside_calls.py" > "$BROKEN"
check "the control's mutation applied" test "$(grep -c 'if False:' "$BROKEN")" -eq 1
python3 -I "$BROKEN" --calls "${CASE_TMP}/ctl_calls.tsv" --consensus "${CASE_TMP}/ctl_cons.tsv" --hla "$GT" --optitype "$OT" >/dev/null
check "control: without the agreement test the disagreeing HLA-B reaches the outside-call file" grep -q $'^HLA-B\t' "${CASE_TMP}/ctl_calls.tsv"
check "control: the check above rejects that output" \
  not_held_back "${CASE_TMP}/ctl_calls.tsv" "${CASE_TMP}/ctl_cons.tsv" HLA-B 'T1K and OptiType disagree at two fields'

# Asked for, but no OptiType result: neither gene is passed on.
rm -f "$OT"
HLA_QC=true run_step 36-pgx-consensus.sh "$S"
check_step_exit 36-pgx-consensus.sh
for g in HLA-A HLA-B; do
  check "${g} held back: HLA_QC asked for, no OptiType result" held_back "$CALLS" "$CONS" "$g" 'the second caller was asked for \(HLA_QC\)'
done

# Step 08 without HLA_QC leaves no OptiType result behind.
cp "${CASE_TMP}/optitype_real.tsv" "$OT"
run_step 08-hla-typing.sh "$S"
check_step_exit 08-hla-typing.sh
check "step 08 without HLA_QC removed the old OptiType result" test ! -e "$OT"
check "and the old concordance table" test ! -e "$CONC"

# --- HG002sp: the positive control (spike-in of B*57:01 and B*07:02) ------------
S=${SRC}sp
SP="${G}/${S}/spike"
rm -rf "${G:?}/${S}"
mkdir -p "$SP"
cat > "${CASE_TMP}/spike.py" <<'PY'
"""spike.py HLA_DAT OUTDIR ALLELE...: read pairs (150 bp, 350 bp fragments,
one every 5 bp, the strands alternating: about 50x) from the longest genomic
sequence in hla.dat of each allele named."""
import sys
dat, out, wanted = sys.argv[1], sys.argv[2], sys.argv[3:]
best, name, seq, in_sq = {}, None, [], False
with open(dat) as f:
    for line in f:
        if line.startswith("DE   "):
            name = line[5:].split(",")[0].strip()
        elif line.startswith("SQ   "):
            in_sq, seq = True, []
        elif line.startswith("//"):
            s = "".join(seq).upper()
            for w in wanted:
                if name and (name == w or name.startswith(w + ":")) and len(s) > len(best.get(w, ("", ""))[1]):
                    best[w] = (name, s)
            in_sq, name = False, None
        elif in_sq:
            seq.append("".join(ch for ch in line if ch.isalpha()))
comp = str.maketrans("ACGTN", "TGCAN")
rc = lambda s: s.translate(comp)[::-1]
n = 0
with open(f"{out}/R1.fq", "w") as r1, open(f"{out}/R2.fq", "w") as r2:
    for w in wanted:
        if w not in best or len(best[w][1]) < 2500:
            sys.exit(f"no genomic sequence of {w} in {dat}")
        nm, s = best[w]
        print(f"{w}: {nm}, {len(s)} bp")
        for start in range(0, len(s) - 350 + 1, 5):
            frag = s[start:start + 350]
            a, b = frag[:150], rc(frag[-150:])
            if (start // 5) % 2:
                a, b = b, a
            n += 1
            r1.write(f"@spike{n}/1\n{a}\n+\n{'I' * 150}\n")
            r2.write(f"@spike{n}/2\n{b}\n+\n{'I' * 150}\n")
print(f"read pairs: {n}")
PY
python3 -I "${CASE_TMP}/spike.py" "$HLA_DAT" "$SP" 'HLA-B*57:01:01' 'HLA-B*07:02:01' | tee "${CASE_TMP}/spike.log"
check_ge "simulated read pairs" "$(sed -n 's/^read pairs: //p' "${CASE_TMP}/spike.log")" 500
MMI="${G}/reference/GRCh38_no_alt_analysis_set.sr.mmi"
[ -s "$MMI" ] || docker run --rm -u "$(id -u):$(id -g)" -v "${G}:/genome" -w /genome "$MINIMAP2_IMAGE" \
  minimap2 -x sr -d reference/GRCh38_no_alt_analysis_set.sr.mmi reference/GRCh38_no_alt_analysis_set.fasta
printf 'chr6\t31352872\t31358188\n' > "${SP}/hla_b.bed"
docker run --rm -u "$(id -u):$(id -g)" -v "${G}:/genome" -w /genome "$MINIMAP2_SAMTOOLS_IMAGE" sh -c '
  set -e
  minimap2 -t 4 -a -x sr reference/GRCh38_no_alt_analysis_set.sr.mmi "$1/spike/R1.fq" "$1/spike/R2.fq" \
    | samtools sort -o "$1/spike/spike.bam" -
  samtools view -b -L "$1/spike/hla_b.bed" -U "$1/spike/rest.bam" -o /dev/null "$2"
  samtools merge -f -o "$1/spike/merged.bam" "$1/spike/rest.bam" "$1/spike/spike.bam"
  samtools index "$1/spike/merged.bam"
  samtools index "$1/spike/spike.bam"
  echo "spike reads mapped near HLA-B: $(samtools view -c -F 0x904 "$1/spike/spike.bam" chr6:31352872-31358188)"' sh \
  "$S" "${SRC}/aligned/${SRC}_sorted.bam" 2>&1 | tee "${CASE_TMP}/align.log"
check_ge "spike reads mapped over HLA-B" "$(sed -n 's/^spike reads mapped near HLA-B: //p' "${CASE_TMP}/align.log")" 500
mkdir -p "${G}/${S}/aligned"
mv "${SP}/merged.bam" "${G}/${S}/aligned/${S}_sorted.bam"
mv "${SP}/merged.bam.bai" "${G}/${S}/aligned/${S}_sorted.bam.bai"
rm -f "${SP}"/*.bam "${SP}"/*.fq

HLA_QC=true run_step 08-hla-typing.sh "$S"
check_step_exit 08-hla-typing.sh
OT="${G}/${S}/hla_optitype/${S}_optitype_result.tsv"
GT="${G}/${S}/hla_t1k/${S}_hla_genotype.tsv"
cat "$OT" "${G}/${S}/hla_optitype/${S}_hla_concordance.tsv" 2>/dev/null
WANT=$(sorted_dip '*57:01/*07:02')
check_eq "T1K calls the spiked-in HLA-B pair" "$(t1k_pair "$GT" HLA-B)" "$WANT"
check_eq "OptiType calls the spiked-in HLA-B pair" "$(ot_pair "$OT" B)" "$WANT"
run_step 36-pgx-consensus.sh "$S"
check_step_exit 36-pgx-consensus.sh
CALLS="${G}/${S}/pgx_consensus/${S}_outside_calls.tsv"
cat "$CALLS" 2>/dev/null
check_eq "step 36 passes the spiked-in HLA-B pair to PharmCAT" \
  "$(sorted_dip "$(awk -F'\t' '$1 == "HLA-B" {print $2}' "$CALLS" 2>/dev/null)")" "$WANT"
rm -rf "${G:?}/${S}" "${G:?}/${SRC}ot"

finish
