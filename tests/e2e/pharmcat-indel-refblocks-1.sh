#!/usr/bin/env bash
# PharmCAT positions longer than one base (indels such as UGT1A1 *28, and
# two-base substitutions) are called reference when the gVCF's reference
# blocks cover them: step 07 runs bin/pgx_refblock_fill.py after expanding the
# blocks. The sample is the invented gVCF tests/fixtures/pharmcat/refblocks.g.vcf
# on the fixture reference: blocks over five such positions, and four that
# must stay missing (a ./. hole inside the span, a block with GQ 10, a variant
# record inside the span, no block at all). It lives in a GENOME_DIR of its
# own (the fixture reference hard-linked into it), so no later case sees it.
#
# Control, in the same run: the step 07 recipe of before the fix (blocks
# expanded, nothing filled) run on the same gVCF with PharmCAT's preprocessor
# and matcher. Its preprocessor must list the five positions as missing; its
# matcher must have no call at the four whose PharmCAT REF is longer than one
# base (it discards the one-base record there), and UGT1A1 *28 must be
# uncallable. At the TPMT insertion (one-base REF) the matcher reads the
# one-base record as C/C even so. Step 07 must list none of the five as
# missing and call all five reference in its match.json, while the four
# planted to stay missing do.
. "$(dirname "$0")/lib.sh"

P=REFBLOCKS
G="${E2E_WORK}/genome-refblocks"
D="${G}/${P}/vcf"
FX="${REPO}/tests/fixtures/pharmcat/refblocks.g.vcf"
REF="reference/GRCh38_no_alt_analysis_set.fasta"
# chrom pos REF gene, and what the control's matcher has there: filled from the blocks
FILL="chr1 97450065 TG DPYD missing
chr2 233760233 CAT UGT1A1 missing
chr6 18149032 C TPMT call:C/C
chr19 38572266 TC RYR1 missing
chrX 154532055 CTCT G6PD missing"
# chrom pos why: must stay missing
LEFT="chr10 94942205 a_./._hole_inside_the_span
chr10 94949161 a_block_with_GQ_10
chr10 94949281 a_variant_record_inside_the_span
chr6 18149103 no_block_over_the_span"

check "the fixture gVCF exists" test -s "$FX"
rm -rf "$G"
mkdir -p "${G}/reference" "$D"
for f in "$REF" "${REF}.fai"; do
  ln -f "${GENOME_DIR}/${f}" "${G}/${f}" 2>/dev/null || cp "${GENOME_DIR}/${f}" "${G}/${f}"
done
check "the fixture reference is where step 07 looks" test -s "${G}/${REF}"

# PharmCAT's own files, read by the step from the pinned image.
IMG_FILES=$(docker run --rm --network none "$PHARMCAT_IMAGE" ls -l /pharmcat/pharmcat_positions.vcf.bgz /pharmcat/pharmcat_regions.bed 2>&1)
echo "$IMG_FILES"
check "the image has /pharmcat/pharmcat_positions.vcf.bgz" has 'pharmcat_positions\.vcf\.bgz$' "$IMG_FILES"
check "the image has /pharmcat/pharmcat_regions.bed" has 'pharmcat_regions\.bed$' "$IMG_FILES"
MULTI=$(docker run --rm --network none "$PHARMCAT_IMAGE" gzip -dc /pharmcat/pharmcat_positions.vcf.bgz \
  | awk -F'\t' '!/^#/ {n = split($5, a, ","); m = length($4) > 1; for (i = 1; i <= n; i++) if (length(a[i]) > 1) m = 1; k += m} END {print k + 0}')
check_eq "PharmCAT positions longer than one base in the image" "$MULTI" 51

# The sample: the gVCF and a variants-only VCF beside it, as step 03 leaves them.
cp "$FX" "${D}/${P}.g.vcf"
GENOME_DIR="$G" bcf view -Oz -o "${P}/vcf/${P}.g.vcf.gz" --write-index=tbi "${P}/vcf/${P}.g.vcf"
GENOME_DIR="$G" bcf view -i 'GT="alt"' -Oz -o "${P}/vcf/${P}.vcf.gz" --write-index=tbi "${P}/vcf/${P}.g.vcf"
rm -f "${D}/${P}.g.vcf"
check_eq "the variants-only VCF holds the one variant record" "$(GENOME_DIR="$G" vcf_count "${P}/vcf/${P}.vcf.gz")" 1

# match_table FILE: chrom, position, gene and state (missing, call:<vcfCall>,
# or for position 0 uncallable:<haplotypes>) of PharmCAT's match.json.
match_table() {
  python3 - "$1" <<'PY'
import json, sys
try:
    d = json.load(open(sys.argv[1]))
except Exception as e:
    print(f"unreadable match.json: {e}", file=sys.stderr)
    sys.exit(0)
for g in d.get("results") or []:
    ch, gene = g.get("chromosome"), g.get("gene")
    for m in (g.get("matchData") or {}).get("missingPositions") or []:
        print(f"{m.get('chromosome') or ch}\t{m.get('position')}\t{gene}\tmissing")
    for v in g.get("variants") or []:
        print(f"{ch}\t{v.get('position')}\t{gene}\tcall:{v.get('vcfCall')}")
    print(f"{ch}\t0\t{gene}\tuncallable:{','.join(sorted(g.get('uncallableHaplotypes') or []))}")
PY
}
listed() { [ -f "$1" ] || { echo "no file"; return; }; awk -F'\t' -v c="$2" -v p="$3" '!/^#/ && $1 == c && $2 == p {k++} END {print k + 0}' "$1"; }
state() { awk -F'\t' -v c="$2" -v p="$3" '$1 == c && $2 == p {print $4}' <<< "$1" | LC_ALL=C sort | paste -sd' ' -; }
ugt1a1_uncallable() { awk -F'\t' '$2 == 0 && $3 == "UGT1A1" {sub(/^uncallable:/, "", $4); print $4}' <<< "$1"; }

# --- control: the recipe without the filler -----------------------------------
CTL="${CASE_TMP}/control"
rm -rf "$CTL"
mkdir -p "$CTL"
# shellcheck disable=SC2016  # $1 to $3 belong to the inner bash
docker run --rm --network none --user "$(id -u):$(id -g)" -e HOME=/tmp -w /ctl \
  -v "${G}:/genome:ro" -v "${CTL}:/ctl" "$PHARMCAT_IMAGE" \
  bash -euo pipefail -c '
    bcftools convert --gvcf2vcf -f "$1" -R /pharmcat/pharmcat_regions.bed -Ou "$2" \
      | bcftools view --trim-alt-alleles -i "GT!=\"mis\"" -Oz -o "$3/before.pgx_regions.vcf.gz" --write-index=tbi
    python3 /pharmcat/pharmcat_vcf_preprocessor -vcf "$3/before.pgx_regions.vcf.gz" -refFna "$1" -o "$3" -bf before
    gzip -dc "$3/before.preprocessed.vcf.bgz" > "$3/before.vcf"
    java -jar /pharmcat/pharmcat.jar -vcf "$3/before.vcf" -o "$3" -bf before -matcher' \
  _ "/genome/${REF}" "/genome/${P}/vcf/${P}.g.vcf.gz" /ctl > "${CTL}/log.txt" 2>&1
CTL_RC=$?
cat "${CTL}/log.txt"
check_eq "control: the recipe without the filler runs" "$CTL_RC" 0
CTL_MISSING="${CTL}/before.missing_pgx_var.vcf"
CTL_MATCH=$(match_table "${CTL}/before.match.json")
check "control: PharmCAT's preprocessor wrote a missing-positions file" test -s "$CTL_MISSING"
CTL_LOG=$(cat "${CTL}/log.txt")
check "control: the preprocessor sets aside the one-base record at chr2:233760233" \
  has 'chr2:233760233 REF=C ALT=\. ?" which is not a valid GT format for INDELs' "$CTL_LOG"
check "control: the matcher discards it (REF C, PharmCAT expects CAT)" \
  has 'REF in VCF \(C\) does not match expected reference \(CAT\)' "$CTL_LOG"
while read -r c p _ g before; do
  check_ge "control: ${g} ${c}:${p} is listed missing" "$(listed "$CTL_MISSING" "$c" "$p")" 1
  check_eq "control: ${g} ${c}:${p} in match.json" "$(state "$CTL_MATCH" "$c" "$p")" "$before"
done <<< "$FILL"
check "control: UGT1A1 *28 is uncallable" has '(^|,)\*28(,|$)' "$(ugt1a1_uncallable "$CTL_MATCH")"

# --- step 07 ----------------------------------------------------------------------
GENOME_DIR="$G" run_step 07-pharmacogenomics.sh "$P"
check_step_exit 07-pharmacogenomics.sh
LOG=$(cat "$STEP_LOG")
check "step 07 reads the gVCF" has "^Input: .*${P}\.g\.vcf\.gz \(reference blocks expanded" "$LOG"
check "step 07 fills 5 of the 51 positions" \
  has 'pgx_refblock_fill: 5 of 51 PharmCAT positions longer than one base called 0/0 from reference blocks with GQ >= 20' "$LOG"
check "the fill names UGT1A1 chr2:233760233" has 'UGT1A1 chr2:233760233' "$LOG"
check "no one-base record at chr2:233760233 reaches the preprocessor" lacks 'chr2:233760233 REF=C ALT=\.' "$LOG"
check "the matcher discards nothing at chr2:233760233" lacks 'does not match expected reference \(CAT\)' "$LOG"
MISSING="${D}/${P}.missing_pgx_var.vcf"
check "step 07 wrote the missing-positions file" test -s "$MISSING"
MATCH=$(match_table "${D}/${P}.match.json")
check "step 07 wrote match.json" test -s "${D}/${P}.match.json"
echo "match.json at the planted positions:"
while read -r c p _; do awk -F'\t' -v c="$c" -v p="$p" '$1 == c && $2 == p' <<< "$MATCH"; done <<< "$FILL
$LEFT"
while read -r c p r g _; do
  check_eq "${g} ${c}:${p} is no longer listed missing" "$(listed "$MISSING" "$c" "$p")" 0
  call=$(state "$MATCH" "$c" "$p")
  check "${g} ${c}:${p} is called reference in match.json (${call:-nothing})" \
    python3 -c 'import re, sys; a = re.split(r"[/|]", sys.argv[1].removeprefix("call:")); sys.exit(not (sys.argv[1].startswith("call:") and a and all(x == sys.argv[2] for x in a)))' \
    "$call" "$r"
done <<< "$FILL"
while read -r c p why; do
  check_ge "${c}:${p} stays missing (${why//_/ })" "$(listed "$MISSING" "$c" "$p")" 1
  check_eq "${c}:${p} is missing in match.json" "$(state "$MATCH" "$c" "$p")" missing
done <<< "$LEFT"
U=$(ugt1a1_uncallable "$MATCH")
check "UGT1A1 *28 is callable now (uncallable: ${U:-none})" lacks '(^|,)\*28(,|$)' "$U"

finish
