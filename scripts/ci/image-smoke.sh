#!/usr/bin/env bash
# image-smoke.sh: run pinned images on the e2e fixture and check what they
# wrote, or check that every pinned image still exists.
#
# Usage:
#   scripts/ci/image-smoke.sh [--full] VAR...  run the rows of tests/smoke/commands.tsv
#                                              for these *_IMAGE variables of versions.env
#   scripts/ci/image-smoke.sh --manifest       `docker manifest inspect` every image in
#                                              versions.env and print its platforms
#   scripts/ci/image-smoke.sh --manifest-self-test
#                                              prove --manifest passes a tag rebuilt in
#                                              place and fails a digest that does not exist
#   scripts/ci/image-smoke.sh --check          parse the whole table (fields, options, input
#                                              names) and check that every *_IMAGE of
#                                              versions.env has a row (no docker)
#   scripts/ci/image-smoke.sh --self-test      prove --check fails on a bad table, and that a
#                                              host check never imports a module the image
#                                              planted in its row directory
#
# --full also runs the rows marked `full`: they need data too big for every
# pull request (the monthly run and a run started by hand pass it).
#
# Env: SMOKE_WORK   work area (default ${RUNNER_TEMP:-/tmp}/image-smoke)
#      GH_TOKEN     for `gh release download` of the fixture
#      ROW_TIMEOUT  seconds per row (default 3600)
#
# The table format and the inputs a row can ask for: tests/smoke/README.md.
# Each row runs `sh -c COMMAND` in its image with the prepared inputs at /in
# (read-only), tests/smoke at /smoke (read-only) and an empty /out as the
# working directory, as the calling user (root for a `root` row), without
# network (unless `net`), and with every versions.env variable in its
# environment. Then EXPECT runs on the host in
# that /out directory. A row passes when COMMAND exits 0 and every check of
# EXPECT passes; the run exits 1 when any row fails. The helper binaries the
# scripts call in an image (scripts/ci/check-container-helpers.sh --list) are
# probed in that image too.
set -uo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
SMOKE_DIR="${REPO}/tests/smoke"
# SMOKE_TABLE and SMOKE_VERSIONS exist for --self-test.
TABLE="${SMOKE_TABLE:-${SMOKE_DIR}/commands.tsv}"
VERSIONS="${SMOKE_VERSIONS:-${REPO}/versions.env}"
WORK="${SMOKE_WORK:-${RUNNER_TEMP:-/tmp}/image-smoke}"
FX="${WORK}/fixture"
IN="${WORK}/in"
LOGS="${WORK}/logs"
ROWS_DIR="${WORK}/rows"
ROW_TIMEOUT=${ROW_TIMEOUT:-3600}
GH_REPO=${GITHUB_REPOSITORY:-GeiserX/Personal-Genome-Pipeline}
TAG=$(tr -d '[:space:]' < "${REPO}/tests/fixtures/VERSION")
SUMMARY=${GITHUB_STEP_SUMMARY:-/dev/null}
NEEDS_KNOWN="ref reads bam longreads longbam truth fullref slice vcf vcf50 pgxvcf revel sv dels bundle cyrius mito hlareads qc ehcatalog chain annotsv"
OPTS_KNOWN="root net full also= rw="

# The mini reference: chr20:10,000,001-10,500,000 of the fixture reference as
# contig chr20 (local coordinates: subtract 10,000,000), with 3,000 random
# bases inserted after local position 450,000. Reads of the sample have no
# such bases, so every caller sees a 3 kb deletion at 450,001-453,000.
# Truth, regions and checks of small variants stop at 449,000.
MINI_LEN=503000
PLANT_POS=450000
PLANT_LEN=3000
SMALL_END=449000

# shellcheck source=../../versions.env
. "$VERSIONS"

# ------------------------------------------------------------------- table
R_VAR=() R_OPTS=() R_NEEDS=() R_CMD=() R_EXPECT=() R_LINE=()
TABLE_ERRORS=()
load_table() {
  local line n=0 tabs a b c d e
  while IFS= read -r line || [ -n "$line" ]; do
    n=$((n + 1))
    [ -z "$line" ] && continue
    [[ "$line" == \#* ]] && continue
    tabs=${line//[^$'\t']/}
    if [ "${#tabs}" -ne 4 ]; then
      TABLE_ERRORS+=("commands.tsv:${n}: ${#tabs} tabs, want 4 (five columns)")
      continue
    fi
    IFS=$'\t' read -r a b c d e <<< "$line"
    R_VAR+=("$a") R_OPTS+=("$b") R_NEEDS+=("$c") R_CMD+=("$d") R_EXPECT+=("$e") R_LINE+=("$n")
  done < "$TABLE"
  local i o v var list
  for i in "${!R_VAR[@]}"; do
    n=${R_LINE[$i]}
    var=${R_VAR[$i]}
    [[ "$var" =~ ^[A-Z0-9_]+_IMAGE$ ]] || TABLE_ERRORS+=("commands.tsv:${n}: '${var}' is not an *_IMAGE name")
    [ -n "${!var:-}" ] || TABLE_ERRORS+=("commands.tsv:${n}: ${var} is not set by versions.env")
    if [ "${R_OPTS[$i]}" != "-" ]; then
      IFS=, read -r -a list <<< "${R_OPTS[$i]}"
      for o in "${list[@]}"; do
        case "$o" in
          root|net|full) ;;
          also=*) v=${o#also=}; [ -n "${!v:-}" ] || TABLE_ERRORS+=("commands.tsv:${n}: also=${v} is not set by versions.env") ;;
          rw=*) grep -qw -- "${o#rw=}" <<< "$NEEDS_KNOWN" || TABLE_ERRORS+=("commands.tsv:${n}: rw=${o#rw=} is not an input name") ;;
          *) TABLE_ERRORS+=("commands.tsv:${n}: unknown option '${o}' (known: ${OPTS_KNOWN})") ;;
        esac
      done
    fi
    if [ "${R_NEEDS[$i]}" != "-" ]; then
      IFS=, read -r -a list <<< "${R_NEEDS[$i]}"
      for o in "${list[@]}"; do
        grep -qw -- "$o" <<< "$NEEDS_KNOWN" || TABLE_ERRORS+=("commands.tsv:${n}: unknown input '${o}' (known: ${NEEDS_KNOWN})")
      done
    fi
  done
}

# Every *_IMAGE of versions.env needs a row: a new image cannot be added untested.
check_coverage() {
  local var
  while read -r var; do
    printf '%s\n' "${R_VAR[@]}" | grep -qx -- "$var" ||
      TABLE_ERRORS+=("${var} (versions.env) has no row in commands.tsv: add one that runs the tool on the fixture")
  done < <(grep -oE '^[A-Z][A-Z0-9_]*_IMAGE=' "$VERSIONS" | tr -d = | sort -u)
}

# Plants one fault per case in a scratch table and versions.env and expects
# --check to name it; the clean pair must pass.
self_test() {
  local tmp fails=0 out rc
  tmp=$(mktemp -d)
  printf '%s\n' 'A_IMAGE="a:1"' 'B_IMAGE="b:1"' 'A_DATA="1"' > "${tmp}/versions.env"
  printf 'A_IMAGE\talso=A_DATA\tref\ttrue\ttrue\nB_IMAGE\t-\t-\ttrue\ttrue\n' > "${tmp}/good.tsv"
  t() {  # t DESCRIPTION WANT_REGEX TABLE_LINES...
    local desc=$1 want=$2; shift 2
    printf '%s\n' "$@" > "${tmp}/t.tsv"
    out=$(env -i PATH="$PATH" SMOKE_TABLE="${tmp}/t.tsv" SMOKE_VERSIONS="${tmp}/versions.env" bash "$0" --check 2>&1); rc=$?
    if [ -z "$want" ]; then
      if [ "$rc" -eq 0 ]; then echo "[PASS] ${desc}"; else echo "[FAIL] ${desc}: ${out}"; fails=$((fails + 1)); fi
    elif [ "$rc" -ne 0 ] && grep -Eq -- "$want" <<< "$out"; then echo "[PASS] ${desc}"
    else echo "[FAIL] ${desc} (exit ${rc}, want /${want}/): ${out}"; fails=$((fails + 1)); fi
  }
  local A=$'A_IMAGE\talso=A_DATA\tref\ttrue\ttrue' B=$'B_IMAGE\t-\t-\ttrue\ttrue'
  t "the clean table passes" "" "$A" "$B"
  t "an image without a row fails" "B_IMAGE .*has no row" "$A"
  t "a row for an unset variable fails" "C_IMAGE is not set" "$A" "$B" $'C_IMAGE\t-\t-\ttrue\ttrue'
  t "an unknown option fails" "unknown option 'fast'" "$A" $'B_IMAGE\tfast\t-\ttrue\ttrue'
  t "an unknown input fails" "unknown input 'genome'" "$A" $'B_IMAGE\t-\tgenome\ttrue\ttrue'
  t "also= on an unset variable fails" "also=B_DATA is not set" "$A" $'B_IMAGE\talso=B_DATA\t-\ttrue\ttrue'
  t "a row with four columns fails" "3 tabs, want 4" "$A" $'B_IMAGE\t-\t-\ttrue'
  # A row's checks run in its row directory, which the image wrote. A json.py
  # or csv.py the image planted there must not be imported by a host-side
  # check (python3 -I keeps the working directory off sys.path).
  mkdir "${tmp}/row"
  for m in json csv; do printf 'raise SystemExit("planted %s.py was imported")\n' "$m" > "${tmp}/row/${m}.py"; done
  printf '{"genes": {"CYP2C19": {"sourceDiplotypes": [{"allele1": {"name": "*1"}, "allele2": {"name": "*2"}}]}}}\n' > "${tmp}/row/r.json"
  printf 'a,b\n1,2\n' > "${tmp}/row/t.csv"
  out=$(cd "${tmp}/row" && CHECK_FAILS=0 && py "json loads" "import json; print(json.dumps(1))" \
          && pharmcat_called r.json 2>/dev/null && csv_cell t.csv a=1 b && echo "fails=${CHECK_FAILS}")
  if [ "$out" = $'[PASS] json loads (1)\n1\n2\nfails=0' ]; then echo "[PASS] a module the image planted in the row directory is not imported"
  else echo "[FAIL] a module the image planted in the row directory is not imported: ${out//$'\n'/ | }"; fails=$((fails + 1)); fi
  rm -rf "$tmp"
  if [ "$fails" -gt 0 ]; then echo "self-test: ${fails} case(s) failed" >&2; return 1; fi
  echo "self-test: all cases passed"
}

has_opt() { [[ ",${1}," == *",${2},"* ]]; }

# --------------------------------------------------------------- manifest
# manifest_ref IMAGE: the reference the manifest check resolves for a
# versions.env image. A name:tag@digest pin resolves as name@digest: with the
# tag in the reference Docker fetches the tag and fails with "manifest
# verification failed" as soon as the publisher rebuilds it in place, though
# the pinned digest still exists and is what every step pulls. Tag drift is
# Renovate's to report, as a digest update. An image without a digest is
# resolved as written.
manifest_ref() {
  local img=$1 name last
  case "$img" in *@*) ;; *) printf '%s\n' "$img"; return ;; esac
  name=${img%@*}
  last=${name##*/}
  # A colon in the last path component is the tag; one before a slash is a
  # registry port.
  [[ "$last" == *:* ]] && name=${name%:*}
  printf '%s@%s\n' "$name" "${img##*@}"
}

# inspect_ref REF: `docker manifest inspect -v`, three tries; the JSON on
# stdout, or "FAILED: <docker's message>" and exit 1.
inspect_ref() {
  local out tries
  for tries in 1 2 3; do
    if out=$(docker manifest inspect -v "$1" 2>&1); then printf '%s\n' "$out"; return 0; fi
    [ "$tries" -lt 3 ] && sleep 10
  done
  printf 'FAILED: %s\n' "$out"
  return 1
}

manifest() {
  local var img ref out plats fail=0 rows=()
  while read -r var; do
    img=${!var}
    ref=$(manifest_ref "$img")
    if ! out=$(inspect_ref "$ref"); then
      echo "FAIL ${var}: ${ref} cannot be resolved: $(head -c 300 <<< "${out#FAILED: }")"
      rows+=("| ${var} | \`${img}\` | **missing** |")
      fail=1
      continue
    fi
    plats=$(jq -r 'if type == "array" then .[] else . end | .Descriptor.platform // empty
                   | select(.os != "unknown") | "\(.os)/\(.architecture)\(if .variant then "/" + .variant else "" end)"' \
              <<< "$out" 2>/dev/null | sort -u | paste -sd ' ' -)
    echo "OK   ${var}: ${img} (${plats:-platform not stated})"
    rows+=("| ${var} | \`${img}\` | ${plats:-not stated} |")
  done < <(grep -oE '^[A-Z0-9_]+_IMAGE=' "$VERSIONS" | tr -d =)
  {
    echo "### Pinned images in versions.env"
    echo
    echo "| Variable | Image | Platforms |"
    echo "|---|---|---|"
    printf '%s\n' "${rows[@]}"
  } >> "$SUMMARY"
  if [ "$fail" -ne 0 ]; then
    echo "ERROR: at least one pinned image no longer exists (see above)." >&2
    return 1
  fi
  echo "All ${#rows[@]} pinned images exist."
}

# manifest_self_test: the manifest check against the real registry. A pin
# whose tag no longer serves its digest must still pass, and a digest that
# does not exist must fail. The moved case pairs PYTHON_IMAGE's own digest with
# a tag that does not exist, so it needs nothing from the registry that the
# manifest check does not need anyway. With the tag kept in the reference,
# Docker fetches the tag and the case fails.
manifest_self_test() {
  local fails=0 got img want moved wrong name
  while IFS='|' read -r img want; do
    got=$(manifest_ref "$img")
    if [ "$got" = "$want" ]; then echo "[PASS] manifest_ref ${img} -> ${got}"
    else echo "[FAIL] manifest_ref ${img}: got ${got}, want ${want}"; fails=$((fails + 1)); fi
  done <<'CASES'
example/tool:1.0--h1@sha256:abc|example/tool@sha256:abc
localhost:5000/example/tool:2@sha256:abc|localhost:5000/example/tool@sha256:abc
localhost:5000/example/tool@sha256:abc|localhost:5000/example/tool@sha256:abc
example/tool:1.0|example/tool:1.0
CASES
  name=${PYTHON_IMAGE%@*}
  name=${name%:*}
  moved="$(manifest_ref "${name}:no-such-tag@${PYTHON_IMAGE##*@}")"
  if inspect_ref "$moved" >/dev/null; then echo "[PASS] a pin whose tag moved resolves by its digest (${name}:no-such-tag -> ${moved})"
  else echo "[FAIL] a pin whose tag moved resolves by its digest (${name}:no-such-tag -> ${moved})"; fails=$((fails + 1)); fi
  wrong="$(manifest_ref "${name}@sha256:$(printf '0%.0s' $(seq 64))")"
  if inspect_ref "$wrong" >/dev/null; then echo "[FAIL] a digest that does not exist fails (${wrong})"; fails=$((fails + 1))
  else echo "[PASS] a digest that does not exist fails (${wrong})"; fi
  if [ "$fails" -gt 0 ]; then echo "manifest self-test: ${fails} case(s) failed" >&2; return 1; fi
  echo "manifest self-test: all cases passed"
}

# ----------------------------------------------------------- docker helpers
ME="$(id -u):$(id -g)"
# hrun IMAGE ARGS...: a helper container over the work area, as the caller.
hrun() {
  local img=$1; shift
  docker run --rm -i --network none --user "$ME" -e HOME=/tmp \
    -v "${FX}:/fx:ro" -v "${IN}:/in" -v "${ROW_DIR:-$IN}:/out" -w /out "$img" "$@"
}
hsam() { hrun "$SAMTOOLS_IMAGE" samtools "$@"; }
hbcf() { hrun "$BCFTOOLS_IMAGE" bcftools "$@"; }

# ---------------------------------------------------------------- fixture
fx_get() {
  local f
  mkdir -p "$FX"
  if [ ! -s "${FX}/SHA256SUMS" ]; then
    gh release download "$TAG" -R "$GH_REPO" -p SHA256SUMS -D "$FX" --clobber || return 1
  fi
  for f in "$@"; do
    if [ ! -s "${FX}/${f}" ]; then
      gh release download "$TAG" -R "$GH_REPO" -p "$f" -D "$FX" --clobber || return 1
    fi
    (cd "$FX" && grep -E "  ${f//./\\.}\$" SHA256SUMS | sha256sum --quiet -c -) || {
      echo "fixture file ${f} does not match SHA256SUMS" >&2; rm -f "${FX}/${f}"; return 1; }
  done
}
# fx_link FILE...: the fixture files under /in too.
fx_link() {
  local f
  fx_get "$@" || return 1
  for f in "$@"; do ln -f "${FX}/${f}" "${IN}/${f}" 2>/dev/null || cp "${FX}/${f}" "${IN}/${f}"; done
}

# ------------------------------------------------------------------ inputs
# need NAME: prepares input NAME once; later calls return its first result.
declare -A NEED_STATE=()
need() {
  local n=$1 rc
  case "${NEED_STATE[$n]:-}" in
    ok) return 0 ;;
    failed) return 1 ;;
  esac
  echo "--- preparing input: ${n}"
  "need_${n}"; rc=$?
  if [ "$rc" -eq 0 ]; then NEED_STATE[$n]=ok; else NEED_STATE[$n]=failed; echo "input ${n} could not be prepared (exit ${rc})" >&2; fi
  return "$rc"
}

need_ref() {
  fx_get fixture_ref.fa.gz fixture_ref.fa.gz.fai fixture_ref.fa.gz.gzi || return 1
  hsam faidx /fx/fixture_ref.fa.gz chr20:10000001-10500000 > "${IN}/slice.fa.tmp" || return 1
  python3 -I - "${IN}/slice.fa.tmp" "${IN}/orig.fa" "${IN}/mini.fa" "$PLANT_POS" "$PLANT_LEN" <<'PY' || return 1
import random, sys
src, orig, mini, pos, n = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
seq = "".join(l.strip() for l in open(src) if not l.startswith(">")).upper()
assert len(seq) == 500000, len(seq)
def write(path, s):
    with open(path, "w") as f:
        f.write(">chr20\n")
        for i in range(0, len(s), 60):
            f.write(s[i:i + 60] + "\n")
write(orig, seq)
random.seed(20)
write(mini, seq[:pos] + "".join(random.choice("ACGT") for _ in range(n)) + seq[pos:])
PY
  rm -f "${IN}/slice.fa.tmp"
  hsam faidx /in/orig.fa && hsam faidx /in/mini.fa && hsam dict -o /in/mini.dict /in/mini.fa || return 1
  [ "$(cut -f2 "${IN}/mini.fa.fai")" = "$MINI_LEN" ]
}

# Read pairs of the chr20 slice (secondary and supplementary records dropped).
slice_reads() {  # slice_reads REGION PREFIX
  fx_get HG002_slice.bam HG002_slice.bam.bai || return 1
  hrun "$SAMTOOLS_IMAGE" sh -c "samtools view -u -F 0x900 /fx/HG002_slice.bam $1 \
    | samtools collate -O -u - /tmp/collate \
    | samtools fastq -n -1 /in/$2_R1.fq.gz -2 /in/$2_R2.fq.gz -0 /dev/null -s /dev/null -"
}
need_reads()    { slice_reads chr20:10000001-10500000 reads; }
need_hlareads() { slice_reads chr6:29900001-33100000 hla; }

need_bam() {
  need ref && need reads || return 1
  docker run --rm --network none --user "$ME" -v "${IN}:/in" "$MINIMAP2_IMAGE" \
      minimap2 -t 4 -a -x sr -R '@RG\tID:HG002\tSM:HG002\tPL:ILLUMINA\tLB:HG002' \
      /in/mini.fa /in/reads_R1.fq.gz /in/reads_R2.fq.gz \
    | hsam sort -@ 2 -m 500M -o /in/mini.bam - && hsam index /in/mini.bam
}

need_truth() {
  need ref || return 1
  fx_get HG002_truth_chr20.vcf.gz HG002_truth_chr20.vcf.gz.tbi HG002_truth_chr20.bed || return 1
  gzip -dc "${FX}/HG002_truth_chr20.vcf.gz" | awk -v off=10000000 -v max="$SMALL_END" -v len="$MINI_LEN" '
    BEGIN {OFS = "\t"}
    /^##contig/ {next}
    /^#CHROM/ {print "##contig=<ID=chr20,length=" len ">"; print; next}
    /^#/ {print; next}
    $1 == "chr20" {p = $2 - off; if (p >= 1 && p + length($4) - 1 < max) {$2 = p; print}}' \
    > "${IN}/truth_local.vcf" || return 1
  awk -v off=10000000 -v max="$SMALL_END" 'BEGIN {OFS = "\t"}
    $1 == "chr20" {s = $2 - off; e = $3 - off; if (s < 0) s = 0; if (e > max) e = max; if (e > s) print "chr20", s, e}' \
    "${FX}/HG002_truth_chr20.bed" > "${IN}/truth.bed" || return 1
  hbcf view -T /in/truth.bed -Oz -o /in/truth.vcf.gz /in/truth_local.vcf && hbcf index -f -t /in/truth.vcf.gz || return 1
  # The truth without every tenth SNV: hap.py must report a recall near 0.9.
  hbcf view /in/truth.vcf.gz | awk '/^#/ {print; next}
    length($4) == 1 && length($5) == 1 {n++; if (n % 10 == 0) next} {print}' > "${IN}/query90.vcf" || return 1
  hbcf view -Oz -o /in/query90.vcf.gz /in/query90.vcf && hbcf index -f -t /in/query90.vcf.gz || return 1
  echo "truth: $(grep -vc '^#' "${IN}/truth_local.vcf") records before the bed, $(hbcf view -H /in/truth.vcf.gz | wc -l) after"
}

# Synthetic HiFi reads: both haplotypes of the truth applied to the original
# slice (without the inserted bases), cut into 10 kb reads every 1 kb, every
# other read reverse-complemented, base quality 30.
need_longreads() {
  need ref && need truth || return 1
  hbcf consensus -H 1 -f /in/orig.fa /in/truth.vcf.gz > "${IN}/hap1.fa" &&
    hbcf consensus -H 2 -f /in/orig.fa /in/truth.vcf.gz > "${IN}/hap2.fa" || return 1
  python3 -I - "${IN}/hap1.fa" "${IN}/hap2.fa" "${IN}/long.fq.gz" <<'PY'
import gzip, sys
comp = str.maketrans("ACGTN", "TGCAN")
with gzip.open(sys.argv[3], "wt") as out:
    k = 0
    for h, path in enumerate(sys.argv[1:3], 1):
        seq = "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()
        for start in range(0, len(seq) - 10000 + 1, 1000):
            r = seq[start:start + 10000]
            if k % 2:
                r = r.translate(comp)[::-1]
            out.write(f"@hap{h}_{start}\n{r}\n+\n{'?' * len(r)}\n")
            k += 1
print(k, "reads")
PY
}
need_longbam() {
  need longreads || return 1
  docker run --rm --network none --user "$ME" -v "${IN}:/in" "$MINIMAP2_IMAGE" \
      minimap2 -t 4 -a -x map-hifi -R '@RG\tID:HG002\tSM:HG002\tPL:PACBIO' /in/mini.fa /in/long.fq.gz \
    | hsam sort -o /in/long.bam - && hsam index /in/long.bam
}

need_fullref() {
  fx_get fixture_ref.fa.gz fixture_ref.fa.gz.fai fixture_ref.dict || return 1
  gzip -dc "${FX}/fixture_ref.fa.gz" > "${IN}/ref.fa.tmp" && mv "${IN}/ref.fa.tmp" "${IN}/ref.fa" || return 1
  cp "${FX}/fixture_ref.fa.gz.fai" "${IN}/ref.fa.fai" && cp "${FX}/fixture_ref.dict" "${IN}/ref.dict"
}
# The fixture's GIAB BAM has no read group; tools that name the sample after
# it (pypgx, bcftools mpileup) would take the file path instead. Add one.
need_slice() {
  fx_get HG002_slice.bam HG002_slice.bam.bai || return 1
  hsam addreplacerg -r '@RG\tID:HG002\tSM:HG002\tPL:ILLUMINA' -o /in/HG002_slice.bam /fx/HG002_slice.bam &&
    hsam index /in/HG002_slice.bam
}
need_vcf() {
  fx_get HG002_vep.vcf || return 1
  hbcf view -Oz -o /in/sample.vcf.gz /fx/HG002_vep.vcf && hbcf index -f -t /in/sample.vcf.gz
}
need_vcf50() {
  need vcf || return 1
  hbcf annotate -x INFO/CSQ /in/sample.vcf.gz | awk '/^#/ || n++ < 50' > "${IN}/sample50.vcf" &&
    [ "$(grep -vc '^#' "${IN}/sample50.vcf")" -eq 50 ]
}
# Every site of the CYP2C19 and CYP2C9 slice called from the GIAB reads,
# reference sites included: PharmCAT counts a position the VCF lacks as
# missing, not as reference, and the fixture's VCF holds variants only.
# Only FORMAT/GT is kept: PharmCAT's normalisation cannot merge FORMAT/AD
# where an indel and a reference site share a position.
need_pgxvcf() {
  need fullref && need slice || return 1
  hrun "$BCFTOOLS_IMAGE" sh -c 'bcftools mpileup -r chr10:94700001-95000000 -f /in/ref.fa -Ou /in/HG002_slice.bam \
    | bcftools call -m -Ou | bcftools annotate -x "^FORMAT/GT" -Oz -o /in/pgx.vcf.gz && bcftools index -f -t /in/pgx.vcf.gz' || return 1
  echo "pgx.vcf.gz: $(hbcf view -H /in/pgx.vcf.gz | wc -l) sites, $(hbcf view -H -i 'GT="alt"' /in/pgx.vcf.gz | wc -l) with an ALT allele"
}
need_revel() { fx_link revel_synthetic.tsv.gz revel_synthetic.tsv.gz.tbi; }
need_sv() { fx_link HG002_sv_manta_style.vcf.gz HG002_sv_manta_style.vcf.gz.tbi; }
# Two deletions for duphold: the planted one (no reads) and a control of the
# same size where the depth is normal.
need_dels() {
  need ref || return 1
  {
    printf '##fileformat=VCFv4.2\n##contig=<ID=chr20,length=%s>\n' "$MINI_LEN"
    printf '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
    printf '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
    printf '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length">\n'
    printf '##ALT=<ID=DEL,Description="Deletion">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n'
    printf 'chr20\t200000\tcontrol\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=%s;SVLEN=-%s\tGT\t0/1\n' $((200000 + PLANT_LEN)) "$PLANT_LEN"
    printf 'chr20\t%s\tplanted\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=%s;SVLEN=-%s\tGT\t1/1\n' "$PLANT_POS" $((PLANT_POS + PLANT_LEN)) "$PLANT_LEN"
  } > "${IN}/dels.vcf"
}
need_bundle() {
  rm -rf "${IN}/pypgx-bundle"
  git clone -q --branch "$PYPGX_BUNDLE_VERSION" --depth 1 https://github.com/sbslee/pypgx-bundle.git "${IN}/pypgx-bundle"
}
need_cyrius() {
  fx_link HG002_cyrius.bam HG002_cyrius.bam.bai || return 1
  cp "${REPO}/scripts/cyrius-constraints.txt" "${IN}/cyrius-constraints.txt"
}
need_mito() {
  need fullref && need slice || return 1
  hrun "$BCFTOOLS_IMAGE" sh -c 'bcftools mpileup -r chrM -f /in/ref.fa -d 5000 -a AD,DP -Ou /in/HG002_slice.bam \
    | bcftools call -mv --ploidy 1 -Oz -o /in/mito.vcf.gz && bcftools index -f -t /in/mito.vcf.gz' || return 1
  echo "chrM variants: $(hbcf view -H /in/mito.vcf.gz | wc -l)"
}
need_qc() {
  need bam || return 1
  mkdir -p "${IN}/qc"
  hsam flagstat /in/mini.bam > "${IN}/qc/HG002.flagstat" && hsam stats /in/mini.bam > "${IN}/qc/HG002.stats"
}
# A one-locus ExpansionHunter catalog: the first (CA)n run of 12 or more units
# in the mini reference, away from the slice ends and the planted bases.
need_ehcatalog() {
  need ref || return 1
  python3 -I - "${IN}/mini.fa" "${IN}/eh_catalog.json" "$SMALL_END" <<'PY'
import json, re, sys
seq = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">"))
end = int(sys.argv[3])
m = next(m for m in re.finditer(r"(?:CA){12,}", seq) if 20000 < m.start() and m.end() < end - 20000)
json.dump([{"LocusId": "SMOKE1", "LocusStructure": "(CA)*",
            "ReferenceRegion": f"chr20:{m.start()}-{m.end()}", "VariantType": "Repeat"}],
          open(sys.argv[2], "w"), indent=1)
print("repeat", m.start(), m.end(), len(m.group()) // 2, "units")
PY
}
# A chain that maps chr20 onto itself over the region of the small variants.
need_chain() {
  need ref || return 1
  printf 'chain 1000 chr20 %s + 0 %s chr20 %s + 0 %s 1\n%s\n\n' \
    "$MINI_LEN" "$SMALL_END" "$MINI_LEN" "$SMALL_END" "$SMALL_END" > "${IN}/identity.chain"
}
# AnnotSV's annotation data, as setup.sh fetches it (5.3 GB; the server is slow).
need_annotsv() {
  local url="https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_${ANNOTSV_ANNOTATIONS_VERSION}.tar.gz"
  [ -d "${IN}/annotsv/Annotations_Human/Genes/GRCh38" ] && return 0
  rm -rf "${IN}/annotsv.part" && mkdir -p "${IN}/annotsv.part" || return 1
  echo "downloading ${url}"
  curl -fsSL --retry 5 --retry-delay 30 "$url" | tar -xz -C "${IN}/annotsv.part" || return 1
  mv "${IN}/annotsv.part" "${IN}/annotsv"
  [ -d "${IN}/annotsv/Annotations_Human/Genes/GRCh38" ]
}

# ------------------------------------------------------------ the checks
# EXPECT is evaluated in the row's output directory with these helpers. Each
# check prints one [PASS] or [FAIL] line; value helpers print a value.
CHECK_FAILS=0
pass() { echo "[PASS] $*"; }
fail() { echo "[FAIL] $*"; CHECK_FAILS=$((CHECK_FAILS + 1)); }
num() { [[ "$1" =~ ^-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$ ]]; }
at_least() {
  if num "$2" && awk -v a="$2" -v b="$3" 'BEGIN {exit !(a + 0 >= b + 0)}'; then pass "$1 (${2} >= ${3})"
  else fail "$1 (got '${2}', want >= ${3})"; fi
}
between() {
  if num "$2" && awk -v a="$2" -v lo="$3" -v hi="$4" 'BEGIN {exit !(a + 0 >= lo + 0 && a + 0 <= hi + 0)}'; then
    pass "$1 (${2} in ${3}..${4})"
  else fail "$1 (got '${2}', want ${3}..${4})"; fi
}
same() { if [ -n "$2" ] && [ "$2" = "$3" ]; then pass "$1 (${2})"; else fail "$1 (got '${2}', want '${3}')"; fi; }
has() { if [ -n "$2" ] && grep -Eq -- "$3" <<< "$2"; then pass "$1 (${2})"; else fail "$1 (got '${2}', want /${3}/)"; fi; }
contains() { if [ -f "$2" ] && grep -Eq -- "$3" "$2"; then pass "$1"; else fail "$1 (/${3}/ not in ${2})"; fi; }
nonempty() { if [ -s "$2" ]; then pass "$1"; else fail "$1 (${2} is missing or empty)"; fi; }
py() {
  local out
  if out=$(python3 -I -c "$2" 2>&1); then pass "$1${out:+ (${out//$'\n'/ })}"; else fail "$1: $(tail -n 3 <<< "$out" | tr '\n' ' ')"; fi
}
bcf() { hbcf "$@" 2>/dev/null; }
# vcf_records FILE: data records (paths relative to the row directory or under /in).
vcf_records() { bcf view -H "$1" | wc -l | tr -d ' '; }
# mapped_pct SAM|BAM: mapped primary reads, percent.
mapped_pct() {
  hsam flagstat "$1" 2>/dev/null | awk '/ primary mapped \(/ {gsub(/[(%]/, "", $6); print $6; exit}'
}
# snv_recall VCF: share of the truth SNVs that VCF calls (PASS or no filter,
# a non-reference genotype), multi-allelic records split first.
snv_recall() {
  local calls truth
  calls=$(bcf norm -m- "$1" | bcf view -v snps -f PASS,. - | bcf query -f '%POS:%REF:%ALT\t[%GT]\n' - \
            | awk -F'\t' '$2 ~ /[1-9]/ {print $1}' | sort -u)
  truth=$(bcf norm -m- /in/truth.vcf.gz | bcf view -v snps - | bcf query -f '%POS:%REF:%ALT\n' - | sort -u)
  awk 'NR == FNR {c[$1]; next} {n++; if ($1 in c) k++} END {if (n) printf "%.4f\n", k / n}' \
    <(printf '%s\n' "$calls") <(printf '%s\n' "$truth")
}
# sv_near VCF|BCF: records whose position is within 500 bp of either end of
# the planted deletion.
sv_near() {
  bcf view -H "$1" | awk -v a="$PLANT_POS" -v b=$((PLANT_POS + PLANT_LEN)) \
    '{d1 = $2 - a; d2 = $2 - b; if ((d1 < 0 ? -d1 : d1) <= 500 || (d2 < 0 ? -d2 : d2) <= 500) n++} END {print n + 0}'
}
# vcf_field VCF ID KEY: the INFO value, or else the first sample's FORMAT value,
# of KEY in the record whose ID column is ID.
vcf_field() {
  awk -F'\t' -v id="$2" -v k="$3" '!/^#/ && $3 == id {
      n = split($8, info, ";"); for (i = 1; i <= n; i++) if (index(info[i], k "=") == 1) {print substr(info[i], length(k) + 2); exit}
      m = split($9, f, ":"); split($10, v, ":"); for (i = 1; i <= m; i++) if (f[i] == k) {print v[i]; exit}
    }' "$1"
}
# tsv_cell FILE COLUMN: the first data row's value in COLUMN (a header name;
# quotes and a leading # are ignored).
tsv_cell() {
  awk -F'\t' -v c="$2" 'NR == 1 {for (i = 1; i <= NF; i++) {h = $i; gsub(/^#|"/, "", h); if (h == c) k = i}; next}
                        k {v = $k; gsub(/"/, "", v); print v; exit}' "$1"
}
# tsv_count FILE COLUMN ERE: data rows whose COLUMN matches ERE.
tsv_count() {
  awk -F'\t' -v c="$2" -v re="$3" 'NR == 1 {for (i = 1; i <= NF; i++) {h = $i; gsub(/^#|"/, "", h); if (h == c) k = i}; next}
                                   k && $k ~ re {n++} END {print n + 0}' "$1"
}
# csv_cell FILE KEY=VALUE[,KEY=VALUE] COLUMN: COLUMN of the first row that matches.
csv_cell() {
  python3 -I - "$1" "$2" "$3" <<'PY'
import csv, sys
want = dict(kv.split("=", 1) for kv in sys.argv[2].split(","))
for row in csv.DictReader(open(sys.argv[1])):
    if all(row.get(k) == v for k, v in want.items()):
        print(row[sys.argv[3]])
        break
PY
}
# pharmcat_called report.json: genes with a named diplotype.
pharmcat_called() {
  python3 -I - "$1" <<'PY'
import json, sys
try:
    data = json.load(open(sys.argv[1]))
except Exception as e:
    print(f"unreadable report: {e}", file=sys.stderr)
    sys.exit(0)
def called(g):
    for d in (g.get("sourceDiplotypes") or g.get("recommendationDiplotypes") or []):
        names = [(d.get(a) or {}).get("name", "") for a in ("allele1", "allele2")]
        if any(n and n.lower() not in ("unknown", "none", "?") for n in names):
            return True
    return False
out = set()
for key, val in (data.get("genes") or {}).items():
    if not isinstance(val, dict):
        continue
    if "sourceDiplotypes" in val or "recommendationDiplotypes" in val:
        if called(val):
            out.add(key)
    else:
        out.update(n for n, g in val.items() if isinstance(g, dict) and called(g))
print(" ".join(sorted(out)), file=sys.stderr)
print(len(out))
PY
}

# ------------------------------------------------------------- the runner
HELPER_TABLE=""
helpers_of() {  # helpers_of VAR: helper binaries the scripts call in VAR's image
  if [ -z "$HELPER_TABLE" ]; then
    HELPER_TABLE=$("${REPO}/scripts/ci/check-container-helpers.sh" --list 2>/dev/null || echo "?")
  fi
  awk -v v="$1" '$1 == v {print $2}' <<< "$HELPER_TABLE" | sort -u | tr '\n' ' '
}

pull() {
  local t
  for t in 1 2 3; do
    docker pull -q "$1" && return 0
    sleep 15
  done
  return 1
}

RES_VAR=() RES_ROW=() RES_RESULT=() RES_TIME=()
run_row() {  # run_row INDEX K
  local i=$1 k=$2 var=${R_VAR[$1]} opts=${R_OPTS[$1]} needs=${R_NEEDS[$1]} img start rc result
  local name="${var}-${k}" log o ok=true missing hb
  local -a list
  img=${!var}
  ROW_DIR="${ROWS_DIR}/${name}"
  log="${LOGS}/${name}.log"
  rm -rf "$ROW_DIR" && mkdir -p "$ROW_DIR"
  start=$(date +%s)
  local -a args=(run --rm --env-file "$ENV_FILE" -e HOME=/tmp -v "${IN}:/in:ro" -v "${SMOKE_DIR}:/smoke:ro"
                 -v "${ROW_DIR}:/out" -w /out --entrypoint sh)
  # An image may set a non-root USER, so root is asked for, not assumed.
  if has_opt "$opts" root; then args+=(--user 0:0); else args+=(--user "$ME"); fi
  has_opt "$opts" net || args+=(--network none)
  IFS=, read -r -a list <<< "$opts"
  for o in "${list[@]}"; do
    [[ "$o" == rw=* ]] && args+=(-v "${IN}/${o#rw=}:/in/${o#rw=}")
  done
  echo "::group::${name}: ${img}"
  {
    echo "row: tests/smoke/commands.tsv line ${R_LINE[$i]}"
    echo "image: ${img}"
    if [ "$needs" != "-" ]; then
      IFS=, read -r -a list <<< "$needs"
      for o in "${list[@]}"; do need "$o" || ok=false; done
    fi
    if ! $ok; then
      echo "ERROR: an input of this row could not be prepared"
      rc=98
    elif ! pull "$img"; then
      echo "ERROR: cannot pull ${img}"
      rc=99
    else
      echo "+ ${R_CMD[$i]}"
      timeout -k 30 "$ROW_TIMEOUT" docker "${args[@]}" "$img" -c "${R_CMD[$i]}"
      rc=$?
      echo "+ exit ${rc}"
    fi
    CHECK_FAILS=0
    if [ "$rc" -eq 0 ]; then
      pass "command exits 0"
      pushd "$ROW_DIR" >/dev/null || exit 2
      eval "${R_EXPECT[$i]}"
      popd >/dev/null || exit 2
      hb=$(helpers_of "$var")
      if [ "$k" -eq 1 ] && [ -n "${hb// /}" ]; then
        # shellcheck disable=SC2016,SC2086  # the probe runs in the image; one argument per binary
        missing=$(docker run --rm --entrypoint sh "$img" -c \
          'for b in "$@"; do command -v "$b" >/dev/null 2>&1 || printf "%s " "$b"; done' sh $hb 2>&1)
        if [ -z "$missing" ]; then pass "helper binaries the scripts call in this image: ${hb}"
        else fail "helper binaries missing from the image: ${missing}"; fi
      fi
    else
      fail "command exits 0 (exit ${rc})"
    fi
    echo "${name}: ${CHECK_FAILS} check(s) failed"
  } > >(tee "$log") 2>&1
  # The tee in the process substitution may still be writing.
  sleep 1
  echo "::endgroup::"
  if [ "$CHECK_FAILS" -eq 0 ]; then result=pass; else result=FAIL; FAILED+=("$var"); fi
  RES_VAR+=("$var") RES_ROW+=("$k (line ${R_LINE[$i]})") RES_RESULT+=("$result") RES_TIME+=("$(( $(date +%s) - start ))s")
  echo "${name}: ${result}"
  ROW_DIR=""
}

KEEP_IMAGES=" ${SAMTOOLS_IMAGE:-} ${BCFTOOLS_IMAGE:-} ${MINIMAP2_IMAGE:-} "
smoke() {
  local full=$1; shift
  local var i k ran img
  mkdir -p "$FX" "$IN" "$LOGS" "$ROWS_DIR"
  ENV_FILE="${WORK}/versions.envfile"
  # shellcheck disable=SC2016  # expanded by the inner bash
  env -i bash --noprofile --norc -c '. "$1"; for n in $(grep -oE "^[A-Z][A-Z0-9_]*=" "$1" | tr -d =); do printf "%s=%s\n" "$n" "${!n}"; done' \
    _ "${REPO}/versions.env" > "$ENV_FILE"
  FAILED=()
  declare -A seen=()
  for var in "$@"; do
    [ -n "${seen[$var]:-}" ] && continue
    seen[$var]=1
    k=0 ran=0
    for i in "${!R_VAR[@]}"; do
      [ "${R_VAR[$i]}" = "$var" ] || continue
      k=$((k + 1))
      if has_opt "${R_OPTS[$i]}" full && [ "$full" != true ]; then
        echo "::warning::${var} row ${k} (commands.tsv line ${R_LINE[$i]}) runs only with --full: start the Container Test workflow by hand on this branch (Run workflow, images: ${var}) before merging a change to it."
        RES_VAR+=("$var") RES_ROW+=("$k (line ${R_LINE[$i]})") RES_RESULT+=("not run (full only)") RES_TIME+=("-")
        continue
      fi
      run_row "$i" "$k"
      ran=$((ran + 1))
    done
    if [ "$k" -eq 0 ]; then
      echo "ERROR: ${var} has no row in tests/smoke/commands.tsv" >&2
      RES_VAR+=("$var") RES_ROW+=("-") RES_RESULT+=("FAIL (no row)") RES_TIME+=("-")
      FAILED+=("$var")
    fi
    img=${!var:-}
    if [ -n "$img" ] && [[ "$KEEP_IMAGES" != *" ${img} "* ]]; then docker rmi -f "$img" >/dev/null 2>&1 || true; fi
    df -h "$WORK" | awk 'NR == 2 {print "disk: " $4 " free"}'
  done
  {
    echo "### Images on fixture ${TAG}: ${#FAILED[@]} failed"
    echo
    echo "| Variable | Row | Result | Time |"
    echo "|---|---|---|---|"
    for i in "${!RES_VAR[@]}"; do
      echo "| ${RES_VAR[$i]} | ${RES_ROW[$i]} | ${RES_RESULT[$i]} | ${RES_TIME[$i]} |"
    done
    echo
    echo "Row logs are in the image-smoke-logs-* artifacts of this run."
    for f in "${LOGS}"/*.log; do
      [ -f "$f" ] || continue
      grep -q '^\[FAIL\]' "$f" || continue
      echo
      echo "#### $(basename "$f" .log)"
      echo '```'
      grep -E '^\[FAIL\]|^ERROR' "$f" | head -n 20
      echo '```'
    done
  } | tee -a "$SUMMARY"
  printf '%s\n' "${FAILED[@]}" | sort -u | grep . > "${WORK}/failed.txt" || true
  [ "${#FAILED[@]}" -eq 0 ]
}

# ------------------------------------------------------------------ main
case "${1:-}" in
  --manifest)
    manifest
    ;;
  --manifest-self-test)
    manifest_self_test
    ;;
  --check)
    load_table
    check_coverage
    if [ "${#TABLE_ERRORS[@]}" -gt 0 ]; then printf 'ERROR: %s\n' "${TABLE_ERRORS[@]}" >&2; exit 1; fi
    echo "OK: ${#R_VAR[@]} rows in tests/smoke/commands.tsv parse."
    ;;
  --self-test)
    self_test
    ;;
  ""|-h|--help)
    sed -n '2,/^set -uo/p' "$0" | sed '$d'
    exit 2
    ;;
  *)
    FULL=false
    if [ "$1" = "--full" ]; then FULL=true; shift; fi
    load_table
    if [ "${#TABLE_ERRORS[@]}" -gt 0 ]; then printf 'ERROR: %s\n' "${TABLE_ERRORS[@]}" >&2; exit 1; fi
    [ $# -gt 0 ] || { echo "ERROR: name at least one *_IMAGE variable" >&2; exit 2; }
    smoke "$FULL" "$@"
    ;;
esac
