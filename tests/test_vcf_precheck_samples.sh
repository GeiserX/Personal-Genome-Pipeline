#!/usr/bin/env bash
# One samplesheet row is one sample: a VCF with more than one sample column
# is refused at intake, with the count and the names, and a single-sample VCF
# goes on.
#
# Two places apply the rule, and both are run here on the synthetic VCFs in
# tests/fixtures/vcf/:
#   - VCF_PRECHECK (modules/local/vcf_precheck/main.nf): its script: block is
#     rendered as Nextflow would (Groovy escapes, the few ${...} it uses) and
#     run under bash -euo pipefail, with a fake bcftools on PATH that answers
#     the calls the block makes from the plain-text fixture;
#   - vcf_header_samples in scripts/lib/common.sh, which validate-setup.sh
#     uses (its own wiring is the validate-multisample-vcf fake-docker case).
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
FIX="${REPO}/tests/fixtures/vcf"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

# --- the script: block of VCF_PRECHECK, as a bash file ---------------------------
python3 - "${REPO}/modules/local/vcf_precheck/main.nf" "${WORK}/precheck.sh" <<'PY'
import re, sys, textwrap
src, out = sys.argv[1], sys.argv[2]
text = open(src).read()
m = re.search(r'\n    script:\n(.*?)\n    """\n(.*?)\n    """\n\n    stub:', text, re.S)
if not m:
    sys.exit("could not find the script: block of VCF_PRECHECK")
body = m.group(2)
# The Groovy values the block interpolates. A new one fails here, so this
# test is updated instead of running a script with a hole in it.
values = {
    'vcf': 'input.vcf.gz',
    'vcf.name': 'input.vcf.gz',
    'meta.id': 'S1',
    'allow_unfiltered': 'false',
    'task.process': 'VCF_PRECHECK',
}
def interp(mo):
    expr = mo.group(1)
    if expr.startswith('task.container'):
        return 'test'
    if expr not in values:
        sys.exit("unknown Groovy value in the script: block: ${%s}" % expr)
    return values[expr]
# Groovy escapes in a triple-quoted GString, then ${...} that is not escaped.
escapes = {'\\': '\\', '$': '\0DOLLAR\0', 'n': '\n', 't': '\t', '"': '"', "'": "'"}
def unescape(mo):
    c = mo.group(1)
    if c not in escapes:
        sys.exit("unexpected Groovy escape \\%s in the script: block" % c)
    return escapes[c]
body = re.sub(r'\\(.)', unescape, body)
body = re.sub(r'\$\{([^}]*)\}', interp, body)
body = body.replace('\0DOLLAR\0', '$')
open(out, 'w').write(textwrap.dedent(body) + '\n')
PY

# --- a fake bcftools, from the plain-text VCF beside the input -----------------
mkdir -p "${WORK}/bin"
cat > "${WORK}/bin/bcftools" <<'FAKE'
#!/usr/bin/env bash
# Answers the calls VCF_PRECHECK makes. The input is gzip-compressed text.
case "$1 ${2:-}" in
  "index -s")
    gzip -dc "$3" | awk -F'\t' '!/^#/ { n[$1]++ } END { for (c in n) printf "%s\t.\t%d\n", c, n[c] }' ;;
  "view -h")
    gzip -dc "$3" | grep '^#' ;;
  "query -f")
    # %FILTER, or %FILTER\t%ALT\t%INFO/END: the reference-block columns read '.'
    case "$3" in
      *ALT*) gzip -dc "$4" | awk -F'\t' -v OFS='\t' '!/^#/ { print $7, $5, "." }' ;;
      *) gzip -dc "$4" | awk -F'\t' '!/^#/ { print $7 }' ;;
    esac ;;
  *) echo "fake bcftools: unexpected call: $*" >&2; exit 2 ;;
esac
FAKE
chmod +x "${WORK}/bin/bcftools"

# run_precheck VCF_TEXT_FILE: run the block in a fresh task directory; sets OUT, RC
run_precheck() {
  local dir
  dir=$(mktemp -d "${WORK}/task.XXXXXX")
  gzip -c "$1" > "${dir}/input.vcf.gz"
  : > "${dir}/input.vcf.gz.tbi"
  set +e
  OUT=$(cd "$dir" && PATH="${WORK}/bin:${PATH}" bash -euo pipefail "${WORK}/precheck.sh" 2>&1)
  RC=$?
  set -e
}

run_precheck "${FIX}/two_samples.vcf"
if [ "$RC" -ne 0 ] && grep -q "Sample 'S1': input.vcf.gz holds 2 samples (SAMPLE_A, SAMPLE_B)" <<<"$OUT" \
   && grep -q 'One samplesheet row is one sample' <<<"$OUT" \
   && grep -q 'bcftools view -s SAMPLE_A -a -c 1 -Oz -o SAMPLE_A.vcf.gz input.vcf.gz' <<<"$OUT"; then
  pass "VCF_PRECHECK refuses a two-sample VCF, naming the count, the samples and the fix"
else
  fail "VCF_PRECHECK on a two-sample VCF: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

awk -F'\t' -v OFS='\t' '/^#CHROM/ { $10 = "S1"; $11 = "S2\tS3\tS4\tS5\tS6\tS7" } { print }' \
  "${FIX}/two_samples.vcf" | grep '^#' > "${WORK}/seven.vcf"
printf 'chr1\t10000100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\n' >> "${WORK}/seven.vcf"
run_precheck "${WORK}/seven.vcf"
if [ "$RC" -ne 0 ] && grep -q 'holds 7 samples (S1, S2, S3, S4, S5, ...)' <<<"$OUT"; then
  pass "VCF_PRECHECK lists the first five of seven samples, then ..."
else
  fail "VCF_PRECHECK on a seven-sample VCF: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

run_precheck "${FIX}/one_sample.vcf"
if [ "$RC" -eq 0 ] && ! grep -q 'samples (' <<<"$OUT" && grep -q 'FILTER counts PASS=3' <<<"$OUT"; then
  pass "VCF_PRECHECK accepts a one-sample VCF and goes on to count FILTER values"
else
  fail "VCF_PRECHECK on a one-sample VCF: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# --- vcf_header_samples in scripts/lib/common.sh ------------------------------------
# shellcheck source=../scripts/lib/common.sh
. "${REPO}/scripts/lib/common.sh"
got=$(vcf_header_samples < "${FIX}/two_samples.vcf" | paste -sd, -)
if [ "$got" = "SAMPLE_A,SAMPLE_B" ]; then pass "vcf_header_samples: two samples"; else fail "vcf_header_samples two: '${got}'"; fi
got=$(vcf_header_samples < "${FIX}/one_sample.vcf" | paste -sd, -)
if [ "$got" = "SAMPLE_A" ]; then pass "vcf_header_samples: one sample"; else fail "vcf_header_samples one: '${got}'"; fi
got=$(cut -f1-8 "${FIX}/one_sample.vcf" | vcf_header_samples | paste -sd, -)
if [ -z "$got" ]; then pass "vcf_header_samples: a sites-only VCF has none"; else fail "vcf_header_samples sites-only: '${got}'"; fi

echo ""
if [ "$FAILS" -gt 0 ]; then
  echo "${FAILS} check(s) failed"
  exit 1
fi
echo "All checks passed"
