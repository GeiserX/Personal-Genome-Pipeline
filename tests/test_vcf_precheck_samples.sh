#!/usr/bin/env bash
# One samplesheet row is one sample: a VCF with more than one sample column
# is refused at intake, with the count and the names, and a single-sample VCF
# goes on. The same block refuses a VCF whose ##contig lengths are not
# GRCh38's (a GRCh37 file named the chr way or the Ensembl way), and lets a
# GRCh38 header through. A header without contig lengths gets its REF bases
# compared with the reference instead: more than 5% of the first 1000
# records different stops the run, 5% or less goes on.
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
    'reference': 'ref.fa',
    'reference_fai': 'ref.fa.fai',
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
  "view -H")
    gzip -dc "$3" | grep -v '^#' ;;
  "norm -c")
    # norm -c x -f REF -o OUT IN: keep the records whose REF is the reference's
    # bases there. Like the real one it fails on a record whose contig the
    # header does not declare, or the reference does not have.
    [ "$3 $4 $6" = "x -f -o" ] || { echo "fake bcftools: unexpected norm call: $*" >&2; exit 2; }
    awk -F'\t' -v OFS='\t' -v out="$7" '
      FNR == 1 { f++ }
      f == 1 { if (/^>/) { name = substr($1, 2); sub(/[ \t].*/, "", name) } else seq[name] = seq[name] $0; next }
      /^##contig=<ID=/ { id = $0; sub(/^##contig=<ID=/, "", id); sub(/[,>].*/, "", id); declared[id] = 1 }
      /^#/ { print > out; next }
      !($1 in declared) { print "[E::vcf_format] Invalid BCF, CONTIG id=0 not present in the header" > "/dev/stderr"; exit 255 }
      !($1 in seq) { print "faidx_fetch_seq failed at " $1 ":" $2 > "/dev/stderr"; exit 255 }
      toupper(substr(seq[$1], $2, length($4))) == toupper($4) { print > out }' "$5" "$8" ;;
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

# --- a small reference -------------------------------------------------------
# chr1 and chrM, every base G except the REF bases of the fixture's three
# records: chr1:10000100 A, chr1:10000200 C, chrM:100 T.
python3 - "${WORK}/ref.fa" <<'PY2'
import sys
seqs = {'chr1': (10001500, {10000100: 'A', 10000200: 'C'}), 'chrM': (16569, {100: 'T'})}
with open(sys.argv[1], 'w') as fa, open(sys.argv[1] + '.fai', 'w') as fai:
    for name, (n, bases) in seqs.items():
        s = bytearray(b'G' * n)
        for pos, b in bases.items():
            s[pos - 1] = ord(b)
        fa.write('>%s\n' % name)
        fai.write('%s\t%d\t%d\t%d\t%d\n' % (name, n, fa.tell(), n, n + 1))
        fa.write(s.decode() + '\n')
PY2

# run_precheck VCF_TEXT_FILE: run the block in a fresh task directory; sets OUT, RC
run_precheck() {
  local dir
  dir=$(mktemp -d "${WORK}/task.XXXXXX")
  gzip -c "$1" > "${dir}/input.vcf.gz"
  : > "${dir}/input.vcf.gz.tbi"
  ln -s "${WORK}/ref.fa" "${dir}/ref.fa"
  ln -s "${WORK}/ref.fa.fai" "${dir}/ref.fa.fai"
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

# A sample name is free text: the printed command quotes it, so copying it
# cannot run what the name holds.
# shellcheck disable=SC2016  # the $(...) is the literal sample name under test
awk -F'\t' -v OFS='\t' '/^#CHROM/ { $10 = "x$(touch pwned) y" } { print }' \
  "${FIX}/two_samples.vcf" > "${WORK}/odd.vcf"
run_precheck "${WORK}/odd.vcf"
# shellcheck disable=SC2016
if [ "$RC" -ne 0 ] && grep -qF 'bcftools view -s x\$\(touch\ pwned\)\ y -a -c 1 -Oz -o x\$\(touch\ pwned\)\ y.vcf.gz input.vcf.gz' <<<"$OUT"; then
  pass "VCF_PRECHECK shell-quotes the sample name in the command it prints"
else
  fail "VCF_PRECHECK on an odd sample name: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

run_precheck "${FIX}/one_sample.vcf"
if [ "$RC" -eq 0 ] && ! grep -q 'samples (' <<<"$OUT" && grep -q 'FILTER counts PASS=3' <<<"$OUT"; then
  pass "VCF_PRECHECK accepts a one-sample VCF and goes on to count FILTER values"
else
  fail "VCF_PRECHECK on a one-sample VCF: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# --- build: ##contig lengths that are not GRCh38's --------------------------------------
# GRCh37's chr1 under the chr name: the stop names the contig, both lengths and GRCh37.
sed 's/^##contig=<ID=chr1,length=248956422>/##contig=<ID=chr1,length=249250621>/' \
  "${FIX}/one_sample.vcf" > "${WORK}/grch37.vcf"
grep -q 'length=249250621' "${WORK}/grch37.vcf" || fail "the GRCh37 test input was not built"
run_precheck "${WORK}/grch37.vcf"
if [ "$RC" -ne 0 ] && grep -q "Sample 'S1': input.vcf.gz is not on GRCh38" <<<"$OUT" \
   && grep -q 'chr1 length 249250621 (GRCh38: 248956422)' <<<"$OUT" \
   && grep -q 'A chr1 length of 249250621 is GRCh37 (hg19)' <<<"$OUT"; then
  pass "VCF_PRECHECK stops a chr-named VCF with GRCh37's chr1 length, naming the build"
else
  fail "VCF_PRECHECK on a GRCh37 header: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# The same under Ensembl names: the build stop comes before the rename advice.
awk -F'\t' -v OFS='\t' '/^##contig=<ID=chr1,/ { $0 = "##contig=<ID=1,length=249250621>" } $1 == "chr1" { $1 = "1" } { print }' \
  "${FIX}/one_sample.vcf" > "${WORK}/grch37-ensembl.vcf"
run_precheck "${WORK}/grch37-ensembl.vcf"
if [ "$RC" -ne 0 ] && grep -q '    1 length 249250621 (GRCh38: 248956422)' <<<"$OUT"; then
  pass "VCF_PRECHECK stops an Ensembl-named VCF with GRCh37's chr1 length"
else
  fail "VCF_PRECHECK on an Ensembl GRCh37 header: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# One wrong length on another contig is enough, and the GRCh37 line is left out.
awk '/^#CHROM/ { print "##contig=<ID=chr2,length=243199373,assembly=hg19>" } { print }' \
  "${FIX}/one_sample.vcf" > "${WORK}/chr2.vcf"
run_precheck "${WORK}/chr2.vcf"
if [ "$RC" -ne 0 ] && grep -q 'chr2 length 243199373 (GRCh38: 242193529)' <<<"$OUT" \
   && ! grep -q 'is GRCh37' <<<"$OUT"; then
  pass "VCF_PRECHECK stops on a wrong chr2 length and names it"
else
  fail "VCF_PRECHECK on a wrong chr2 length: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# Every GRCh38 length of chr1-22, X and Y, in any attribute order, goes on;
# so does a header without ##contig lines (validate-setup.sh spot-checks those).
{
  grep '^##fileformat' "${FIX}/one_sample.vcf"
  awk -F'[(", )]+' '/^    \("chr/ { for (i = 2; i < NF; i += 2) print $i, $(i + 1) }' \
    "${REPO}/tests/demo/make_demo_sample.py" \
    | awk '{ printf "##contig=<ID=%s,assembly=GRCh38,length=%s>\n", $1, $2 }'
  grep -v '^##fileformat' "${FIX}/one_sample.vcf" | grep -v '^##contig'
} > "${WORK}/grch38-all.vcf"
[ "$(grep -c '^##contig' "${WORK}/grch38-all.vcf")" -eq 24 ] || fail "the GRCh38 test header does not have 24 contigs"
run_precheck "${WORK}/grch38-all.vcf"
if [ "$RC" -eq 0 ] && grep -q 'FILTER counts PASS=3' <<<"$OUT" && ! grep -q 'REF bases' <<<"$OUT"; then
  pass "VCF_PRECHECK accepts GRCh38's 24 lengths (the demo sample's NCBI table), with no REF spot-check"
else
  fail "VCF_PRECHECK on a full GRCh38 header: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# --- build: REF bases, when the header gives no contig length ------------------------
# The fixture's records match the small reference: the file goes on.
grep -v '^##contig' "${FIX}/one_sample.vcf" > "${WORK}/nocontig.vcf"
run_precheck "${WORK}/nocontig.vcf"
if [ "$RC" -eq 0 ] && ! grep -q 'not on GRCh38' <<<"$OUT" \
   && grep -q 'REF bases differ from the reference in 0 of the first 3 records (0%)' <<<"$OUT"; then
  pass "VCF_PRECHECK spot-checks a header without ##contig lines and lets matching REF bases through"
else
  fail "VCF_PRECHECK on a header without ##contig lines: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# spot_vcf FILE N_SAME N_DIFF [CONTIG]: a header without contig lengths and
# N_SAME records whose REF is the reference's G, then N_DIFF whose REF is T.
spot_vcf() {
  local c=${4:-chr1}
  {
    grep -v '^##contig' "${FIX}/one_sample.vcf" | grep '^#'
    awk -v s="$2" -v d="$3" -v c="$c" 'BEGIN { OFS = "\t"
      for (i = 1; i <= s + d; i++) print c, 10000300 + i, ".", (i <= s ? "G" : "T"), "A", 50, "PASS", ".", "GT", "0/1" }'
  } > "$1"
}

# Every REF different: the stop names the sample, the count, the rate and the fix.
spot_vcf "${WORK}/spot_all.vcf" 0 20
run_precheck "${WORK}/spot_all.vcf"
if [ "$RC" -ne 0 ] && grep -q "Sample 'S1': input.vcf.gz is not on GRCh38. Its header gives no contig lengths" <<<"$OUT" \
   && grep -q 'they differ in 20 of the first 20 records (100%)' <<<"$OUT" \
   && grep -q 'start from the FASTQ' <<<"$OUT" && ! grep -q 'FILTER counts' <<<"$OUT"; then
  pass "VCF_PRECHECK stops a header without contig lengths whose REF bases are not the reference's"
else
  fail "VCF_PRECHECK on REF bases all different: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# The line is validate-setup.sh's: 1 of 20 (5%) goes on, 2 of 20 (10%) stops.
spot_vcf "${WORK}/spot_5.vcf" 19 1
run_precheck "${WORK}/spot_5.vcf"
if [ "$RC" -eq 0 ] && grep -q 'REF bases differ from the reference in 1 of the first 20 records (5%)' <<<"$OUT"; then
  pass "VCF_PRECHECK lets 5% different REF bases through"
else
  fail "VCF_PRECHECK on 1 of 20 REF bases different: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi
spot_vcf "${WORK}/spot_10.vcf" 18 2
run_precheck "${WORK}/spot_10.vcf"
if [ "$RC" -ne 0 ] && grep -q 'they differ in 2 of the first 20 records (10%)' <<<"$OUT"; then
  pass "VCF_PRECHECK stops at 10% different REF bases"
else
  fail "VCF_PRECHECK on 2 of 20 REF bases different: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# Only the first 1000 records count.
spot_vcf "${WORK}/spot_1000.vcf" 1000 200
run_precheck "${WORK}/spot_1000.vcf"
if [ "$RC" -eq 0 ] && grep -q 'differ from the reference in 0 of the first 1000 records' <<<"$OUT"; then
  pass "VCF_PRECHECK compares the first 1000 records only"
else
  fail "VCF_PRECHECK on 1200 records: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# Records on a contig the reference lacks are left out, not compared.
spot_vcf "${WORK}/spot_chr2.vcf" 0 20 chr2
run_precheck "${WORK}/spot_chr2.vcf"
if [ "$RC" -eq 0 ] && grep -q 'no record is on a reference contig: REF bases not compared' <<<"$OUT"; then
  pass "VCF_PRECHECK leaves out records on contigs the reference lacks"
else
  fail "VCF_PRECHECK on records off the reference: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
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
