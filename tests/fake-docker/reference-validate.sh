#!/usr/bin/env bash
# validate-setup.sh and the reference a run uses:
#   - a reference .fai with ALT or HLA contigs fails, unless
#     ALLOW_ALT_REFERENCE=true turns the failure into a warning;
#   - a BAM whose @SQ lines are not the .fai's (another length, an extra
#     contig, another order) fails and names the first difference;
#   - a BAM sorted by name, one with no sort order in its header and no .bai
#     or a .bai older than itself,
#     one samtools quickcheck rejects, or one whose header samtools cannot
#     read, fails (quickcheck runs in every case);
#   - the default reference with a BAM aligned to it passes every one.
# The fake docker answers `samtools view -H` with HEADER and `samtools
# quickcheck` with QUICKCHECK_RC through a run hook.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
FAI="${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta.fai"
[ -f "$FAI" ] || fail "seed_reference wrote no ${FAI}: the default reference name and the seed disagree"

# The first sequences of the no-ALT analysis set, and of a reference that
# also has ALT and HLA contigs.
NOALT=$'chr1\t248956422\t112\t70\t71\nchr2\t242193529\t252513167\t70\t71\nchrM\t16569\t3099750718\t70\t71'
WITHALT="${NOALT}"$'\nchr22_KI270879v1_alt\t304135\t3100000000\t70\t71\nHLA-A*01:01:01:01\t3503\t3100400000\t70\t71'

# header FAI_TEXT [SO]: a BAM header with one @SQ line per .fai line.
header() {
  printf '@HD\tVN:1.6\tSO:%s\n' "${2:-coordinate}"
  awk -F'\t' '{printf "@SQ\tSN:%s\tLN:%s\n", $1, $2}' <<< "$1"
  printf '@RG\tID:sample1\tSM:sample1\n'
}

cat > "${CASE_WORK}/bam-hook" <<'HOOK'
#!/usr/bin/env bash
shift   # the image
case "$*" in
  *"samtools view -H"*)
    # view.rc, when present, makes the header unreadable.
    if [ -f "${CASE_WORK}/view.rc" ]; then exit "$(cat "${CASE_WORK}/view.rc")"; fi
    cat "${CASE_WORK}/header.sam" ;;
  *"samtools quickcheck"*) exit "$(cat "${CASE_WORK}/quickcheck.rc")" ;;
esac
exit 0
HOOK
chmod +x "${CASE_WORK}/bam-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/bam-hook"

# scenario FAI_TEXT HEADER_TEXT QUICKCHECK_RC: lay out one case.
scenario() {
  printf '%s\n' "$1" > "$FAI"
  printf '%s\n' "$2" > "${CASE_WORK}/header.sam"
  echo "$3" > "${CASE_WORK}/quickcheck.rc"
}

# --- the default reference and a BAM aligned to it ------------------------------
scenario "$NOALT" "$(header "$NOALT")" 0
run_expect 0 clean "${SCRIPTS}/validate-setup.sh" sample1
output_has clean '\[OK\].*Reference has no ALT or HLA contigs \(3 sequences\)'
output_has clean '\[OK\].*BAM header matches the reference: the same 3 sequences in the same order'
output_has clean '\[OK\].*BAM is coordinate-sorted'
output_has clean '\[OK\].*BAM passes samtools quickcheck'
output_lacks clean '^ *\[FAIL\]'

# --- a reference with ALT and HLA contigs ---------------------------------------
scenario "$WITHALT" "$(header "$WITHALT")" 0
run_expect 1 with-alt "${SCRIPTS}/validate-setup.sh" sample1
output_has with-alt '\[FAIL\].*Reference has 2 ALT/HLA contigs \(first: chr22_KI270879v1_alt\)'
output_has with-alt 'ALLOW_ALT_REFERENCE=true'
output_has with-alt 'docs/realignment\.md'
[ "$(grep -c '^ *\[FAIL\]' "${CASE_WORK}/with-alt.out")" -eq 1 ] || fail "with-alt: expected exactly one [FAIL], the ALT check"

ALLOW_ALT_REFERENCE=true run_expect 0 with-alt-allowed "${SCRIPTS}/validate-setup.sh" sample1
output_has with-alt-allowed '\[WARN\].*Reference has 2 ALT/HLA contigs \(allowed by ALLOW_ALT_REFERENCE=true\)'
output_lacks with-alt-allowed '^ *\[FAIL\]'

# --- a BAM aligned to another reference --------------------------------------------
# Another length for chr2: the first difference is sequence 2.
OTHER_LEN=$'chr1\t248956422\t112\t70\t71\nchr2\t111\t252513167\t70\t71\nchrM\t16569\t3099750718\t70\t71'
scenario "$NOALT" "$(header "$OTHER_LEN")" 0
run_expect 1 other-length "${SCRIPTS}/validate-setup.sh" sample1
output_has other-length '\[FAIL\].*this BAM was aligned to a different reference: realign \(docs/realignment\.md\)'
output_has other-length 'First difference: sequence 2 is chr2 \(111 bp\) in the BAM and chr2 \(242193529 bp\) in the reference\.'

# The BAM of the with-ALT reference against the no-ALT one: two more sequences.
scenario "$NOALT" "$(header "$WITHALT")" 0
run_expect 1 extra-contigs "${SCRIPTS}/validate-setup.sh" sample1
output_has extra-contigs 'First difference: sequence 4 is chr22_KI270879v1_alt \(304135 bp\) in the BAM and missing in the reference\.'

# The same sequences in another order.
SWAPPED=$'chr2\t242193529\t112\t70\t71\nchr1\t248956422\t252513167\t70\t71\nchrM\t16569\t3099750718\t70\t71'
scenario "$NOALT" "$(header "$SWAPPED")" 0
run_expect 1 swapped "${SCRIPTS}/validate-setup.sh" sample1
output_has swapped 'First difference: sequence 1 is chr2 \(242193529 bp\) in the BAM and chr1 \(248956422 bp\) in the reference\.'

# --- a BAM sorted by name, and one samtools quickcheck rejects ------------------------
scenario "$NOALT" "$(header "$NOALT" queryname)" 0
run_expect 1 by-name "${SCRIPTS}/validate-setup.sh" sample1
output_has by-name '\[FAIL\].*BAM is sorted by queryname, not by coordinate'

scenario "$NOALT" "$(header "$NOALT")" 1
run_expect 1 quickcheck "${SCRIPTS}/validate-setup.sh" sample1
output_has quickcheck '\[FAIL\].*BAM fails samtools quickcheck'

# No @HD SO: tag: a .bai not older than the BAM shows it is sorted (samtools
# index refuses an unsorted one); an older .bai or none is a failure.
scenario "$NOALT" "$(header "$NOALT" | grep -v '^@HD')" 0
run_expect 0 no-so "${SCRIPTS}/validate-setup.sh" sample1
output_has no-so '\[OK\].*BAM is coordinate-sorted \(no @HD SO: tag, but samtools only indexes a sorted BAM and its \.bai is not older than it\)'
# An index older than the BAM may belong to a BAM it replaced.
touch -t 202001010000 "${GENOME_DIR}/sample1/aligned/sample1_sorted.bam.bai"
run_expect 1 no-so-old-bai "${SCRIPTS}/validate-setup.sh" sample1
output_has no-so-old-bai '\[FAIL\].*no @HD SO: tag\) and its \.bai is older than the BAM'
touch "${GENOME_DIR}/sample1/aligned/sample1_sorted.bam.bai"
mv "${GENOME_DIR}/sample1/aligned/sample1_sorted.bam.bai" "${CASE_WORK}/bai.aside"
run_expect 1 no-so-no-bai "${SCRIPTS}/validate-setup.sh" sample1
output_has no-so-no-bai '\[FAIL\].*no @HD SO: tag\) and it has no \.bai'
mv "${CASE_WORK}/bai.aside" "${GENOME_DIR}/sample1/aligned/sample1_sorted.bam.bai"

# A header samtools cannot read: a failure, and quickcheck still runs.
scenario "$NOALT" "$(header "$NOALT")" 1
echo 1 > "${CASE_WORK}/view.rc"
run_expect 1 no-header "${SCRIPTS}/validate-setup.sh" sample1
output_has no-header '\[FAIL\].*Could not read the BAM header'
output_has no-header '\[FAIL\].*BAM fails samtools quickcheck'
rm -f "${CASE_WORK}/view.rc"

echo "reference-validate: ALT contigs, a BAM from another reference, a name-sorted BAM and a broken BAM fail; the default reference passes"
