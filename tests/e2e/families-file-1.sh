#!/usr/bin/env bash
# The Nextflow pipeline given a families file (--pedigree with a leading
# family column, tests/fixtures/families/stub_families.tsv), in stub runs of
# the stub samplesheet (no containers):
#   - two families without --family_name stop the run before any task and
#     name both families;
#   - --family_name family_b runs family_b and publishes it under
#     <outdir>/family_b/: members.tsv, the family's own pedigree.tsv,
#     family_plan.tsv and the couple tables, and nothing under <outdir>/family/;
#   - red control: the same run on a copy whose modules publish to the
#     literal family/ folder leaves no family_plan.tsv under family_b/;
#   - the stub pedigree without the column still publishes under
#     <outdir>/family/ (couple-carrier-1.sh checks that on real tools).
# Cost: about 5 minutes, three stub runs and a failed start.
. "$(dirname "$0")/lib.sh"

FAMS="${REPO}/tests/fixtures/families/stub_families.tsv"
# The profile turns KIR typing on, which needs hla_typing: a params file turns it off.
echo '{"kir": false}' > "${CASE_TMP}/no-kir.json"
stub() {  # stub TREE NAME [nextflow options]: a stub run of TREE's main.nf, output in CASE_TMP/NAME
  local tree=$1 name=$2; shift 2
  rm -rf "${CASE_TMP:?}/${name}"
  ( cd "$tree" && nextflow run main.nf -profile test_all -stub -ansi-log false --outdir "${CASE_TMP}/${name}/out" \
      -work-dir "${CASE_TMP}/${name}/work" --tools clinvar,html_report -params-file "${CASE_TMP}/no-kir.json" "$@" ) \
    > "${CASE_TMP}/${name}.log" 2>&1
}
tasks() { find "$1" -name .command.run -exec sed -n "s/^### name: [\"']\(.*\)[\"']\$/\1/p" {} + | grep -c "$2" || true; }

# --- 1. two families, no --family_name: the run stops and lists them --------------------
stub "$REPO" both --pedigree "$FAMS"
RC=$?
tail -n 5 "${CASE_TMP}/both.log"
check "two families without --family_name: the run fails (exit ${RC})" test "$RC" -ne 0
check "and names both families" grep -q 'holds 2 families (family_a, family_b): name the one to run with --family_name' "${CASE_TMP}/both.log"
check_eq "before any family task" "$(tasks "${CASE_TMP}/both/work" FAMILY_PREPARE)" 0

# --- 2. --family_name family_b: published under family_b/ -----------------------------------
stub "$REPO" famb --pedigree "$FAMS" --family_name family_b
RC=$?
tail -n 5 "${CASE_TMP}/famb.log"
O="${CASE_TMP}/famb/out"
check_eq "--family_name family_b: exit" "$RC" 0
check_eq "one FAMILY_PREPARE task" "$(tasks "${CASE_TMP}/famb/work" FAMILY_PREPARE)" 1
for f in family_plan.tsv members.tsv pedigree.tsv family.vcf.gz relationships.tsv; do
  check "family_b/${f} is published" test -e "${O}/family_b/${f}"
done
check "the couple table of the two partners is under family_b/" sh -c "ls '${O}/family_b/' | grep -q '^couple_stub_bam_stub_fastq\.tsv$'"
check "each partner's copy stays in its own folder" test -e "${O}/stub_bam/couple/stub_bam_couple_carrier.tsv"
check "nothing is published under family/" test ! -e "${O}/family"

# --- 3. red control: modules that publish to the literal family/ ----------------------------
RED="${CASE_TMP}/red-tree"
rm -rf "$RED"; mkdir -p "$RED"
git -C "$REPO" ls-files -z | (cd "$REPO" && tar --null -T - -cf -) | tar -xf - -C "$RED"
# shellcheck disable=SC2016  # the Groovy expression, literally
sed -i 's/\${params.family_name ?: '"'"'family'"'"'}/family/g' "$RED"/modules/local/*/main.nf
check "red control: no publishDir of the copy's modules names --family_name" sh -c "! grep -q 'params.family_name ?:' '${RED}'/modules/local/*/main.nf"
stub "$RED" red --pedigree "$FAMS" --family_name family_b
check "red: without the change family_plan.tsv is not under family_b/" test ! -e "${CASE_TMP}/red/out/family_b/family_plan.tsv"
check "red: it went to family/ instead" test -e "${CASE_TMP}/red/out/family/family_plan.tsv"
rm -rf "$RED" "${CASE_TMP}/red"

# --- 4. a pedigree without the column: family/ as before ------------------------------------
stub "$REPO" plain --pedigree "${REPO}/assets/stub/pedigree.tsv"
RC=$?
check_eq "the stub pedigree without the family column: exit" "$RC" 0
check "it publishes under family/" test -e "${CASE_TMP}/plain/out/family/family_plan.tsv" -a -e "${CASE_TMP}/plain/out/family/members.tsv"
check "and has no family_b/ folder" test ! -e "${CASE_TMP}/plain/out/family_b"
rm -rf "${CASE_TMP}/both" "${CASE_TMP}/famb" "${CASE_TMP}/plain"
finish
