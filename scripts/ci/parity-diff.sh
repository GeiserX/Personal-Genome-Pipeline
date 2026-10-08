#!/usr/bin/env bash
# parity-diff.sh: compare what the bash single-step scripts and the Nextflow
# pipeline wrote for one sample, item by item, so the two cannot drift apart
# without a red check.
#
# Items, each read the same way on both sides:
#   alignment     samtools flagstat of the BAM (reads, mapped, duplicates)
#   variants      bcftools isec of the two VCFs: records only in one side, and
#                 shared records whose FILTER or GT differ
#   gvcf          every gVCF record (CHROM POS REF ALT FILTER END GT)
#   clinvar       the rows of <sample>_clinvar_hits.tsv
#   roh           the RG segments of bcftools roh, and the >= 5 Mb summary
#   pharmcat      the diplotype PharmCAT reports for each gene
#   prs           each score's sum, matched count and input (gvcf or vcf)
#   sv_merged     each SV consensus record (CHROM POS END SVTYPE SVLEN SUPP
#                 SUPP_VEC and the callers' GT)
#   y_haplogroup  Yleaf's prediction table
#
# PARITY_ITEMS (space separated) picks the items; the default is the first
# seven, which the E2E workflow compares on the run from FASTQ (docs/nextflow.md
# has their table). It compares sv_merged and y_haplogroup on the Nextflow
# hardening run (docs/testing.md).
#
# A difference listed in KNOWN below (or in PARITY_KNOWN, same format) is
# reported with its reason and does not fail the check; a listed difference
# that is no longer seen fails it, so the list stays true.
#
# Usage:
#   scripts/ci/parity-diff.sh BASH_DIR NF_DIR SAMPLE
#       BASH_DIR  the sample directory the scripts wrote ($GENOME_DIR/<sample>)
#       NF_DIR    the pipeline's <outdir>/<sample>
#   scripts/ci/parity-diff.sh --self-test BASH_DIR SAMPLE
#       copies the bash outputs into the Nextflow layout, requires that copy
#       to compare equal, then plants one difference per item and requires
#       each to be reported (and a known one to pass, a stale one to fail,
#       a missing file to fail), on the items PARITY_ITEMS picks.
# Needs docker (bcftools and samtools from versions.env) and python3. Writes
# a markdown table to stdout and to $GITHUB_STEP_SUMMARY when it is set.
set -euo pipefail
export LC_ALL=C

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
# shellcheck source=../../versions.env
. "${ROOT}/versions.env"

# Differences that stay, one per line: ITEM<TAB>reason.
KNOWN=''

ALL_ITEMS="alignment variants gvcf clinvar roh pharmcat prs sv_merged y_haplogroup"
ITEMS=${PARITY_ITEMS:-alignment variants gvcf clinvar roh pharmcat prs}
for item in $ITEMS; do
  case " ${ALL_ITEMS} " in
    *" ${item} "*) ;;
    *) echo "parity-diff: unknown item '${item}' in PARITY_ITEMS (known: ${ALL_ITEMS})" >&2; exit 2 ;;
  esac
done

# path SIDE ITEM SAMPLE: where a side keeps an item's file(s), relative to its
# sample directory (roh lists two files).
path() {
  local s=$3
  case "$1:$2" in
    *:alignment) echo "aligned/${s}_sorted.bam" ;;
    *:variants)  echo "vcf/${s}.vcf.gz" ;;
    *:gvcf)      echo "vcf/${s}.g.vcf.gz" ;;
    *:clinvar)   echo "clinvar/${s}_clinvar_hits.tsv" ;;
    bash:roh)    echo "vcf/${s}_roh.txt vcf/${s}_roh_summary.txt" ;;
    nf:roh)      echo "roh/${s}_roh.txt roh/${s}_roh_summary.txt" ;;
    bash:pharmcat) echo "vcf/${s}.report.json" ;;
    nf:pharmcat)   echo "pharmcat/${s}.report.json" ;;
    *:prs)       echo "prs/${s}_prs_summary.tsv" ;;
    *:sv_merged) echo "sv_merged/${s}_sv_consensus.vcf.gz" ;;
    *:y_haplogroup) echo "y_haplogroup/${s}_y_haplogroup.txt" ;;
  esac
}

# in_image IMAGE A B CMD...: run CMD in IMAGE with side A at /a and B at /b.
in_image() {
  local image=$1 a=$2 b=$3
  shift 3
  docker run --rm -i -u "$(id -u):$(id -g)" -v "${a}:/a:ro" -v "${b}:/b:ro" "$image" "$@"
}

# canon ITEM DIR SIDE SAMPLE: the item as plain sorted text, for diff.
canon() {
  local item=$1 dir=$2 side=$3 s=$4 f
  f="${dir}/$(path "$side" "$item" "$s" | awk '{print $1}')"
  case "$item" in
    alignment)
      docker run --rm -u "$(id -u):$(id -g)" -v "$(dirname "$f"):/d:ro" "$SAMTOOLS_IMAGE" \
        samtools flagstat "/d/$(basename "$f")" ;;
    gvcf)
      docker run --rm -u "$(id -u):$(id -g)" -v "$(dirname "$f"):/d:ro" "$BCFTOOLS_IMAGE" \
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/END\t[%GT]\n' "/d/$(basename "$f")" ;;
    clinvar)
      awk 'NR > 1' "$f" | sort ;;
    roh)
      { { grep '^RG' "$f" || true; } | cut -f2- | sort
        echo "summary:"
        cat "${dir}/$(path "$side" roh "$s" | awk '{print $2}')"; } ;;
    pharmcat)
      python3 - "$f" <<'PY'
import json, sys
data = json.load(open(sys.argv[1]))
rows = set()
def dips(gene, g):
    for d in (g.get("sourceDiplotypes") or g.get("recommendationDiplotypes") or []):
        names = [(d.get(a) or {}).get("name", "") or "" for a in ("allele1", "allele2")]
        rows.add(f"{gene}\t{'/'.join(names)}")
for key, val in (data.get("genes") or {}).items():
    if not isinstance(val, dict):
        continue
    if "sourceDiplotypes" in val or "recommendationDiplotypes" in val:
        dips(key, val)
    else:
        for name, g in val.items():
            if isinstance(g, dict):
                dips(name, g)
print("\n".join(sorted(rows)))
PY
      ;;
    prs)
      awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
        {print $c["PGS_ID"] "\t" $c["Score_SUM"] "\t" $c["Variants_Matched"] "\t" (("Input" in c) ? $c["Input"] : "none")}' "$f" | sort ;;
    sv_merged)
      docker run --rm -u "$(id -u):$(id -g)" -v "$(dirname "$f"):/d:ro" "$BCFTOOLS_IMAGE" \
        bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/SUPP\t%INFO/SUPP_VEC[\t%GT]\n' \
        "/d/$(basename "$f")" | sort ;;
    y_haplogroup)
      cat "$f" ;;
  esac
}

# variants A B SAMPLE: "only_a only_b shared mismatched" from bcftools isec.
variants() {
  local fa fb
  fa=$(path bash variants "$3"); fb=$(path nf variants "$3")
  # shellcheck disable=SC2016  # $1 and $2 belong to the inner bash
  in_image "$BCFTOOLS_IMAGE" "$1" "$2" bash -euo pipefail -c '
    d=$(mktemp -d)
    bcftools isec -p "$d" "/a/$1" "/b/$2" >/dev/null
    n() { grep -vc "^#" "$1" || true; }
    q() { bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT]\n" "$1"; }
    # Each side sorted, so records at one position pair by key whatever
    # order each VCF holds them in; a pair whose keys differ counts too.
    mis=$(paste <(q "$d/0002.vcf" | LC_ALL=C sort) <(q "$d/0003.vcf" | LC_ALL=C sort) \
      | awk -F"\t" "\$1 != \$7 || \$2 != \$8 || \$3 != \$9 || \$4 != \$10 || \$5 != \$11 || \$6 != \$12" | wc -l)
    echo "$(n "$d/0000.vcf") $(n "$d/0001.vcf") $(n "$d/0002.vcf") $mis"' _ "$fa" "$fb"
}

# compare BASH_DIR NF_DIR SAMPLE: print the table, return 1 on a difference
# that is not listed, or a listed one that is not seen.
compare() {
  local a=$1 b=$2 s=$3 item rc=0 bad result detail reason known tmp f side dir
  known=$(printf '%s\n%s\n' "$KNOWN" "${PARITY_KNOWN:-}" | grep -v '^[[:space:]]*$' || true)
  tmp=$(mktemp -d)
  echo "### Parity: bash scripts against the Nextflow pipeline (${s})"
  echo
  echo "| Item | Result | Detail |"
  echo "|---|---|---|"
  for item in $ITEMS; do
    bad=""
    for side in bash nf; do
      dir=$a; [ "$side" = nf ] && dir=$b
      for f in $(path "$side" "$item" "$s"); do
        [ -s "${dir}/${f}" ] || bad="${bad} ${side}:${f}"
      done
    done
    reason=$(awk -F'\t' -v i="$item" '$1 == i {print $2; exit}' <<<"$known")
    if [ -n "$bad" ]; then
      result=MISSING detail="no file:${bad}"
    elif [ "$item" = variants ]; then
      read -r oa ob sh mm < <(variants "$a" "$b" "$s" 2>"${tmp}/err" || true) || true
      if ! [[ "${oa:-}${ob:-}${sh:-}${mm:-}" =~ ^[0-9]+$ ]]; then
        result=ERROR detail="bcftools isec failed: $(head -c 200 "${tmp}/err" | tr '\n' ' ')"
      else
        detail="bash only ${oa}, Nextflow only ${ob}, shared ${sh}, shared with another FILTER or GT ${mm}"
        if [ "$oa" -eq 0 ] && [ "$ob" -eq 0 ] && [ "$mm" -eq 0 ] && [ "$sh" -gt 0 ]; then result=equal; else result=differs; fi
      fi
    else
      canon "$item" "$a" bash "$s" > "${tmp}/a" 2>"${tmp}/err" || { result=ERROR; detail=$(head -c 200 "${tmp}/err" | tr '\n' ' '); }
      canon "$item" "$b" nf "$s" > "${tmp}/b" 2>>"${tmp}/err" || { result=ERROR; detail=$(head -c 200 "${tmp}/err" | tr '\n' ' '); }
      if [ "${result:-}" != ERROR ]; then
        if cmp -s "${tmp}/a" "${tmp}/b"; then
          result=equal
          detail="$(grep -c . "${tmp}/a" || true) lines on each side"
        else
          result=differs
          detail="bash $(grep -c . "${tmp}/a" || true) lines, Nextflow $(grep -c . "${tmp}/b" || true); $(diff "${tmp}/a" "${tmp}/b" | grep -c '^[<>]' || true) lines differ"
          diff "${tmp}/a" "${tmp}/b" | head -n 20 | sed "s/^/  ${item}: /" >&2 || true
        fi
      fi
    fi
    case "$result" in
      equal)
        if [ -n "$reason" ]; then result="FAIL: listed as a known difference but equal"; rc=1; fi ;;
      differs)
        if [ -n "$reason" ]; then result="differs (known: ${reason})"; else result="FAIL: differs"; rc=1; fi ;;
      *) result="FAIL: ${result}"; rc=1 ;;
    esac
    echo "| ${item} | ${result} | ${detail} |"
    result=""
  done
  rm -rf "$tmp"
  return "$rc"
}

# self_test BASH_DIR SAMPLE
self_test() {
  local a=$1 s=$2 t fail=0 item out rc
  t=$(mktemp -d)
  # shellcheck disable=SC2064
  trap "rm -rf '$t'" EXIT
  # fresh: the bash outputs under the Nextflow layout, in $t/nf
  fresh() {
    rm -rf "${t}/nf"
    local item i bf nf
    for item in $ITEMS; do
      read -r -a bf <<<"$(path bash "$item" "$s")"
      read -r -a nf <<<"$(path nf "$item" "$s")"
      for i in "${!bf[@]}"; do
        mkdir -p "$(dirname "${t}/nf/${nf[$i]}")"
        cp "${a}/${bf[$i]}" "${t}/nf/${nf[$i]}"
        for x in .tbi .bai; do [ -f "${a}/${bf[$i]}${x}" ] && cp "${a}/${bf[$i]}${x}" "${t}/nf/${nf[$i]}${x}"; done
      done
    done
    return 0
  }
  # in_t CMD...: run CMD in BCFTOOLS_IMAGE (or SAMTOOLS_IMAGE with --sam) with $t/nf writable at /n
  in_t() {
    local image=$BCFTOOLS_IMAGE
    if [ "$1" = --sam ]; then image=$SAMTOOLS_IMAGE; shift; fi
    docker run --rm -u "$(id -u):$(id -g)" -v "${t}/nf:/n" "$image" bash -euo pipefail -c "$1"
  }
  # drop_first FILE: remove the first record of a bgzipped VCF in $t/nf
  drop_first() {
    # awk reads to the end: head would close the pipe and pipefail would
    # turn bcftools' SIGPIPE into a failed step.
    in_t "f=/n/$1; p=\$(bcftools query -f '%CHROM:%POS\n' \$f | awk 'NR == 1')
          bcftools view -t \"^\$p\" -Oz -o /n/x.vcf.gz \$f; mv /n/x.vcf.gz \$f; bcftools index -f -t \$f"
  }
  # expect NAME PATTERN [KNOWN]: compare must exit 1 (0 with KNOWN) and print PATTERN
  expect() {
    rc=0
    out=$(PARITY_KNOWN="${3:-}" compare "$a" "${t}/nf" "$s" 2>/dev/null) || rc=$?
    local want=1
    [ "${4:-}" = pass ] && want=0
    if [ "$rc" -ne "$want" ] || ! grep -qE -- "$2" <<<"$out"; then
      echo "self-test: '$1' exited ${rc} (want ${want}) and did not report /$2/:"; printf '%s\n' "$out"; fail=1
    else
      echo "self-test: '$1' caught: $(grep -E -- "$2" <<<"$out" | head -n 1)"
    fi
  }

  # The control: the same files on both sides compare equal.
  fresh
  rc=0; out=$(compare "$a" "${t}/nf" "$s" 2>&1) || rc=$?
  if [ "$rc" -ne 0 ] || [ "$(grep -cE '^\| [a-z_]+ \| equal \|' <<<"$out" || true)" -ne "$(wc -w <<<"$ITEMS")" ]; then
    echo "self-test: the bash outputs against a copy of themselves did not compare equal:"; printf '%s\n' "$out"; fail=1
  else
    echo "self-test: a copy compares equal on every item"
  fi

  # plant ITEM: one difference in the Nextflow copy of ITEM
  plant() {
    case "$1" in
      alignment)
        in_t --sam "samtools view -b -F 1024 -o /n/x.bam /n/$(path nf alignment "$s"); mv /n/x.bam /n/$(path nf alignment "$s"); samtools index /n/$(path nf alignment "$s")" ;;
      variants) drop_first "$(path nf variants "$s")" ;;
      variants-gt)
        in_t "f=/n/$(path nf variants "$s")
              bcftools view \$f | awk -F'\t' -v OFS='\t' '/^#/ {print; next} !d && \$10 ~ /^0\/1/ {sub(/^0\/1/, \"1/1\", \$10); d=1} {print}' \
                | bcftools view -Oz -o /n/x.vcf.gz
              mv /n/x.vcf.gz \$f; bcftools index -f -t \$f" ;;
      gvcf) drop_first "$(path nf gvcf "$s")" ;;
      clinvar) printf 'chr1\t1\tA\tT\t0/1\t1\tGENE:1\tPathogenic\tplanted\n' >> "${t}/nf/$(path nf clinvar "$s")" ;;
      roh) printf 'RG\t%s\tchr1\t1\t6000000\t6000000\t1\t99.0\n' "$s" >> "${t}/nf/$(path nf roh "$s" | awk '{print $1}')" ;;
      pharmcat)
        python3 - "${t}/nf/$(path nf pharmcat "$s")" <<'PY'
import json, sys
p = sys.argv[1]
d = json.load(open(p))
d["genes"] = {}
json.dump(d, open(p, "w"))
PY
        ;;
      prs)
        awk -F'\t' -v OFS='\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; print; next}
            NR == 2 {$c["Score_SUM"] = $c["Score_SUM"] + 1} {print}' "${a}/$(path bash prs "$s")" > "${t}/nf/$(path nf prs "$s")" ;;
      sv_merged) drop_first "$(path nf sv_merged "$s")" ;;
      y_haplogroup)
        awk -F'\t' -v OFS='\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; print; next}
            NR == 2 {$c["Hg"] = "PLANTED"} {print}' "${a}/$(path bash y_haplogroup "$s")" > "${t}/nf/$(path nf y_haplogroup "$s")" ;;
    esac
  }
  # The name of each item's planted difference.
  declare -A what=(
    [alignment]="duplicates dropped from the BAM" [variants]="a VCF record dropped" [gvcf]="a gVCF record dropped"
    [clinvar]="a ClinVar hit added" [roh]="a ROH segment added" [pharmcat]="PharmCAT's genes emptied"
    [prs]="a PRS sum changed" [sv_merged]="an SV consensus record dropped" [y_haplogroup]="the Y haplogroup changed"
  )
  for item in $ITEMS; do
    fresh; plant "$item"
    if [ "$item" = variants ]; then
      expect "${what[$item]}" '^\| variants \| FAIL: differs \| bash only 1, '
      fresh; plant variants-gt
      expect "a genotype changed" '^\| variants \| FAIL: differs \| bash only 0, Nextflow only 0, shared [0-9]+, shared with another FILTER or GT [1-9]'
    else
      expect "${what[$item]}" "^\\| ${item} \\| FAIL: differs"
    fi
  done
  # The list of known differences, and a missing file, on the last item.
  local k=${ITEMS##* }
  fresh; plant "$k"
  expect "the same difference, listed as known" "^\\| ${k} \\| differs \\(known: planted\\)" "${k}"$'\tplanted' pass
  fresh
  expect "a known difference that is not seen" "^\\| ${k} \\| FAIL: listed as a known difference but equal" "${k}"$'\tplanted'
  fresh; rm -f "${t}/nf/$(path nf "$k" "$s" | awk '{print $1}')"
  expect "a missing file" "^\\| ${k} \\| FAIL: MISSING"

  [ "$fail" -eq 0 ] && echo "self-test: OK"
  return "$fail"
}

case "${1:-}" in
  --self-test)
    [ $# -eq 3 ] || { echo "usage: $0 --self-test BASH_DIR SAMPLE" >&2; exit 2; }
    self_test "$2" "$3" ;;
  -*|'') echo "usage: $0 BASH_DIR NF_DIR SAMPLE | --self-test BASH_DIR SAMPLE" >&2; exit 2 ;;
  *)
    [ $# -eq 3 ] || { echo "usage: $0 BASH_DIR NF_DIR SAMPLE" >&2; exit 2; }
    rc=0
    out=$(compare "$1" "$2" "$3") || rc=$?
    printf '%s\n' "$out"
    [ -n "${GITHUB_STEP_SUMMARY:-}" ] && printf '%s\n\n' "$out" >> "$GITHUB_STEP_SUMMARY"
    exit "$rc" ;;
esac
