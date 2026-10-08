#!/usr/bin/env bash
# Step 23's ClinVar tier follows the ClinVar file step 13 annotates with
# --custom (CSQ field ClinVar_CLNSIG), not the cache release's CLIN_SIG. The
# fixture's VEP output (built with --database, no --custom) gets the field
# added: record A is Pathogenic in ClinVar_CLNSIG and not in CLIN_SIG (newly
# classified), record B is pathogenic in CLIN_SIG only (reclassified since the
# cache release). On a copy of the sample with no step 06 hits, the tier must
# hold A and not B.
. "$(dirname "$0")/lib.sh"

T="${SAMPLE}cv"
rm -rf "${GENOME_DIR:?}/${T}"
mkdir -p "${GENOME_DIR}/${T}/vep"
SRC="${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_vep.vcf"
OUT="${GENOME_DIR}/${T}/vep/${T}_vep.vcf"
python3 - "$SRC" "$OUT" "${CASE_TMP}/ab.tsv" <<'PY'
import re, sys
src, out, ab = sys.argv[1:]
fmt, a, b, lines = None, None, None, []
for line in open(src):
    if line.startswith("##INFO=<ID=CSQ"):
        fmt = re.search(r"Format: ([^\"]+)", line).group(1).split("|")
        old = "|".join(fmt)
        missing = "CLIN_SIG" not in fmt
        fmt = fmt + (["CLIN_SIG"] if missing else [])
        line = line.replace(old, "|".join(fmt + ["ClinVar_CLNSIG"]))
        lines.append(line)
        continue
    if line.startswith("#"):
        lines.append(line)
        continue
    r = line.rstrip("\n").split("\t")
    if r[6] != "PASS":
        lines.append(line)
        continue
    info = r[7].split(";")
    k = next(i for i, x in enumerate(info) if x.startswith("CSQ="))
    trs = info[k][4:].split(",")
    ci = fmt.index("CLIN_SIG")
    if a is None:
        a, add, cache = (r[0], r[1]), "Pathogenic", ""
    elif b is None:
        b, add, cache = (r[0], r[1]), "", "pathogenic"
    else:
        add, cache = "", None
    new = []
    for tr in trs:
        v = tr.split("|") + ([""] if missing else [])
        if cache is not None:
            v[ci] = cache
        new.append("|".join(v + [add]))
    info[k] = "CSQ=" + ",".join(new)
    r[7] = ";".join(info)
    lines.append("\t".join(r) + "\n")
open(out, "w").writelines(lines)
open(ab, "w").write(f"{a[0]}\t{a[1]}\n{b[0]}\t{b[1]}\n")
PY
read -r A_CHR A_POS < <(sed -n 1p "${CASE_TMP}/ab.tsv")
read -r B_CHR B_POS < <(sed -n 2p "${CASE_TMP}/ab.tsv")
echo "A ${A_CHR}:${A_POS} (ClinVar_CLNSIG Pathogenic), B ${B_CHR}:${B_POS} (CLIN_SIG pathogenic only)"

run_step 23-clinical-filter.sh "$T"
check_step_exit 23-clinical-filter.sh
check "the log names ClinVar_CLNSIG as the tier's source" has 'ClinVar tier: ClinVar_CLNSIG' "$(cat "$STEP_LOG")"
TIER="${T}/clinical/${T}_clinvar_pathogenic.vcf.gz"
check "the ClinVar tier VCF is readable" vcf_ok "$TIER"
check_ge "record A (Pathogenic in the current ClinVar file) is in the tier" "$(vcf_count -r "${A_CHR}:${A_POS}" "$TIER")" 1
check_eq "record B (pathogenic in the cache only) is not" "$(vcf_count -r "${B_CHR}:${B_POS}" "$TIER")" 0

finish
