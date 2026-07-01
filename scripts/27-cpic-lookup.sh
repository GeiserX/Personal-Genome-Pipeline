#!/usr/bin/env bash
# 27-cpic-lookup.sh — Look up CPIC drug-gene recommendations from PharmCAT results
# Usage: ./scripts/27-cpic-lookup.sh <sample_name>
#
# Parses PharmCAT output to extract gene phenotypes, then looks up CPIC
# drug-gene recommendations using a built-in static table. Produces a
# plain-text report of medications that may require dosing adjustments
# based on your pharmacogenomic profile.
#
# Requires: PharmCAT output from step 7. No internet connection needed.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=../versions.env
. "${SCRIPT_DIR}/../versions.env"

# Find PharmCAT JSON output
PHARMCAT_JSON=""
for DIR in "${GENOME_DIR}/${SAMPLE}/pharmcat" "${GENOME_DIR}/${SAMPLE}/vcf"; do
  for FILE in "${DIR}"/*.report.json "${DIR}"/*_pharmcat.json; do
    if [ -f "$FILE" ] 2>/dev/null; then
      PHARMCAT_JSON="$FILE"
      break 2
    fi
  done
done

if [ -z "$PHARMCAT_JSON" ]; then
  echo "ERROR: PharmCAT output not found. Run step 7 first."
  echo "  Expected in: ${GENOME_DIR}/${SAMPLE}/pharmcat/ or ${GENOME_DIR}/${SAMPLE}/vcf/"
  exit 1
fi

OUTDIR="${GENOME_DIR}/${SAMPLE}/cpic"
mkdir -p "$OUTDIR"
OUTPUT="${OUTDIR}/${SAMPLE}_cpic_recommendations.txt"

echo "============================================"
echo "  Step 27: CPIC Drug-Gene Recommendations"
echo "  Sample: ${SAMPLE}"
echo "  Input:  ${PHARMCAT_JSON}"
echo "  Output: ${OUTPUT}"
echo "============================================"
echo ""

# Known CPIC gene-drug pairs with direct clinical recommendations
# Source: https://cpicpgx.org/guidelines/
declare -A CPIC_DRUGS
CPIC_DRUGS=(
  ["CYP2C19"]="clopidogrel,voriconazole,escitalopram,citalopram,sertraline,amitriptyline,clomipramine,doxepin,imipramine,trimipramine,lansoprazole,omeprazole,pantoprazole,dexlansoprazole"
  ["CYP2C9"]="warfarin,phenytoin,flurbiprofen,celecoxib,ibuprofen,lornoxicam,meloxicam,piroxicam,tenoxicam,siponimod"
  ["CYP2D6"]="codeine,tramadol,hydrocodone,oxycodone,amitriptyline,clomipramine,desipramine,doxepin,imipramine,nortriptyline,trimipramine,fluvoxamine,paroxetine,atomoxetine,ondansetron,tropisetron,tamoxifen,eliglustat"
  ["CYP3A5"]="tacrolimus"
  ["DPYD"]="fluorouracil,capecitabine,tegafur"
  ["TPMT"]="azathioprine,mercaptopurine,thioguanine"
  ["NUDT15"]="azathioprine,mercaptopurine,thioguanine"
  ["UGT1A1"]="atazanavir,belinostat,irinotecan"
  # CPIC 2022 statin guideline covers all statins, but clinical impact varies:
  # simvastatin has strongest evidence, atorvastatin/rosuvastatin moderate, others weaker.
  ["SLCO1B1"]="simvastatin,atorvastatin,rosuvastatin,pravastatin,pitavastatin,fluvastatin,lovastatin"
  ["VKORC1"]="warfarin"
  ["HLA-A"]="carbamazepine,oxcarbazepine"
  ["HLA-B"]="abacavir,carbamazepine,oxcarbazepine,phenytoin,allopurinol"
  # IFNL3/peginterferon is largely historical — DAAs (sofosbuvir, etc.) have replaced
  # interferon-based HCV therapy. Retained for completeness.
  ["IFNL3"]="peginterferon-alfa-2a,peginterferon-alfa-2b"
  # CYP2B6-methadone was evaluated by CPIC (2024) but classified as optional, not
  # a direct prescribing recommendation. Only efavirenz has actionable CPIC guidance.
  ["CYP2B6"]="efavirenz"
  ["RYR1"]="desflurane,enflurane,halothane,isoflurane,methoxyflurane,sevoflurane,succinylcholine"
  ["CACNA1S"]="desflurane,enflurane,halothane,isoflurane,methoxyflurane,sevoflurane,succinylcholine"
  ["G6PD"]="rasburicase,dapsone,chloroquine,primaquine,nitrofurantoin,methylene-blue"
  ["MT-RNR1"]="aminoglycosides"
)

echo "Pharmacogenomic Drug Recommendations" > "$OUTPUT"
echo "=====================================" >> "$OUTPUT"
echo "Sample: ${SAMPLE}" >> "$OUTPUT"
echo "Date: $(date '+%Y-%m-%d')" >> "$OUTPUT"
echo "Source: PharmCAT + CPIC Guidelines" >> "$OUTPUT"
echo "" >> "$OUTPUT"

# Parse PharmCAT phenotype data
# Try to extract from the JSON or fall back to the HTML report
echo "Parsing PharmCAT results..."

# Use Python in Docker to parse JSON properly
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  "${PYTHON_IMAGE}" \
  python3 -c "
import json, sys, glob

candidate_paths = [
  '/genome/${SAMPLE}/pharmcat/${SAMPLE}.report.json',
  '/genome/${SAMPLE}/vcf/${SAMPLE}.report.json',
  '/genome/${SAMPLE}/pharmcat/report.json',
]

data = None
for path in candidate_paths:
    try:
        with open(path) as f:
            data = json.load(f)
        break
    except (FileNotFoundError, json.JSONDecodeError):
        continue

if data is None:
    for pattern in ['/genome/${SAMPLE}/pharmcat/*.json', '/genome/${SAMPLE}/vcf/*.json']:
        for fp in glob.glob(pattern):
            try:
                with open(fp) as fh:
                    data = json.load(fh)
                break
            except Exception:
                continue
        if data is not None:
            break

if data is None:
    print('NO_JSON_FOUND')
    sys.exit(0)

# NOTE: PharmCAT 3.x renamed some JSON fields (e.g. wildtypeAllele -> referenceAllele)
# and may emit either a flat {gene -> data} or a nested {source -> {gene -> data}} 'genes'
# map. Auto-detect both so this parser works across versions, and accept either
# sourceDiplotypes or recommendationDiplotypes.
def parse_gene(name, g):
    if not isinstance(g, dict):
        return None
    dips = g.get('sourceDiplotypes') or g.get('recommendationDiplotypes') or []
    if not dips:
        return None
    dip = dips[0]
    a1 = (dip.get('allele1') or {}).get('name', '?')
    a2 = (dip.get('allele2') or {}).get('name', '?')
    diplotype = dip.get('label') or (a1 + '/' + a2)
    phenos = dip.get('phenotypes') or []
    phenotype = phenos[0] if phenos else dip.get('phenotype', 'N/A')
    return (name, diplotype, phenotype)

results = []
genes = data.get('genes')
if isinstance(genes, dict):
    for key, val in genes.items():
        if not isinstance(val, dict):
            continue
        if 'sourceDiplotypes' in val or 'recommendationDiplotypes' in val:
            r = parse_gene(key, val)            # flat: key is the gene
            if r:
                results.append(r)
        else:
            for gene_name, g in val.items():    # nested: key is the source
                r = parse_gene(gene_name, g)
                if r:
                    results.append(r)
elif isinstance(genes, list):
    for entry in genes:
        name = entry.get('geneSymbol', entry.get('gene', 'Unknown'))
        r = parse_gene(name, entry)
        if r:
            results.append(r)
        else:
            dl = entry.get('diplotype', 'N/A')
            ph = entry.get('phenotype', 'N/A')
            if dl != 'N/A' or ph != 'N/A':
                results.append((name, dl, ph))
elif isinstance(data.get('geneResults'), list):
    for gr in data['geneResults']:
        results.append((gr.get('gene', 'Unknown'), gr.get('diplotype', 'N/A'), gr.get('phenotype', 'N/A')))
else:
    print('UNKNOWN_FORMAT')
    sys.exit(0)

seen = set()
deduped = []
for g, d, p in results:
    if g in seen:
        continue
    seen.add(g)
    deduped.append((g, d, p))

# FAIL LOUD: a recognized PharmCAT report that yields zero genes is a format break,
# never a clean result — emit a sentinel so the report cannot say 'all clear'.
if not deduped:
    print('PARSE_EMPTY')
    sys.exit(0)

for g, d, p in deduped:
    print(g + '\t' + d + '\t' + p)
" 2>/dev/null > "${OUTDIR}/${SAMPLE}_phenotypes.tsv" || true

# Generate recommendations
echo "" >> "$OUTPUT"
echo "Gene Results:" >> "$OUTPUT"
echo "─────────────" >> "$OUTPUT"
printf "%-12s %-25s %-35s\n" "Gene" "Diplotype" "Phenotype" >> "$OUTPUT"
printf "%-12s %-25s %-35s\n" "────" "─────────" "─────────" >> "$OUTPUT"

if [ -f "${OUTDIR}/${SAMPLE}_phenotypes.tsv" ] && [ -s "${OUTDIR}/${SAMPLE}_phenotypes.tsv" ]; then
  while IFS=$'\t' read -r GENE DIPLOTYPE PHENOTYPE; do
    if [ "$GENE" = "NO_JSON_FOUND" ] || [ "$GENE" = "UNKNOWN_FORMAT" ] || [ "$GENE" = "PARSE_EMPTY" ]; then
      echo "  WARNING: Could not parse PharmCAT JSON output." >> "$OUTPUT"
      break
    fi
    printf "%-12s %-25s %-35s\n" "$GENE" "$DIPLOTYPE" "$PHENOTYPE" >> "$OUTPUT"
  done < "${OUTDIR}/${SAMPLE}_phenotypes.tsv"
fi

# Add drug recommendations based on known CPIC pairs
echo "" >> "$OUTPUT"
echo "Affected Medications:" >> "$OUTPUT"
echo "─────────────────────" >> "$OUTPUT"
echo "" >> "$OUTPUT"

while IFS=$'\t' read -r GENE DIPLOTYPE PHENOTYPE; do
  if [ -z "$GENE" ] || [ "$GENE" = "NO_JSON_FOUND" ] || [ "$GENE" = "UNKNOWN_FORMAT" ] || [ "$GENE" = "PARSE_EMPTY" ]; then
    continue
  fi

  # Skip normal metabolizers — exact word match to avoid substring false positives
  if echo "$PHENOTYPE" | grep -qiw "normal\|typical\|extensive"; then
    continue
  fi

  # Flag uncallable genes — absence from report ≠ normal
  if echo "$PHENOTYPE" | grep -qi "no result\|N/A\|indeterminate"; then
    continue
  fi

  DRUGS="${CPIC_DRUGS[$GENE]:-}"
  if [ -n "$DRUGS" ]; then
    echo "  ${GENE} — ${PHENOTYPE}:" >> "$OUTPUT"
    echo "    Diplotype: ${DIPLOTYPE}" >> "$OUTPUT"
    echo "    Affected drugs: ${DRUGS}" >> "$OUTPUT"
    echo "    Action: Consult CPIC guidelines at https://cpicpgx.org/guidelines/" >> "$OUTPUT"
    echo "" >> "$OUTPUT"
  fi
done < "${OUTDIR}/${SAMPLE}_phenotypes.tsv" 2>/dev/null || true

# List genes that could not be called — absence from report ≠ normal
echo "Uncallable Genes:" >> "$OUTPUT"
echo "─────────────────" >> "$OUTPUT"
echo "" >> "$OUTPUT"
# Check if PharmCAT parsing failed entirely
PARSE_FAILED=false
if [ ! -f "${OUTDIR}/${SAMPLE}_phenotypes.tsv" ] || [ ! -s "${OUTDIR}/${SAMPLE}_phenotypes.tsv" ]; then
  # No TSV, or an empty one, means nothing was parsed — never report that as "all clear".
  PARSE_FAILED=true
else
  FIRST_LINE=$(head -1 "${OUTDIR}/${SAMPLE}_phenotypes.tsv")
  case "$FIRST_LINE" in
    NO_JSON_FOUND*|UNKNOWN_FORMAT*|PARSE_EMPTY*)
      PARSE_FAILED=true
      ;;
  esac
fi

UNCALLED=0
if [ "$PARSE_FAILED" = true ]; then
  echo "  WARNING: PharmCAT output could not be parsed (no gene phenotypes extracted)." >> "$OUTPUT"
  echo "  This is NOT a clean result — likely a PharmCAT output-format change." >> "$OUTPUT"
  echo "  Consult the authoritative PharmCAT HTML report; no gene results can be trusted." >> "$OUTPUT"
  UNCALLED=1
else
  while IFS=$'\t' read -r GENE DIPLOTYPE PHENOTYPE; do
    if echo "$PHENOTYPE" | grep -qi "no result\|N/A\|indeterminate"; then
      echo "  ${GENE} — ${PHENOTYPE} (not callable from available data)" >> "$OUTPUT"
      UNCALLED=$((UNCALLED + 1))
    fi
  done < "${OUTDIR}/${SAMPLE}_phenotypes.tsv" 2>/dev/null || true
  if [ "$UNCALLED" -eq 0 ]; then
    echo "  None — all genes were successfully called." >> "$OUTPUT"
  fi
fi
echo "" >> "$OUTPUT"
echo "NOTE: Uncallable genes may lack coverage, have complex structural" >> "$OUTPUT"
echo "variants, or require data not present in your VCF. Their absence" >> "$OUTPUT"
echo "from the recommendations section does NOT mean normal function." >> "$OUTPUT"
echo "" >> "$OUTPUT"

echo "─────────────────────────────────────────────────────────────" >> "$OUTPUT"
echo "DISCLAIMER: These recommendations are based on CPIC clinical" >> "$OUTPUT"
echo "guidelines. Always consult your healthcare provider before" >> "$OUTPUT"
echo "making any medication changes. This is NOT medical advice." >> "$OUTPUT"
echo "─────────────────────────────────────────────────────────────" >> "$OUTPUT"

echo ""
echo "============================================"
echo "  CPIC recommendations complete: ${SAMPLE}"
echo "  Output: ${OUTPUT}"
echo "============================================"
echo ""
cat "$OUTPUT"
