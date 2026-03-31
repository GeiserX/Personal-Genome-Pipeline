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
  python:3.11-slim \
  python3 -c "
import json, sys

json_path = '/genome/${SAMPLE}/pharmcat/${SAMPLE}.report.json'
alt_paths = [
  '/genome/${SAMPLE}/vcf/${SAMPLE}.report.json',
  '/genome/${SAMPLE}/pharmcat/report.json',
]

data = None
for path in [json_path] + alt_paths:
    try:
        with open(path) as f:
            data = json.load(f)
        break
    except (FileNotFoundError, json.JSONDecodeError):
        continue

if data is None:
    # Try to find any JSON file
    import glob
    for pattern in ['/genome/${SAMPLE}/pharmcat/*.json', '/genome/${SAMPLE}/vcf/*.json']:
        files = glob.glob(pattern)
        for f in files:
            try:
                with open(f) as fh:
                    data = json.load(fh)
                break
            except:
                continue
        if data:
            break

if data is None:
    print('NO_JSON_FOUND')
    sys.exit(0)

# Extract gene results — handle multiple PharmCAT JSON versions
if 'genes' in data and isinstance(data['genes'], dict):
    # PharmCAT 2.15+ format: genes -> {source -> {gene_name -> data}}
    for source, gene_dict in data['genes'].items():
        if not isinstance(gene_dict, dict):
            continue
        for gene_name, g in gene_dict.items():
            dips = g.get('sourceDiplotypes', [])
            if not dips:
                continue
            dip = dips[0]
            a1_obj = dip.get('allele1')
            a2_obj = dip.get('allele2')
            a1 = a1_obj.get('name', '?') if a1_obj else '?'
            a2 = a2_obj.get('name', '?') if a2_obj else '?'
            diplotype = f'{a1}/{a2}'
            phenotypes = dip.get('phenotypes', [])
            phenotype = phenotypes[0] if phenotypes else 'N/A'
            if phenotype != 'No Result':
                print(f'{gene_name}\t{diplotype}\t{phenotype}')
        break  # Only use first source (CPIC)
elif 'genes' in data and isinstance(data['genes'], list):
    # PharmCAT older list format
    for gene_entry in data['genes']:
        gene = gene_entry.get('geneSymbol', gene_entry.get('gene', 'Unknown'))
        diplotype = 'N/A'
        phenotype = 'N/A'
        src_dips = gene_entry.get('sourceDiplotypes', [])
        if src_dips:
            diplotype = src_dips[0].get('label', 'N/A')
            phenotype = src_dips[0].get('phenotype', 'N/A')
        if diplotype != 'N/A' or phenotype != 'N/A':
            print(f'{gene}\t{diplotype}\t{phenotype}')
elif 'geneResults' in data:
    for gene_result in data['geneResults']:
        gene = gene_result.get('gene', 'Unknown')
        diplotype = gene_result.get('diplotype', 'N/A')
        phenotype = gene_result.get('phenotype', 'N/A')
        if diplotype or phenotype:
            print(f'{gene}\t{diplotype}\t{phenotype}')
else:
    print('UNKNOWN_FORMAT')
" 2>/dev/null > "${OUTDIR}/${SAMPLE}_phenotypes.tsv" || true

# Generate recommendations
echo "" >> "$OUTPUT"
echo "Gene Results:" >> "$OUTPUT"
echo "─────────────" >> "$OUTPUT"
printf "%-12s %-25s %-35s\n" "Gene" "Diplotype" "Phenotype" >> "$OUTPUT"
printf "%-12s %-25s %-35s\n" "────" "─────────" "─────────" >> "$OUTPUT"

if [ -f "${OUTDIR}/${SAMPLE}_phenotypes.tsv" ] && [ -s "${OUTDIR}/${SAMPLE}_phenotypes.tsv" ]; then
  while IFS=$'\t' read -r GENE DIPLOTYPE PHENOTYPE; do
    if [ "$GENE" = "NO_JSON_FOUND" ] || [ "$GENE" = "UNKNOWN_FORMAT" ]; then
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
  if [ -z "$GENE" ] || [ "$GENE" = "NO_JSON_FOUND" ] || [ "$GENE" = "UNKNOWN_FORMAT" ]; then
    continue
  fi

  # Skip normal metabolizers — no dosing changes needed
  if echo "$PHENOTYPE" | grep -qi "normal\|typical\|extensive"; then
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
