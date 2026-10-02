#!/usr/bin/env bash
# Step 09 (ExpansionHunter) and step 09b (Stranger with its GRCh38 catalog).
# Stranger writes its output through a temporary name and annotates every
# ExpansionHunter locus it finds in the catalog.
. "$(dirname "$0")/lib.sh"

# The fixture reference holds 14 chromosomes, and ExpansionHunter stops on a
# catalog locus whose contig the reference lacks. The case keeps the bundled
# GRCh38 catalog's loci on the fixture's contigs and passes them as EH_CATALOG.
CAT="${GENOME_DIR}/reference/eh_catalog_fixture.json"
docker run --rm "$EXPANSIONHUNTER_IMAGE" cat /usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json \
  > "${CASE_TMP}/catalog.json"
python3 - "${CASE_TMP}/catalog.json" "${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta.fai" "$CAT" <<'PY2'
import json, sys
loci = json.load(open(sys.argv[1]))
contigs = {line.split("\t")[0] for line in open(sys.argv[2])}
def regions(locus):
    r = locus["ReferenceRegion"]
    return r if isinstance(r, list) else [r]
kept = [l for l in loci if all(r.split(":")[0] in contigs for r in regions(l))]
print(f"catalog loci: {len(loci)}, on the fixture's contigs: {len(kept)}")
json.dump(kept, open(sys.argv[3], "w"), indent=1)
PY2
N_CAT=$(python3 -c 'import json, sys; print(len(json.load(open(sys.argv[1]))))' "$CAT" 2>/dev/null || echo 0)
check_ge "catalog loci on the fixture's contigs" "$N_CAT" 5

EH_CATALOG="$CAT" run_step 09-expansion-hunter.sh "$SAMPLE" male
check_step_exit 09-expansion-hunter.sh
check "the log names EH_CATALOG" has 'Variant catalog: .*eh_catalog_fixture\.json \(EH_CATALOG\)' "$(cat "$STEP_LOG")"
EH="${SAMPLE}/expansion_hunter/${SAMPLE}_eh.vcf"
N_EH=$(grep -vc '^#' "${GENOME_DIR}/${EH}" 2>/dev/null || true)
check_ge "ExpansionHunter records" "$N_EH" 1

run_step 09b-stranger.sh "$SAMPLE"
check_step_exit 09b-stranger.sh
check "the log names the GRCh38 catalog" has 'variant_catalog_grch38' "$(cat "$STEP_LOG")"
OUT="${GENOME_DIR}/${SAMPLE}/expansion_hunter/${SAMPLE}_eh_stranger.vcf"
check "no temporary file is left" test ! -e "${OUT}.tmp"
check "the output declares STR_STATUS" grep -q '^##INFO=<ID=STR_STATUS' "$OUT"
check_eq "records in = records out" "$(grep -vc '^#' "$OUT" 2>/dev/null || true)" "$N_EH"
check_ge "records with the catalog's pathologic threshold" "$(grep -v '^#' "$OUT" 2>/dev/null | grep -c 'STR_PATHOLOGIC_MIN=' || true)" 3

finish
