/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CPIC Lookup — Drug-gene recommendations from PharmCAT results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Parses PharmCAT JSON output to extract gene phenotypes, then looks up CPIC
    drug-gene recommendations using a built-in static table. Produces a plain-text
    report of medications that may require dosing adjustments based on the
    pharmacogenomic profile.

    Equivalent to: scripts/27-cpic-lookup.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CPIC_LOOKUP {
    tag "$meta.id"
    label 'process_single'

    container 'python:3.11'

    publishDir "${params.outdir}/${meta.id}/cpic", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(pharmcat_json)

    output:
    tuple val(meta), path("${meta.id}_cpic_recommendations.txt"), emit: recommendations
    tuple val(meta), path("${meta.id}_phenotypes.tsv"),           emit: phenotypes
    path "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 << 'PYEOF'
import json, sys, os
from datetime import date

sample = "${meta.id}"
json_path = "${pharmcat_json}"
recommendations_path = f"{sample}_cpic_recommendations.txt"
phenotypes_path = f"{sample}_phenotypes.tsv"

# --- Known CPIC gene-drug pairs with direct clinical recommendations ---
# Source: https://cpicpgx.org/guidelines/
CPIC_DRUGS = {
    "CYP2C19": "clopidogrel,voriconazole,escitalopram,citalopram,sertraline,amitriptyline,clomipramine,doxepin,imipramine,trimipramine,lansoprazole,omeprazole,pantoprazole,dexlansoprazole",
    "CYP2C9": "warfarin,phenytoin,flurbiprofen,celecoxib,ibuprofen,lornoxicam,meloxicam,piroxicam,tenoxicam,siponimod",
    "CYP2D6": "codeine,tramadol,hydrocodone,oxycodone,amitriptyline,clomipramine,desipramine,doxepin,imipramine,nortriptyline,trimipramine,fluvoxamine,paroxetine,atomoxetine,ondansetron,tropisetron,tamoxifen,eliglustat",
    "CYP3A5": "tacrolimus",
    "DPYD": "fluorouracil,capecitabine,tegafur",
    "TPMT": "azathioprine,mercaptopurine,thioguanine",
    "NUDT15": "azathioprine,mercaptopurine,thioguanine",
    "UGT1A1": "atazanavir,belinostat,irinotecan",
    "SLCO1B1": "simvastatin,atorvastatin,rosuvastatin,pravastatin,pitavastatin,fluvastatin,lovastatin",
    "VKORC1": "warfarin",
    "HLA-A": "carbamazepine,oxcarbazepine",
    "HLA-B": "abacavir,carbamazepine,oxcarbazepine,phenytoin,allopurinol",
    "IFNL3": "peginterferon-alfa-2a,peginterferon-alfa-2b",
    "CYP2B6": "efavirenz",
    "RYR1": "desflurane,enflurane,halothane,isoflurane,methoxyflurane,sevoflurane,succinylcholine",
    "CACNA1S": "desflurane,enflurane,halothane,isoflurane,methoxyflurane,sevoflurane,succinylcholine",
    "G6PD": "rasburicase,dapsone,chloroquine,primaquine,nitrofurantoin,methylene-blue",
    "MT-RNR1": "aminoglycosides",
}

# --- Parse PharmCAT JSON ---
data = None
try:
    with open(json_path) as f:
        data = json.load(f)
except (FileNotFoundError, json.JSONDecodeError) as e:
    print(f"ERROR: Could not read PharmCAT JSON: {e}", file=sys.stderr)
    sys.exit(1)

gene_results = []  # list of (gene, diplotype, phenotype)

if "genes" in data and isinstance(data["genes"], dict):
    # PharmCAT 3.x format: genes -> {gene_name -> gene_data}
    for gene_name, g in data["genes"].items():
        if not isinstance(g, dict):
            continue
        dips = g.get("sourceDiplotypes", [])
        if not dips:
            continue
        dip = dips[0]
        a1_obj = dip.get("allele1")
        a2_obj = dip.get("allele2")
        a1 = a1_obj.get("name", "?") if a1_obj else "?"
        a2 = a2_obj.get("name", "?") if a2_obj else "?"
        diplotype = f"{a1}/{a2}"
        phenotypes = dip.get("phenotypes", [])
        phenotype = phenotypes[0] if phenotypes else "N/A"
        gene_results.append((gene_name, diplotype, phenotype))
elif "genes" in data and isinstance(data["genes"], list):
    # PharmCAT older list format
    for gene_entry in data["genes"]:
        gene = gene_entry.get("geneSymbol", gene_entry.get("gene", "Unknown"))
        diplotype = "N/A"
        phenotype = "N/A"
        src_dips = gene_entry.get("sourceDiplotypes", [])
        if src_dips:
            diplotype = src_dips[0].get("label", "N/A")
            phenotype = src_dips[0].get("phenotype", "N/A")
        if diplotype != "N/A" or phenotype != "N/A":
            gene_results.append((gene, diplotype, phenotype))
elif "geneResults" in data:
    for gene_result in data["geneResults"]:
        gene = gene_result.get("gene", "Unknown")
        diplotype = gene_result.get("diplotype", "N/A")
        phenotype = gene_result.get("phenotype", "N/A")
        if diplotype or phenotype:
            gene_results.append((gene, diplotype, phenotype))
else:
    print("ERROR: Unknown PharmCAT JSON format — no gene results extracted", file=sys.stderr)
    print("  Available top-level keys: " + ", ".join(data.keys() if data else []), file=sys.stderr)
    sys.exit(1)

# --- Write phenotypes TSV ---
with open(phenotypes_path, "w") as f:
    f.write("Gene\\tDiplotype\\tPhenotype\\n")
    for gene, diplotype, phenotype in gene_results:
        f.write(f"{gene}\\t{diplotype}\\t{phenotype}\\n")

# --- Write recommendations report ---
with open(recommendations_path, "w") as out:
    out.write("Pharmacogenomic Drug Recommendations\\n")
    out.write("=====================================\\n")
    out.write(f"Sample: {sample}\\n")
    out.write(f"Date: {date.today().isoformat()}\\n")
    out.write("Source: PharmCAT + CPIC Guidelines\\n\\n")

    out.write("Gene Results:\\n")
    out.write("-" * 72 + "\\n")
    out.write(f"{'Gene':<12} {'Diplotype':<25} {'Phenotype':<35}\\n")
    out.write(f"{'----':<12} {'---------':<25} {'---------':<35}\\n")
    for gene, diplotype, phenotype in gene_results:
        out.write(f"{gene:<12} {diplotype:<25} {phenotype:<35}\\n")

    out.write("\\nAffected Medications:\\n")
    out.write("-" * 72 + "\\n\\n")

    affected_count = 0
    for gene, diplotype, phenotype in gene_results:
        # Skip normal metabolizers (exact word match to avoid substring false positives)
        phenotype_lower = phenotype.lower()
        phenotype_words = set(phenotype_lower.replace('-', ' ').split())
        if phenotype_words & {"normal", "typical", "extensive"}:
            continue
        # Skip uncallable
        if any(kw in phenotype_lower for kw in ["no result", "n/a", "indeterminate"]):
            continue
        drugs = CPIC_DRUGS.get(gene, "")
        if drugs:
            out.write(f"  {gene} -- {phenotype}:\\n")
            out.write(f"    Diplotype: {diplotype}\\n")
            out.write(f"    Affected drugs: {drugs}\\n")
            out.write(f"    Action: Consult CPIC guidelines at https://cpicpgx.org/guidelines/\\n\\n")
            affected_count += 1

    if affected_count == 0:
        out.write("  No non-normal phenotypes with CPIC drug recommendations detected.\\n\\n")

    # Uncallable genes
    out.write("Uncallable Genes:\\n")
    out.write("-" * 72 + "\\n\\n")
    uncalled = 0
    for gene, diplotype, phenotype in gene_results:
        phenotype_lower = phenotype.lower()
        if any(kw in phenotype_lower for kw in ["no result", "n/a", "indeterminate"]):
            out.write(f"  {gene} -- {phenotype} (not callable from available data)\\n")
            uncalled += 1
    if uncalled == 0:
        out.write("  None -- all genes were successfully called.\\n")

    out.write("\\n")
    out.write("NOTE: Uncallable genes may lack coverage, have complex structural\\n")
    out.write("variants, or require data not present in your VCF. Their absence\\n")
    out.write("from the recommendations section does NOT mean normal function.\\n\\n")
    out.write("-" * 72 + "\\n")
    out.write("DISCLAIMER: These recommendations are based on CPIC clinical\\n")
    out.write("guidelines. Always consult your healthcare provider before\\n")
    out.write("making any medication changes. This is NOT medical advice.\\n")
    out.write("-" * 72 + "\\n")

print(f"CPIC recommendations written: {recommendations_path}")
print(f"Phenotypes written: {phenotypes_path}")
print(f"Gene results: {len(gene_results)}, affected medications: {affected_count}")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_cpic_recommendations.txt
    printf 'Gene\\tDiplotype\\tPhenotype\\n' > ${meta.id}_phenotypes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
