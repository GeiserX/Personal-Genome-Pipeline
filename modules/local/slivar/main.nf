/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SLIVAR — Variant prioritization and compound heterozygote detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Prioritizes clinically interesting variants using tiered filters (rare HIGH,
    rare MODERATE + deleterious predictors, ClinVar pathogenic) and detects
    compound heterozygote candidates. Optionally annotates results with gnomAD
    gene constraint metrics (LOEUF, pLI).

    Equivalent to: scripts/31-slivar.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SLIVAR {
    tag "$meta.id"
    label 'process_medium'

    container 'staphb/bcftools:1.21'

    publishDir "${params.outdir}/${meta.id}/slivar", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(gnomad_constraint)
    path(slivar_bin)

    output:
    tuple val(meta), path("*_prioritized.vcf.gz"),     emit: vcf
    tuple val(meta), path("*_prioritized.vcf.gz.tbi"), emit: vcf_index
    tuple val(meta), path("*_compound_hets.vcf.gz"),   emit: compound_het_vcf, optional: true
    tuple val(meta), path("*_slivar_summary.tsv"),     emit: summary_tsv
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def has_constraint = gnomad_constraint ? true : false
    """
    # --- Make slivar binary executable (staged into workdir by Nextflow) ---
    chmod +x ${slivar_bin}
    export PATH="\${PWD}:\${PATH}"

    # --- Generate PED file for single sample ---
    SAMPLE_NAME=\$(bcftools query -l ${vcf} | head -1)
    echo -e "\${SAMPLE_NAME}\\t\${SAMPLE_NAME}\\t0\\t0\\t0\\t-9" > ${meta.id}.ped

    # --- Filter 1: rare_high (PASS + HIGH impact + gnomAD AF < 1%) ---
    bcftools view -f PASS ${vcf} | \\
        bcftools +split-vep - -c IMPACT,gnomADe_AF -s worst \\
            -i 'IMPACT="HIGH" && (gnomADe_AF<0.01 || gnomADe_AF=".")' \\
            -Oz -o ${meta.id}_rare_high.vcf.gz || \\
        bcftools view -f PASS ${vcf} | \\
            bcftools +split-vep - -c IMPACT -s worst \\
                -i 'IMPACT="HIGH"' \\
                -Oz -o ${meta.id}_rare_high.vcf.gz
    bcftools index -t ${meta.id}_rare_high.vcf.gz

    # --- Filter 2: rare_moderate_deleterious (two-pass: CSQ then INFO predictors) ---
    # Pass 1: extract rare MODERATE via split-vep (CSQ fields)
    bcftools view -f PASS ${vcf} | \\
        bcftools +split-vep - -c IMPACT,gnomADe_AF -s worst \\
            -i 'IMPACT="MODERATE" && (gnomADe_AF<0.01 || gnomADe_AF=".")' \\
            -Oz -o ${meta.id}_rare_moderate_all.vcf.gz || \\
        bcftools view -f PASS ${vcf} | \\
            bcftools +split-vep - -c IMPACT -s worst \\
                -i 'IMPACT="MODERATE"' \\
                -Oz -o ${meta.id}_rare_moderate_all.vcf.gz
    bcftools index -t ${meta.id}_rare_moderate_all.vcf.gz

    # Pass 2: gate on deleteriousness predictors in INFO fields (added by vcfanno)
    # Build predictor expression from available INFO fields in the VCF header
    PREDICTOR_PARTS=""
    INFO_HEADER=\$(bcftools view -h ${meta.id}_rare_moderate_all.vcf.gz | grep '^##INFO' || true)

    echo "\$INFO_HEADER" | grep -q 'ID=CADD_PHRED,' && \\
        PREDICTOR_PARTS="\${PREDICTOR_PARTS:+\${PREDICTOR_PARTS} || }INFO/CADD_PHRED>=20"
    echo "\$INFO_HEADER" | grep -q 'ID=CADD_PHRED_indel,' && \\
        PREDICTOR_PARTS="\${PREDICTOR_PARTS:+\${PREDICTOR_PARTS} || }INFO/CADD_PHRED_indel>=20"
    echo "\$INFO_HEADER" | grep -q 'ID=REVEL' && \\
        PREDICTOR_PARTS="\${PREDICTOR_PARTS:+\${PREDICTOR_PARTS} || }INFO/REVEL>=0.5"
    echo "\$INFO_HEADER" | grep -q 'ID=AM_class' && \\
        PREDICTOR_PARTS="\${PREDICTOR_PARTS:+\${PREDICTOR_PARTS} || }INFO/AM_class=\"likely_pathogenic\""
    echo "\$INFO_HEADER" | grep -q 'ID=SpliceAI,' && \\
        PREDICTOR_PARTS="\${PREDICTOR_PARTS:+\${PREDICTOR_PARTS} || }INFO/SpliceAI!=\".\""
    echo "\$INFO_HEADER" | grep -q 'ID=SpliceAI_indel,' && \\
        PREDICTOR_PARTS="\${PREDICTOR_PARTS:+\${PREDICTOR_PARTS} || }INFO/SpliceAI_indel!=\".\""

    if [ -n "\${PREDICTOR_PARTS}" ]; then
        echo "  Filtering MODERATE variants with predictors: \${PREDICTOR_PARTS}"
        bcftools view -i "\${PREDICTOR_PARTS}" \\
            ${meta.id}_rare_moderate_all.vcf.gz \\
            -Oz -o ${meta.id}_rare_moderate_del.vcf.gz
    else
        echo "  WARNING: No vcfanno predictor annotations found — including all rare MODERATE variants."
        echo "  Run vcfanno for CADD/REVEL/AlphaMissense/SpliceAI filtering."
        cp ${meta.id}_rare_moderate_all.vcf.gz ${meta.id}_rare_moderate_del.vcf.gz
    fi
    bcftools index -t ${meta.id}_rare_moderate_del.vcf.gz

    # --- Filter 3: clinvar_pathogenic ---
    bcftools view -f PASS ${vcf} | \\
        bcftools +split-vep - -c CLIN_SIG \\
            -i 'CLIN_SIG~"pathogenic" && CLIN_SIG!~"conflicting"' \\
            -Oz -o ${meta.id}_clinvar_path.vcf.gz && \\
        bcftools index -t ${meta.id}_clinvar_path.vcf.gz || \\
        echo "ClinVar filter skipped (CLIN_SIG not in annotations)"

    # --- Merge tiers into prioritized VCF ---
    MERGE_FILES="${meta.id}_rare_high.vcf.gz ${meta.id}_rare_moderate_del.vcf.gz"
    if [ -f "${meta.id}_clinvar_path.vcf.gz" ] && [ -s "${meta.id}_clinvar_path.vcf.gz" ]; then
        MERGE_FILES="\${MERGE_FILES} ${meta.id}_clinvar_path.vcf.gz"
    fi
    bcftools concat -a -D \${MERGE_FILES} | \\
        bcftools sort -Oz -o ${meta.id}_prioritized.vcf.gz
    bcftools index -t ${meta.id}_prioritized.vcf.gz

    # --- Compound heterozygote detection ---
    slivar compound-hets \\
        --allow-non-trios \\
        --vcf ${meta.id}_prioritized.vcf.gz \\
        --ped ${meta.id}.ped \\
        2>${meta.id}_compound_hets.log | \\
        bcftools view -Oz -o ${meta.id}_compound_hets.vcf.gz || \\
        echo "Compound het detection returned no results"

    # --- Generate summary TSV with optional gnomAD constraint enrichment ---
    # Extract raw variant info
    bcftools +split-vep \\
        ${meta.id}_prioritized.vcf.gz \\
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%IMPACT\\t%SYMBOL\\t%Consequence\\t%Existing_variation[\\t%GT]\\n' \\
        -s worst -d > ${meta.id}_variants_raw.tsv 2>/dev/null || true

    if [ "${has_constraint}" = "true" ] && [ -f "${gnomad_constraint}" ] && command -v python3 &>/dev/null; then
        # Join with gnomAD constraint metrics (LOEUF, pLI, mis_z)
        python3 -c "
import csv, sys

constraint = {}
try:
    with open('${gnomad_constraint}') as f:
        reader = csv.DictReader(f, delimiter='\\t')
        for row in reader:
            gene = row.get('gene', row.get('gene_symbol', row.get('symbol', '')))
            if not gene:
                continue
            loeuf = row.get('oe_lof_upper', row.get('lof.oe_ci.upper', '.'))
            pli = row.get('pLI', row.get('lof.pLI', '.'))
            mis_z = row.get('mis_z', row.get('missense.z_score', '.'))
            constraint[gene] = (loeuf, pli, mis_z)
except Exception as e:
    print(f'WARNING: Could not load constraint file: {e}', file=sys.stderr)

header = 'CHROM\\tPOS\\tREF\\tALT\\tIMPACT\\tSYMBOL\\tConsequence\\tExisting_variation\\tGT\\tLOEUF\\tpLI\\tmis_z\\tCONSTRAINED'
print(header)

try:
    with open('${meta.id}_variants_raw.tsv') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\\t')
            if len(fields) < 6:
                continue
            gene = fields[5]
            loeuf, pli, mis_z = constraint.get(gene, ('.', '.', '.'))
            constrained = 'NO'
            try:
                if loeuf != '.' and float(loeuf) < 0.35:
                    constrained = 'YES'
                elif pli != '.' and float(pli) > 0.9:
                    constrained = 'YES'
            except ValueError:
                pass
            print(f'{line}\\t{loeuf}\\t{pli}\\t{mis_z}\\t{constrained}')
except FileNotFoundError:
    pass
" > ${meta.id}_slivar_summary.tsv
    else
        {
            echo -e "CHROM\\tPOS\\tREF\\tALT\\tIMPACT\\tSYMBOL\\tConsequence\\tExisting_variation\\tGT"
            cat ${meta.id}_variants_raw.tsv 2>/dev/null || true
        } > ${meta.id}_slivar_summary.tsv
    fi

    # Clean up intermediate files
    rm -f ${meta.id}_variants_raw.tsv ${meta.id}_rare_high.vcf.gz* ${meta.id}_rare_moderate_del.vcf.gz* \\
          ${meta.id}_rare_moderate_all.vcf.gz* ${meta.id}_clinvar_path.vcf.gz*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slivar: \$(slivar expr 2>&1 | head -1 | grep -oP '[\\d.]+' || echo '0.3.3')
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_prioritized.vcf.gz
    touch ${meta.id}_prioritized.vcf.gz.tbi
    touch ${meta.id}_compound_hets.vcf.gz
    printf 'CHROM\\tPOS\\tREF\\tALT\\tIMPACT\\tSYMBOL\\tConsequence\\tExisting_variation\\tGT\\n' > ${meta.id}_slivar_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slivar: 0.3.3
        bcftools: 1.21
    END_VERSIONS
    """
}
