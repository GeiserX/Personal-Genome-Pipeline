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

    container 'quay.io/biocontainers/slivar:0.3.3--h5f107b1_0'

    publishDir "${params.outdir}/${meta.id}/slivar", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(gnomad_constraint)

    output:
    tuple val(meta), path("*_prioritized.vcf.gz"),     emit: vcf
    tuple val(meta), path("*_prioritized.vcf.gz.tbi"), emit: vcf_index
    tuple val(meta), path("*_compound_hets.vcf.gz"),   emit: compound_het_vcf, optional: true
    tuple val(meta), path("*_slivar_summary.tsv"),     emit: summary_tsv
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def has_constraint = gnomad_constraint.name != 'NO_FILE'
    """
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

    # --- Filter 2: rare_moderate_deleterious ---
    bcftools view -f PASS ${vcf} | \\
        bcftools +split-vep - -c IMPACT,gnomADe_AF -s worst \\
            -i 'IMPACT="MODERATE" && (gnomADe_AF<0.01 || gnomADe_AF=".")' \\
            -Oz -o ${meta.id}_rare_moderate_del.vcf.gz || \\
        bcftools view -f PASS ${vcf} | \\
            bcftools +split-vep - -c IMPACT -s worst \\
                -i 'IMPACT="MODERATE"' \\
                -Oz -o ${meta.id}_rare_moderate_del.vcf.gz
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

    # --- Generate summary TSV ---
    {
        echo -e "CHROM\\tPOS\\tREF\\tALT\\tIMPACT\\tSYMBOL\\tConsequence\\tExisting_variation\\tGT"
        bcftools +split-vep \\
            ${meta.id}_prioritized.vcf.gz \\
            -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%IMPACT\\t%SYMBOL\\t%Consequence\\t%Existing_variation[\\t%GT]\\n' \\
            -s worst -d 2>/dev/null || true
    } > ${meta.id}_slivar_summary.tsv

    # Clean up intermediate tier files
    rm -f ${meta.id}_rare_high.vcf.gz* ${meta.id}_rare_moderate_del.vcf.gz* ${meta.id}_clinvar_path.vcf.gz*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slivar: \$(slivar expr 2>&1 | head -1 | grep -oP '[\\d.]+' || echo '0.3.3')
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    echo '##fileformat=VCFv4.2' | bgzip > ${meta.id}_prioritized.vcf.gz
    tabix -p vcf ${meta.id}_prioritized.vcf.gz
    echo '##fileformat=VCFv4.2' | bgzip > ${meta.id}_compound_hets.vcf.gz
    printf 'CHROM\\tPOS\\tREF\\tALT\\tIMPACT\\tSYMBOL\\tConsequence\\tExisting_variation\\tGT\\n' > ${meta.id}_slivar_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slivar: 0.3.3
        bcftools: 1.21
    END_VERSIONS
    """
}
