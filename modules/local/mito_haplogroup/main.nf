/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Mitochondrial Haplogroup — Determine maternal lineage from mtDNA variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Two-step process:
    1. Extract chrM variants from full-genome VCF (bcftools)
    2. Classify haplogroup with haplogrep3

    Equivalent to: scripts/12-mito-haplogroup.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MITO_EXTRACT_CHRM {
    tag "$meta.id"
    label 'process_single'

    container 'staphb/bcftools:1.21'

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("*_chrM.vcf.gz"), path("*_chrM.vcf.gz.tbi"), emit: chrm_vcf
    path "versions.yml",                                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools view -r chrM ${vcf} -Oz -o ${meta.id}_chrM.vcf.gz
    bcftools index -t ${meta.id}_chrM.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_chrM.vcf.gz
    touch ${meta.id}_chrM.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}

process MITO_HAPLOGROUP {
    tag "$meta.id"
    label 'process_single'

    container 'jtb114/haplogrep3:latest'

    publishDir "${params.outdir}/${meta.id}/mito", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(chrm_vcf), path(chrm_vcf_index)

    output:
    tuple val(meta), path("*_haplogroup.txt"), emit: haplogroup
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    haplogrep3 classify \\
        --tree phylotree-fu-rcrs@1.2 \\
        --input ${chrm_vcf} \\
        --output ${meta.id}_haplogroup.txt \\
        --extend-report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplogrep3: \$(haplogrep3 --version 2>&1 | grep -oP '[\\d.]+' | head -1 || echo 'latest')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_haplogroup.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplogrep3: latest
    END_VERSIONS
    """
}
