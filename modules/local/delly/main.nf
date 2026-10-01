/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DELLY — Structural variant caller (paired-end + split-read + read-depth)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Calls structural variants (DEL, DUP, INV, BND, INS) into a BCF.

    Two processes, because the delly image carries htslib but not bcftools:
      1. DELLY          — delly call (optionally with an exclude map, -x)
      2. DELLY_BCF2VCF  — bcftools image, BCF -> bgzipped VCF + tabix index

    Equivalent to: scripts/19-delly.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DELLY {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/delly:2.1.0--h3752d28_0'

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_fai)
    path(exclude)

    output:
    tuple val(meta), path("${meta.id}_sv.bcf"), emit: bcf
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def exclude_arg = exclude ? "-x ${exclude}" : ""
    """
    delly call \\
        -g ${reference} \\
        ${exclude_arg} \\
        -o ${meta.id}_sv.bcf \\
        ${bam}

    DELLY_VERSION=\$(delly --version 2>&1 | sed -n 's/.*Delly version: v\\{0,1\\}\\([0-9][0-9.]*\\).*/\\1/p' | head -1)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \${DELLY_VERSION:-unknown}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_sv.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: 2.1.0
    END_VERSIONS
    """
}

process DELLY_BCF2VCF {
    tag "$meta.id"
    label 'process_single'

    container 'staphb/bcftools:1.21'

    publishDir { "${params.outdir}/${meta.id}/delly" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bcf)

    output:
    tuple val(meta), path("${meta.id}_sv.vcf.gz"),     emit: sv_vcf
    tuple val(meta), path("${meta.id}_sv.vcf.gz.tbi"), emit: sv_vcf_index
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools view ${bcf} -Oz -o ${meta.id}_sv.vcf.gz
    bcftools index -t ${meta.id}_sv.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_sv.vcf.gz
    touch ${meta.id}_sv.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.21
    END_VERSIONS
    """
}
