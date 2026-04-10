/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DELLY — Structural variant caller (paired-end + split-read + read-depth)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Calls structural variants (DEL, DUP, INV, BND, INS) and converts
    the native BCF output to an indexed VCF.

    Equivalent to: scripts/19-delly.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DELLY {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/delly:1.7.3--hd6466ae_0'

    publishDir "${params.outdir}/${meta.id}/delly", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_fai)

    output:
    tuple val(meta), path("${meta.id}_sv.vcf.gz"),     emit: sv_vcf
    tuple val(meta), path("${meta.id}_sv.vcf.gz.tbi"), emit: sv_vcf_index
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Call structural variants (outputs BCF)
    delly call \\
        -g ${reference} \\
        -o ${meta.id}_sv.bcf \\
        ${bam}

    # Convert BCF to VCF
    bcftools view \\
        ${meta.id}_sv.bcf \\
        -Oz -o ${meta.id}_sv.vcf.gz

    # Index VCF
    bcftools index -t ${meta.id}_sv.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly --version 2>&1 | grep -oP 'Delly version: \\K[\\d.]+' || echo '1.7.3')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_sv.vcf.gz
    touch ${meta.id}_sv.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: 1.7.3
    END_VERSIONS
    """
}
