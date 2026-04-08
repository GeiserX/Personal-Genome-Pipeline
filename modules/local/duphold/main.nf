/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DUPHOLD — Annotate structural variants with depth-based quality metrics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Adds DHBFC (depth fold-change at breakpoints) and DHFFC (depth fold-change
    at flanks) annotations to SV VCFs for filtering false positives.

    Equivalent to: scripts/15-duphold.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DUPHOLD {
    tag "$meta.id"
    label 'process_medium'

    container 'brentp/duphold:v0.2.3'

    publishDir "${params.outdir}/${meta.id}/sv_filtered", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(sv_vcf), path(bam), path(bai)
    path(reference)
    path(reference_fai)

    output:
    tuple val(meta), path("${meta.id}_sv_duphold.vcf"), emit: annotated_vcf
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    duphold \\
        -v ${sv_vcf} \\
        -b ${bam} \\
        -f ${reference} \\
        -o ${meta.id}_sv_duphold.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        duphold: \$(duphold --version 2>&1 | grep -oP '[\\d.]+' || echo '0.2.3')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_sv_duphold.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        duphold: 0.2.3
    END_VERSIONS
    """
}
