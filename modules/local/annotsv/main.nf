/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANNOTSV — Annotate structural variants with ACMG pathogenicity classification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Classifies each SV into ACMG class 1-5 and adds gene/disease annotations.

    Equivalent to: scripts/05-annotsv.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ANNOTSV {
    tag "$meta.id"
    label 'process_medium'

    container 'getwilds/annotsv:3.4.4'

    publishDir "${params.outdir}/${meta.id}/annotsv", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(sv_vcf)

    output:
    tuple val(meta), path("${meta.id}_sv_annotated.tsv"), emit: annotated_tsv
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    AnnotSV \\
        -SVinputFile ${sv_vcf} \\
        -outputFile ${meta.id}_sv_annotated.tsv \\
        -genomeBuild GRCh38 \\
        -annotationMode both

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(AnnotSV -help 2>&1 | grep -oP 'AnnotSV \\K[\\d.]+' || echo '3.4.4')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_sv_annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: 3.4.4
    END_VERSIONS
    """
}
