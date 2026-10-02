/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANNOTSV — Annotate structural variants with ACMG pathogenicity classification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Classifies each SV into ACMG class 1-5 and adds gene/disease annotations.

    The biocontainer holds AnnotSV's code but none of its annotation data, so
    the annotation directory (built once with AnnotSV's INSTALL_annotations.sh,
    it contains Annotations_Human/) comes in through --annotsv_annotations.

    Equivalent to: scripts/05-annotsv.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ANNOTSV {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/annotsv" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(sv_vcf)
    path(annotations_dir)

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
        -outputDir . \\
        -genomeBuild GRCh38 \\
        -annotationsDir ${annotations_dir} \\
        -annotationMode both

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_sv_annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
