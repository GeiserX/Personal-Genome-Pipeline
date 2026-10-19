/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANNOTSV — Annotate structural variants with ACMG pathogenicity classification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Classifies each SV into ACMG class 1-5 and adds gene/disease annotations.

    The biocontainer holds AnnotSV's code but none of its annotation data, so
    the annotation directory (setup.sh's annotsv_annotations/, which holds
    Annotations_Human/) comes in through --annotsv_annotations. setup.sh adds
    gnomAD-SV's common SVs to its benign set and mane_select_refseq.txt, the
    MANE Select transcripts passed as -txFile; the task stops without them.

    Equivalent to: scripts/05-annotsv.sh, which reads the same duphold-filtered
    calls when step 15 wrote them, and Manta's calls (with a note) otherwise.
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
    for f in ${annotations_dir}/mane_select_refseq.txt \\
             ${annotations_dir}/Annotations_Human/SVincludedInFt/BenignSV/GRCh38/gnomad-sv.installed; do
        if [ ! -s "\$f" ]; then
            echo "ERROR: \$f not found: AnnotSV's data lacks the MANE Select list or gnomAD-SV's common SVs." >&2
            echo "  Run ./scripts/setup.sh --annotsv-data <genome_dir> to add them, or see docs/05-annotsv.md." >&2
            exit 1
        fi
    done

    AnnotSV \\
        -SVinputFile ${sv_vcf} \\
        -outputFile ${meta.id}_sv_annotated.tsv \\
        -outputDir . \\
        -genomeBuild GRCh38 \\
        -annotationsDir ${annotations_dir} \\
        -annotationMode both \\
        -txFile ${annotations_dir}/mane_select_refseq.txt

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
