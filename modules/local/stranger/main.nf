/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STRANGER — Clinical STR pathogenicity annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Annotates an ExpansionHunter VCF with STR_STATUS (normal / pre_mutation /
    full_mutation), disease name, OMIM number, and repeat-size thresholds from
    the bundled ClinGen/OMIM catalog.

    Equivalent to: scripts/09b-stranger.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process STRANGER {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/stranger:0.10.2--pyhdfd78af_0'

    publishDir { "${params.outdir}/${meta.id}/expansion_hunter" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_eh_stranger.vcf"), emit: vcf
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    stranger ${vcf} > ${prefix}_eh_stranger.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stranger: \$(stranger --version 2>&1 | grep -oP '[\\d.]+' | head -1 || echo '0.10.2')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_eh_stranger.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stranger: 0.10.2
    END_VERSIONS
    """
}
