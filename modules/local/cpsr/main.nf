/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CPSR — Cancer Predisposition Sequencing Reporter (ACMG SF v3.2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Screens germline VCF for cancer predisposition variants using PCGR 2.x.
    Requires the VEP 113 cache (separate from the VEP step's release 116 cache).

    Equivalent to: scripts/17-cpsr.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CPSR {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/cpsr" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(pcgr_data)
    path(vep_cache_cpsr)

    output:
    tuple val(meta), path("*.cpsr.grch38.html"),                  emit: html_report
    tuple val(meta), path("*.cpsr.grch38.classification.tsv.gz"), emit: tsv_report
    path "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cpsr \\
        --input_vcf ${vcf} \\
        --vep_dir ${vep_cache_cpsr} \\
        --refdata_dir ${pcgr_data} \\
        --output_dir ./ \\
        --genome_assembly grch38 \\
        --sample_id ${meta.id} \\
        --panel_id 0 \\
        --classify_all \\
        --secondary_findings \\
        --force_overwrite

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cpsr: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.cpsr.grch38.html
    touch ${meta.id}.cpsr.grch38.classification.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cpsr: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
