/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CPSR — Cancer Predisposition Sequencing Reporter (ACMG SF v3.2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Screens germline VCF for cancer predisposition variants using PCGR 2.x.
    Requires VEP 113 cache (separate from the VEP step's release_112).

    Equivalent to: scripts/17-cpsr.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CPSR {
    tag "$meta.id"
    label 'process_high'

    container 'sigven/pcgr:2.2.5'

    publishDir "${params.outdir}/${meta.id}/cpsr", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(pcgr_data)
    path(vep_cache_cpsr)

    output:
    tuple val(meta), path("*.cpsr.grch38.html"),                   emit: html_report
    tuple val(meta), path("*.cpsr.grch38.snvs_indels.tiers.tsv"), emit: tsv_report
    path "versions.yml",                                           emit: versions

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
        --force_overwrite

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cpsr: \$(cpsr --version 2>&1 | grep -oP '[\\d.]+' | head -1 || echo '2.2.5')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.cpsr.grch38.html
    touch ${meta.id}.cpsr.grch38.snvs_indels.tiers.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cpsr: 2.2.5
    END_VERSIONS
    """
}
