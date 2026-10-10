/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CPSR — Cancer Predisposition Sequencing Reporter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Screens germline VCF for cancer predisposition variants using PCGR 2.x.
    Requires the VEP cache of PCGR_VEP_CACHE_RELEASE in versions.env (115 for
    PCGR 2.3.2), separate from the VEP step's VEP_CACHE_RELEASE cache. CPSR 2.3
    has no --classify_all: it classifies every panel variant by itself.

    CPSR stops on a --sample_id outside 3 to 40 characters (PCGR 2.3.2), so
    it gets bin/cpsr_sample_id's id (S1 -> S1_cpsr, a longer id cut to 40),
    and its files are renamed back to <meta.id>.cpsr.*. CPSR's own HTML
    title then shows that id.

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
    cpsr_id=\$(cpsr_sample_id '${meta.id}')
    cpsr \\
        --input_vcf ${vcf} \\
        --vep_dir ${vep_cache_cpsr} \\
        --refdata_dir ${pcgr_data} \\
        --output_dir ./ \\
        --genome_assembly grch38 \\
        --sample_id "\$cpsr_id" \\
        --panel_id 0 \\
        --secondary_findings \\
        --force_overwrite
    cpsr_sample_id --rename . '${meta.id}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cpsr: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    cpsr_id=\$(cpsr_sample_id '${meta.id}')
    touch "\$cpsr_id".cpsr.grch38.html "\$cpsr_id".cpsr.grch38.classification.tsv.gz
    cpsr_sample_id --rename . '${meta.id}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cpsr: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
