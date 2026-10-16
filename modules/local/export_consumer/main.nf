/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONSUMER_EXPORT — a consumer raw-data file (23andMe v5 layout) from the gVCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Opt-in ('consumer_export' in --tools, with --consumer_sites, the site
    table `scripts/setup.sh --consumer-sites` builds). bin/export_consumer.py
    genotypes every site of the 23andMe v5 chip list from the gVCF: a site
    inside a reference block is homozygous reference, a site with no call is
    '--'. It writes the 23andMe layout and the same sites in the AncestryDNA
    layout, both on GRCh37 positions and zipped, reads both back with its
    strict reader, and writes a summary table for the report. The samplesheet
    sex sets the ploidy rules of chrX and chrY (a male's non-pseudoautosomal
    sites get one letter; a female's chrY sites are '--').

    It writes files only: uploading one to a service is the user's decision
    (docs/40-export-consumer.md).

    Equivalent to: scripts/40-export-consumer.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CONSUMER_EXPORT {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/consumer_export" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(gvcf), path(gvcf_index)
    path sites
    path par_bed

    output:
    tuple val(meta), path("${meta.id}_23andme_v5.zip"),         emit: file_23andme
    tuple val(meta), path("${meta.id}_ancestrydna_layout.zip"), emit: file_ancestrydna
    tuple val(meta), path("${meta.id}_consumer_export.tsv"),    emit: summary
    path "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export_consumer.py export \\
        --gvcf ${gvcf} \\
        --sites ${sites} \\
        --sample ${meta.id} \\
        --sex "${meta.sex ?: ''}" \\
        --par ${par_bed} \\
        --out-dir .

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    printf 'key\\tvalue\\nsites\\t0\\ncalled\\t0\\nno_call\\t0\\n' > ${meta.id}_consumer_export.tsv
    touch ${meta.id}_23andme_v5.zip ${meta.id}_ancestrydna_layout.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
