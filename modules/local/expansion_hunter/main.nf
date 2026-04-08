/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXPANSION_HUNTER — Short tandem repeat (STR) expansion detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Screens for pathogenic repeat expansions: Huntington's, Fragile X,
    Friedreich's ataxia, ALS/FTD, SCAs, myotonic dystrophy, etc.
    Uses the bundled GRCh38 variant catalog (31 pathogenic loci).

    Equivalent to: scripts/09-expansion-hunter.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process EXPANSION_HUNTER {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/expansionhunter:5.0.0--hc26b3af_5'

    publishDir "${params.outdir}/${meta.id}/expansion_hunter", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_fai)
    path(variant_catalog)

    output:
    tuple val(meta), path("*.vcf"),    emit: vcf
    tuple val(meta), path("*.json"),   emit: json
    tuple val(meta), path("*_realigned.bam"), emit: bamlet, optional: true
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sex_flag = meta.sex ? "--sex ${meta.sex}" : ""
    """
    ExpansionHunter \\
        --reads ${bam} \\
        --reference ${reference} \\
        --variant-catalog ${variant_catalog} \\
        --output-prefix ${prefix}_eh \\
        --threads ${task.cpus} \\
        ${sex_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunter: \$(ExpansionHunter --version 2>&1 | grep -oP '[\\d.]+' | head -1 || echo '5.0.0')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_eh.vcf ${prefix}_eh.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunter: 5.0.0
    END_VERSIONS
    """
}
