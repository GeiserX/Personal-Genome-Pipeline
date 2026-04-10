/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MOSDEPTH — Fast per-base and per-region coverage statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Computes precise per-base depth from actual alignments (not just the index).
    Outputs coverage distributions, threshold summaries, and region BED files.

    Equivalent to: scripts/16b-mosdepth.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MOSDEPTH {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/mosdepth:0.3.13--hba6dcaf_0'

    publishDir "${params.outdir}/${meta.id}/coverage", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.mosdepth.summary.txt"),     emit: summary
    tuple val(meta), path("*.mosdepth.global.dist.txt"), emit: global_dist
    tuple val(meta), path("*.regions.bed.gz"),           emit: regions_bed, optional: true
    tuple val(meta), path("*.thresholds.bed.gz"),        emit: thresholds_bed, optional: true
    path "versions.yml",                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mosdepth \\
        --by 500 \\
        --fast-mode \\
        --no-per-base \\
        --threads ${task.cpus} \\
        --thresholds 1,5,10,15,20,30,50 \\
        ${prefix} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mosdepth.summary.txt
    touch ${prefix}.mosdepth.global.dist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: 0.3.13
    END_VERSIONS
    """
}
