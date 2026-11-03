/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MOSDEPTH — Fast per-base and per-region coverage statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Computes precise per-base depth from actual alignments (not just the index).
    Outputs coverage distributions, threshold summaries, and region BED files.

    A second pass over the genes of assets/callability_genes.tsv counts, per
    gene on chr1-22, X and Y, the bases with 10 or more reads of MAPQ 20 or
    more: ${id}_callability.bed.gz, which the report reads to list the genes
    this genome could not read reliably (bin/collect_summary.py
    sec_callability).

    Equivalent to: scripts/16b-mosdepth.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MOSDEPTH {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/coverage" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(callability_genes)   // assets/callability_genes.tsv

    output:
    tuple val(meta), path("*.mosdepth.summary.txt"),     emit: summary
    tuple val(meta), path("*.mosdepth.global.dist.txt"), emit: global_dist
    tuple val(meta), path("*.regions.bed.gz"),           emit: regions_bed, optional: true
    tuple val(meta), path("*.thresholds.bed.gz"),        emit: thresholds_bed, optional: true
    tuple val(meta), path("*_callability.bed.gz"),       emit: callability
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

    # Gene callability, as scripts/16b-mosdepth.sh: the genes with a position
    # (not the ALT-only ones), reads of MAPQ 20 or more, bases at 10x or more.
    # It writes into callability/, so its own summary and region files never
    # match the globs of the outputs above.
    mkdir -p callability
    awk -F'\\t' -v OFS='\\t' '!/^#/ && \$1 != "gene" && \$4 ~ /^[0-9]+\$/ {print \$3, \$4, \$5, \$1}' \\
        ${callability_genes} > callability/genes.bed
    mosdepth \\
        --by callability/genes.bed \\
        --mapq 20 \\
        --fast-mode \\
        --no-per-base \\
        --threads ${task.cpus} \\
        --thresholds 10 \\
        callability/c \\
        ${bam}
    mv callability/c.thresholds.bed.gz ${prefix}_callability.bed.gz
    rm -rf callability

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mosdepth.summary.txt
    touch ${prefix}.mosdepth.global.dist.txt
    echo | gzip > ${prefix}_callability.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
