/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TELOMERE_HUNTER — Telomere length estimation from WGS BAM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Estimates telomere length via GC-corrected telomeric reads per million.
    Higher tel_content values = longer telomeres. Provides biological age baseline.

    WARNING: Reads the entire BAM (~30-40 GB). Takes 30-60 minutes.

    Equivalent to: scripts/10-telomere-hunter.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process TELOMERE_HUNTER {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/telomere" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(cytoband)  // UCSC GRCh38 chromosome bands (--cytoband) or []

    output:
    tuple val(meta), path("${meta.id}"), emit: telomere_results
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Without -b TelomereHunter classifies reads by its own hg19 bands;
    // workflows/bam_analysis.nf warns when --cytoband is not set.
    // --plotNone: the Bioconda build's plots fail (its R has no dplyr, and its
    // PyPDF2 is written for Python 3); the summary does not need them.
    def band_arg = cytoband ? "-b ${cytoband}" : ''
    """
    telomerehunter \\
        -ibt ${bam} \\
        -o ./ \\
        -p ${meta.id} \\
        --plotNone \\
        ${band_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomerehunter: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
        telomerehunter_reported: \$( { telomerehunter --version 2>&1 || true; } | awk '!v && match(\$0, /[0-9]+\\.[0-9]+(\\.[0-9]+)*/) { v = substr(\$0, RSTART, RLENGTH) } END { print (v != "" ? v : "unknown") }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}
    touch ${meta.id}/${meta.id}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomerehunter: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
