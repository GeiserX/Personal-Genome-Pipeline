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
    label 'process_high'

    container 'lgalarno/telomerehunter:latest'

    publishDir "${params.outdir}/${meta.id}/telomere", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}"), emit: telomere_results
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    telomerehunter \\
        -ibt ${bam} \\
        -o ./ \\
        -p ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomerehunter: \$(telomerehunter --version 2>&1 | grep -oP '[\\d.]+' | head -1 || echo 'unknown')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}
    touch ${meta.id}/${meta.id}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomerehunter: unknown
    END_VERSIONS
    """
}
