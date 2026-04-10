/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MULTIQC — Aggregate QC reports into a single HTML dashboard
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MultiQC auto-discovers supported tool outputs (fastp JSON, mosdepth,
    samtools flagstat/stats, etc.) and produces a combined interactive HTML report.
    Runs once across ALL samples — does not use per-sample meta maps.

    Equivalent to: scripts/28-multiqc.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MULTIQC {
    label 'process_low'

    container 'quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0'

    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode

    input:
    path(multiqc_files, stageAs: "?/*")

    output:
    path "multiqc_report.html",  emit: report
    path "multiqc_report_data",  emit: data
    path "versions.yml",         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    multiqc \\
        . \\
        -f \\
        -o . \\
        -n "multiqc_report.html" \\
        --title "Personal Genome Pipeline QC"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version 2>&1 | sed 's/.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p multiqc_report_data
    touch multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version 2>&1 | sed 's/.*version //')
    END_VERSIONS
    """
}
