/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CYRIUS — CYP2D6 star allele calling from WGS BAM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CYP2D6 is the hardest pharmacogene to call because of its pseudogene (CYP2D7)
    and complex structural variants (deletions, duplications, hybrids).
    Cyrius uses depth-based analysis specifically designed for CYP2D6.

    NOTE: Cyrius is installed at runtime via pip (pinned to 1.1.1) because no
    pre-built container image exists. This requires network access on first run.
    The tool may return "None" for complex CYP2D6 arrangements. Verify results
    against PharmCAT or clinical lab calls before acting on them.

    Equivalent to: scripts/21-cyrius.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CYRIUS {
    tag "$meta.id"
    label 'process_low'

    container 'python:3.11'

    publishDir "${params.outdir}/${meta.id}/cyrius", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*_cyp2d6.tsv"), emit: cyp2d6_results
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cyrius_version = '1.1.1'
    """
    pip install -q 'cyrius==${cyrius_version}'

    echo "${bam}" > manifest.txt

    cyrius \\
        --manifest manifest.txt \\
        --genome 38 \\
        --prefix ${prefix}_cyp2d6 \\
        --outDir ./ \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cyrius: ${cyrius_version}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'Sample\\tGenotype\\tFilter\\n${prefix}\\t*1/*1\\tPASS\\n' > ${prefix}_cyp2d6.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cyrius: unknown
    END_VERSIONS
    """
}
