/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CYRIUS — CYP2D6 star allele calling from WGS BAM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CYP2D6 is the hardest pharmacogene to call because of its pseudogene (CYP2D7)
    and complex structural variants (deletions, duplications, hybrids).
    Cyrius uses depth-based analysis specifically designed for CYP2D6.

    EXPERIMENTAL: Cyrius is installed at runtime via pip and may return "None"
    for complex CYP2D6 arrangements. Verify results against PharmCAT or
    clinical lab calls before acting on them.

    Equivalent to: scripts/21-cyrius.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CYRIUS {
    tag "$meta.id"
    label 'process_low'

    container 'python:3.11-slim'

    publishDir "${params.outdir}/${meta.id}/cyrius", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)

    output:
    tuple val(meta), path("*_cyp2d6.tsv"), emit: cyp2d6_results
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pip install -q cyrius 2>/dev/null

    echo '${bam}' > manifest.txt

    star_caller \\
        --manifest manifest.txt \\
        --genome 38 \\
        --prefix ${prefix}_cyp2d6 \\
        --outDir ./ \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cyrius: \$(pip show cyrius 2>/dev/null | grep Version | sed 's/Version: //' || echo 'unknown')
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
