/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CPIC Lookup — Drug-gene recommendations from PharmCAT results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Reads PharmCAT's report.json with bin/pgx_parse.py (on the task PATH; the
    same code scripts/27-cpic-lookup.sh and tests/test_cpic_parser.py run) and
    writes the gene phenotypes and the medications PharmCAT's own report
    matches to them. With the PYPGX summary it also writes the PharmCAT/pypgx
    comparison (published to pypgx/, as the bash steps do) and warns, in the
    recommendations, about a gene PharmCAT could not call but pypgx did.
    Fails when the report yields no gene.

    Equivalent to: scripts/27-cpic-lookup.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CPIC_LOOKUP {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/cpic" }, mode: params.publish_dir_mode,
        pattern: "*_{cpic_recommendations.txt,phenotypes.tsv}"
    publishDir { "${params.outdir}/${meta.id}/pypgx" }, mode: params.publish_dir_mode,
        pattern: "*_pharmcat_comparison.tsv"

    input:
    // pypgx_summary is [] when pypgx is not in --tools
    tuple val(meta), path(pharmcat_json), path(pypgx_summary)

    output:
    tuple val(meta), path("${meta.id}_cpic_recommendations.txt"), emit: recommendations
    tuple val(meta), path("${meta.id}_phenotypes.tsv"),           emit: phenotypes
    tuple val(meta), path("${meta.id}_pharmcat_comparison.tsv"),  emit: comparison, optional: true
    path "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def pypgx_args = pypgx_summary ? "--pypgx ${pypgx_summary} --comparison ${meta.id}_pharmcat_comparison.tsv" : ''
    """
    pgx_parse.py cpic-report \\
        --sample ${meta.id} \\
        --report ${pharmcat_json} \\
        --outdir . \\
        ${pypgx_args}

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    touch ${meta.id}_cpic_recommendations.txt
    printf 'Gene\\tDiplotype\\tPhenotype\\tStatus\\n' > ${meta.id}_phenotypes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
