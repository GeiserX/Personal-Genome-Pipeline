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
    With PGX_CONSENSUS's table it lists the calls other tools gave PharmCAT
    (HLA from T1K, an agreed CYP2D6) and why a gene was held back. When
    PharmCAT has no single result for CYP4F2 or VKORC1, whose warfarin rule
    rests on one position, it adds the genotype PharmCAT read there, labelled
    as a single-position fallback, not a PharmCAT call.
    With --medications (opt-in, a list the user gives) and
    assets/pgx_inhibitors.tsv, a strong or moderate CYP2D6 inhibitor on the
    list gives CYP2D6 a phenoconverted phenotype by CPIC's activity-score
    rule, beside the genotype one; an inhibitor of CYP2C19, CYP2C9, CYP3A4
    or CYP3A5 gets a note.
    With CLINVAR_SCREEN's hits table and assets/pgx_contraindications.tsv, a
    pathogenic or likely pathogenic hit in one of its disease genes (POLG
    and the urea-cycle genes with valproate, HPRT1, SCN1A) adds a WARNING
    naming the drug; without the hits the report says they were not checked.
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
    // pypgx_summary is [] when pypgx is not in --tools; consensus is [] when
    // PGX_CONSENSUS did not run for the sample; clinvar_hits is [] when
    // CLINVAR_SCREEN did not run for it
    tuple val(meta), path(pharmcat_json), path(pypgx_summary), path(consensus), path(clinvar_hits)
    // medications is [] without --medications; inhibitors is assets/pgx_inhibitors.tsv;
    // contraindications is assets/pgx_contraindications.tsv
    path(medications)
    path(inhibitors)
    path(contraindications)

    output:
    tuple val(meta), path("${meta.id}_cpic_recommendations.txt"), emit: recommendations
    tuple val(meta), path("${meta.id}_phenotypes.tsv"),           emit: phenotypes
    tuple val(meta), path("${meta.id}_pharmcat_comparison.tsv"),  emit: comparison, optional: true
    path "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def pypgx_args = pypgx_summary ? "--pypgx ${pypgx_summary} --comparison ${meta.id}_pharmcat_comparison.tsv" : ''
    def consensus_args = consensus ? "--consensus ${consensus}" : ''
    def medication_args = medications ? "--medications ${medications} --inhibitors ${inhibitors}" : ''
    def contra_args = contraindications ? "--contraindications ${contraindications}" : ''
    def hits_args = contraindications && clinvar_hits ? "--clinvar-hits ${clinvar_hits}" : ''
    """
    pgx_parse.py cpic-report \\
        --sample ${meta.id} \\
        --report ${pharmcat_json} \\
        --outdir . \\
        ${pypgx_args} \\
        ${consensus_args} \\
        ${medication_args} \\
        ${contra_args} \\
        ${hits_args}

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
