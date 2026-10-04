/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PharmCAT — Clinical pharmacogenomics (star alleles + drug recommendations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Two-step process:
    1. Preprocess VCF (normalize, filter to PGx positions)
    2. Run PharmCAT (star allele calling + drug recommendation reports)

    PharmCAT reads a PGx position missing from its input as "not covered",
    and a variants-only VCF lists only where the sample differs from the
    reference. When the sample has a gVCF (DEEPVARIANT wrote one), step 1
    first expands its reference blocks over PharmCAT's gene regions into a
    plain VCF, as scripts/07-pharmacogenomics.sh does: a covered position
    becomes a 0/0 call, an uncovered one (./.) stays missing. The expanded
    file is named without .g.vcf, which PharmCAT refuses.

    Equivalent to: scripts/07-pharmacogenomics.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PHARMCAT_PREPROCESS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(vcf_index), path(gvcf), path(gvcf_index)  // gVCF pair or []
    path(reference)
    path(reference_fai)  // staged beside the FASTA, so no task builds its own

    output:
    tuple val(meta), path("*.preprocessed.vcf.bgz"), emit: preprocessed_vcf
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def pgx_input = gvcf ? "${meta.id}.pgx_regions.vcf.gz" : "${vcf}"
    def expand = gvcf ? """
    echo "Input: ${gvcf} (reference blocks expanded over PharmCAT's gene regions)"
    bcftools convert --gvcf2vcf -f ${reference} -R /pharmcat/pharmcat_regions.bed -Ou ${gvcf} \\
        | bcftools view --trim-alt-alleles -i 'GT!="mis"' -Oz -o ${meta.id}.pgx_regions.vcf.gz --write-index=tbi
    """ : """
    echo "Input: ${vcf} (no gVCF: PGx positions where the sample matches the reference read as missing)"
    """
    """
    ${expand}
    python3 /pharmcat/pharmcat_vcf_preprocessor \\
        -vcf ${pgx_input} \\
        -refFna ${reference} \\
        -o ./ \\
        -bf ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.preprocessed.vcf.bgz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process PHARMCAT {
    tag "$meta.id"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/pharmcat" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(preprocessed_vcf)

    output:
    tuple val(meta), path("*.report.html"),  emit: html_report
    tuple val(meta), path("*.report.json"),  emit: json_report
    tuple val(meta), path("*.match.json"),   emit: match_json, optional: true
    tuple val(meta), path("*.phenotype.json"), emit: phenotype_json, optional: true
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // PharmCAT up to 3.4.0 bundles vcf-parser 0.3.1, which stops with "Error
    // parsing metadata: character to be escaped is missing" on a backslash in
    // a ## header line. Such lines are valid VCF: bcftools writes one for a
    // soft filter with a quoted string (-s LowDP -e 'GT!="0/0"'). The awk
    // below rewrites PharmCAT's own copy only, on ## lines: \" becomes ' and
    // any other \ becomes /. Remove it once a PharmCAT release bundles
    // vcf-parser newer than 0.3.1 (scripts/07-pharmacogenomics.sh does the same).
    """
    gzip -dc ${preprocessed_vcf} \\
        | awk '/^##/ { gsub(/\\\\"/, "\\047"); gsub(/\\\\/, "/") } { print }' \\
        > ${meta.id}.pharmcat_input.vcf

    java -jar /pharmcat/pharmcat.jar \\
        -vcf ${meta.id}.pharmcat_input.vcf \\
        -o ./ \\
        -bf ${meta.id} \\
        -reporterJson \\
        -reporterHtml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.report.html ${meta.id}.report.json
    touch ${meta.id}.match.json ${meta.id}.phenotype.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
