/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PharmCAT — Clinical pharmacogenomics (star alleles + drug recommendations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Two-step process:
    1. Preprocess VCF (normalize, filter to PGx positions)
    2. Run PharmCAT (star allele calling + drug recommendation reports)

    Equivalent to: scripts/07-pharmacogenomics.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PHARMCAT_PREPROCESS {
    tag "$meta.id"
    label 'process_low'

    container 'pgkb/pharmcat:3.2.0'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(reference)

    output:
    tuple val(meta), path("*.preprocessed.vcf.bgz"), emit: preprocessed_vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 /pharmcat/pharmcat_vcf_preprocessor \\
        -vcf ${vcf} \\
        -refFna ${reference} \\
        -o ./ \\
        -bf ${meta.id}
    """

    stub:
    """
    touch ${meta.id}.preprocessed.vcf.bgz
    """
}

process PHARMCAT {
    tag "$meta.id"
    label 'process_low'

    container 'pgkb/pharmcat:3.2.0'

    publishDir "${params.outdir}/${meta.id}/pharmcat", mode: params.publish_dir_mode

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
    """
    java -jar /pharmcat/pharmcat.jar \\
        -vcf ${preprocessed_vcf} \\
        -o ./ \\
        -bf ${meta.id} \\
        -reporterJson \\
        -reporterHtml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: \$(java -jar /pharmcat/pharmcat.jar -version 2>&1 | grep -oP '[\\d.]+' | head -1 || echo '3.2.0')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.report.html ${meta.id}.report.json
    touch ${meta.id}.match.json ${meta.id}.phenotype.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: 3.2.0
    END_VERSIONS
    """
}
