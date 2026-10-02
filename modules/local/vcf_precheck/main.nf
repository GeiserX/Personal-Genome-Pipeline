/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VCF_PRECHECK — Count FILTER values once per sample before any analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ClinVar screen, clinical filter and slivar keep only FILTER=PASS records.
    A VCF from a caller that leaves FILTER as '.' (unfiltered GATK
    HaplotypeCaller, FreeBayes) would give zero records everywhere and a
    report with zero hits.

    Status emitted per sample:
      pass        — at least one PASS record; the VCF is used as given
      no_pass     — no PASS record; main.nf stops the run naming the sample
                    (also with --allow_unfiltered when no record has FILTER '.',
                    since the relaxed copy would still have no PASS record)
      unfiltered  — no PASS record, at least one FILTER '.' record, and
                    --allow_unfiltered is set: a copy with
                    FILTER '.' rewritten to PASS is emitted, which is what
                    `bcftools view -f .,PASS` would select (records with any
                    other FILTER value stay excluded)

    The fallback applies only when the file has no PASS record at all, and
    the log says so. A file with some PASS records is never relaxed.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process VCF_PRECHECK {
    tag "$meta.id"
    label 'process_single'

    container 'staphb/bcftools:1.21'

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), env(FILTER_STATUS), env(FILTER_COUNTS),          emit: status
    tuple val(meta), path("${meta.id}.unfiltered_as_pass.vcf.gz"),
                     path("${meta.id}.unfiltered_as_pass.vcf.gz.tbi"), emit: relaxed, optional: true
    path "versions.yml",                                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def allow_unfiltered = params.allow_unfiltered ? 'true' : 'false'
    """
    bcftools query -f '%FILTER\\n' ${vcf} \\
        | awk '{ if (\$1 == "PASS") p++; else if (\$1 == ".") d++; else o++ } END { printf "%d %d %d\\n", p, d, o }' \\
        > filter_counts.txt
    read -r N_PASS N_DOT N_OTHER < filter_counts.txt
    FILTER_COUNTS="PASS=\${N_PASS} .=\${N_DOT} other=\${N_OTHER}"
    echo "${meta.id}: FILTER counts \${FILTER_COUNTS}"

    if [ "\${N_PASS}" -gt 0 ]; then
        FILTER_STATUS=pass
    elif [ "${allow_unfiltered}" = "true" ] && [ "\${N_DOT}" -gt 0 ]; then
        FILTER_STATUS=unfiltered
        echo "WARNING: ${meta.id} has no PASS record; --allow_unfiltered set, treating FILTER '.' as PASS" >&2
        bcftools view ${vcf} \\
            | awk -F'\\t' -v OFS='\\t' '/^#/ { print; next } \$7 == "." { \$7 = "PASS" } { print }' \\
            | bcftools view -Oz -o ${meta.id}.unfiltered_as_pass.vcf.gz
        bcftools index -t ${meta.id}.unfiltered_as_pass.vcf.gz
    else
        FILTER_STATUS=no_pass
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    FILTER_STATUS=pass
    FILTER_COUNTS="stub"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.21
    END_VERSIONS
    """
}
