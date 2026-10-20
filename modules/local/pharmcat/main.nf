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
    becomes a 0/0 call, an uncovered one (./.) stays missing. The expansion
    writes one-base records, which do not count at a PharmCAT position whose
    REF or ALT is longer than one base (an indel such as UGT1A1 *28), so
    pgx_refblock_fill.py (bin/) writes a 0/0 record with PharmCAT's own
    alleles there when reference blocks with GQ >= 20 cover the whole REF
    span and no other record overlaps it. The expanded file is named without
    .g.vcf, which PharmCAT refuses.

    PharmCAT types neither HLA nor CYP2D6 from a VCF. When PGX_CONSENSUS
    wrote outside calls for the sample (HLA-A and HLA-B from T1K, a CYP2D6
    call pypgx and Cyrius agree on), step 2 reads them with -po; an empty
    file is not passed.

    PharmCAT 3.4.0 reads MT-RNR1 only as an outside call: its matcher has
    no MT-RNR1 definition, so the gene is No Result from any VCF. When
    MITO_VARIANTS called the sample's chrM, pgx_mito_merge.py (bin/) adds
    an MT-RNR1 line to a copy of the outside calls: a homoplasmic (AF >=
    0.95) PASS allele PharmCAT names, or Reference when there is none, and
    no line, with the reason in the task log, for a heteroplasmic or
    filtered named allele.

    The jar's built-in logback.xml writes pharmcat.log into the working
    directory. A task runs in its own writable work directory, so the log
    lands there and logback prints nothing; scripts/07-pharmacogenomics.sh
    runs the jar in /tmp for the same reason.

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
    bcftools view -R /pharmcat/pharmcat_regions.bed -Ov -o gvcf_records.vcf ${gvcf}
    bcftools convert --gvcf2vcf -f ${reference} -R /pharmcat/pharmcat_regions.bed -Ou ${gvcf} \\
        | bcftools view --trim-alt-alleles -i 'GT!="mis"' -Ov \\
        | pgx_refblock_fill.py /pharmcat/pharmcat_positions.vcf.bgz gvcf_records.vcf \\
        | bcftools view -Oz -o ${meta.id}.pgx_regions.vcf.gz --write-index=tbi
    rm -f gvcf_records.vcf
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
    // outside_calls is [] when PGX_CONSENSUS did not run for the sample,
    // mito_vcf is [] when MITO_VARIANTS did not
    tuple val(meta), path(preprocessed_vcf), path(outside_calls), path(mito_vcf)

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
    def outside = outside_calls ? "${outside_calls}" : ''
    def mito = mito_vcf ? "${mito_vcf}" : ''
    """
    gzip -dc ${preprocessed_vcf} \\
        | awk '/^##/ { gsub(/\\\\"/, "\\047"); gsub(/\\\\/, "/") } { print }' \\
        > ${meta.id}.pharmcat_input.vcf

    OUTSIDE="${outside}" FROM=PGX_CONSENSUS
    if [ -n "${mito}" ] && [ -s "${mito}" ]; then
        echo "chrM calls (MITO_VARIANTS): ${mito}"
        MERGE=(--mito "${mito}")
        if [ -n "\${OUTSIDE}" ] && [ -s "\${OUTSIDE}" ]; then MERGE+=(--outside "\${OUTSIDE}"); fi
        pgx_mito_merge.py --jar /pharmcat/pharmcat.jar "\${MERGE[@]}" --out ${meta.id}.pharmcat_outside_calls.tsv
        OUTSIDE=${meta.id}.pharmcat_outside_calls.tsv FROM="PGX_CONSENSUS, and MT-RNR1 from MITO_VARIANTS"
    else
        echo "No chrM calls (MITO_VARIANTS did not run): MT-RNR1 stays No Result."
    fi
    PO=()
    if [ -n "\${OUTSIDE}" ] && [ -s "\${OUTSIDE}" ]; then
        echo "Outside calls (\${FROM}):"
        cat "\${OUTSIDE}"
        PO=(-po "\${OUTSIDE}")
    else
        echo "No outside calls: HLA and CYP2D6 stay uncalled."
    fi
    java -jar /pharmcat/pharmcat.jar \\
        -vcf ${meta.id}.pharmcat_input.vcf \\
        \${PO[@]+"\${PO[@]}"} \\
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
