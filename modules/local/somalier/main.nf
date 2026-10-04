/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOMALIER, SOMALIER_RELATE, SAMPLE_QC — Is the sample the person you think?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOMALIER reads each BAM's genotypes and depth at somalier's known sites
    (--somalier_sites). It names the sample after the BAM's @RG SM tag, which
    <id>.somalier_id holds.

    SOMALIER_RELATE runs once over every sample of the run: per sample, the
    sex somalier infers from heterozygosity at chrX sites and chrY depth; per
    pair, relatedness, so two rows that are the same person show up.

    SAMPLE_QC writes <id>_sample_qc.tsv with bin/collect_summary.py's
    sample-qc command (the rule scripts/33-sample-qc.sh runs too): somalier's
    sex against the samplesheet's, VERIFYBAMID2's FREEMIX against
    --freemix_warn, and the other samples somalier finds to be the same
    person. workflows/bam_analysis.nf stops the run on a sex mismatch
    (--sex_check warn logs it), the check INDEXCOV makes from the BAM index,
    made again from the reads.

    Equivalent to: scripts/33-sample-qc.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SOMALIER {
    tag "$meta.id"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/qc/somalier" }, mode: params.publish_dir_mode,
        pattern: "*.somalier"

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_fai)
    path(sites)

    output:
    tuple val(meta), path("${meta.id}.somalier"), path("${meta.id}.somalier_id"), emit: extract
    path "versions.yml",                                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    somalier extract -d extract --sites ${sites} -f ${reference} ${bam}
    n=\$(find extract -name '*.somalier' | wc -l)
    if [ "\$n" -ne 1 ]; then
        echo "ERROR: somalier extract wrote \$n files, expected one" >&2
        exit 1
    fi
    f=\$(find extract -name '*.somalier')
    # relate reports the sample by the name inside the file (the @RG SM)
    basename "\$f" .somalier > ${meta.id}.somalier_id
    mv "\$f" ${meta.id}.somalier

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.somalier
    echo ${meta.id} > ${meta.id}.somalier_id

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process SOMALIER_RELATE {
    label 'process_single'

    publishDir { "${params.outdir}/somalier" }, mode: params.publish_dir_mode

    input:
    path(extracted, stageAs: 'extract/*')
    path(sites)

    output:
    path "somalier.samples.tsv", emit: samples
    path "somalier.pairs.tsv",   emit: pairs
    path "somalier.html",        emit: html, optional: true
    path "versions.yml",         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # --infer: without it somalier leaves the sex column as the pedigree's (-9)
    somalier relate --infer --sites ${sites} -o somalier extract/*.somalier

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    // One row per sample, every one male; SAMPLE_QC reads the sex column.
    """
    printf '#family_id\\tsample_id\\tpaternal_id\\tmaternal_id\\tsex\\tphenotype\\toriginal_pedigree_sex\\tgt_depth_mean\\tn_hom_ref\\tn_het\\tn_hom_alt\\tX_depth_mean\\tX_n\\tX_het\\tX_hom_alt\\tY_depth_mean\\tY_n\\n' > somalier.samples.tsv
    for f in extract/*.somalier; do
        s=\$(basename "\$f" .somalier)
        printf '%s\\t%s\\t-9\\t-9\\t-9\\t-9\\t-9\\t30\\t1\\t1\\t1\\t15\\t20\\t0\\t20\\t15\\t5\\n' "\$s" "\$s" >> somalier.samples.tsv
    done
    printf '#sample_a\\tsample_b\\trelatedness\\n' > somalier.pairs.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process SAMPLE_QC {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/qc" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(somalier_id), path(selfsm), path(marker_check)
    path(samples_tsv)
    path(pairs_tsv)

    output:
    tuple val(meta), path("${meta.id}_sample_qc.tsv"), emit: table
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    collect_summary.py sample-qc \\
        --sample ${meta.id} \\
        --somalier-id "\$(cat ${somalier_id})" \\
        --somalier-samples ${samples_tsv} \\
        --somalier-pairs ${pairs_tsv} \\
        --selfsm ${selfsm} \\
        --declared-sex "${meta.sex ?: ''}" \\
        --freemix-warn ${params.freemix_warn} \\
        --marker-check "\$(cat ${marker_check})" \\
        --out ${meta.id}_sample_qc.tsv

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    // somalier "cannot tell" here, so the stub run never stops on the sex.
    """
    printf 'key\\tvalue\\nsample\\t%s\\ninferred_sex\\tunknown\\nsex_check\\tnot_checked\\nsex_check_reason\\tstub\\nfreemix\\t0.0010\\nfreemix_warn_above\\t%s\\ncontamination\\tok\\nsame_person_as\\t\\n' \\
        ${meta.id} ${params.freemix_warn} > ${meta.id}_sample_qc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
