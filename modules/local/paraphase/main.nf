/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAPHASE_CHECK, PARAPHASE, PARAPHASE_REPORT — paralog genes from PacBio
    HiFi reads (opt-in: --tools paraphase)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Near-identical gene copies (SMN1/SMN2, HBA1/HBA2, CYP21A2, GBA1, PMS2,
    STRC, NCF1, OTOA and the rest of Paraphase's GRCh38 regions) hide
    copy-number changes and variants from a normal alignment. Paraphase
    phases the reads of each gene family into haplotypes and calls each
    gene's copy number, the HBA 3.7 kb and 4.2 kb deletions by name, and
    phased variants.

    PacBio HiFi only: Paraphase is written for HiFi reads, and its licence
    (PacBio Software License Agreement) allows it only on data from a PacBio
    instrument. PARAPHASE_CHECK reads the BAM header with
    bin/paraphase_report.py check-bam and fails unless every read group says
    PL:PACBIO, so Paraphase never reads another BAM. All three are
    report-only (nextflow.config): a failed check is listed, the step is
    marked failed, and the rest of the run goes on.

    PARAPHASE_REPORT turns Paraphase's JSON and VCFs into one row per gene
    and one per variant (bin/paraphase_report.py table).

    Equivalent to: scripts/35b-paraphase.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PARAPHASE_CHECK {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_paraphase_check.txt"), emit: ok
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    paraphase_report.py check-bam ${bam} > ${meta.id}_paraphase_check.txt

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    echo "stub: ${bam}" > ${meta.id}_paraphase_check.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process PARAPHASE {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/paralogs" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_fai)

    output:
    tuple val(meta), path("${meta.id}_paraphase"), emit: results
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def regions = params.paraphase_regions ? "-g ${params.paraphase_regions}" : ''
    """
    paraphase -b ${bam} -r ${reference} -o ${meta.id}_paraphase -p ${meta.id} -t ${task.cpus} ${regions}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paraphase: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}_paraphase/${meta.id}_paraphase_vcfs
    echo '{"smn1": {"total_cn": 4, "region_specific_info": {"smn1_cn": 2, "smn2_cn": 2, "smn_del78_cn": 0}}}' \\
        > ${meta.id}_paraphase/${meta.id}.paraphase.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paraphase: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process PARAPHASE_REPORT {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/paralogs" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(results)

    output:
    tuple val(meta), path("${meta.id}_paraphase_copy_number.tsv"), emit: copy_number
    tuple val(meta), path("${meta.id}_paraphase_variants.tsv"),    emit: variants
    path "versions.yml",                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    paraphase_report.py table --sample ${meta.id} \\
        --json ${results}/${meta.id}.paraphase.json \\
        --vcf-dir ${results}/${meta.id}_paraphase_vcfs \\
        --out-cn ${meta.id}_paraphase_copy_number.tsv \\
        --out-variants ${meta.id}_paraphase_variants.tsv

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    printf 'region\\tgene\\tcopy_number\\tstatus\\tdetail\\ttotal_cn\\tregion_depth\\tgenome_depth\\n' > ${meta.id}_paraphase_copy_number.tsv
    printf 'region\\tchrom\\tpos\\tref\\talt\\tfilter\\thaplotypes\\tsupport\\n' > ${meta.id}_paraphase_variants.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
