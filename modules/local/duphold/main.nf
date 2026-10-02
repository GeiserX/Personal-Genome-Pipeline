/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DUPHOLD — Annotate structural variants with depth-based quality metrics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Adds per-sample FORMAT fields DHFC (depth fold-change against the
    chromosome), DHBFC (against GC-matched bins) and DHFFC (against the
    flanking regions) to an SV VCF.

    Two processes:
      1. DUPHOLD         — duphold image, annotation only    -> sv_duphold/
      2. DUPHOLD_FILTER  — bcftools image, drops DELs with DHFFC >= 0.7 and
                           DUPs with DHBFC <= 1.3 (duphold's recommended
                           cut-offs); other types and records without a value
                           are kept                           -> sv_filtered/

    Equivalent to: scripts/15-duphold.sh (which prints the filter commands
    instead of applying them)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DUPHOLD {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/sv_duphold" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(sv_vcf), path(bam), path(bai)
    path(reference)
    path(reference_fai)

    output:
    tuple val(meta), path("${meta.id}_sv_duphold.vcf"), emit: annotated_vcf
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    duphold \\
        -v ${sv_vcf} \\
        -b ${bam} \\
        -f ${reference} \\
        -o ${meta.id}_sv_duphold.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        duphold: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_sv_duphold.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        duphold: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process DUPHOLD_FILTER {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/sv_filtered" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(annotated_vcf)

    output:
    tuple val(meta), path("${meta.id}_sv_filtered.vcf.gz"),     emit: filtered_vcf
    tuple val(meta), path("${meta.id}_sv_filtered.vcf.gz.tbi"), emit: filtered_vcf_index
    tuple val(meta), path("${meta.id}_sv_filtered.log"),        emit: log
    path "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Deletions keep only when the depth drop against the flanks is real
    # (DHFFC < 0.7); duplications only when the gain against GC-matched bins
    # is real (DHBFC > 1.3). Missing values compare false, so those records stay.
    bcftools view \\
        -e '(INFO/SVTYPE="DEL" && FMT/DHFFC[0] >= 0.7) || (INFO/SVTYPE="DUP" && FMT/DHBFC[0] <= 1.3)' \\
        ${annotated_vcf} -Oz -o ${meta.id}_sv_filtered.vcf.gz
    bcftools index -t ${meta.id}_sv_filtered.vcf.gz

    N_IN=\$(bcftools view -H ${annotated_vcf} | wc -l)
    N_OUT=\$(bcftools view -H ${meta.id}_sv_filtered.vcf.gz | wc -l)
    N_DEL=\$(bcftools view -H -i 'INFO/SVTYPE="DEL" && FMT/DHFFC[0] >= 0.7' ${annotated_vcf} | wc -l)
    N_DUP=\$(bcftools view -H -i 'INFO/SVTYPE="DUP" && FMT/DHBFC[0] <= 1.3' ${annotated_vcf} | wc -l)
    {
        echo "duphold filter for ${meta.id}: kept \${N_OUT} of \${N_IN} records"
        echo "  removed \${N_DEL} DEL with DHFFC >= 0.7 and \${N_DUP} DUP with DHBFC <= 1.3"
        if [ "\${N_OUT}" -eq "\${N_IN}" ]; then
            echo "  nothing removed: no DEL or DUP failed the depth check"
        fi
    } | tee ${meta.id}_sv_filtered.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_sv_filtered.vcf.gz ${meta.id}_sv_filtered.vcf.gz.tbi ${meta.id}_sv_filtered.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
