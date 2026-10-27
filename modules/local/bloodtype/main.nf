/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BLOODTYPE — Blood group antigens (red cell, platelet, neutrophil) predicted from the genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Opt-in ('bloodtype' in --tools). Two processes, as in scripts/41-bloodtype.sh:
      BLOODTYPE_DEPTH  samtools bedcov -Q 20 -j over assets/bloodtype_markers.tsv,
                       for each sample with a BAM: the read depth at every marker
                       site and over RHD and its flanks
      BLOODTYPE        bin/bloodtype_report.py on the VCF and that depth; without
                       a BAM ([] in its place) RhD is not determined and sites
                       with no VCF record are read as reference (confidence medium)

    Equivalent to: scripts/41-bloodtype.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BLOODTYPE_DEPTH {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/bloodtype" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(markers)   // assets/bloodtype_markers.tsv

    output:
    tuple val(meta), path("${meta.id}_bloodtype_depth.tsv"), emit: depth
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The marker table is a BED with extra columns; bedcov reads it once the
    // comment and header lines are gone. bedcov stops on a contig the BAM has
    // no @SQ line for, so those rows are left out (their sites are then not
    // determined). The same commands as scripts/41-bloodtype.sh.
    """
    samtools view -H ${bam} | awk -F '\\t' '\$1 == "@SQ" {for (i = 2; i <= NF; i++) if (\$i ~ /^SN:/) print substr(\$i, 4)}' > contigs.txt
    grep -v '^#' ${markers} | tail -n +2 | awk -F '\\t' 'NR == FNR {c[\$1] = 1; next} \$1 in c' contigs.txt - > targets.bed
    if [ -s targets.bed ]; then samtools bedcov -Q 20 -j targets.bed ${bam} > ${meta.id}_bloodtype_depth.tsv; else : > ${meta.id}_bloodtype_depth.tsv; fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_bloodtype_depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process BLOODTYPE {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/bloodtype" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(depth)   // depth: BLOODTYPE_DEPTH's, or [] without a BAM
    path(markers)

    output:
    tuple val(meta), path("${meta.id}_bloodtype.tsv"),          emit: table
    tuple val(meta), path("${meta.id}_bloodtype_evidence.tsv"), emit: evidence
    path "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def depth_arg = depth ? "--depth ${depth}" : ''
    """
    bloodtype_report.py --sample ${meta.id} --vcf ${vcf} --markers ${markers} ${depth_arg} \\
        --out ${meta.id}_bloodtype.tsv --evidence ${meta.id}_bloodtype_evidence.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf 'system\\tphenotype\\tgenotype\\tconfidence\\tnote\\nABO\\tnot determined\\t.\\tnot determined\\tstub\\n' > ${meta.id}_bloodtype.tsv
    touch ${meta.id}_bloodtype_evidence.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
