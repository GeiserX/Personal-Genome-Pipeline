/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Y_STR — Y-STR repeat counts (FamilyTreeDNA's Y-111 markers) with HipSTR
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Opt-in ('y_str' in --tools). Runs on male samples only:
    workflows/bam_analysis.nf passes the BAMs whose sex is male, the sex
    INDEXCOV has checked against the reads.

    One task: bin/ystr_report.py writes HipSTR's region file from
    assets/ystr_panel.tsv, HipSTR genotypes those chrY regions, and the
    same script, with the image's Python, writes the table in
    FamilyTreeDNA's order and its CSV. Multi-copy markers and repeats
    longer than the reads can span are listed as not callable.

    Equivalent to: scripts/45-y-str.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process Y_STR {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/y_str" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_fai)
    path(panel)        // assets/ystr_panel.tsv

    output:
    tuple val(meta), path("${meta.id}_y_str.tsv"), emit: table
    tuple val(meta), path("${meta.id}_y_str.csv"), emit: csv
    tuple val(meta), path("${meta.id}_y_str.vcf.gz"), path("${meta.id}_y_str.log"), emit: calls
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ystr_report.py bed --panel ${panel} --out regions.bed
    HipSTR --bams ${bam} --fasta ${reference} --regions regions.bed --str-vcf ${meta.id}_y_str.vcf.gz \\
        --bam-samps ${meta.id} --bam-libs ${meta.id} --haploid-chrs chrY --def-stutter-model --min-reads 5 \\
        --max-str-len 150 --output-filters --log ${meta.id}_y_str.log
    ystr_report.py report --panel ${panel} --vcf ${meta.id}_y_str.vcf.gz \\
        --tsv ${meta.id}_y_str.tsv --csv ${meta.id}_y_str.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hipstr: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf 'order\\tpanel\\tmarker\\tcount\\tconfidence\\treads\\tstutter_reads\\tflank_indel_reads\\tq\\tnote\\n1\\tY-12\\tDYS393\\t\\tno call\\t\\t\\t\\t\\tstub\\n' > ${meta.id}_y_str.tsv
    printf 'marker,value,confidence\\nDYS393,,no call\\n' > ${meta.id}_y_str.csv
    touch ${meta.id}_y_str.vcf.gz ${meta.id}_y_str.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hipstr: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
