/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEEPVARIANT — Small variants from the BAM: a VCF and a gVCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    run_deepvariant with the flags of scripts/03-deepvariant.sh: the WGS
    model, the sample name from the samplesheet, one shard per CPU, and a
    gVCF beside the VCF. The gVCF also records where the sample matches the
    reference; PRS and PharmCAT read their hom-ref genotypes from it.

    Sex: for a male sample chrX and chrY are called haploid outside the
    pseudoautosomal regions (assets/par_grch38.bed), so no impossible
    heterozygous call is made there. A female sample is diploid everywhere.
    main.nf requires the sex on every row this process calls, and INDEXCOV
    has checked it against the BAM before this runs.

    --intervals (the script's INTERVALS) limits calling to space-separated
    regions, passed to --regions.

    Equivalent to: scripts/03-deepvariant.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'

    publishDir { "${params.outdir}/${meta.id}/vcf" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_fai)  // staged beside the FASTA
    path(par_bed)

    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), path("${meta.id}.vcf.gz.tbi"),     emit: vcf
    tuple val(meta), path("${meta.id}.g.vcf.gz"), path("${meta.id}.g.vcf.gz.tbi"), emit: gvcf
    tuple val(meta), path("${meta.id}.visual_report.html"),                       emit: report, optional: true
    path "versions.yml",                                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def haploid = meta.sex == 'male' ? "--haploid_contigs=chrX,chrY --par_regions_bed=${par_bed}" : ''
    def regions = params.intervals ? "--regions \"${params.intervals}\"" : ''
    """
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=WGS \\
        --ref=${reference} \\
        --reads=${bam} \\
        --output_vcf=${meta.id}.vcf.gz \\
        --output_gvcf=${meta.id}.g.vcf.gz \\
        --intermediate_results_dir=deepvariant_tmp \\
        --sample_name=${meta.id} \\
        --num_shards=${task.cpus} \\
        ${haploid} \\
        ${regions}
    rm -rf deepvariant_tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > ${meta.id}.vcf.gz
    echo | gzip > ${meta.id}.g.vcf.gz
    touch ${meta.id}.vcf.gz.tbi ${meta.id}.g.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
