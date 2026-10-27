/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHIP — Clonal haematopoiesis of step 49, from the BAM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHIP_PREPARE   the MANE Select transcripts of assets/chip_genes.tsv from
                   Ensembl's GFF3, and their coding exons (10 bases wider each
                   side) on the reference's contigs as a BED (once per run)
    CHIP_MUTECT2   Mutect2 tumour-only on that BED, the read orientation
                   model and FilterMutectCalls; gnomAD's germline resource
                   and a panel of normals when given
    CHIP_ANNOTATE  bcftools csq names each change on its transcript, and
                   bcftools mpileup counts the depth over the BED
    CHIP_REPORT    <id>_chip.tsv (each protein change with its status) and
                   <id>_chip_depth.tsv (depth and detection limit per gene),
                   bin/chip_report.py

    Runs for the samples with a BAM only; main.nf logs the samples it skips.

    Equivalent to: scripts/49-chip.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CHIP_PREPARE {
    label 'process_single'

    input:
    path(gff)
    path(genes)
    path(reference_fai)

    output:
    path("chip.gff3"),   emit: gff
    path("chip.bed"),    emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    chip_report.py prepare --gff ${gff} --genes ${genes} --fai ${reference_fai} --out-gff chip.gff3 --out-bed chip.bed
    if [ ! -s chip.bed ]; then
        echo "ERROR: none of the CHIP genes is on a contig of the reference." >&2
        exit 1
    fi

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    touch chip.gff3
    printf 'chr2\\t25234268\\t25234430\\tDNMT3A\\n' > chip.bed
    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}

process CHIP_MUTECT2 {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/chip" }, mode: params.publish_dir_mode, pattern: '*_chip_filtered.vcf.gz*'

    input:
    tuple val(meta), path(bam), path(bai)
    path(bed)
    path(reference)
    path(reference_fai)
    path(reference_dict)
    path(germline)        // [] without --chip_germline_resource
    path(germline_index)
    path(pon)             // [] without --chip_pon
    path(pon_index)

    output:
    tuple val(meta), path("${meta.id}_chip_filtered.vcf.gz"), emit: vcf
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The same commands as scripts/49-chip.sh
    def resources = (germline ? "--germline-resource ${germline} " : '') + (pon ? "-pon ${pon}" : '')
    """
    gatk Mutect2 -R ${reference} -I ${bam} -L ${bed} -O unfiltered.vcf.gz --f1r2-tar-gz f1r2.tar.gz \\
        --native-pair-hmm-threads ${task.cpus} --max-mnp-distance 0 ${resources}
    gatk LearnReadOrientationModel -I f1r2.tar.gz -O orientation.tar.gz
    gatk FilterMutectCalls -R ${reference} -V unfiltered.vcf.gz --ob-priors orientation.tar.gz \\
        -O ${meta.id}_chip_filtered.vcf.gz

    printf '"%s":\\n    gatk: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    touch ${meta.id}_chip_filtered.vcf.gz
    printf '"%s":\\n    gatk: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}

process CHIP_ANNOTATE {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(vcf), path(bam), path(bai)
    path(gff)
    path(bed)
    path(reference)
    path(reference_fai)

    output:
    tuple val(meta), path("${meta.id}_chip_calls.tsv"), path("${meta.id}_chip_pileup.tsv"), emit: tables
    path "versions.yml",                                                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The same commands as scripts/49-chip.sh
    """
    bcftools norm -m -any -f ${reference} -Ou ${vcf} \\
        | bcftools csq -v 0 -f ${reference} -g ${gff} -s - --local-csq -Ou \\
        | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%INFO/BCSQ\\t[%AD]\\t[%AF]\\t[%DP]\\n' > ${meta.id}_chip_calls.tsv
    bcftools mpileup -f ${reference} -R ${bed} -q 20 -Q 10 -B -d 100000 -a FORMAT/DP -Ou ${bam} \\
        | bcftools query -f '%CHROM\\t%POS\\t[%DP]\\n' > ${meta.id}_chip_pileup.tsv

    printf '"%s":\\n    bcftools: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    touch ${meta.id}_chip_calls.tsv ${meta.id}_chip_pileup.tsv
    printf '"%s":\\n    bcftools: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}

process CHIP_REPORT {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/chip" }, mode: params.publish_dir_mode, pattern: '*_chip{,_depth}.tsv'

    input:
    tuple val(meta), path(calls), path(pileup)
    path(genes)
    path(gff)
    path(bed)

    output:
    tuple val(meta), path("${meta.id}_chip.tsv"),       emit: chip
    tuple val(meta), path("${meta.id}_chip_depth.tsv"), emit: depth
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    chip_report.py report --genes ${genes} --gff ${gff} --bed ${bed} --calls ${calls} --depth ${pileup} \\
        --sample ${meta.id} --out ${meta.id}_chip.tsv --depth-out ${meta.id}_chip_depth.tsv

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    printf 'gene\\tchrom\\tpos\\tref\\talt\\tconsequence\\tprotein_change\\texon\\tvaf\\talt_reads\\tdepth\\tfilter\\tstatus\\tclass\\tnote\\n' > ${meta.id}_chip.tsv
    printf 'gene\\ttranscript\\tchrom\\tcoding_bases\\tmean_depth\\tdetection_limit_vaf\\tstatus\\n' > ${meta.id}_chip_depth.tsv
    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}
