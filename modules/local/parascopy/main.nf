/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARASCOPY — copy number of SMN1 and SMN2 (opt-in: --tools parascopy)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SMN1 and SMN2 are near-identical copies: short reads cannot be placed on
    one of them, and no VCF-based step sees an SMN1 copy-number loss, the
    usual cause of SMA carrier status. Parascopy (MIT licence) estimates the
    aggregate copy number of the pair (agCN) and, where paralogous sequence
    variants allow, the copy number of each (psCN), each with a Phred
    quality. It reads a BAM aligned to a reference without ALT contigs (the
    default) and uses the precomputed GRCh38 homology table and the 1000
    Genomes model parameters `setup.sh --parascopy-data` installs
    (--parascopy_data; --parascopy_population picks the model, EUR by
    default).

    Background depth comes from Parascopy's own GRCh38 windows, or from
    --parascopy_depth_bed (windows of one size) for a BAM that covers only
    part of the genome.

    Equivalent to: scripts/35-paralogs.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PARASCOPY {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/paralogs" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_fai)
    path(parascopy_data)  // setup.sh --parascopy-data: homology_table/ and models_GRCh38_1KGP/
    path(depth_bed)       // background windows, or [] for Parascopy's GRCh38 windows

    output:
    tuple val(meta), path("${meta.id}_smn_copy_number.tsv"), emit: copy_number
    tuple val(meta), path("${meta.id}_parascopy"),           emit: results
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def background = depth_bed ? "-b ${depth_bed}" : '-g GRCh38'
    def model = "${parascopy_data}/models_GRCh38_1KGP/${params.parascopy_population}/SMN1.gz"
    """
    [ -f "${model}" ] || { echo "ERROR: no Parascopy model ${model} (populations: AFR AMR EAS EUR SAS)" >&2; exit 1; }
    parascopy depth -i ${bam}::${meta.id} -f ${reference} ${background} -o depth -@ ${task.cpus}
    parascopy cn-using "${model}" \\
        -i ${bam}::${meta.id} \\
        -f ${reference} \\
        -t ${parascopy_data}/homology_table/GRCh38.bed.gz \\
        -d depth \\
        -o ${meta.id}_parascopy \\
        -@ ${task.cpus}

    # One row per region of the SMN1/SMN2 profile: aggregate and
    # paralog-specific copy number, each with its filter and quality.
    {
        printf 'chrom\\tstart\\tend\\tlocus\\tagCN_filter\\tagCN\\tagCN_qual\\tpsCN_filter\\tpsCN\\tpsCN_qual\\thomologous_regions\\n'
        gzip -dc ${meta.id}_parascopy/res.samples.bed.gz \\
            | awk -F'\\t' -v OFS='\\t' '!/^#/ {print \$1, \$2, \$3, \$4, \$6, \$7, \$8, \$9, \$10, \$11, \$13}'
    } > ${meta.id}_smn_copy_number.tsv
    rm -rf depth

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parascopy: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}_parascopy
    printf 'chrom\\tstart\\tend\\tlocus\\tagCN_filter\\tagCN\\tagCN_qual\\tpsCN_filter\\tpsCN\\tpsCN_qual\\thomologous_regions\\n' \\
        > ${meta.id}_smn_copy_number.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parascopy: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
