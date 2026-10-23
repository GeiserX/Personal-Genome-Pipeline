/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARASCOPY — copy number of SMN1 and SMN2, HBA1/HBA2 from depth (opt-in: --tools parascopy)
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
    part of the genome, without GC stratification (--no-gc): a few regions
    rarely span the GC range Parascopy's GC model needs, and it stops on them.

    HBA1 and HBA2 have no Parascopy model: samtools bedcov (in the same
    image) sums the depth over both gene bodies and two single-copy flanks,
    over all reads and over reads of MAPQ 1 or more, and
    bin/hba_copy_number.py turns it into a count of alpha-globin genes (not
    assessed when no read over the genes has MAPQ 0: the aligner then left
    out the reads that fit both genes).
    A failure there writes an HBA table that says not assessed (so a
    published table from an earlier run is replaced) and does not fail the
    task.

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
    tuple val(meta), path("${meta.id}_hba_copy_number.tsv"), emit: hba
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def background = depth_bed ? "-b ${depth_bed} --no-gc" : '-g GRCh38'
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

    # HBA1/HBA2 from depth (as scripts/35-paralogs.sh): all MAPQ and MAPQ 1 or more, then the verdict
    if hba_copy_number.py bed > hba_regions.bed \\
        && samtools bedcov hba_regions.bed ${bam} > hba_bedcov.tsv \\
        && samtools bedcov -Q 1 hba_regions.bed ${bam} > hba_bedcov_q1.tsv \\
        && hba_copy_number.py check --bedcov hba_bedcov.tsv --bedcov-mapq1 hba_bedcov_q1.tsv \\
            --out ${meta.id}_hba_copy_number.tsv.tmp; then
        mv ${meta.id}_hba_copy_number.tsv.tmp ${meta.id}_hba_copy_number.tsv
    else
        rm -f ${meta.id}_hba_copy_number.tsv.tmp
        echo "WARNING: the HBA1/HBA2 depth check failed; the report lists HBA1/HBA2 as not assessed." >&2
        printf 'metric\\tvalue\\nstatus\\tnot assessed\\nalpha_globin_genes\\t.\\nestimate\\t.\\nhba2_depth\\t.\\nhba1_depth\\t.\\nflank_left_depth\\t.\\nflank_right_depth\\t.\\nmapq0_share\\t.\\nmessage\\t%s\\n' \\
            'not assessed: the HBA1/HBA2 depth check failed; see the log of step 35' > ${meta.id}_hba_copy_number.tsv
    fi
    rm -f hba_regions.bed hba_bedcov.tsv hba_bedcov_q1.tsv

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
    printf 'metric\\tvalue\\nstatus\\tnot assessed\\nalpha_globin_genes\\t.\\nestimate\\t.\\nhba2_depth\\t0\\nhba1_depth\\t0\\nflank_left_depth\\t0\\nflank_right_depth\\t0\\nmapq0_share\\t0\\nmessage\\tstub\\n' \\
        > ${meta.id}_hba_copy_number.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parascopy: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
