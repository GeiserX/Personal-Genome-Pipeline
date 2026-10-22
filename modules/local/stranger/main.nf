/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STRANGER — Clinical STR pathogenicity annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Annotates an ExpansionHunter VCF with STR_STATUS (normal / pre_mutation /
    full_mutation), disease name, OMIM number, and repeat-size thresholds from
    Stranger's bundled GRCh38 catalog.

    STRANGER_RFC1, after it, counts the repeat motif of the reads at RFC1
    (bin/rfc1_motif.py) when Stranger flags RFC1, from the sample's BAM.

    Equivalent to: scripts/09b-stranger.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process STRANGER {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/expansion_hunter" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_eh_stranger.vcf"), emit: vcf
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Stranger's default catalog is its bundled GRCh37 one; the loci come from
    # ExpansionHunter's GRCh38 catalog, so pass the bundled GRCh38 catalog.
    CATALOG=\$(python -c 'import os, stranger; print(os.path.join(os.path.dirname(stranger.__file__), "resources", "variant_catalog_grch38.json"))')
    if [ ! -f "\$CATALOG" ]; then
        echo "ERROR: no GRCh38 catalog in the Stranger image (\${CATALOG:-not found})" >&2
        exit 1
    fi
    stranger --repeats-file "\$CATALOG" ${vcf} > ${prefix}_eh_stranger.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stranger: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_eh_stranger.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stranger: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}


/*
    STRANGER_RFC1: Stranger's catalog sizes RFC1's degenerate AARRG unit and
    flags a low count as full_mutation, so the flag says nothing until the
    motif is known. When Stranger flags RFC1 (pre or full mutation), this
    counts the pentamer motif of the reads at the locus with
    bin/rfc1_motif.py and writes <id>_rfc1_motif.tsv; without the flag it
    writes nothing, as step 9b. pysam is in PYPGX_IMAGE.

    Equivalent to: the RFC1 part of scripts/09b-stranger.sh
*/
process STRANGER_RFC1 {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/expansion_hunter" }, mode: params.publish_dir_mode

    input:
    // Fixed names: a samplesheet path may hold any character, and these go on a command line.
    tuple val(meta), path(vcf, stageAs: 'stranger.vcf'), path(bam, stageAs: 'sample.bam'), path(bai, stageAs: 'sample.bam.bai')

    output:
    tuple val(meta), path("${meta.id}_rfc1_motif.tsv"), emit: motif, optional: true
    path "versions.yml",                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if awk -F'\\t' '!/^#/ && \$8 ~ /(^|;)REPID=RFC1(;|\$)/ && \$8 ~ /(^|;)STR_STATUS=[^;]*(pre|full)_mutation/ {f = 1; exit} END {exit !f}' ${vcf}; then
        echo "RFC1 is flagged: counting the repeat motif of its reads (bin/rfc1_motif.py)"
        rfc1_motif.py --vcf ${vcf} --bam ${bam} --out ${meta.id}_rfc1_motif.tsv
    else
        echo "Stranger does not flag RFC1: no motif count"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf 'locus\\tchrom\\tstart\\tend\\trepeat_reads\\tdominant_motif\\tdominant_class\\tdominant_share\\tmotif_shares\\tverdict\\tnote\\n' > ${meta.id}_rfc1_motif.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
