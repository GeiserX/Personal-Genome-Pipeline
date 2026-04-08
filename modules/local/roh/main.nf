/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ROH — Runs of Homozygosity (consanguinity screening)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Detects long stretches of homozygous DNA that indicate shared ancestry.
    Centromeric regions (chr1 125-143MB, chr9 42-60MB, chr18 15-20MB) are
    known artifacts and should be filtered during interpretation.

    Equivalent to: scripts/11-roh-analysis.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ROH {
    tag "$meta.id"
    label 'process_low'

    container 'staphb/bcftools:1.21'

    publishDir "${params.outdir}/${meta.id}/roh", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("*_roh.txt"),         emit: roh_regions
    tuple val(meta), path("*_roh_summary.txt"), emit: roh_summary
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Auto-detect chip data: if FORMAT/PL is absent, use -G30 (genotype-only mode)
    HAS_PL=\$(bcftools view -h ${vcf} | grep -c '##FORMAT=<ID=PL' || true)

    ROH_FLAGS="--AF-dflt 0.4"
    if [ "\${HAS_PL}" -eq 0 ]; then
        ROH_FLAGS="\${ROH_FLAGS} -G30"
    fi

    bcftools roh \${ROH_FLAGS} \\
        -o ${meta.id}_roh.txt \\
        ${vcf}

    # Generate summary: autosomal ROH segments >1MB
    echo "# ROH Summary for ${meta.id}" > ${meta.id}_roh_summary.txt
    echo "# Segments >1MB on autosomes (excludes chrX/chrY)" >> ${meta.id}_roh_summary.txt
    echo -e "chrom\\tstart\\tend\\tlength_bp\\tlength_mb" >> ${meta.id}_roh_summary.txt
    grep '^RG' ${meta.id}_roh.txt 2>/dev/null | \\
        awk '\$3 !~ /chrX|chrY/ && \$6 > 1000000 {printf "%s\\t%s\\t%s\\t%s\\t%.1f\\n", \$3,\$4,\$5,\$6,\$6/1e6}' \\
        >> ${meta.id}_roh_summary.txt || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_roh.txt ${meta.id}_roh_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
