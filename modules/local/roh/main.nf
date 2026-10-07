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

    publishDir { "${params.outdir}/${meta.id}/roh" }, mode: params.publish_dir_mode

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

    # Only PASS calls (and records with no filter, as chip VCFs have), as in
    # scripts/11-roh-analysis.sh. pipefail: a failed view must fail the task.
    set -o pipefail
    bcftools view -f PASS,. -Ou ${vcf} \\
        | bcftools roh \${ROH_FLAGS} -o ${meta.id}_roh.txt -

    # Summary: autosomal segments of 5 Mb or more, the threshold of
    # docs/11-roh-analysis.md; the bash step writes the same file.
    echo "# ROH Summary for ${meta.id}" > ${meta.id}_roh_summary.txt
    echo "# Segments >=5MB on autosomes (excludes chrX/chrY)" >> ${meta.id}_roh_summary.txt
    printf 'chrom\\tstart\\tend\\tlength_bp\\tlength_mb\\n' >> ${meta.id}_roh_summary.txt
    awk '\$1 == "RG" && \$3 !~ /chrX|chrY/ && \$6 >= 5000000 {printf "%s\\t%s\\t%s\\t%s\\t%.1f\\n", \$3, \$4, \$5, \$6, \$6 / 1e6}' \\
        ${meta.id}_roh.txt >> ${meta.id}_roh_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_roh.txt ${meta.id}_roh_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
