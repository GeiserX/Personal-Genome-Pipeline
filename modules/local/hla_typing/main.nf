/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HLA_TYPING — HLA allele typing from WGS BAM using T1K
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Types HLA-A, B, C (Class I) and DRB1, DQB1, DPB1 (Class II)
    at 4-digit resolution using IPD-IMGT/HLA database against GRCh38.

    Two-step process:
    1. Build HLA coordinate reference from genome (one-time, cached)
    2. Run T1K genotyping against the BAM

    Equivalent to: scripts/08-hla-typing.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process HLA_TYPING {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0'

    publishDir "${params.outdir}/${meta.id}/hla", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(hla_dat)

    output:
    tuple val(meta), path("*_hla_genotype.tsv"), emit: hla_alleles
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Step 1: Build coordinate file from reference genome using pre-downloaded hla.dat
    t1k-build.pl \\
        -d ${hla_dat} \\
        -g ${reference} \\
        -o hlaidx_grch38

    # Step 3: Run HLA typing
    # Locate build output (file naming varies across T1K versions)
    SEQ_FA=\$(ls hlaidx_grch38/*dna_seq.fa 2>/dev/null | head -1)
    COORD_FA=\$(ls hlaidx_grch38/*dna_coord.fa 2>/dev/null | head -1)
    if [ -z "\$SEQ_FA" ] || [ -z "\$COORD_FA" ]; then
        echo "ERROR: t1k-build did not produce expected output files in hlaidx_grch38/"
        ls -la hlaidx_grch38/ 2>/dev/null
        exit 1
    fi

    run-t1k \\
        -b ${bam} \\
        -f "\$SEQ_FA" \\
        -c "\$COORD_FA" \\
        --preset hla-wgs \\
        -t ${task.cpus} \\
        --od ./ \\
        -o ${prefix}_hla

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: \$(run-t1k --version 2>&1 | grep -oP '[\\d.]+' | head -1 || echo '1.0.9')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_hla_genotype.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: 1.0.9
    END_VERSIONS
    """
}
