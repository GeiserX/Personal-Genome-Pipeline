/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CNVPYTOR — Depth-based CNV detection (maintained Python successor to CNVnator)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Same read-depth method and lab (Abyzov) as CNVnator, but HDF5 .pytor storage
    instead of the CERN ROOT dependency. Orthogonal to Manta/Delly (paired-end).

    Two processes:
      1. CNVPYTOR      — rd -> his -> partition -> call (TSV) -> view (raw VCF)
      2. CNVPYTOR_VCF  — reheader from reference .fai, sort, bgzip, index

    The 1.3.2 biocontainer ships WITHOUT the GC/mask resource files and its
    built-in `-download` is broken, so the pinned resource dir is staged in and
    copied into cnvpytor's package data dir at runtime (see docs/00-reference-setup.md).

    Equivalent to: scripts/18-cnvpytor.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CNVPYTOR {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/cnvpytor:1.3.2--pyhdfd78af_0'

    input:
    tuple val(meta), path(bam), path(bai)
    path(cnvpytor_resources)

    output:
    tuple val(meta), path("${meta.id}_cnvs.txt"),     emit: cnv_calls
    tuple val(meta), path("${meta.id}_cnvs.raw.vcf"), emit: raw_vcf
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bin_size = task.ext.bin_size ?: 1000
    """
    # Stage the pinned GC/mask resources into cnvpytor's data dir (the container
    # ships without them and its -download is broken in 1.3.2).
    DATADIR=\$(python -c 'import cnvpytor, os; print(os.path.dirname(cnvpytor.__file__) + "/data")')
    cp -f ${cnvpytor_resources}/*.pytor "\$DATADIR"/

    cnvpytor -root ${meta.id}.pytor -rd ${bam} -j ${task.cpus}
    cnvpytor -root ${meta.id}.pytor -his ${bin_size}
    cnvpytor -root ${meta.id}.pytor -partition ${bin_size}
    cnvpytor -root ${meta.id}.pytor -call ${bin_size} > ${meta.id}_cnvs.txt

    printf 'set print_filename %s_cnvs.raw.vcf\\nprint calls\\n' ${meta.id} \\
        | cnvpytor -root ${meta.id}.pytor -view ${bin_size} > /dev/null

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(python -c 'import cnvpytor; print(cnvpytor.__version__)' 2>/dev/null || echo '1.3.2')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_cnvs.txt ${meta.id}_cnvs.raw.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: 1.3.2
    END_VERSIONS
    """
}

process CNVPYTOR_VCF {
    tag "$meta.id"
    label 'process_single'

    container 'staphb/bcftools:1.21'

    publishDir { "${params.outdir}/${meta.id}/cnvpytor" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(raw_vcf)
    path(reference_fai)

    output:
    tuple val(meta), path("${meta.id}_cnvs.vcf.gz"),     emit: cnv_vcf
    tuple val(meta), path("${meta.id}_cnvs.vcf.gz.tbi"), emit: cnv_vcf_tbi
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # CNVpytor's VCF only carries ##contig lines for processed chromosomes;
    # reheader from the reference .fai so headers match the other SV callers for
    # consensus merging. Emit a valid header-only VCF when there are no calls.
    if [ -s ${raw_vcf} ] && grep -qv '^#' ${raw_vcf}; then
        bcftools reheader --fai ${reference_fai} ${raw_vcf} | bcftools sort -Oz -o ${meta.id}_cnvs.vcf.gz -
    else
        { [ -s ${raw_vcf} ] && bcftools view -h ${raw_vcf} \\
            || printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n'; } \\
          | bcftools view -Oz -o ${meta.id}_cnvs.vcf.gz -
    fi
    bcftools index -t ${meta.id}_cnvs.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_cnvs.vcf.gz ${meta.id}_cnvs.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.21
    END_VERSIONS
    """
}
