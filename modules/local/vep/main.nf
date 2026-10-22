/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VEP_NORM — the sample's PASS records, split and left-aligned, for VEP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VEP matches an allele to gnomAD and ClinVar as those files write it. A
    deletion written right-shifted, or inside a multiallelic record, finds no
    known variant and gets no frequency, so the clinical filter would show a
    common allele as "frequency unknown". This keeps the PASS records ('.'
    too), splits multiallelic records and left-aligns indels against the
    reference (bcftools norm -m -any -f), keeps the AD sum of a split when
    the header defines AD (--keep-sum AD), drops the half whose genotype
    carries no ALT (0/0), and leaves out the records of a contig the
    reference lacks, saying how many contigs that is.

    Equivalent to: the normalisation in scripts/13-vep-annotation.sh, the
    same commands in the same order
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process VEP_NORM {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(reference)
    path(reference_fai)  // staged beside the FASTA

    output:
    tuple val(meta), path("${meta.id}.norm.vcf.gz"), path("${meta.id}.norm.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Temporary names after the sample id: the given VCF may have any name.
    P=${meta.id}.vep_norm_pass.tmp.vcf.gz N=${meta.id}.vep_norm_split.tmp.vcf.gz T=${meta.id}.vep_norm_contigs.tmp.tsv
    bcftools view -f PASS,. -Oz -o \$P ${vcf}
    # -T keeps the contigs the reference holds: bcftools norm stops at the
    # first record whose contig the FASTA lacks. -c w warns on a REF base
    # that differs.
    awk -v OFS='\\t' '{print \$1, 1, 2147483647}' ${reference_fai} > \$T
    OFF_REF=\$(bcftools index -s ${vcf} \\
        | awk -F'\\t' 'NR == FNR {ref[\$1] = 1; next} \$3 > 0 && !(\$1 in ref) {n++} END {print n + 0}' ${reference_fai} -)
    if [ "\$OFF_REF" -gt 0 ]; then
        echo "NOTICE: \$OFF_REF contig(s) of ${vcf} are not in ${reference}; their records are not annotated."
    fi
    # --keep-sum AD stops with an error when the header defines no FORMAT/AD.
    KEEP_SUM=""
    # Matched in a variable, as step 13 does: a pipe into grep -q could stop
    # bcftools early and, under pipefail, fail the test.
    HDR=\$(bcftools view -h \$P)
    case "\$HDR" in *'##FORMAT=<ID=AD,'*) KEEP_SUM="--keep-sum AD" ;; esac
    bcftools norm -m -any \$KEEP_SUM -c w -f ${reference} -T \$T -Oz -o \$N \$P
    bcftools view -i 'GT="alt"' -Oz -o ${meta.id}.norm.vcf.gz \$N
    bcftools index -t ${meta.id}.norm.vcf.gz
    rm -f \$P \$N \$T

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.norm.vcf.gz ${meta.id}.norm.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VEP — Ensembl Variant Effect Predictor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Full functional annotation: consequence, SIFT, PolyPhen, gnomAD AF, ClinVar,
    regulatory features, etc. Uses offline cache for reproducibility. With
    --clinvar, the same ClinVar file the screen reads is added with --custom:
    its CLNSIG, CLNREVSTAT and CLNDN as ClinVar_CLNSIG, ClinVar_CLNREVSTAT and
    ClinVar_CLNDN, so the clinical filter's ClinVar tier follows a refresh of
    that file instead of the cache release.

    VEP reads the PASS records only ('.' too, for a caller that writes no
    FILTER): the clinical filter and slivar, after it, keep PASS records
    alone. VEP_NORM, before it, selects them and normalises them in the
    bcftools image, because the VEP image has no bcftools; the awk below
    keeps the same records again, so a VCF given straight to VEP reads the
    same way.

    Equivalent to: scripts/13-vep-annotation.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process VEP {
    tag "$meta.id"
    label 'process_high'

    publishDir { "${params.outdir}/${meta.id}/vep" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(reference)
    path(reference_fai)  // staged beside the FASTA, so no task builds its own
    path(vep_cache)
    path(clinvar)        // [] without --clinvar
    path(clinvar_index)  // staged beside it

    output:
    tuple val(meta), path("*_vep.vcf.gz"),     emit: vcf
    tuple val(meta), path("*_vep.vcf.gz.tbi"), emit: vcf_index
    tuple val(meta), path("*_vep.vcf_summary.html"), emit: stats, optional: true
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def custom = clinvar ? "--custom file=${clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN" : ''
    """
    bgzip -dc ${vcf} | awk -F '\\t' '/^#/ || \$7 == "." || (";" \$7 ";") ~ /;PASS;/' | bgzip -c > ${meta.id}.pass.vcf.gz

    vep \\
        --input_file ${meta.id}.pass.vcf.gz \\
        --output_file ${meta.id}_vep.vcf \\
        --vcf \\
        --cache \\
        --dir_cache ${vep_cache} \\
        --offline \\
        --assembly GRCh38 \\
        --everything \\
        --af_gnomade \\
        --force_overwrite \\
        --fork ${task.cpus} \\
        --fasta ${reference} \\
        ${custom} \\
        ${args}

    bgzip -c ${meta.id}_vep.vcf > ${meta.id}_vep.vcf.gz
    tabix -p vcf ${meta.id}_vep.vcf.gz
    rm -f ${meta.id}_vep.vcf ${meta.id}.pass.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_vep.vcf.gz
    touch ${meta.id}_vep.vcf.gz.tbi
    touch ${meta.id}_vep.vcf_summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
