/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VCFANNO — Annotate VCF with CADD, SpliceAI, REVEL, and AlphaMissense scores
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Enriches a VEP-annotated VCF with pathogenicity scores from external databases.

    Two processes, because the vcfanno image holds only the vcfanno binary
    (its conda package depends on glibc alone: no bcftools, bgzip or tabix):
      1. VCFANNO        — vcfanno image, one pass over every score file, plain VCF out
      2. VCFANNO_INDEX  — bcftools image, bgzip + tabix index of that VCF

    Chromosome names: CADD uses bare names (1, 2, 3) while the VCF and the other
    score files use chr1, chr2, chr3. vcfanno handles this itself: its tabix
    reader retries a region with the "chr" prefix added or removed (brentp/bix,
    ChunkedReader) and its sweep compares positions within one chromosome only,
    so a single pass annotates both kinds of file.

    Score files are optional: only those given are used. SpliceAI accepts the
    masked or the raw precomputed files under any name.

    Equivalent to: scripts/30-vcfanno.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process VCFANNO {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(cadd_snv)
    path(cadd_snv_index)
    path(cadd_indel)
    path(cadd_indel_index)
    path(spliceai_snv)
    path(spliceai_snv_index)
    path(spliceai_indel)
    path(spliceai_indel_index)
    path(revel)
    path(revel_index)
    path(alphamissense)
    path(alphamissense_index)

    output:
    tuple val(meta), path("${meta.id}_vcfanno.vcf"), emit: vcf
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def toml = ""
    if (cadd_snv) {
        toml += """
[[annotation]]
file="${cadd_snv}"
columns=[6]
names=["CADD_PHRED"]
ops=["self"]
"""
    }
    if (cadd_indel) {
        toml += """
[[annotation]]
file="${cadd_indel}"
columns=[6]
names=["CADD_PHRED_indel"]
ops=["self"]
"""
    }
    if (spliceai_snv) {
        toml += """
[[annotation]]
file="${spliceai_snv}"
fields=["SpliceAI"]
names=["SpliceAI"]
ops=["self"]
"""
    }
    if (spliceai_indel) {
        toml += """
[[annotation]]
file="${spliceai_indel}"
fields=["SpliceAI"]
names=["SpliceAI_indel"]
ops=["self"]
"""
    }
    if (revel) {
        toml += """
[[annotation]]
file="${revel}"
columns=[5]
names=["REVEL"]
ops=["self"]
"""
    }
    if (alphamissense) {
        toml += """
[[annotation]]
file="${alphamissense}"
columns=[9,10]
names=["AM_pathogenicity","AM_class"]
ops=["self","self"]
"""
    }
    if (!toml) {
        error "VCFANNO needs at least one score file (--cadd_snv, --cadd_indel, --spliceai_snv, --spliceai_indel, --revel or --alphamissense)"
    }
    """
    cat > vcfanno.toml <<'TOML_END'
${toml}
TOML_END

    vcfanno -p ${task.cpus} vcfanno.toml ${vcf} > ${meta.id}_vcfanno.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfanno: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_vcfanno.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfanno: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process VCFANNO_INDEX {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/vep" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}_annotated.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}_annotated.vcf.gz.tbi"), emit: vcf_index
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools view ${vcf} -Oz -o ${meta.id}_annotated.vcf.gz
    bcftools index -t ${meta.id}_annotated.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_annotated.vcf.gz
    touch ${meta.id}_annotated.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
