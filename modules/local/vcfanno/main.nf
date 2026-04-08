/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VCFANNO — Annotate VCF with CADD, SpliceAI, REVEL, and AlphaMissense scores
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Enriches a VEP-annotated VCF with pathogenicity scores from external databases.
    Handles chromosome naming mismatch: CADD uses bare names (1, 2, 3) while
    the VCF and other databases use chr-prefixed names (chr1, chr2, chr3).
    Solved via two-pass annotation with chromosome renaming between passes.

    Annotation files are optional — only present files are annotated.

    Equivalent to: scripts/30-vcfanno.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process VCFANNO {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/vcfanno:0.3.7--he881be0_0'

    publishDir "${params.outdir}/${meta.id}/vep", mode: params.publish_dir_mode

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
    tuple val(meta), path("*_annotated.vcf.gz"),     emit: vcf
    tuple val(meta), path("*_annotated.vcf.gz.tbi"), emit: vcf_index
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Build TOML for chr-prefixed tracks (SpliceAI, REVEL, AlphaMissense)
    def chr_toml = ""
    if (spliceai_snv) {
        chr_toml += """
[[annotation]]
file="${spliceai_snv}"
fields=["SpliceAI"]
names=["SpliceAI"]
ops=["self"]
"""
    }
    if (spliceai_indel) {
        chr_toml += """
[[annotation]]
file="${spliceai_indel}"
fields=["SpliceAI"]
names=["SpliceAI_indel"]
ops=["self"]
"""
    }
    if (revel) {
        chr_toml += """
[[annotation]]
file="${revel}"
columns=[5]
names=["REVEL"]
ops=["self"]
"""
    }
    if (alphamissense) {
        chr_toml += """
[[annotation]]
file="${alphamissense}"
columns=[9,10]
names=["AM_pathogenicity","AM_class"]
ops=["self","self"]
"""
    }

    // Build TOML for no-chr tracks (CADD)
    def nochr_toml = ""
    if (cadd_snv) {
        nochr_toml += """
[[annotation]]
file="${cadd_snv}"
columns=[6]
names=["CADD_PHRED"]
ops=["self"]
"""
    }
    if (cadd_indel) {
        nochr_toml += """
[[annotation]]
file="${cadd_indel}"
columns=[6]
names=["CADD_PHRED_indel"]
ops=["self"]
"""
    }

    def has_nochr = (cadd_snv || cadd_indel) ? true : false
    def has_chr   = (spliceai_snv || spliceai_indel || revel || alphamissense) ? true : false
    """
    CURRENT_VCF="${vcf}"

    # --- Pass 1: CADD annotation (no-chr tracks) ---
    if [ "${has_nochr}" = "true" ]; then
        # Generate chromosome rename maps
        for i in \$(seq 1 22) X Y M; do
            echo "chr\${i} \${i}" >> strip_chr.txt
            echo "\${i} chr\${i}" >> add_chr.txt
        done

        cat > nochr.toml <<'TOML_END'
${nochr_toml}
TOML_END

        # Strip chr prefix -> annotate -> re-add chr prefix
        bcftools annotate --rename-chrs strip_chr.txt \${CURRENT_VCF} -Oz -o nochr_input.vcf.gz
        tabix -p vcf nochr_input.vcf.gz

        vcfanno -p ${task.cpus} nochr.toml nochr_input.vcf.gz > nochr_annotated.vcf

        bgzip -c nochr_annotated.vcf > nochr_annotated.vcf.gz
        bcftools annotate --rename-chrs add_chr.txt nochr_annotated.vcf.gz -Oz -o pass1_output.vcf.gz
        tabix -p vcf pass1_output.vcf.gz

        CURRENT_VCF="pass1_output.vcf.gz"
    fi

    # --- Pass 2: chr-prefixed tracks (SpliceAI, REVEL, AlphaMissense) ---
    if [ "${has_chr}" = "true" ]; then
        cat > chr.toml <<'TOML_END'
${chr_toml}
TOML_END

        vcfanno -p ${task.cpus} chr.toml \${CURRENT_VCF} > pass2_output.vcf
        CURRENT_VCF="pass2_output.vcf"
    fi

    # --- If nothing to annotate, just copy input ---
    if [ "${has_nochr}" = "false" ] && [ "${has_chr}" = "false" ]; then
        cp ${vcf} ${meta.id}_annotated.vcf.gz
        cp ${vcf_index} ${meta.id}_annotated.vcf.gz.tbi
    else
        # Compress and index final output
        if [ "\${CURRENT_VCF}" = "pass1_output.vcf.gz" ]; then
            cp pass1_output.vcf.gz ${meta.id}_annotated.vcf.gz
            tabix -p vcf ${meta.id}_annotated.vcf.gz
        else
            bgzip -c \${CURRENT_VCF} > ${meta.id}_annotated.vcf.gz
            tabix -p vcf ${meta.id}_annotated.vcf.gz
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfanno: \$(vcfanno 2>&1 | grep -oP 'version\\s+\\K[\\d.]+' || echo '0.3.7')
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_annotated.vcf.gz
    touch ${meta.id}_annotated.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfanno: 0.3.7
        bcftools: 1.21
    END_VERSIONS
    """
}
