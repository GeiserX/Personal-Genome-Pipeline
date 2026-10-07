/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SURVIVOR_MERGE — SV consensus: the calls two or more callers agree on
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SURVIVOR_PREP   each caller's PASS (or unfiltered) records as plain VCF,
                    one sample column named after the caller (bcftools)
    SURVIVOR_MERGE  `SURVIVOR merge LIST 1000 2 1 1 0 50`: both breakpoints
                    within 1,000 bp, same type and strands, >= 50 bp, kept
                    when 2+ callers support it; SUPP and SUPP_VEC say which
    SURVIVOR_SORT   sorted, bgzipped and indexed (bcftools)

    Equivalent to: scripts/22-survivor-merge.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SURVIVOR_PREP {
    tag "$meta.id:$caller"
    label 'process_single'

    input:
    tuple val(meta), val(caller), path(vcf)

    output:
    tuple val(meta), val(caller), path("${caller}.vcf"), emit: vcf
    path "versions.yml",                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // SURVIVOR names its output columns after the input samples: three
    // columns all named after the sample would be one name three times.
    """
    bcftools view -f PASS,. ${vcf} | awk -v s=${caller} 'BEGIN { FS = OFS = "\\t" }
        /^##/ { print; next }
        /^#CHROM/ {
            if (NF < 10) { print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\\"Genotype\\">"; sites = 1; print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, "FORMAT", s }
            else print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, s
            next
        }
        { if (sites) print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, "GT", "./."; else print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10 }' > ${caller}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${caller}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process SURVIVOR_MERGE {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), val(callers), path(vcfs)

    output:
    tuple val(meta), path("${meta.id}_sv_merged.vcf"), emit: vcf
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The list in a fixed caller order, so SUPP_VEC reads the same on every run.
    def order = ['manta', 'delly', 'cnvpytor']
    def listed = callers.sort(false) { c -> order.indexOf(c) }.collect { c -> "${c}.vcf" }.join(' ')
    """
    printf '%s\\n' ${listed} > sv_files.txt
    SURVIVOR merge sv_files.txt 1000 2 1 1 0 50 ${meta.id}_sv_merged.vcf
    # SURVIVOR exits 0 when it cannot open an input, so its output is the check.
    head -n 1 ${meta.id}_sv_merged.vcf | grep -q '^##fileformat=VCF' \\
        || { echo "ERROR: SURVIVOR wrote no VCF" >&2; exit 1; }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_sv_merged.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process SURVIVOR_SORT {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/sv_merged" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(merged)

    output:
    tuple val(meta), path("${meta.id}_sv_consensus.vcf.gz"),     emit: merged_vcf
    tuple val(meta), path("${meta.id}_sv_consensus.vcf.gz.tbi"), emit: merged_vcf_index
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools sort ${merged} -Oz -o ${meta.id}_sv_consensus.vcf.gz
    bcftools index -t ${meta.id}_sv_consensus.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_sv_consensus.vcf.gz
    touch ${meta.id}_sv_consensus.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
