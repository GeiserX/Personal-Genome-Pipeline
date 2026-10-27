/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TRAITS — The traits panel of step 42, from the gVCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TRAITS_SITES  genotypes the positions of assets/traits.tsv and of the
                  rule sites of assets/genotype_rules.tsv from the gVCF:
                  each reference block over a position becomes a record with
                  the reference base and the block's GQ, and a position
                  without reads has no line (it is "not covered", never the
                  reference)
    TRAITS_PHASE  for a sample with a BAM: WhatsHap phases the heterozygous
                  sites of the rules that need phase (TAS2R38) from the reads;
                  an empty file when no rule has two such sites
    TRAITS        <id>_traits.tsv: the reading of each genotype from the
                  table, and the APOE line; <id>_genotype_rules.tsv: one
                  line per rule (bin/traits_report.py)

    Runs for the samples with a gVCF only (DEEPVARIANT's, or one the
    samplesheet gives); main.nf logs the samples it skips. A sample without
    a BAM gets [] for the phase: a rule that needs it reads only when at most
    one of its sites is heterozygous.

    Equivalent to: scripts/42-traits.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process TRAITS_SITES {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/traits" }, mode: params.publish_dir_mode, pattern: '*_traits_genotypes.tsv'

    input:
    tuple val(meta), path(gvcf), path(gvcf_index)
    path(table)
    path(rules)     // assets/genotype_rules.tsv
    path(reference)
    path(reference_fai)  // staged beside the FASTA

    output:
    tuple val(meta), path("${meta.id}_traits_genotypes.tsv"), emit: genotypes
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The same commands as scripts/42-traits.sh
    """
    n=\$(bcftools query -l ${gvcf} | wc -l | tr -d ' ')
    if [ "\$n" -ne 1 ]; then
        echo "ERROR: Sample '${meta.id}': ${gvcf} holds \$n samples; step 42 reads a one-sample gVCF (keep one with bcftools view -s <sample>)." >&2
        exit 1
    fi
    { grep -v '^#' ${table} | awk -F'\\t' 'NR > 1 { print \$1 "\\t" \$2 }'
      awk -F'\\t' '\$1 == "site" { print \$4 "\\t" \$5 }' ${rules}; } | LC_ALL=C sort -k1,1 -k2,2n -u > sites.tsv
    bcftools convert --gvcf2vcf -f ${reference} -R sites.tsv -Ou ${gvcf} \\
        | bcftools view -T sites.tsv --trim-alt-alleles -i 'GT!="mis"' -Ou \\
        | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%TGT]\\t[%GQ]\\n' > ${meta.id}_traits_genotypes.tsv
    rm -f sites.tsv

    printf '"%s":\\n    bcftools: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    touch ${meta.id}_traits_genotypes.tsv
    printf '"%s":\\n    bcftools: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}

process TRAITS_PHASE {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(genotypes), path(bam), path(bai)
    path(rules)
    path(reference)
    path(reference_fai)  // staged beside the FASTA

    output:
    tuple val(meta), path("${meta.id}_traits_phased.vcf"), emit: phased
    path "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The same commands as scripts/42-traits.sh
    """
    traits_report.py --rules ${rules} --genotypes ${genotypes} --sample ${meta.id} --phase-sites phase_sites.vcf
    if grep -qv '^#' phase_sites.vcf; then
        whatshap phase --reference ${reference} --ignore-read-groups -o ${meta.id}_traits_phased.vcf phase_sites.vcf ${bam}
    else
        : > ${meta.id}_traits_phased.vcf
    fi
    rm -f phase_sites.vcf

    printf '"%s":\\n    whatshap: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    touch ${meta.id}_traits_phased.vcf
    printf '"%s":\\n    whatshap: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}

process TRAITS {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/traits" }, mode: params.publish_dir_mode, pattern: '{*_traits.tsv,*_genotype_rules.tsv}'

    input:
    tuple val(meta), path(genotypes), path(phased)   // phased: TRAITS_PHASE's, or [] without a BAM
    path(table)
    path(rules)

    output:
    tuple val(meta), path("${meta.id}_traits.tsv"),         emit: traits
    tuple val(meta), path("${meta.id}_genotype_rules.tsv"), emit: rules
    path "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def phased_arg = phased ? "--phased ${phased}" : ''
    """
    traits_report.py --table ${table} --rules ${rules} --genotypes ${genotypes} ${phased_arg} --sample ${meta.id} \\
        --out ${meta.id}_traits.tsv --rules-out ${meta.id}_genotype_rules.tsv

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    printf 'rsid\\tgene\\tchrom\\tpos\\ttrait\\tgenotype\\tstatus\\treading\\tancestry\\tnote\\tpmids\\n' > ${meta.id}_traits.tsv
    printf 'rule\\tgene\\ttrait\\trsids\\tchrom\\tpos\\tgenotypes\\tgq\\tphase\\tresult\\tstatus\\treading\\tancestry\\tnote\\tpmids\\n' > ${meta.id}_genotype_rules.tsv
    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}
