/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ClinVar Pathogenic Screen — intersect sample VCF with ClinVar pathogenic/LP variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Finds known disease-causing variants the person carries by intersecting
    PASS variants against the ClinVar pathogenic subset.

    Equivalent to: scripts/06-clinvar-screen.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CLINVAR_SCREEN {
    tag "$meta.id"
    label 'process_low'

    container 'staphb/bcftools:1.21'

    publishDir "${params.outdir}/${meta.id}/clinvar", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(clinvar)
    path(clinvar_index)

    output:
    tuple val(meta), path("isec/"),              emit: isec_dir
    tuple val(meta), path("*_pass.vcf.gz"),      emit: pass_vcf
    tuple val(meta), path("*_pass.vcf.gz.tbi"),  emit: pass_vcf_index
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Step 1: Extract PASS variants
    bcftools view -f PASS ${vcf} -Oz -o ${meta.id}_pass.vcf.gz
    bcftools index -t ${meta.id}_pass.vcf.gz

    # Step 2: Intersect with ClinVar pathogenic
    bcftools isec -p isec ${meta.id}_pass.vcf.gz ${clinvar}

    # Log hit count
    HITS=\$(grep -c -v '^#' isec/0002.vcf 2>/dev/null || echo 0)
    echo "ClinVar pathogenic hits: \${HITS}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
