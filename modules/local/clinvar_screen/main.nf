/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ClinVar Pathogenic Screen — intersect sample VCF with ClinVar pathogenic/LP variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Finds known disease-causing variants the person carries by intersecting
    PASS variants against the ClinVar pathogenic subset. Writes
    <id>_clinvar_hits.vcf: the sample's matching records with ClinVar's ID,
    GENEINFO, CLNSIG and CLNREVSTAT copied on, and prints the hits per ClinVar
    review-status star tier.

    Equivalent to: scripts/06-clinvar-screen.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CLINVAR_SCREEN {
    tag "$meta.id"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/clinvar" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(clinvar)
    path(clinvar_index)
    path(reference)

    output:
    // The emit keeps its old name (isec_dir) so workflows/pgx.nf needs no change;
    // it now carries the annotated hits VCF, not the isec directory.
    tuple val(meta), path("${meta.id}_clinvar_hits.vcf"), emit: isec_dir
    tuple val(meta), path("*_pass.vcf.gz"),      emit: pass_vcf
    tuple val(meta), path("*_pass.vcf.gz.tbi"),  emit: pass_vcf_index
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    # Step 0: the two files must share contig names, or zero hits would look clean
    SHARED=\$(comm -12 <(bcftools index -s ${vcf} | cut -f1 | sort -u) <(bcftools index -s ${clinvar} | cut -f1 | sort -u) | grep -c . || true)
    if [ "\${SHARED}" -eq 0 ]; then
        echo "ERROR: ${vcf} and ${clinvar} have no contig name in common (chr1 vs 1?). Use the chr-prefixed ClinVar file." >&2
        exit 1
    fi

    # Step 1: ClinVar split and left-aligned. Reuse the copy step 06 (bash) keeps beside
    # the source file when it is newer than the source; otherwise normalise here.
    SRC=\$(readlink -f ${clinvar})
    NORM_BESIDE="\${SRC%.vcf.gz}.norm.vcf.gz"
    if [ -f "\${NORM_BESIDE}" ] && [ -f "\${NORM_BESIDE}.tbi" ] && [ "\${NORM_BESIDE}" -nt "\${SRC}" ]; then
        CLINVAR_NORM="\${NORM_BESIDE}"
    else
        bcftools norm -m -any -c w -f ${reference} ${clinvar} -Oz -o clinvar_norm.vcf.gz
        bcftools index -t clinvar_norm.vcf.gz
        CLINVAR_NORM=clinvar_norm.vcf.gz
    fi

    # Step 2: PASS records, or '.,PASS' when the caller never writes PASS
    HAS_PASS=\$( (bcftools view -H -f PASS ${vcf} || true) | head -n 1 | wc -l)
    if [ "\${HAS_PASS}" -gt 0 ]; then
        FILTER=PASS
    else
        FILTER=.,PASS
        echo "NOTICE: ${vcf} has no PASS record; using records with FILTER '.' or PASS"
    fi
    bcftools view -f "\${FILTER}" -Ou ${vcf} | \\
        bcftools norm -m -any -c w -f ${reference} -Oz -o ${meta.id}_pass.vcf.gz -
    bcftools index -t ${meta.id}_pass.vcf.gz
    if [ "\$(bcftools index -n ${meta.id}_pass.vcf.gz)" -eq 0 ]; then
        echo "ERROR: no records left after the '\${FILTER}' filter" >&2
        exit 1
    fi

    # Step 3: sample records matching a ClinVar allele, with ClinVar's ID, gene,
    # significance and review status copied on
    # (annotate -a needs an indexed target, so the shared records go to a file first)
    bcftools isec -n=2 -w1 -Oz -o shared.vcf.gz ${meta.id}_pass.vcf.gz "\${CLINVAR_NORM}"
    bcftools index -t shared.vcf.gz
    bcftools annotate -a "\${CLINVAR_NORM}" --pair-logic exact \\
        -c ID,INFO/GENEINFO,INFO/CLNSIG,INFO/CLNREVSTAT -Ov -o ${meta.id}_clinvar_hits.vcf shared.vcf.gz

    HITS=\$(grep -c -v '^#' ${meta.id}_clinvar_hits.vcf || true)
    echo "ClinVar pathogenic hits: \${HITS}"
    # Grouped by ClinVar review stars, as scripts/06 prints them
    # (bin/clinvar_hits.awk is on the task PATH).
    if [ "\${HITS}" -gt 0 ]; then
        echo "By review status (ClinVar stars):"
        awk -f "\$(command -v clinvar_hits.awk)" ${meta.id}_clinvar_hits.vcf | cut -f1 | sort -nr | uniq -c \\
            | awk '{printf "  %d star%s: %d\\n", \$2, (\$2 == 1 ? "" : "s"), \$1}'
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_clinvar_hits.vcf
    touch ${meta.id}_pass.vcf.gz
    touch ${meta.id}_pass.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
