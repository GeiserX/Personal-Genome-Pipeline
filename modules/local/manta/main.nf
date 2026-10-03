/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MANTA — Structural variant calling (deletions, duplications, inversions, translocations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Three steps:
    1. Configure Manta (configManta.py), on --manta_call_regions only when set
    2. Run workflow (runWorkflow.py)
    3. Inversions: Manta writes each one as two breakend (BND) records, and
       duphold, AnnotSV and the SV consensus read SVTYPE. Manta's
       libexec/convertInversion.py (with the samtools and bgzip beside it in
       the image) turns each pair into one SVTYPE=INV record. Manta's own file
       is kept as diploidSV.raw.vcf.gz; diploidSV.vcf.gz is the converted one.

    Equivalent to: scripts/04-manta.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MANTA {
    tag "$meta.id"
    label 'process_high'

    publishDir { "${params.outdir}/${meta.id}/manta" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_fai)
    path(call_regions)        // bgzipped BED or []
    path(call_regions_index)  // its .tbi or []

    output:
    tuple val(meta), path("results/variants/diploidSV.vcf.gz"),     emit: diploid_sv
    tuple val(meta), path("results/variants/diploidSV.raw.vcf.gz"), emit: diploid_sv_raw
    tuple val(meta), path("results/variants/candidateSV.vcf.gz"),   emit: candidate_sv
    path "versions.yml",                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def regions_arg = call_regions ? "--callRegions ${call_regions}" : ''
    """
    # Step 1: Configure Manta
    configManta.py \\
        --bam ${bam} \\
        --referenceFasta ${reference} \\
        ${regions_arg} \\
        --runDir manta_run

    # Step 2: Run Manta workflow
    manta_run/runWorkflow.py -j ${task.cpus}

    # Move results to expected output location
    mv manta_run/results .

    # Step 3: inversion breakend pairs become SVTYPE=INV records
    LIBEXEC="\$(dirname "\$(readlink -f "\$(command -v configManta.py)")")/../libexec"
    V=results/variants
    mv \$V/diploidSV.vcf.gz \$V/diploidSV.raw.vcf.gz
    mv \$V/diploidSV.vcf.gz.tbi \$V/diploidSV.raw.vcf.gz.tbi
    "\$LIBEXEC/convertInversion.py" "\$LIBEXEC/samtools" ${reference} \$V/diploidSV.raw.vcf.gz \\
        | "\$LIBEXEC/bgzip" -c > \$V/diploidSV.vcf.gz
    "\$LIBEXEC/tabix" -f -p vcf \$V/diploidSV.vcf.gz
    N_BND=\$(gzip -cd \$V/diploidSV.raw.vcf.gz | awk -F'\\t' '!/^#/ && (\$5 ~ /^\\[/ || \$5 ~ /\\]\$/) {
        m = \$5; sub(/^[^][]*[][]/, "", m); sub(/:.*/, "", m); if (m == \$1) n++ } END { print n + 0 }')
    N_INV=\$(gzip -cd \$V/diploidSV.vcf.gz | awk -F'\\t' '!/^#/ && \$8 ~ /(^|;)SVTYPE=INV(;|\$)/ { n++ } END { print n + 0 }')
    echo "Inversion conversion: \${N_BND} inversion breakend records in, \${N_INV} SVTYPE=INV records out"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p results/variants
    touch results/variants/diploidSV.vcf.gz
    touch results/variants/diploidSV.raw.vcf.gz
    touch results/variants/candidateSV.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
