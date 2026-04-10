/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CNVNATOR — Depth-based CNV detection (orthogonal to Manta/Delly)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Five-step process:
    1. Extract read mapping (tree)
    2. Generate read-depth histogram
    3. Compute statistics
    4. Partition
    5. Call CNVs

    Equivalent to: scripts/18-cnvnator.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CNVNATOR {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/cnvnator:0.4.1--py312h99c8fb2_11'

    publishDir "${params.outdir}/${meta.id}/cnvnator", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_fai)

    output:
    tuple val(meta), path("${meta.id}_cnvs.txt"),    emit: cnv_calls
    tuple val(meta), path("${meta.id}_cnvs.vcf.gz"), emit: cnv_vcf, optional: true
    path "versions.yml",                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bin_size = task.ext.bin_size ?: 1000
    """
    # Step 1: Extract read mapping
    cnvnator \\
        -root ${meta.id}.root \\
        -tree ${bam}

    # Step 2: Generate read-depth histogram
    cnvnator \\
        -root ${meta.id}.root \\
        -his ${bin_size} \\
        -fasta ${reference}

    # Step 3: Compute statistics
    cnvnator \\
        -root ${meta.id}.root \\
        -stat ${bin_size}

    # Step 4: Partition
    cnvnator \\
        -root ${meta.id}.root \\
        -partition ${bin_size}

    # Step 5: Call CNVs
    cnvnator \\
        -root ${meta.id}.root \\
        -call ${bin_size} \\
        > ${meta.id}_cnvs.txt

    # Convert CNVnator TXT to VCF for downstream merging
    {
        echo '##fileformat=VCFv4.2'
        echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
        echo '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">'
        echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">'
        awk '{printf "##contig=<ID=%s,length=%s>\\n", \$1, \$2}' ${reference_fai} || true
        printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n'
        awk '{
            split(\$2,a,":");
            split(a[2],b,"-");
            svtype="DEL"; alt="<DEL>"; svlen=-(b[2]-b[1]);
            if(\$1=="duplication") { svtype="DUP"; alt="<DUP>"; svlen=b[2]-b[1]; }
            printf "%s\\t%s\\t.\\tN\\t%s\\t.\\t.\\tSVTYPE=%s;END=%s;SVLEN=%d\\n",
                a[1], b[1], alt, svtype, b[2], svlen;
        }' ${meta.id}_cnvs.txt
    } > ${meta.id}_cnvs_unsorted.vcf

    # Sort, compress, and index — fail on real errors, emit header-only VCF if no calls
    if grep -q -v '^#' ${meta.id}_cnvs_unsorted.vcf; then
        bcftools sort ${meta.id}_cnvs_unsorted.vcf -Oz -o ${meta.id}_cnvs.vcf.gz
        bcftools index -t ${meta.id}_cnvs.vcf.gz
    else
        bcftools view -h ${meta.id}_cnvs_unsorted.vcf | bgzip -c > ${meta.id}_cnvs.vcf.gz
        bcftools index -t ${meta.id}_cnvs.vcf.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvnator: \$(cnvnator 2>&1 | grep -oP 'CNVnator v\\K[\\d.]+' || echo '0.4.1')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_cnvs.txt
    touch ${meta.id}_cnvs.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvnator: 0.4.1
    END_VERSIONS
    """
}
