/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MEI — Mobile element insertions: Alu, L1, SVA with MEGAnE (step 50)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Opt-in ('mei' in --tools, with --megane_kmers: setup.sh --megane-kmers).
    Three processes, as in scripts/50-mei.sh:
      MEI_DEPTH   samtools: the chromosomes MEGAnE reads (chr1 to chr22, chrX,
                  chrY as the BAM names them) and samtools bedcov over the
                  windows of assets/mei_sites.tsv on the BAM's contigs
      MEI         MEGAnE call_genotype_38 -only_ins on the BAM
      MEI_REPORT  bin/mei_report.py: the calls, and the catalog as seen / not
                  seen / not callable

    MEGAnE's image has no ps, which Nextflow's task wrapper needs:
    resources/usr/bin/ps here gives the MEI task one (a module binary).

    Equivalent to: scripts/50-mei.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MEI_DEPTH {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/mei" }, mode: params.publish_dir_mode,
        saveAs: { f -> f.endsWith('_mei_depth.tsv') ? f : null }

    input:
    tuple val(meta), path(bam), path(bai)
    path(sites)   // assets/mei_sites.tsv

    output:
    tuple val(meta), path("${meta.id}_mei_mainchr.txt"), path("${meta.id}_mei_depth.tsv"), emit: depth
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The same commands as scripts/50-mei.sh, step 1.
    def p = meta.id
    """
    samtools view -H ${bam} | awk -F '\\t' '\$1 == "@SQ" {sub(/^SN:/, "", \$2); print \$2}' > ${p}_bam_contigs.txt
    grep -E '^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)\$' ${p}_bam_contigs.txt > ${p}_mei_mainchr.txt || true
    awk -F '\\t' 'NR == FNR {c[\$1]; next} !/^#/ && \$1 != "gene" && (\$2 in c) {print \$2 "\\t" \$3 - 1 "\\t" \$4 "\\t" \$1}' \\
        ${p}_bam_contigs.txt ${sites} > ${p}_mei_sites.bed
    if [ -s ${p}_mei_sites.bed ]; then samtools bedcov -c ${p}_mei_sites.bed ${bam} > ${p}_mei_depth.tsv; else : > ${p}_mei_depth.tsv; fi
    rm -f ${p}_bam_contigs.txt ${p}_mei_sites.bed
    [ -s ${p}_mei_mainchr.txt ] || { echo "ERROR: the BAM has none of chr1-chr22, chrX, chrY: nothing for MEGAnE to read." >&2; exit 1; }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    echo chr20 > ${meta.id}_mei_mainchr.txt
    touch ${meta.id}_mei_depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process MEI {
    tag "$meta.id"
    label 'process_high'

    publishDir { "${params.outdir}/${meta.id}/mei" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai), path(mainchr)
    path(fasta)
    path(fai)
    tuple path(kmers), path(kmers_count)   // reference/megane/<name>.mk and .mi

    output:
    tuple val(meta), path("${meta.id}_megane.vcf"), emit: vcf
    path "${meta.id}_megane.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // MEGAnE refuses a .fai that is not newer than its FASTA, and a staged
    // link carries the FASTA's own time: a fresh copy of the .fai sits beside
    // a link to the FASTA. MEGAnE's .mi count file must sit beside the .mk
    // under the same name, which the staging keeps. Its -mainchr comes after
    // the wrapper's own, so it wins. MEGAnE's single-sample VCF is
    // MEI_final_percentile_genotyped.vcf; when too few insertions are found
    // to set its filter thresholds (a BAM of one chromosome), it writes
    // MEI_final_failed_genotyped.vcf instead, with fixed thresholds.
    def sex = meta.sex in ['male', 'female'] ? meta.sex : 'unknown'
    """
    # MEGAnE stops when chrX is missing from the BAM unless told the sample has
    # no sex chromosomes; a BAM cut to some autosomes needs -no_sex_chr.
    nosex=""
    grep -qE '^(chr)?X\$' ${mainchr} || nosex="-no_sex_chr"
    ln -s \$(readlink -f ${fasta}) mei_ref.fa
    cp ${fai} mei_ref.fa.fai
    call_genotype_38 -i ${bam} -fa mei_ref.fa -mk ${kmers} -mainchr ${mainchr} -sex ${sex} \$nosex \\
        -sample_name ${meta.id} -p ${task.cpus} -only_ins -outdir megane
    vcf=megane/MEI_final_percentile_genotyped.vcf
    if [ ! -s \$vcf ] && [ -s megane/MEI_final_failed_genotyped.vcf ]; then
        vcf=megane/MEI_final_failed_genotyped.vcf
        echo "WARNING: MEGAnE found too few insertions to set its filter thresholds and used fixed ones (\$vcf)." >&2
    fi
    [ -s \$vcf ] || { echo "ERROR: MEGAnE wrote no MEI_final_percentile_genotyped.vcf" >&2; exit 1; }
    cp \$vcf ${meta.id}_megane.vcf
    cat megane/*.log > ${meta.id}_megane.log 2>/dev/null || true
    rm -rf megane

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megane: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf '##fileformat=VCFv4.1\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${meta.id}\\n' > ${meta.id}_megane.vcf
    touch ${meta.id}_megane.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megane: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process MEI_REPORT {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/mei" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(depth)
    path(sites)

    output:
    tuple val(meta), path("${meta.id}_mei.tsv"), path("${meta.id}_mei_sites.tsv"), emit: tables
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mei_report.py --sample ${meta.id} --vcf ${vcf} --sites ${sites} --depth ${depth} \\
        --out-calls ${meta.id}_mei.tsv --out-sites ${meta.id}_mei_sites.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf 'chrom\\tpos\\tend\\telement\\tsubfamily\\tcopies\\tfilter\\tconfidence\\n' > ${meta.id}_mei.tsv
    printf 'gene\\tchrom\\tstart\\tend\\tsite\\telement\\tinheritance\\tcondition\\tclinvar\\tpmid\\tstatus\\tmean_depth\\tcall\\n' > ${meta.id}_mei_sites.tsv
    printf 'FKTN\\tchr9\\t1\\t2\\t1\\tSVA\\tautosomal recessive\\tstub\\t.\\t.\\tnot callable\\t.\\t.\\n' >> ${meta.id}_mei_sites.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
