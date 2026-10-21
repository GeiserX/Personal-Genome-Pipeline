/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HLA_TYPING — HLA allele typing from WGS BAM using T1K
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Types HLA-A, B, C (Class I) and DRB1, DQB1, DPB1 (Class II)
    at 4-digit resolution using IPD-IMGT/HLA database against GRCh38.

    The index (allele sequences and their GRCh38 coordinates) comes from
    T1K_BUILD, built once per run for every sample; this process types one
    BAM against it.

    Equivalent to: step 2 of scripts/08-hla-typing.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process HLA_TYPING {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/hla" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(seq_fa)    // T1K_BUILD: allele sequences
    path(coord_fa)  // T1K_BUILD: their GRCh38 coordinates

    output:
    tuple val(meta), path("*_hla_genotype.tsv"), emit: hla_alleles
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    run-t1k \\
        -b ${bam} \\
        -f ${seq_fa} \\
        -c ${coord_fa} \\
        --preset hla-wgs \\
        -t ${task.cpus} \\
        --od ./ \\
        -o ${prefix}_hla

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_hla_genotype.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KIR_BUILD and KIR_TYPING — KIR genes, a second T1K pass (opt-in: --kir)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KIR_BUILD turns IPD-KIR's kir.dat (--kir_dat, release KIR_DB_RELEASE that
    `setup.sh --kir-data` installs) into T1K's allele sequences and GRCh38
    coordinates, taken from the same GENCODE gene lines as HLA. Several KIR
    genes are not on the GRCh38 primary assembly and get no coordinates; T1K
    still types them from the reads of the genes that are. Built once per run.
    KIR_TYPING runs T1K with its kir-wgs preset and writes the release beside
    the genotypes; a BAM with too few KIR reads gets a genotype file that says
    so instead of a failed task.

    Equivalent to: the KIR=true pass of scripts/08-hla-typing.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KIR_BUILD {
    label 'process_low'

    input:
    path(kir_dat)
    path(genes_gtf)

    output:
    path("kir_idx/*dna_seq.fa"),   emit: seq_fa
    path("kir_idx/*dna_coord.fa"), emit: coord_fa
    path("kir_release.txt"),       emit: release
    path "versions.yml",           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    grep -m 1 'IPD-KIR Release Version' ${kir_dat} | sed 's/^CC *//' > kir_release.txt
    if [ ! -s kir_release.txt ]; then
        echo "ERROR: ${kir_dat} names no IPD-KIR release (no 'IPD-KIR Release Version' line); is it IPD-KIR's kir.dat?" >&2
        exit 1
    fi
    t1k-build.pl \\
        -d ${kir_dat} \\
        -g ${genes_gtf} \\
        --prefix kir \\
        -o kir_idx
    for f in kir_idx/*dna_seq.fa kir_idx/*dna_coord.fa; do
        [ -f "\$f" ] || { echo "ERROR: t1k-build.pl wrote no *dna_seq.fa and *dna_coord.fa in kir_idx/" >&2; exit 1; }
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p kir_idx
    touch kir_idx/kir_dna_seq.fa kir_idx/kir_dna_coord.fa
    echo 'IPD-KIR Release Version stub' > kir_release.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process KIR_TYPING {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/kir" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(seq_fa)    // KIR_BUILD: allele sequences
    path(coord_fa)  // KIR_BUILD: their GRCh38 coordinates
    path(release)   // KIR_BUILD: the IPD-KIR release line of kir.dat

    output:
    tuple val(meta), path("${meta.id}_kir_genotype.tsv"), path("database_release.txt"), emit: kir_genotype
    path "versions.yml",                                                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # T1K may stop when it extracts no KIR read at all: that is "too few
    # reads" when its candidate reads file exists and is empty; any other
    # failure (no candidate file, or candidate reads it failed to type) is an error.
    rc=0
    run-t1k \\
        -b ${bam} \\
        -f ${seq_fa} \\
        -c ${coord_fa} \\
        --preset kir-wgs \\
        -t ${task.cpus} \\
        --od ./ \\
        -o ${prefix}_kir_t1k || rc=\$?
    if [ "\$rc" -ne 0 ]; then
        cand=""
        for c in ${prefix}_kir_t1k_candidate_1.fq ${prefix}_kir_t1k_candidate.fq; do
            if [ -e "\$c" ]; then cand=\$c; break; fi
        done
        if [ -z "\$cand" ] || [ -s "\$cand" ]; then
            echo "ERROR: run-t1k failed (exit \$rc); candidate reads: \${cand:-none written}" >&2
            exit "\$rc"
        fi
        echo "run-t1k exited \$rc and extracted no KIR read (\$cand is empty)"
    fi
    { echo "database: \$(cat ${release})"; echo "t1k: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}"; } > database_release.txt
    if [ -s ${prefix}_kir_t1k_genotype.tsv ] && awk -F'\\t' '\$5 > 0 || \$8 > 0 {found = 1} END {exit !found}' ${prefix}_kir_t1k_genotype.tsv; then
        cp ${prefix}_kir_t1k_genotype.tsv ${meta.id}_kir_genotype.tsv
    else
        printf '# KIR not typed: T1K found too few reads at the KIR genes (chr19 leukocyte receptor complex) in this BAM\\n' \\
            > ${meta.id}_kir_genotype.tsv
        echo "KIR: too few reads at the KIR genes to type any of them."
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_kir_genotype.tsv database_release.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HLA_QC_READS and HLA_QC — OptiType, a second HLA caller (opt-in: --hla_qc)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HLA_QC_READS takes the read pairs of HLA-A, HLA-B and HLA-C from the BAM:
    each gene's span in T1K_BUILD's coordinate file, 1 kb either side.
    HLA_QC types them with OptiType against the HLA sequences bundled in its
    image and writes the two-field concordance of T1K and OptiType
    (pgx_outside_calls.py hla-concordance). PGX_CONSENSUS passes HLA-A and
    HLA-B to PharmCAT only when the two agree. A BAM with no read pairs
    there gets a result file that says so, which PGX_CONSENSUS reads as no
    OptiType call.

    Equivalent to: the HLA_QC=true pass of scripts/08-hla-typing.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process HLA_QC_READS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)
    path(coord_fa)  // T1K_BUILD: the GRCh38 coordinates of the allele sequences

    output:
    tuple val(meta), path("${meta.id}_hla_R1.fq"), path("${meta.id}_hla_R2.fq"), emit: reads
    path "versions.yml",                                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk '/^>HLA-[ABC]\\*/ && \$3 > 0 {g = substr(\$1, 2); sub(/\\*.*/, "", g)
        if (!(g in lo) || \$3 < lo[g]) lo[g] = \$3; if (\$4 > hi[g]) hi[g] = \$4; c[g] = \$2}
      END {for (g in lo) printf "%s:%d-%d\\n", c[g], (lo[g] > 1000 ? lo[g] - 1000 : 1), hi[g] + 1000}' ${coord_fa} | sort > regions.txt
    if [ "\$(wc -l < regions.txt)" -ne 3 ]; then
        echo "ERROR: ${coord_fa} gives no GRCh38 span for each of HLA-A, HLA-B and HLA-C: \$(tr '\\n' ' ' < regions.txt)" >&2
        exit 1
    fi
    mkdir -p tmp
    samtools view -u -F 0x900 ${bam} \$(cat regions.txt) \\
        | samtools collate -u -O - tmp/collate \\
        | samtools fastq -n -1 ${meta.id}_hla_R1.fq -2 ${meta.id}_hla_R2.fq -0 /dev/null -s /dev/null -
    echo "read pairs at HLA-A, HLA-B and HLA-C: \$(( \$(wc -l < ${meta.id}_hla_R1.fq) / 4 ))"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_hla_R1.fq ${meta.id}_hla_R2.fq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process HLA_QC {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/hla_optitype" }, mode: params.publish_dir_mode

    input:
    // HLA_QC_READS' read pairs joined with HLA_TYPING's genotypes
    tuple val(meta), path(r1), path(r2), path(t1k_genotype)

    output:
    tuple val(meta), path("${meta.id}_optitype_result.tsv"), emit: result
    tuple val(meta), path("${meta.id}_hla_concordance.tsv"), emit: concordance, optional: true
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if [ ! -s ${r1} ]; then
        printf '# OptiType not run: no read pairs at HLA-A, HLA-B or HLA-C in this BAM\\n' > ${meta.id}_optitype_result.tsv
        echo "OptiType: no reads to type; PGX_CONSENSUS will not pass HLA-A or HLA-B to PharmCAT."
    else
        optitype run -i ${r1} -i ${r2} --dna --solver glpk --threads ${task.cpus} -o out -p ${meta.id}
        out=\$(find out -name '*_result.tsv' | head -n 1)
        [ -n "\$out" ] || { echo "ERROR: OptiType wrote no *_result.tsv" >&2; exit 1; }
        cp "\$out" ${meta.id}_optitype_result.tsv
        pgx_outside_calls.py hla-concordance --t1k ${t1k_genotype} \\
            --optitype ${meta.id}_optitype_result.tsv --out ${meta.id}_hla_concordance.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        optitype: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf '\\tA1\\tA2\\tB1\\tB2\\tC1\\tC2\\tReads\\tObjective\\n' > ${meta.id}_optitype_result.tsv
    printf 'Gene\\tT1K\\tOptiType\\tAgree\\tNote\\n' > ${meta.id}_hla_concordance.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        optitype: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
