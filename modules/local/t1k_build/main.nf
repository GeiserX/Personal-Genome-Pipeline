/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    T1K_BUILD — The T1K HLA index, built once per run for every sample
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    t1k-build.pl turns IPD-IMGT/HLA's hla.dat into the allele sequences
    (*dna_seq.fa) and their GRCh38 coordinates (*dna_coord.fa), taking each
    gene's position from a gene annotation (--hla_genes: the gene lines of
    GENCODE's basic GTF, which scripts/setup.sh installs). A FASTA there gives
    every gene "-1 -1" coordinates and T1K then extracts no reads, so the
    build stops when a typed gene has no coordinates.

    Its inputs are value channels, so the two files it writes are value
    channels too: one build serves every HLA_TYPING task of the run. The task
    hash covers the T1K image, hla.dat and the annotation, the key
    scripts/08-hla-typing.sh puts in its index directory name (T1K version,
    database release, GENCODE release), so a change of any builds a new index
    and -resume reuses an unchanged one.

    Equivalent to: step 1 of scripts/08-hla-typing.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process T1K_BUILD {
    label 'process_low'

    input:
    path(hla_dat)
    path(genes_gtf)

    output:
    path("hla_idx/*dna_seq.fa"),   emit: seq_fa
    path("hla_idx/*dna_coord.fa"), emit: coord_fa
    path "versions.yml",           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The genes the HLA report lists; each must have coordinates.
    def typed_genes = 'HLA-A HLA-B HLA-C HLA-DRB1 HLA-DQB1 HLA-DPB1'
    """
    t1k-build.pl \\
        -d ${hla_dat} \\
        -g ${genes_gtf} \\
        --prefix hla \\
        -o hla_idx

    COORD=""
    SEQ=""
    for f in hla_idx/*dna_coord.fa; do if [ -f "\$f" ]; then COORD=\$f; break; fi; done
    for f in hla_idx/*dna_seq.fa; do if [ -f "\$f" ]; then SEQ=\$f; break; fi; done
    if [ -z "\$COORD" ] || [ -z "\$SEQ" ]; then
        echo "ERROR: t1k-build.pl wrote no *dna_seq.fa and *dna_coord.fa in hla_idx/" >&2
        exit 1
    fi
    NO_COORD=""
    for g in ${typed_genes}; do
        if grep -q "^>\${g}\\*[^ ]* [^ ]* -1 -1 " "\$COORD" || ! grep -q "^>\${g}\\*" "\$COORD"; then
            NO_COORD="\${NO_COORD} \${g}"
        fi
    done
    if [ -n "\$NO_COORD" ]; then
        echo "ERROR: genes without GRCh38 coordinates in \$COORD:\${NO_COORD}" >&2
        echo "  T1K would extract no reads for them. Check --hla_genes (${genes_gtf}): it must be a gene annotation (GTF), not a FASTA." >&2
        exit 1
    fi
    echo "\$(grep -c ' -1 -1 ' "\$COORD" || true) alleles of genes outside the GRCh38 primary assembly have no coordinates (not typed)."

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p hla_idx
    touch hla_idx/hla_dna_seq.fa hla_idx/hla_dna_coord.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
