/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARASCOPY — paralog copy number (SMN1 and SMN2 first), HBA1/HBA2 from depth (opt-in: --tools parascopy)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SMN1 and SMN2 are near-identical copies: short reads cannot be placed on
    one of them, and no VCF-based step sees an SMN1 copy-number loss, the
    usual cause of SMA carrier status. Parascopy (MIT licence) estimates the
    aggregate copy number of the pair (agCN) and, where paralogous sequence
    variants allow, the copy number of each (psCN), each with a Phred
    quality. It reads a BAM aligned to a reference without ALT contigs (the
    default) and uses the precomputed GRCh38 homology table and the 1000
    Genomes model parameters `setup.sh --parascopy-data` installs
    (--parascopy_data; --parascopy_population picks the model, EUR by
    default).

    The loci come from assets/paralog_loci.tsv, one `parascopy cn-using` per
    row (paralogSpecs below reads its model, region and table columns). A row
    is not assessed, with the reason in its table, when its model is not
    installed, when the BAM has no read over the row's region, or when
    Parascopy fails on it; the other rows still run.

    Background depth comes from Parascopy's own GRCh38 windows, or from
    --parascopy_depth_bed (windows of one size) for a BAM that covers only
    part of the genome, without GC stratification (--no-gc): a few regions
    rarely span the GC range Parascopy's GC model needs, and it stops on them.

    HBA1 and HBA2 have no Parascopy model: samtools bedcov (in the same
    image) sums the depth over both gene bodies and two single-copy flanks,
    over all reads and over reads of MAPQ 1 or more, and bin/copy_number.py
    (the HBA row of bin/cn_loci.tsv) turns it into a count of alpha-globin
    genes (not assessed when no read over the genes has MAPQ 0: the aligner then left
    out the reads that fit both genes).
    A failure there writes an HBA table that says not assessed (so a
    published table from an earlier run is replaced) and does not fail the
    task.

    The gene-dosage rows of the same catalog (HP, LCE3BC, CCL3L1: loci
    no pinned caller models) go the same way, total copies only, one
    line per locus in <id>_cn_dosage.tsv (docs/cn-loci.md).

    CYP21A2 has no Parascopy model either: bedcov over CYP21A2 and its
    pseudogene CYP21A1P gives the RCCX module total, and samtools mpileup at
    the sites where the two differ gives the CYP21A2 copies at each, so the
    30 kb chimeric deletion and a possible conversion show
    (<id>_cyp21a2_psv.tsv, docs/cn-loci.md).

    The interval rows of that catalog (ClinGen's recurrent CNV regions, the
    carrier deletion alleles and the Y chromosome AZF regions) take one more
    bedcov pass over all their segments, one line per row in
    <id>_cn_intervals.tsv; the chrX and chrY rows read the samplesheet's sex
    (meta.sex; without it they are not assessed).

    Equivalent to: scripts/35-paralogs.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// One spec per row of assets/paralog_loci.tsv, model|region|table, space
// separated (as scripts/35-paralogs.sh reads it). tests/test_paralog_loci.py
// keeps these three columns to characters a shell word can carry.
def paralogSpecs() {
    def lines = file("${projectDir}/assets/paralog_loci.tsv").readLines().findAll { l -> l.trim() && !l.startsWith('#') }
    def head = lines[0].split('\t').toList()
    return lines.drop(1).collect { l ->
        def c = l.split('\t')
        ['model', 'region', 'table'].collect { k -> c[head.indexOf(k)] }.join('|')
    }.join(' ')
}

process PARASCOPY {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/paralogs" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_fai)
    path(parascopy_data)  // setup.sh --parascopy-data: homology_table/ and models_GRCh38_1KGP/
    path(depth_bed)       // background windows, or [] for Parascopy's GRCh38 windows

    output:
    tuple val(meta), path("${meta.id}_smn_copy_number.tsv"), emit: copy_number, optional: true  // the SMN1 row's
    tuple val(meta), path("${meta.id}_parascopy"),           emit: results
    tuple val(meta), path("${meta.id}_hba_copy_number.tsv"), emit: hba
    tuple val(meta), path("${meta.id}_*_copy_number.tsv"),   emit: tables  // every catalog row's table, and HBA's
    tuple val(meta), path("${meta.id}_cn_dosage.tsv"),       emit: dosage  // HP, LCE3BC, CCL3L1 from depth
    tuple val(meta), path("${meta.id}_cyp21a2_psv.tsv"),     emit: cyp21a2  // CYP21A2 from depth and PSVs
    tuple val(meta), path("${meta.id}_cn_intervals.tsv"),    emit: intervals  // the interval rows from depth
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def background = depth_bed ? "-b ${depth_bed} --no-gc" : '-g GRCh38'
    def loci = paralogSpecs()
    def sex = meta.sex in ['male', 'female'] ? meta.sex : 'unknown'
    """
    # Paralog loci (as scripts/35-paralogs.sh): the background depth once, then
    # one cn-using per catalog row, one table per row with one row per region of
    # the locus profile (aggregate and paralog-specific copy number, each with
    # filter and quality).
    # A BAM with ALT or HLA contigs is refused, as scripts/35-paralogs.sh does
    # (bin/copy_number.py check-bam): its depth splits between a gene and its
    # ALT copy, and a depth-based count reads that as a deletion.
    samtools view -H ${bam} > bam_header.txt
    copy_number.py check-bam --header bam_header.txt
    models=${parascopy_data}/models_GRCh38_1KGP/${params.parascopy_population}
    [ -d "\$models" ] || { echo "ERROR: no Parascopy models \$models (populations: AFR AMR EAS EUR SAS)" >&2; exit 1; }
    parascopy depth -i ${bam}::${meta.id} -f ${reference} ${background} -o depth -@ ${task.cpus}
    header() { printf 'chrom\\tstart\\tend\\tlocus\\tagCN_filter\\tagCN\\tagCN_qual\\tpsCN_filter\\tpsCN\\tpsCN_qual\\thomologous_regions\\n'; }
    mkdir -p ${meta.id}_parascopy
    loci='${loci}'  # model|region|table per row; quoted, as | would start a pipe
    for spec in \$loci; do
        model=\${spec%%|*} rest=\${spec#*|}
        region=\${rest%%|*} name=\${rest#*|}
        out=${meta.id}_\${name}_copy_number.tsv
        reason="" reads=x
        if [ -f "\$models/\$model" ]; then
            reads=\$(samtools view -c -F 0x904 ${bam} "\$region") || reads=x
        fi
        if [ ! -f "\$models/\$model" ]; then
            reason="no model \$model in the installed Parascopy data (${params.parascopy_population})"
        elif ! [ "\$reads" -ge 0 ] 2>/dev/null; then
            reason="could not count reads over \$region (samtools failed)"
        elif [ "\$reads" -eq 0 ]; then
            reason="no read over \$region"
        elif parascopy cn-using "\$models/\$model" -i ${bam}::${meta.id} -f ${reference} \\
                -t ${parascopy_data}/homology_table/GRCh38.bed.gz -d depth -o ${meta.id}_parascopy/\$name -@ ${task.cpus}; then
            # no result file although Parascopy exited 0: not assessed, as in step 35
            if {
                header
                gzip -dc ${meta.id}_parascopy/\$name/res.samples.bed.gz \\
                    | awk -F'\\t' -v OFS='\\t' '!/^#/ {print \$1, \$2, \$3, \$4, \$6, \$7, \$8, \$9, \$10, \$11, \$13}'
            } > "\$out"; then
                continue
            fi
            reason="parascopy cn-using wrote no readable result; see the log of step 35"
        else
            reason="parascopy cn-using failed; see the log of step 35"
        fi
        echo "WARNING: \$model: not assessed: \$reason" >&2
        rm -rf ${meta.id}_parascopy/\$name
        chrom=\${region%%:*} span=\${region#*:}
        { header; printf '%s\\t%s\\t%s\\t%s\\t%s\\t.\\t.\\t%s\\t.\\t.\\t.\\n' "\$chrom" "\$((\${span%-*} - 1))" "\${span#*-}" \\
            "\${model%.gz}" "not assessed: \$reason" "not assessed: \$reason"; } > "\$out"
    done
    rm -rf depth
    # end of the paralog loci

    # HBA1/HBA2 from depth (as scripts/35-paralogs.sh): all MAPQ and MAPQ 1 or more, then the verdict
    if copy_number.py bed --locus HBA > hba_regions.bed \\
        && samtools bedcov hba_regions.bed ${bam} > hba_bedcov.tsv \\
        && samtools bedcov -Q 1 hba_regions.bed ${bam} > hba_bedcov_q1.tsv \\
        && copy_number.py call --locus HBA --bedcov hba_bedcov.tsv --bedcov-mapq1 hba_bedcov_q1.tsv \\
            --table ${meta.id}_hba_copy_number.tsv.tmp; then
        mv ${meta.id}_hba_copy_number.tsv.tmp ${meta.id}_hba_copy_number.tsv
    else
        rm -f ${meta.id}_hba_copy_number.tsv.tmp
        echo "WARNING: the HBA1/HBA2 depth check failed; the report lists HBA1/HBA2 as not assessed." >&2
        printf 'metric\\tvalue\\nstatus\\tnot assessed\\nalpha_globin_genes\\t.\\nestimate\\t.\\nhba2_depth\\t.\\nhba1_depth\\t.\\nflank_left_depth\\t.\\nflank_right_depth\\t.\\nmapq0_share\\t.\\nmessage\\t%s\\n' \\
            'not assessed: the HBA1/HBA2 depth check failed; see the log of step 35' > ${meta.id}_hba_copy_number.tsv
    fi
    rm -f hba_regions.bed hba_bedcov.tsv hba_bedcov_q1.tsv

    # Gene-dosage rows of bin/cn_loci.tsv (as scripts/35-paralogs.sh): per locus all MAPQ and MAPQ 1
    # or more, then the verdict; one line per locus, not assessed when its check fails
    dosage_loci="HP LCE3BC CCL3L1"
    printf 'locus\\tstatus\\tcopies\\testimate\\treason\\treading\\n' > cn_dosage.part
    for loc in \$dosage_loci; do
        if copy_number.py bed --locus "\$loc" > dosage_regions.bed \\
            && samtools bedcov dosage_regions.bed ${bam} > dosage_bedcov.tsv \\
            && samtools bedcov -Q 1 dosage_regions.bed ${bam} > dosage_bedcov_q1.tsv \\
            && copy_number.py call --locus "\$loc" --bedcov dosage_bedcov.tsv --bedcov-mapq1 dosage_bedcov_q1.tsv \\
                --table dosage_table.tsv; then
            awk -F'\\t' -v OFS='\\t' -v l="\$loc" 'NR > 1 {v[\$1] = \$2}
                END {print l, v["status"], v["copies"], v["estimate"], v["reason"], v["reading"]}' dosage_table.tsv >> cn_dosage.part
        else
            echo "WARNING: the \$loc depth check failed; ${meta.id}_cn_dosage.tsv lists it as not assessed." >&2
            printf '%s\\tnot assessed\\t.\\t.\\tcheck_failed\\tthe depth check failed; see the log of step 35\\n' "\$loc" >> cn_dosage.part
        fi
        rm -f dosage_regions.bed dosage_bedcov.tsv dosage_bedcov_q1.tsv dosage_table.tsv
    done
    mv cn_dosage.part ${meta.id}_cn_dosage.tsv

    # CYP21A2 from depth and PSVs (as scripts/35-paralogs.sh): bedcov over the row (the RCCX module
    # total), mpileup without a reference at its PSV sites, then the verdict; not assessed when it fails
    if copy_number.py bed --locus CYP21A2 > cyp21a2_regions.bed \\
        && copy_number.py bed --psv --locus CYP21A2 > cyp21a2_psv.bed \\
        && samtools bedcov cyp21a2_regions.bed ${bam} > cyp21a2_bedcov.tsv \\
        && samtools mpileup -A -B -q 0 -Q 13 -d 100000 \\
            -r "\$(awk 'NR == 1 {c = \$1; s = \$2 + 1} {e = \$3} END {print c ":" s "-" e}' cyp21a2_psv.bed)" \\
            -l cyp21a2_psv.bed ${bam} > cyp21a2_pileup.tsv \\
        && copy_number.py call --locus CYP21A2 --bedcov cyp21a2_bedcov.tsv --pileup cyp21a2_pileup.tsv \\
            --table ${meta.id}_cyp21a2_psv.tsv.tmp; then
        mv ${meta.id}_cyp21a2_psv.tsv.tmp ${meta.id}_cyp21a2_psv.tsv
    else
        rm -f ${meta.id}_cyp21a2_psv.tsv.tmp
        echo "WARNING: the CYP21A2 depth and PSV check failed; ${meta.id}_cyp21a2_psv.tsv says not assessed." >&2
        printf 'metric\\tvalue\\nstatus\\tnot assessed\\ncyp21a2_copies\\t.\\nmodule_total\\t.\\nexpected_total\\t.\\nestimate\\t.\\nfinding\\t.\\nconversion_sites\\t.\\ntarget_depth\\t.\\ncontrol_depth\\t.\\nreason\\tcheck_failed\\nreading\\t%s\\n' \\
            'the CYP21A2 depth and PSV check failed; see the log of step 35' > ${meta.id}_cyp21a2_psv.tsv
    fi
    rm -f cyp21a2_regions.bed cyp21a2_psv.bed cyp21a2_bedcov.tsv cyp21a2_pileup.tsv

    # Interval rows of bin/cn_loci.tsv (as scripts/35-paralogs.sh): one bedcov over every row's segments,
    # then the verdict; a failed pass writes every row as not assessed (check_failed)
    if ! { copy_number.py bed --kind interval > interval_regions.bed \
            && samtools bedcov interval_regions.bed ${bam} > interval_bedcov.tsv \
            && copy_number.py intervals --bedcov interval_bedcov.tsv --sex ${sex} --out ${meta.id}_cn_intervals.tsv; }; then
        echo "WARNING: the interval depth check failed; ${meta.id}_cn_intervals.tsv lists every row as not assessed." >&2
        copy_number.py intervals --failed --out ${meta.id}_cn_intervals.tsv
    fi
    rm -f interval_regions.bed interval_bedcov.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parascopy: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    def loci = paralogSpecs()
    """
    mkdir -p ${meta.id}_parascopy
    loci='${loci}'  # model|region|table per row; quoted, as | would start a pipe
    for spec in \$loci; do
        printf 'chrom\\tstart\\tend\\tlocus\\tagCN_filter\\tagCN\\tagCN_qual\\tpsCN_filter\\tpsCN\\tpsCN_qual\\thomologous_regions\\n' \\
            > ${meta.id}_\${spec##*|}_copy_number.tsv
    done
    printf 'metric\\tvalue\\nstatus\\tnot assessed\\nalpha_globin_genes\\t.\\nestimate\\t.\\nhba2_depth\\t0\\nhba1_depth\\t0\\nflank_left_depth\\t0\\nflank_right_depth\\t0\\nmapq0_share\\t0\\nmessage\\tstub\\n' \\
        > ${meta.id}_hba_copy_number.tsv
    printf 'locus\\tstatus\\tcopies\\testimate\\treason\\treading\\n' > ${meta.id}_cn_dosage.tsv
    printf 'metric\\tvalue\\nstatus\\tnot assessed\\nreason\\tstub\\n' > ${meta.id}_cyp21a2_psv.tsv
    printf 'locus\\tname\\tchrom\\tstart\\tend\\tstatus\\tcopies\\texpected\\tchange\\testimate\\treason\\n' > ${meta.id}_cn_intervals.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parascopy: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
