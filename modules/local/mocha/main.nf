/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MOCHA — Mosaic chromosomal alterations: large gains, losses and CN-LOH (step 51)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Opt-in ('mocha' in --tools, with --mocha_data: setup.sh --mocha-data). Four
    processes, the commands of scripts/51-mocha.sh:
      MOCHA_PREP    the VCF's biallelic PASS SNVs of chr1 to chr22 with GT and
                    AD, a genotype under 10 reads or GQ 20 set missing, INFO/GC
                    from the reference (bcftools +mochatools)
      MOCHA_PHASE   Beagle against the panel, one chromosome at a time (those
                    the panel has)
      MOCHA_CALL    the phased genotypes onto the prepared sites, MoChA
                    (bcftools +mocha), and the per-chromosome site counts
      MOCHA_REPORT  bin/mocha_report.py: the calls above the floor

    Equivalent to: scripts/51-mocha.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MOCHA_PREP {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(idx)
    path(reference)
    path(reference_fai)

    output:
    tuple val(meta), path("prep.bcf"), path("prep.bcf.csi"), emit: prep
    tuple val(meta), path("phase_in.vcf.gz"), path("phase_in.vcf.gz.tbi"), path("chroms.txt"), emit: phase_in
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def auto = (1..22).collect { n -> "chr${n}" }.join(',')
    """
    set -euo pipefail
    bcftools view --no-version -f PASS,. -m2 -M2 -v snps -t ${auto} -Ou ${vcf} \\
      | bcftools filter --no-version -e "FMT/DP<10 | FMT/GQ<20" --set-GT . -Ou \\
      | bcftools annotate --no-version -x ID,QUAL,FILTER,INFO,^FMT/GT,^FMT/AD -Ou \\
      | bcftools +mochatools --no-version -Ob -o prep.bcf - -- -t GC -f ${reference}
    bcftools index -f prep.bcf
    bcftools view --no-version -e 'GT="mis"' -Oz -o phase_in.vcf.gz prep.bcf
    bcftools index -f -t phase_in.vcf.gz
    bcftools index -s prep.bcf | cut -f1 > chroms.txt
    [ -s chroms.txt ] || { echo "ERROR: no biallelic PASS SNV on chr1 to chr22 in ${vcf}" >&2; exit 1; }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mocha: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch prep.bcf prep.bcf.csi phase_in.vcf.gz phase_in.vcf.gz.tbi
    echo chr20 > chroms.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mocha: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process MOCHA_PHASE {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(phase_in), path(phase_in_tbi), path(chroms)
    path(mocha_data)   // reference/mocha-<version>: panel/chrN.vcf.gz and maps/chrN.map

    output:
    tuple val(meta), path("phase"), emit: phased
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Beagle's heap: the task's memory less 2 GB for the JVM itself.
    def heap = Math.max(1, (task.memory ? task.memory.toGiga() : 8) - 2)
    """
    set -euo pipefail
    mkdir -p phase
    for c in \$(cat ${chroms}); do
        [ -s ${mocha_data}/panel/\$c.vcf.gz ] || continue
        beagle -Xmx${heap}g gt=${phase_in} chrom=\$c ref=${mocha_data}/panel/\$c.vcf.gz map=${mocha_data}/maps/\$c.map \\
            out=phase/\$c impute=false seed=1 nthreads=${task.cpus} > phase/\$c.out 2>&1 \\
            || { tail -n 20 phase/\$c.out; exit 1; }
        echo \$c >> phase/chroms.txt
    done
    [ -s phase/chroms.txt ] || { echo "ERROR: the panel has none of the VCF's chromosomes (\$(paste -sd' ' ${chroms}))" >&2; exit 1; }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p phase
    echo chr20 > phase/chroms.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process MOCHA_CALL {
    tag "$meta.id"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/mosaic" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(prep), path(prep_csi), path(phase)

    output:
    tuple val(meta), path("${meta.id}_mocha_calls_raw.tsv"), path("${meta.id}_mocha_stats.tsv"),
          path("${meta.id}_mocha_sites.tsv"), emit: tables
    path "${meta.id}_mocha.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // MoChA's GRCh38 MHC and KIR regions (README of freeseek/mocha), as in scripts/51-mocha.sh.
    def p = meta.id
    """
    set -euo pipefail
    chroms=\$(paste -sd, ${phase}/chroms.txt)
    # Beagle writes no ##contig lines and concat refuses records on an
    # undeclared contig: give each Beagle VCF the contig lines of prep.bcf.
    bcftools view -h ${prep} | grep '^##contig=' > contigs.hdr
    f=(); for c in \$(cat ${phase}/chroms.txt); do
        { bcftools view -h ${phase}/\$c.vcf.gz | grep '^##' | grep -v '^##contig='
          cat contigs.hdr; bcftools view -h ${phase}/\$c.vcf.gz | grep '^#CHROM'; } > \$c.hdr
        bcftools reheader -h \$c.hdr -o \$c.ctg.vcf.gz ${phase}/\$c.vcf.gz
        f+=("\$c.ctg.vcf.gz")
    done
    bcftools concat --no-version -Ob -o phased_gt.bcf "\${f[@]}"
    bcftools index -f phased_gt.bcf
    bcftools annotate --no-version -a phased_gt.bcf -c FMT/GT -r "\$chroms" -Ob -o phased.bcf ${prep}
    bcftools index -f phased.bcf
    bcftools +mocha --genome GRCh38 --no-version --mhc chr6:27518932-33480487 --kir chr19:54071493-54992731 \\
      --calls ${p}_mocha_calls_raw.tsv --stats ${p}_mocha_stats.tsv --log ${p}_mocha.log phased.bcf
    bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n' phased.bcf | awk -F '\\t' -v OFS='\\t' '
      \$1 != c { if (c != "") print c, a, b, n, h, ph; c = \$1; a = \$2; n = h = ph = 0 }
      { b = \$2; n++; if (\$3 ~ /^(0[|\\/]1|1[|\\/]0)\$/) { h++; if (\$3 ~ /[|]/) ph++ } }
      END { if (c != "") print c, a, b, n, h, ph }' > ${p}_mocha_sites.tsv
    rm -f phased_gt.bcf phased_gt.bcf.csi phased.bcf phased.bcf.csi contigs.hdr chr*.hdr chr*.ctg.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mocha: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    def p = meta.id
    """
    touch ${p}_mocha_calls_raw.tsv ${p}_mocha_stats.tsv ${p}_mocha_sites.tsv ${p}_mocha.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mocha: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process MOCHA_REPORT {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/mosaic" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(calls), path(stats), path(sites)

    output:
    tuple val(meta), path("${meta.id}_mocha.tsv"), path("${meta.id}_mocha_summary.tsv"), emit: tables
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mocha_report.py --sample ${meta.id} --calls ${calls} --stats ${stats} --sites ${sites} --phasing beagle \\
        --out-calls ${meta.id}_mocha.tsv --out-summary ${meta.id}_mocha_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf 'chrom\\tstart\\tend\\tlength_mb\\tarms\\ttype\\tcell_fraction\\tbdev\\trel_cov\\tn_hets\\tlod_baf_phase\\tlod_lrr_baf\\treported\\treason\\n' > ${meta.id}_mocha.tsv
    printf 'metric\\tvalue\\nstatus\\tnot assessed\\nmessage\\tstub\\n' > ${meta.id}_mocha_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
