/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FAMILY_IBD_SITES, FAMILY_IBD_PHASE, FAMILY_IBD_SEGMENTS, FAMILY_IBD_REPORT
    — The DNA each pair of members shares by descent
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The FAMILY workflow runs these when --ancestry_fine_panel is given: the
    panel folder of step 26b (setup.sh --ancestry-fine-panel) holds the
    phased reference genotypes (panel/<chrom>.vcf.gz) and Beagle's genetic
    maps (maps/<chrom>.map) they need.

      FAMILY_IBD_SITES     on each autosome with the panel's SNVs and a map:
                           the joint call's genotypes at the panel's SNVs
                           (same position and alleles, bcftools isec), GT
                           only, half-missing genotypes made missing; a
                           chromosome with fewer than 100 is left out
      FAMILY_IBD_PHASE     Beagle phases them against the panel
      FAMILY_IBD_SEGMENTS  hap-ibd's segments (Zhou et al. 2020) on the
                           phased genotypes, its default 2 cM minimum
      FAMILY_IBD_REPORT    bin/ibd_report.py: ibd_pairs.tsv, ibd_segments.tsv
                           and ibd_report.md, published under family/

    Equivalent to: scripts/38-family.sh (its ibd_segments function)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FAMILY_IBD_SITES {
    label 'process_low'

    input:
    tuple path(vcf), path(vcf_index)
    path(panel)
    path(relationships)  // FAMILY_CHECK's: the segments wait for the relationship check

    output:
    path "ibd_in",       emit: sites
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p ibd_in
    for c in \$(bcftools index -s ${vcf} | cut -f1); do
        case "\$c" in chr[1-9]|chr1[0-9]|chr2[0-2]) ;; *) continue ;; esac
        [ -s "${panel}/panel/\$c.vcf.gz" ] && [ -s "${panel}/maps/\$c.map" ] || continue
        bcftools isec -n=2 -w1 -r "\$c" ${vcf} "${panel}/panel/\$c.vcf.gz" \\
            | bcftools annotate -x ID,INFO,^FORMAT/GT -Ou \\
            | bcftools +setGT -Oz -o "ibd_in/\$c.in.vcf.gz" -- -t . -n . > /dev/null
        n=\$(bcftools view -H "ibd_in/\$c.in.vcf.gz" | wc -l)
        if [ "\$n" -ge 100 ]; then echo "\$c" >> ibd_in/chroms.txt; else rm -f "ibd_in/\$c.in.vcf.gz"; fi
        echo "\$c: \$n of the panel SNVs in the joint call"
    done
    touch ibd_in/chroms.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ibd_in
    touch ibd_in/chroms.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process FAMILY_IBD_PHASE {
    label 'process_medium'

    input:
    path(sites)
    path(panel)

    output:
    path "phased",       emit: phased
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Beagle's heap: the task's memory less 2 GB for the JVM itself.
    def heap = Math.max(1, (task.memory ? task.memory.toGiga() : 8) - 2)
    """
    mkdir -p phased
    for c in \$(cat ${sites}/chroms.txt); do
        beagle -Xmx${heap}g gt=${sites}/\$c.in.vcf.gz chrom=\$c ref=${panel}/panel/\$c.vcf.gz map=${panel}/maps/\$c.map \\
            out=phased/\$c impute=false seed=1 nthreads=${task.cpus} > phased/\$c.beagle.out 2>&1 \\
            || { tail -n 20 phased/\$c.beagle.out; exit 1; }
        mv -f phased/\$c.log phased/\$c.beagle.log
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p phased
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process FAMILY_IBD_SEGMENTS {
    label 'process_low'

    input:
    path(phased)
    path(sites)
    path(panel)

    output:
    path "segments",     emit: segments
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def heap = Math.max(1, (task.memory ? task.memory.toGiga() : 4) - 1)
    """
    mkdir -p segments
    for c in \$(cat ${sites}/chroms.txt); do
        hap-ibd -Xmx${heap}g gt=${phased}/\$c.vcf.gz map=${panel}/maps/\$c.map out=segments/\$c nthreads=${task.cpus} \\
            > segments/\$c.hapibd.out 2>&1 || { tail -n 20 segments/\$c.hapibd.out; exit 1; }
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hap-ibd: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p segments
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hap-ibd: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process FAMILY_IBD_REPORT {
    label 'process_single'

    publishDir { "${params.outdir}/family" }, mode: params.publish_dir_mode

    input:
    path(pedigree, stageAs: 'pedigree.tsv')
    path(sites)
    path(phased)
    path(segments)
    path(panel)

    output:
    path "ibd_pairs.tsv",    emit: pairs
    path "ibd_segments.tsv", emit: segment_table
    path "ibd_report.md",    emit: report
    path "versions.yml",     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ibd_report.py --pedigree ${pedigree} --chroms ${sites}/chroms.txt --maps ${panel}/maps \\
        --phased ${phased} --ibd ${segments} --out-dir .

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    touch ibd_pairs.tsv ibd_segments.tsv ibd_report.md
    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}
