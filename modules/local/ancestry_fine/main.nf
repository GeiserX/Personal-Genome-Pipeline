/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANCESTRY_FINE — Fine-grained ancestry: admixture fractions and a chromosome painting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The five parts of scripts/26b-ancestry-fine.sh, one process per image,
    against the panel folder setup.sh --ancestry-fine-panel writes
    (--ancestry_fine_panel):

      ANCESTRY_FINE_GENOTYPE  the panel's SNVs in the sample, from the gVCF
                              (or from a VCF that holds reference calls)
      ANCESTRY_FINE_PREPARE   the one-person ADMIXTURE input, the groups'
                              frequencies at its SNVs and the chromosomes to
                              paint (bin/ancestry_fine_report.py)
      ANCESTRY_FINE_ADMIXTURE ADMIXTURE -P: the fractions at each level and
                              their block bootstrap
      ANCESTRY_FINE_PHASE     Beagle phasing against the panel
      ANCESTRY_FINE_PAINT     RFMix painting with the continental groups
      ANCESTRY_FINE_REPORT    the fractions table, the painting (BED, PNG) and
                              the text summary

    ADMIXTURE is free for academic use; RFMix for academic research use only,
    and any other use, a personal one included, needs a licence from Stanford.

    Equivalent to: scripts/26b-ancestry-fine.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ANCESTRY_FINE_GENOTYPE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(vcf_index), val(kind)   // kind: gvcf or vcf
    path panel
    path reference
    path reference_fai

    output:
    tuple val(meta), path("${meta.id}_genotypes.vcf.gz"), path("${meta.id}_genotypes.vcf.gz.tbi"), path("${meta.id}_counts.tsv"), emit: genotypes
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def id = meta.id
    def extract = kind == 'gvcf'
        ? "bcftools convert --gvcf2vcf -f ${reference} -R ${panel}/panel/positions.tsv -Ou ${vcf} | bcftools view -T ${panel}/panel/positions.tsv --trim-alt-alleles -i 'GT!=\"mis\"' -Ov"
        : "bcftools view -T ${panel}/panel/positions.tsv -Ov ${vcf}"
    """
    set -euo pipefail
    # As scripts/26b-ancestry-fine.sh: the awk reads the first sample column
    # only, so a file with more than one sample stops here.
    n=\$(bcftools query -l ${vcf} | wc -l | tr -d ' ')
    if [ "\$n" -ne 1 ]; then
        echo "ERROR: Sample '${id}': ${vcf} holds \$n samples; step 26b reads a one-sample file (keep one with bcftools view -s <sample>)." >&2
        exit 1
    fi
    # As scripts/26b-ancestry-fine.sh: a record is kept only with the panel's
    # REF and ALT (a 0/0 gets the panel's ALT), and the records are counted
    # per chromosome with the first and the last position.
    ${extract} \\
      | awk -F'\\t' -v OFS='\\t' -v sites="${panel}/panel/sites.tsv" -v counts="${id}_counts.tsv" '
          BEGIN { while ((getline l < sites) > 0) { split(l, f, "\\t"); ref[f[1] ":" f[2]] = f[3]; alt[f[1] ":" f[2]] = f[4] } }
          /^##/ { if (\$0 !~ /^##(INFO|FILTER|FORMAT)=/) print; next }
          /^#CHROM/ { print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\\"Genotype\\">"; NF = 10; print; next }
          {
            k = \$1 ":" \$2
            if (!(k in ref) || (k in seen) || \$4 != ref[k]) next
            if (\$5 == ".") \$5 = alt[k]
            else if (\$5 != alt[k]) next
            n = split(\$9, fmt, ":"); split(\$10, val, ":"); gt = "."
            for (i = 1; i <= n; i++) if (fmt[i] == "GT") gt = val[i]
            if (gt !~ /^[01][\\/|][01]\$/) next
            seen[k] = 1; if (!(\$1 in per)) first[\$1] = \$2
            per[\$1]++; last[\$1] = \$2; if (gt ~ /^0[\\/|]0\$/) homref++
            print \$1, \$2, ".", \$4, \$5, ".", ".", ".", "GT", gt
          }
          END { for (c in per) print c "\\t" per[c] "\\t" first[c] "\\t" last[c] > counts; print "homref\\t" homref + 0 > counts }' > ${id}_genotypes.vcf
    n=\$(awk -F'\\t' '\$1 != "homref" {n += \$2} END {print n + 0}' ${id}_counts.tsv)
    h=\$(awk -F'\\t' '\$1 == "homref" {print \$2}' ${id}_counts.tsv)
    echo "\${n} panel SNVs genotyped, \${h} of them as reference (0/0)"
    if [ "${kind}" = vcf ] && [ "\$(( h * 5 ))" -lt "\$n" ]; then
        echo "ERROR: ${vcf} has reference calls (0/0) at only \${h} of the \${n} panel SNVs it genotypes;" \\
             "ancestry_fine needs the gVCF (samplesheet columns gvcf,gvcf_index) or a VCF with reference calls." >&2
        exit 1
    fi
    bcftools view -Oz -o ${id}_genotypes.vcf.gz ${id}_genotypes.vcf
    bcftools index -t ${id}_genotypes.vcf.gz
    rm -f ${id}_genotypes.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_genotypes.vcf.gz ${meta.id}_genotypes.vcf.gz.tbi ${meta.id}_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process ANCESTRY_FINE_PREPARE {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(genotypes), path(genotypes_index), path(counts)
    path panel

    output:
    tuple val(meta), path("admixture"), emit: admixture_input
    tuple val(meta), path("${meta.id}_chroms.txt"), emit: chroms
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def boot = task.ext.bootstrap != null ? task.ext.bootstrap : 100
    """
    set -euo pipefail
    levels=()
    while IFS=\$'\\t' read -r level k _; do levels+=(--level "\${level}:\${k}:${panel}/model/\${level}.P"); done < ${panel}/model/levels.tsv
    ancestry_fine_report.py sample-bed --bim ${panel}/model/ref.bim --vcf ${genotypes} \\
        --sample ${meta.id} "\${levels[@]}" --bootstrap ${boot} --out-dir admixture
    # The chromosomes to phase and paint, as the bash step picks them.
    ancestry_fine_report.py chroms --counts ${counts} --maps ${panel}/maps > ${meta.id}_chroms.txt

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    mkdir -p admixture
    touch ${meta.id}_chroms.txt
    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}

process ANCESTRY_FINE_ADMIXTURE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(admixture_in, stageAs: 'in/*')

    output:
    tuple val(meta), path("admixture"), emit: admixture
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // run_admixture.sh (written by ANCESTRY_FINE_PREPARE): the projection of
    // each level, then its bootstrap replicates.
    """
    set -eu
    mkdir -p admixture
    cp -L in/admixture/* admixture/
    cd admixture
    sh run_admixture.sh ${task.cpus}
    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        admixture: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p admixture
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        admixture: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process ANCESTRY_FINE_PHASE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(genotypes), path(genotypes_index), path(chroms)
    path panel

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
        beagle -Xmx${heap}g gt=${genotypes} chrom=\$c ref=${panel}/panel/\$c.vcf.gz map=${panel}/maps/\$c.map \\
            out=phase/\$c impute=false seed=1 nthreads=${task.cpus} > phase/\$c.out 2>&1 \\
            || { tail -n 20 phase/\$c.out; exit 1; }
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p phase
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process ANCESTRY_FINE_PAINT {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(phased, stageAs: 'in/*')
    path panel

    output:
    tuple val(meta), path("paint"), emit: painting
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail
    mkdir -p paint
    for f in in/phase/chr*.vcf.gz; do
        [ -e "\$f" ] || continue   # no chromosome long enough to paint
        c=\$(basename "\$f" .vcf.gz)
        rfmix -f "\$f" -r ${panel}/panel/\$c.vcf.gz -m ${panel}/model/paint_map.tsv -g ${panel}/maps/\$c.rfmix.map \\
            -o paint/\$c --chromosome=\$c --n-threads=${task.cpus} --random-seed=1 > paint/\$c.log 2>&1 \\
            || { tail -n 20 paint/\$c.log; exit 1; }
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rfmix: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p paint
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rfmix: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}

process ANCESTRY_FINE_REPORT {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}/ancestry_fine" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(admixture, stageAs: 'adm/*'), path(painting, stageAs: 'pnt/*')
    path panel
    val panel_name

    output:
    tuple val(meta), path("${meta.id}_ancestry_fine.tsv"), emit: fractions
    tuple val(meta), path("${meta.id}_painting.bed"),      emit: painting_bed
    tuple val(meta), path("${meta.id}_painting.png"),      emit: painting_png
    tuple val(meta), path("${meta.id}_ancestry_fine.txt"), emit: text
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ancestry_fine_report.py report --sample ${meta.id} --model ${panel}/model \\
        --admixture adm/admixture --paint pnt/paint --panel ${panel_name} --out .

    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """

    stub:
    """
    touch ${meta.id}_ancestry_fine.tsv ${meta.id}_painting.bed ${meta.id}_painting.png ${meta.id}_ancestry_fine.txt
    printf '"%s":\\n    python: %s\\n' "${task.process}" "${task.container.replaceFirst(/^[^:@]+[:@]/, '')}" > versions.yml
    """
}
