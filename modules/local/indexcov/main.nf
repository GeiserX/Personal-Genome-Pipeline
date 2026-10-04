/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INDEXCOV — Coverage QC and the sex check, from the BAM index alone
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    goleft indexcov reads only the .bai (seconds per sample) and writes
    per-chromosome coverage plots and a .ped row with the sex it infers from
    the chrX and chrY copy numbers.

    The process only reports: <sample>_sex_check.tsv holds the inferred sex
    (male, female or unknown) and CNchrX and CNchrY, read from the .ped by
    column name. main.nf compares it with the samplesheet's sex and stops the
    run on a mismatch (--sex_check warn logs it instead), before DeepVariant
    or any other BAM step starts, since they take the declared sex.

    Equivalent to: scripts/16-indexcov.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process INDEXCOV {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.id}" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_sex_check.tsv"), emit: sex
    tuple val(meta), path("indexcov"),                 emit: indexcov
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    goleft indexcov --directory indexcov ${bam}
    PED=indexcov/indexcov-indexcov.ped
    if [ ! -s "\$PED" ]; then
        echo "ERROR: goleft wrote no \$PED" >&2
        exit 1
    fi

    # goleft's .ped columns: #family_id sample_id paternal_id maternal_id sex
    # phenotype CNchrX CNchrY ...; sex in PED coding (1 male, 2 female).
    awk -v OFS='\\t' '
        /^#/ { for (i = 1; i <= NF; i++) { h = \$i; sub(/^#/, "", h); col[h] = i } next }
        { split(\$0, f) }
        END {
            if (!("sex" in col)) { print "ERROR: no sex column in the .ped header" > "/dev/stderr"; exit 1 }
            s = f[col["sex"]]
            sex = (s == "1") ? "male" : (s == "2") ? "female" : "unknown"
            x = ("CNchrX" in col) ? f[col["CNchrX"]] : "NA"
            y = ("CNchrY" in col) ? f[col["CNchrY"]] : "NA"
            print "inferred_sex", "CNchrX", "CNchrY"
            print sex, x, y
        }' "\$PED" > ${meta.id}_sex_check.tsv
    cat ${meta.id}_sex_check.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        goleft: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    // The stub agrees with the samplesheet, so a stub run never stops here.
    """
    mkdir -p indexcov
    printf 'inferred_sex\\tCNchrX\\tCNchrY\\n%s\\tNA\\tNA\\n' ${meta.sex ?: 'unknown'} > ${meta.id}_sex_check.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        goleft: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
