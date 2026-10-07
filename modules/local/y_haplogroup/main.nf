/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Y_HAPLOGROUP — Y-chromosome haplogroup (paternal line) with Yleaf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Opt-in ('y_haplogroup' in --tools). Runs on male samples only:
    workflows/bam_analysis.nf passes the BAMs whose sex is male, the sex
    INDEXCOV has checked against the reads. Writes Yleaf's prediction table;
    Hg NA means too few Y markers had reads to call a haplogroup.

    Equivalent to: scripts/37-y-haplogroup.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process Y_HAPLOGROUP {
    tag "$meta.id"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/y_haplogroup" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)

    output:
    tuple val(meta), path("${meta.id}_y_haplogroup.txt"), emit: haplogroup
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Yleaf downloads the whole hg38 FASTA on its first run unless its config
    // names one, and the image's config is read-only: its constant is pointed
    // at the reference before it starts, as scripts/37-y-haplogroup.sh does.
    """
    python3 -c 'import sys; from pathlib import Path; from yleaf import yleaf_constants; yleaf_constants.HG38_FULL_GENOME = Path(sys.argv[1]); from yleaf import Yleaf; sys.argv = ["Yleaf"] + sys.argv[2:]; Yleaf.main()' "\$(readlink -f ${reference})" -bam ${bam} -o yleaf -rg hg38 -force -t ${task.cpus}
    [ "\$(grep -c . yleaf/hg_prediction.hg)" -ge 2 ] || { echo "ERROR: Yleaf wrote no prediction" >&2; exit 1; }
    cp yleaf/hg_prediction.hg ${meta.id}_y_haplogroup.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yleaf: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """

    stub:
    """
    printf 'Sample_name\\tHg\\tHg_marker\\tTotal_reads\\tValid_markers\\tQC-score\\tQC-1\\tQC-2\\tQC-3\\n' > ${meta.id}_y_haplogroup.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yleaf: ${task.container.replaceFirst(/^[^:@]+[:@]/, '')}
    END_VERSIONS
    """
}
