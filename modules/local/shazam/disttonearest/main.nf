process SHAZAM_DISTTONEAREST {
    tag "$meta.id"
    label 'process_medium'
    label 'scoper'

    input:
    tuple val(meta), path(tabs)   // one or more AIRR TSV files grouped by subject

    output:
    tuple val(meta), path("*_threshold.txt"),     emit: threshold
    tuple val(meta), path("*_distToNearest.pdf"), emit: plot
    path "versions.yml",                          emit: versions

    script:
    def fixed_thr = task.ext.threshold ?: params.clonal_threshold
    def drift_tol = task.ext.drift_tol ?: params.threshold_drift_tol
    def tab_list  = tabs instanceof List ? tabs.join(' ') : tabs
    """
    find_threshold.R \\
        ${meta.id} \\
        ${fixed_thr} \\
        ${drift_tol} \\
        ${task.cpus} \\
        ${tab_list}
    """

    stub:
    """
    printf 'patient\\t${meta.id}\\nstatus\\tOK\\n' > ${meta.id}_threshold.txt
    touch ${meta.id}_distToNearest.pdf
    printf '"SHAZAM_DISTTONEAREST":\\n    shazam: 1.3.2\\n    alakazam: 1.4.3\\n    r-base: 4.4.2\\n' > versions.yml
    """
}
