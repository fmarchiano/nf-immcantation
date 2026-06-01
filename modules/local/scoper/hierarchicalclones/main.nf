process SCOPER_HIERARCHICALCLONES {
    tag "$meta.id"
    label 'process_medium'
    label 'scoper'

    input:
    tuple val(meta), path(tabs)   // one or more AIRR TSV files grouped by subject

    output:
    tuple val(meta), path("*_clone-pass.tsv"), emit: tab
    path "versions.yml",                       emit: versions

    script:
    def method    = task.ext.method    ?: 'nt'
    def linkage   = task.ext.linkage   ?: 'complete'
    def threshold = task.ext.threshold ?: params.clonal_threshold
    def tab_list  = tabs instanceof List ? tabs.join(' ') : tabs
    """
    run_scoper.R \\
        ${meta.id} \\
        ${method} \\
        ${linkage} \\
        ${threshold} \\
        ${task.cpus} \\
        ${tab_list}
    """

    stub:
    """
    touch ${meta.id}_clone-pass.tsv
    printf '"SCOPER_HIERARCHICALCLONES":\\n    scoper: 1.5.0\\n    alakazam: 1.4.3\\n    r-base: 4.4.2\\n' > versions.yml
    """
}
