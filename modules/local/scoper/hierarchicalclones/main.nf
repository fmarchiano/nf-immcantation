process SCOPER_HIERARCHICALCLONES {
    tag "$meta.id"
    label 'process_medium'
    label 'scoper'

    input:
    tuple val(meta), path(tabs)   // one or more AIRR TSV files grouped by subject

    output:
    tuple val(meta), path("*_clone-pass.tsv"), emit: tab
    path "*_scoper_summary.txt", optional: true, emit: summary
    path "*_inter_intra.tsv", optional: true, emit: inter_intra
    path "*_eff_threshold.tsv", optional: true, emit: eff_threshold
    path "*_vjl_groups.tsv", optional: true, emit: vjl_groups
    path "*_spectral_density.pdf", optional: true, emit: spectral_plot
    path "*_clone_summary.pdf", optional: true, emit: clone_plot
    path "*_clonality_measures.tsv", optional: true, emit: clonality_measures
    path "versions.yml",                       emit: versions

    script:
    def method    = task.ext.method    ?: 'nt'
    def linkage   = task.ext.linkage   ?: 'complete'
    def threshold = task.ext.threshold ?: params.clonal_threshold
    def tab_list  = tabs instanceof List ? tabs.join(' ') : tabs
    """
    run_scoper.R ${meta.id} ${method} ${linkage} ${threshold} 1 ${tab_list}
    """

    stub:
    """
    touch ${meta.id}_clone-pass.tsv
    printf '"SCOPER_HIERARCHICALCLONES":\\n    scoper: 1.5.0\\n    alakazam: 1.4.3\\n    r-base: 4.4.2\\n' > versions.yml
    """
}
