process SCOPER_SPECTRALCLONES_FAST {
    tag "$meta.id"
    label 'process_medium'
    label 'scoper_fast'

    input:
    tuple val(meta), path(tabs)

    output:
    tuple val(meta), path("*_clone-pass.tsv"), emit: tab
    path "*_scoper_summary.txt", optional: true, emit: summary
    path "*_inter_intra.tsv", optional: true, emit: inter_intra
    path "*_eff_threshold.tsv", optional: true, emit: eff_threshold
    path "*_vjl_groups.tsv", optional: true, emit: vjl_groups
    path "*_clonality_measures.tsv", optional: true, emit: clonality_measures
    path "versions.yml",                       emit: versions

    script:
    def nproc    = task.cpus ?: 1
    def tab_list = tabs instanceof List ? tabs.join(' ') : tabs
    """
    run_scoper_fast.sh ${meta.id} novj none 0 ${nproc} ${tab_list}
    """

    stub:
    """
    touch ${meta.id}_clone-pass.tsv
    printf '"SCOPER_SPECTRALCLONES_FAST":\\n    scoper-fast: 0.1.0\\n' > versions.yml
    """
}
