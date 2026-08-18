process CHANGEO_CREATEGERMLINES {
    tag "$meta.id"
    label 'process_medium'
    label 'changeo'

    input:
    tuple val(meta), path(tab)
    path germlines

    output:
    tuple val(meta), path("*_germ-pass.tsv"), emit: tab
    path "*.log",                             emit: log     , optional: true
    path "versions.yml",                      emit: versions

    script:
    def vdj_dir = "${germlines}/imgt/${meta.species}/vdj"
    def tab_list = tab instanceof List ? tab.join(' ') : tab
    """
    CreateGermlines.py \\
        -d ${tab_list} \\
        -r ${vdj_dir}/imgt_${meta.species}_IGHV.fasta \\
           ${vdj_dir}/imgt_${meta.species}_IGHD.fasta \\
           ${vdj_dir}/imgt_${meta.species}_IGHJ.fasta \\
        -g full dmask vonly regions \\
        --cloned \\
        --format airr \\
        --outname ${meta.id} \\
        --log ${meta.id}-creategermlines.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        changeo: \$(python3 -c "import changeo; print(changeo.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_germ-pass.tsv
    echo '"${task.process}":' > versions.yml
    echo '    changeo: 1.3.0' >> versions.yml
    """
}
