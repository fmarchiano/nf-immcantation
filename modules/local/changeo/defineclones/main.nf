process CHANGEO_DEFINECLONES {
    tag "$meta.id"
    label 'process_medium'
    label 'immcantation'

    input:
    tuple val(meta), path(tabs)   // one or more TSV files (grouped by subject)

    output:
    tuple val(meta), path("*clone-pass.tsv"), emit: tab
    path "*clone-fail.tsv",                   emit: fail    , optional: true
    path "*.log",                             emit: log     , optional: true
    path "versions.yml",                      emit: versions

    script:
    def args    = task.ext.args ?: ''
    def tab_in  = tabs instanceof List ? tabs.join(' ') : tabs
    def needs_cat = tabs instanceof List && tabs.size() > 1
    """
    ${needs_cat ? "cat ${tab_in} > combined_db.tsv" : ""}

    DefineClones.py \\
        -d ${needs_cat ? 'combined_db.tsv' : tab_in} \\
        ${args} \\
        --nproc ${task.cpus} \\
        --outname ${meta.id} \\
        --log ${meta.id}-defineclones.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        changeo: \$(python3 -c "import changeo; print(changeo.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_clone-pass.tsv
    touch ${meta.id}-defineclones.log
    echo '"${task.process}":' > versions.yml
    echo '    changeo: 1.3.0' >> versions.yml
    """
}
