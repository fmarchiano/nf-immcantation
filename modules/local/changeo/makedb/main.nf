process CHANGEO_MAKEDB {
    tag "$meta.id"
    label 'process_medium'
    label 'changeo'

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(blast)

    output:
    tuple val(meta), path("*db-pass.tsv"), emit: tab
    path "*db-fail.tsv",                   emit: fail    , optional: true
    path "*.log",                          emit: log     , optional: true
    path "versions.yml",                   emit: versions

    script:
    def args = task.ext.args ?: ''
    def vdj_dir = "/usr/local/share/germlines/imgt/${meta.species}/vdj"
    """
    MakeDb.py igblast \\
        -i ${blast} \\
        -s ${fasta} \\
        -r ${vdj_dir} \\
        ${args} \\
        --outname ${meta.id} \\
        --log ${meta.id}-makedb.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        changeo: \$(python3 -c "import changeo; print(changeo.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_db-pass.tsv
    touch ${meta.id}-makedb.log
    echo '"${task.process}":' > versions.yml
    echo '    changeo: 1.3.0' >> versions.yml
    """
}
