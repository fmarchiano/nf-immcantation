process PRESTO_COLLAPSESEQ {
    tag "$meta.id"
    label 'process_medium'
    label 'immcantation'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*collapse-unique.fastq"), emit: reads
    path "*collapse-duplicate.fastq",                emit: dup    , optional: true
    path "*collapse-undetermined.fastq",             emit: und    , optional: true
    path "*.log",                                    emit: log    , optional: true
    path "versions.yml",                             emit: versions

    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    CollapseSeq.py \\
        -s ${reads} \\
        ${args} \\
        --outname ${meta.id}

    ParseLog.py \\
        -l ${meta.id}_collapseseq.log \\
        ${args2} 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$(python3 -c "import presto; print(presto.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_collapse-unique.fastq
    echo '"${task.process}":' > versions.yml
    echo '    presto: 0.7.2' >> versions.yml
    """
}
