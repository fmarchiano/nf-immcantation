process PRESTO_SPLITSEQ {
    tag "$meta.id"
    label 'process_single'
    label 'immcantation'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_atleast-*.fasta"), emit: fasta
    path "versions.yml",                        emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    SplitSeq.py group \\
        -s ${reads} \\
        ${args} \\
        --fasta \\
        --outname ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$(python3 -c "import presto; print(presto.__version__)")
    END_VERSIONS
    """

    stub:
    def min_count = params.splitseq_min_count ?: 2
    """
    touch ${meta.id}_atleast-${min_count}.fasta
    echo '"${task.process}":' > versions.yml
    echo '    presto: 0.7.2' >> versions.yml
    """
}
