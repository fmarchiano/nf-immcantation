process PRESTO_BUILDCONSENSUS {
    tag "$meta.id"
    label 'process_medium'
    label 'presto'

    input:
    tuple val(meta), path(reads)   // [meta, single-end FASTQ annotated with BARCODE]
    val   read_tag                 // suffix for output naming, e.g. 'V' or 'C'

    output:
    tuple val(meta), path("*consensus-pass.fastq"), emit: reads
    path "*.log",                                   emit: log, optional: true
    path "versions.yml",                            emit: versions

    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    BuildConsensus.py \\
        -s ${reads} \\
        ${args} \\
        --nproc ${task.cpus} \\
        --outname ${meta.id}-${read_tag} \\
        --log ${meta.id}-${read_tag}-consensus.log

    ParseLog.py \\
        -l ${meta.id}-${read_tag}-consensus.log \\
        ${args2} 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$(python3 -c "import presto; print(presto.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}-${read_tag}_consensus-pass.fastq
    touch ${meta.id}-${read_tag}-consensus.log
    echo '"${task.process}":' > versions.yml
    echo '    presto: 0.7.2' >> versions.yml
    """
}
