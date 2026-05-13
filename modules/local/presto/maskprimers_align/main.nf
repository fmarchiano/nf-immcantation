process PRESTO_MASKPRIMERS_ALIGN {
    tag "$meta.id"
    label 'process_medium'
    label 'immcantation'

    input:
    tuple val(meta), path(reads)   // [meta, R2.fastq] — C-read only
    path  cprimers                 // C-region primers FASTA

    output:
    tuple val(meta), path("*primers-pass.fastq"), emit: reads
    path "*primers-fail.fastq",                   emit: fail    , optional: true
    path "*.log",                                 emit: log     , optional: true
    path "versions.yml",                          emit: versions

    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    MaskPrimers.py align \\
        -s ${reads} \\
        -p ${cprimers} \\
        ${args} \\
        --nproc ${task.cpus} \\
        --outname ${meta.id}-R2 \\
        --log ${meta.id}-maskprimers.log

    ParseLog.py \\
        -l ${meta.id}-maskprimers.log \\
        ${args2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$(python3 -c "import presto; print(presto.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}-R2_primers-pass.fastq
    touch ${meta.id}-maskprimers.log
    echo '"${task.process}":' > versions.yml
    echo '    presto: 0.7.2' >> versions.yml
    """
}
