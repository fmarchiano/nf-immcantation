process PRESTO_MASKPRIMERS_ALIGN {
    tag "$meta.id"
    label 'process_medium'
    label 'presto'

    input:
    tuple val(meta), path(reads)   // [meta, single-end FASTQ to mask]
    path  primers                  // primers FASTA (C-region for R1, V-region for R2)
    val   read_tag                 // suffix for output naming, e.g. 'C' or 'V'

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
        -p ${primers} \\
        ${args} \\
        --nproc ${task.cpus} \\
        --outname ${meta.id}-${read_tag} \\
        --log ${meta.id}-${read_tag}-maskprimers.log

    ParseLog.py \\
        -l ${meta.id}-${read_tag}-maskprimers.log \\
        ${args2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$(python3 -c "import presto; print(presto.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}-${read_tag}_primers-pass.fastq
    touch ${meta.id}-${read_tag}-maskprimers.log
    echo '"${task.process}":' > versions.yml
    echo '    presto: 0.7.2' >> versions.yml
    """
}
