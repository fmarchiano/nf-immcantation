process ASSEMBLEPAIRS {
    tag "$meta.id"
    label 'process_high'
    label 'presto'

    input:
    tuple val(meta), path(r1_pair), path(r2_pair)   // [meta, R1_pair-pass, R2_pair-pass]

    output:
    tuple val(meta), path("*assemble-pass.fastq"), emit: reads
    path "*assemble-fail.fastq",                   emit: fail    , optional: true
    path "*.log",                                  emit: log     , optional: true
    path "versions.yml",                           emit: versions

    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    def vref  = "/usr/local/share/igblast/fasta/imgt_${meta.species}_ig_v.fasta"
    """
    AssemblePairs.py sequential \\
        -1 ${r2_pair} \\
        -2 ${r1_pair} \\
        -r ${vref} \\
        ${args} \\
        --nproc ${task.cpus} \\
        --outname ${meta.id} \\
        --log ${meta.id}-assemblepairs.log

    ParseLog.py \\
        -l ${meta.id}-assemblepairs.log \\
        ${args2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$(python3 -c "import presto; print(presto.__version__)")
        blastn: \$(blastn -version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_assemble-pass.fastq
    touch ${meta.id}-assemblepairs.log
    echo '"${task.process}":' > versions.yml
    echo '    presto: 0.7.2' >> versions.yml
    echo '    blastn: 2.13.0+' >> versions.yml
    """
}
