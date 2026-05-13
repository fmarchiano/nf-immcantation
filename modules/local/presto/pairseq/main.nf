process PRESTO_PAIRSEQ {
    tag "$meta.id"
    label 'process_medium'
    label 'immcantation'

    input:
    tuple val(meta), path(r1), path(r2)   // R1 (V-read), R2_primers-pass (C-read)

    output:
    tuple val(meta), path("*R1_pair-pass.fastq"), path("*R2_pair-pass.fastq"), emit: reads
    path "*.log",                                                               emit: log, optional: true
    path "versions.yml",                                                        emit: versions

    script:
    def args = task.ext.args ?: ''
    // r1 may be gzipped (from fastp); pRESTO accepts both .gz and plain
    """
    PairSeq.py \\
        -1 ${r1} \\
        -2 ${r2} \\
        ${args} \\
        --outname ${meta.id}

    # Rename outputs to predictable names
    mv ${meta.id}-1_pair-pass.fastq ${meta.id}-R1_pair-pass.fastq || \\
        rename 's/_1_pair-pass/_R1_pair-pass/' *_pair-pass.fastq 2>/dev/null || true
    mv ${meta.id}-2_pair-pass.fastq ${meta.id}-R2_pair-pass.fastq || \\
        rename 's/_2_pair-pass/_R2_pair-pass/' *_pair-pass.fastq 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$(python3 -c "import presto; print(presto.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}-R1_pair-pass.fastq
    touch ${meta.id}-R2_pair-pass.fastq
    touch ${meta.id}-pairseq.log
    echo '"${task.process}":' > versions.yml
    echo '    presto: 0.7.2' >> versions.yml
    """
}
