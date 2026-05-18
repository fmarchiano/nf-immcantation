process PEAR {
    tag "$meta.id"
    label 'process_high'
    label 'pear'

    input:
    tuple val(meta), path(r1_pair), path(r2_pair)   // [meta, R1_pair-pass (V-read), R2_pair-pass (C-read)]

    output:
    tuple val(meta), path("*.assembled.fastq"),      emit: reads
    path "*.unassembled.forward.fastq",              emit: unassembled_fwd, optional: true
    path "*.unassembled.reverse.fastq",              emit: unassembled_rev, optional: true
    path "*.discarded.fastq",                        emit: discarded,       optional: true
    path "versions.yml",                             emit: versions

    script:
    def args = task.ext.args ?: ''
    // -f = V-read (R1_pair, left/5' end of amplicon), -r = C-read (R2_pair, right/3' end — PEAR RCs internally)
    """
    pear \\
        -f ${r1_pair} \\
        -r ${r2_pair} \\
        -o ${meta.id} \\
        -j ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pear: \$(pear 2>&1 | grep -E '^PEAR v' | awk '{print \$2}' | tr -d 'v')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.assembled.fastq
    touch ${meta.id}.unassembled.forward.fastq
    touch ${meta.id}.unassembled.reverse.fastq
    touch ${meta.id}.discarded.fastq
    echo '"${task.process}":' > versions.yml
    echo '    pear: 0.9.2' >> versions.yml
    """
}
