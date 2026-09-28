process PRESTO_CLUSTERSETS {
    tag "$meta.id"
    label 'process_long'
    label 'presto'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("*_R1_cluster-pass.fastq"), path("*_R2_cluster-pass.fastq"), emit: reads
    path "*_command_log.txt", emit: logs
    path "versions.yml",     emit: versions

    script:
    def args = task.ext.args ?: '--exec vsearch --ident 0.9'
    """
    ClusterSets.py set \\
        --nproc ${task.cpus} \\
        -s ${r1} \\
        --outname ${meta.id}_R1 \\
        ${args} \\
        --log ${meta.id}_R1.log \\
        > ${meta.id}_command_log.txt

    ClusterSets.py set \\
        --nproc ${task.cpus} \\
        -s ${r2} \\
        --outname ${meta.id}_R2 \\
        ${args} \\
        --log ${meta.id}_R2.log \\
        >> ${meta.id}_command_log.txt

    ParseHeaders.py copy -s ${meta.id}_R1_cluster-pass.fastq -f BARCODE -k CLUSTER --act cat
    ParseHeaders.py copy -s ${meta.id}_R2_cluster-pass.fastq -f BARCODE -k CLUSTER --act cat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$(python3 -c "import presto; print(presto.__version__)")
        vsearch: \$(vsearch --version 2>&1 | head -1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_R1_cluster-pass.fastq
    touch ${meta.id}_R2_cluster-pass.fastq
    touch ${meta.id}_command_log.txt
    echo '"${task.process}":' > versions.yml
    echo '    presto: 0.7.2' >> versions.yml
    """
}
