process PRESTO_FILTERSEQ {
    tag "$meta.id"
    label 'process_medium'
    label 'presto'

    input:
    tuple val(meta), path(reads)
    val filterseq_q

    output:
    tuple val(meta), path("*_quality-pass.fastq"), emit: reads
    path "*.tab",                                  emit: log_tab
    path "versions.yml",                           emit: versions

    script:
    """
    FilterSeq.py quality -s ${reads} -q ${filterseq_q} --outname ${meta.id} --log ${meta.id}.log --nproc ${task.cpus}
    ParseLog.py -l ${meta.id}.log -f ID QUALITY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$(FilterSeq.py --version 2>&1 | grep -o '[0-9][0-9.]*' | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_quality-pass.fastq
    echo -e "ID\\tQUALITY" > logs.tab
    echo '"${task.process}":' > versions.yml
    echo '    presto: 0.7.2' >> versions.yml
    """
}
