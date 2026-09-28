process PRESTO_FILTERSEQ {
    tag "$meta.id"
    label 'process_single'
    label 'presto'

    input:
    tuple val(meta), path(reads)   // [meta, [R1.fastq, R2.fastq]]
    val filterseq_q

    output:
    tuple val(meta), path("*_quality-pass.fastq"), emit: reads
    path "*.tab",                                  emit: log_tab
    path "versions.yml",                           emit: versions

    script:
    def r1 = reads[0]
    def r2 = reads[1]
    """
    FilterSeq.py quality -s ${r1} -q ${filterseq_q} --outname ${meta.id}_R1 --log ${meta.id}_R1.log --nproc ${task.cpus}
    FilterSeq.py quality -s ${r2} -q ${filterseq_q} --outname ${meta.id}_R2 --log ${meta.id}_R2.log --nproc ${task.cpus}
    ParseLog.py -l ${meta.id}_R1.log ${meta.id}_R2.log -f ID QUALITY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$(FilterSeq.py --version 2>&1 | grep -o '[0-9][0-9.]*' | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_R1_quality-pass.fastq ${meta.id}_R2_quality-pass.fastq
    echo -e "ID\\tQUALITY" > logs.tab
    echo '"${task.process}":' > versions.yml
    echo '    presto: 0.7.2' >> versions.yml
    """
}
