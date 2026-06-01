process GUNZIP {
    tag "$meta.id"
    label 'process_single'
    container 'ubuntu:22.04'

    input:
    tuple val(meta), path(reads)   // [meta, [R1.fastq.gz, R2.fastq.gz]]

    output:
    tuple val(meta), path("*.fastq"), emit: reads
    path "versions.yml",             emit: versions

    script:
    def r1 = reads[0]
    def r2 = reads[1]
    def r1_out = r1.name.replaceAll(/\.gz$/, '')
    def r2_out = r2.name.replaceAll(/\.gz$/, '')
    """
    gunzip -c ${r1} > ${r1_out}
    gunzip -c ${r2} > ${r2_out}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(gunzip --version 2>&1 | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    def r1_out = reads[0].name.replaceAll(/\.gz$/, '')
    def r2_out = reads[1].name.replaceAll(/\.gz$/, '')
    """
    touch ${r1_out} ${r2_out}
    echo '"${task.process}":' > versions.yml
    echo '    gunzip: 1.0' >> versions.yml
    """
}
