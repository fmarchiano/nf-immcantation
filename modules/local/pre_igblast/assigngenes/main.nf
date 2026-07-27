process PRE_IGBLAST_ASSIGNGENES {
    tag "$meta.id"
    label 'process_high'
    label 'pre_igblast'

    input:
    tuple val(meta), path(fasta)
    path igblast
    path germlines

    output:
    tuple val(meta), path("*.fmt7"),           emit: blast
    tuple val(meta), path(fasta),              emit: fasta
    path "versions.yml",                       emit: versions

    script:
    def args = task.ext.args ?: ''
    def germline_v = "${germlines}/imgt/${meta.species}/vdj/imgt_${meta.species}_IGHV.fasta"
    """
    python -m pre_igblast \\
        -s ${fasta} \\
        -b ${igblast} \\
        --germline-v ${germline_v} \\
        --organism ${meta.species} \\
        --nproc ${task.cpus} \\
        ${args} \\
        --outdir . \\
        --outname ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pre_igblast: \$(pip show pre-igblast 2>/dev/null | awk '/^Version:/{print \$2}')
        igblast: \$(igblastn -version 2>&1 | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_igblast.fmt7
    echo '"${task.process}":' > versions.yml
    echo '    pre_igblast: 0.1.0' >> versions.yml
    echo '    igblast: 1.22.0' >> versions.yml
    """
}
