process CHANGEO_ASSIGNGENES {
    tag "$meta.id"
    label 'process_high'
    label 'immcantation'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fmt7"),           emit: blast
    tuple val(meta), path(fasta),              emit: fasta
    path "versions.yml",                       emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    AssignGenes.py igblast \\
        -s ${fasta} \\
        -b /usr/local/share/igblast \\
        --organism ${meta.species} \\
        ${args} \\
        --nproc ${task.cpus} \\
        --outdir . \\
        --outname ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        changeo: \$(python3 -c "import changeo; print(changeo.__version__)")
        igblast: \$(igblastn -version 2>&1 | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_igblast.fmt7
    echo '"${task.process}":' > versions.yml
    echo '    changeo: 1.3.0' >> versions.yml
    echo '    igblast: 1.22.0' >> versions.yml
    """
}
