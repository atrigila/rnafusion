process FUSIONCATCHER_DETECT {
    tag "$meta.id"
    label 'process_high'

    conda "${projectDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d53f36e9e01d14a0ae8e15f8046f52b2883c970c27fe43fdfbd9440a55f5403f/data' :
        'community.wave.seqera.io/library/fusioncatcher:1.33--4733482b637ef92f' }"

    input:
    tuple val(meta), path(fastqs, stageAs: "input/*")
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.fusioncatcher.fusion-genes.txt"), emit: fusions, optional: true
    tuple val(meta), path("*.fusioncatcher.summary.txt")     , emit: summary, optional: true
    tuple val(meta), path("*.fusioncatcher.log")             , emit: log
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def single_end = meta.single_end ? "--single-end" : ""
    """
    fusioncatcher.py \\
        -d ${reference} \\
        -i input \\
        -p ${task.cpus} \\
        -o . \\
        --skip-blat \\
        ${single_end} \\
        ${args}

    mv final-list_candidate-fusion-genes.txt ${prefix}.fusioncatcher.fusion-genes.txt
    mv summary_candidate_fusions.txt ${prefix}.fusioncatcher.summary.txt
    mv fusioncatcher.log ${prefix}.fusioncatcher.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: "\$(fusioncatcher.py --version 2>&1 | awk '{print \$2}')"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.fusioncatcher.fusion-genes.txt
    touch ${prefix}.fusioncatcher.summary.txt
    touch ${prefix}.fusioncatcher.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: "\$(fusioncatcher.py --version 2>&1 | awk '{print \$2}')"
    END_VERSIONS
    """
}
