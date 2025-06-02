process FUSIONCATCHER_DOWNLOAD {
    tag "fusioncatcher_download"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b54fa9135194c72a18d00db6b399c03248103f87e43ca75e4b50d61179994b3/data':
        'community.wave.seqera.io/library/wget:1.21.4--8b0fcde81c17be5e' }"

    input:
    val genome_gencode_version

    output:
    tuple env(meta), path("human_v*"), emit: reference
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    meta = [ id: "human_v${genome_gencode_version}" ]
    """
    wget --no-check-certificate $args ${params.fusioncatcher_download_link}/human_v${genome_gencode_version}.tar.gz.partaa
    wget --no-check-certificate $args ${params.fusioncatcher_download_link}/human_v${genome_gencode_version}.tar.gz.partab
    wget --no-check-certificate $args ${params.fusioncatcher_download_link}/human_v${genome_gencode_version}.tar.gz.partac
    cat human_v${genome_gencode_version}.tar.gz.* | tar xz
    rm human_v${genome_gencode_version}.tar*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -1 | cut -d ' ' -f 3)
        tar: \$(tar --version | head -1 | sed -e 's/tar (GNU tar) //')
    END_VERSIONS
    """

    stub:
    """
    mkdir human_v${genome_gencode_version}
    touch human_v${genome_gencode_version}/ensembl_fully_overlapping_genes.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -1 | cut -d ' ' -f 3)
        tar: \$(tar --version | head -1 | sed -e 's/tar (GNU tar) //')
    END_VERSIONS
    """
}
