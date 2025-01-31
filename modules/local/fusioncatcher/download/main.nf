process FUSIONCATCHER_DOWNLOAD {
    tag "fusioncatcher_download"
    label 'process_medium'

    conda "${projectDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d53f36e9e01d14a0ae8e15f8046f52b2883c970c27fe43fdfbd9440a55f5403f/data' :
        'community.wave.seqera.io/library/fusioncatcher:1.33--4733482b637ef92f' }"

    input:
    val genome_gencode_version

    output:
    tuple env(meta), path("*"), emit: reference
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    meta = [ id: "human_v${genome_gencode_version}" ]
    """
    wget $args http://sourceforge.net/projects/fusioncatcher/files/data/human_${genome_gencode_version}.tar.gz.aa
    wget $args http://sourceforge.net/projects/fusioncatcher/files/data/human_${genome_gencode_version}.tar.gz.ab
    wget $args http://sourceforge.net/projects/fusioncatcher/files/data/human_${genome_gencode_version}.tar.gz.ac
    wget $args http://sourceforge.net/projects/fusioncatcher/files/data/human_${genome_gencode_version}.tar.gz.ad
    cat human_${genome_gencode_version}.tar.gz.* | tar xz
    rm human_${genome_gencode_version}.tar*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: "\$(fusioncatcher --version 2>&1 | awk '{print \$2}')"
    END_VERSIONS
    """

    stub:
    """
    mkdir human_v${genome_gencode_version}
    touch human_v${genome_gencode_version}/ensembl_fully_overlapping_genes.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: "\$(fusioncatcher --version 2>&1 | awk '{print \$2}')"
    END_VERSIONS
    """
}
