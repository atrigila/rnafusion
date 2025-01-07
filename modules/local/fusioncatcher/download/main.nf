process FUSIONCATCHER_DOWNLOAD {
    tag "fusioncatcher_download"
    label 'process_medium'

    container "docker.io/rannickscilifelab/fusioncatcher:1.34"


    input:
    val genome_gencode_version


    output:
    path "human_v${genome_gencode_version}"                , emit: reference
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
        // TODO: move to S3

    // def url =
    """
    wget $args http://sourceforge.net/projects/fusioncatcher/files/data/human_${genome_gencode_version}.tar.gz.aa
    wget $args http://sourceforge.net/projects/fusioncatcher/files/data/human_${genome_gencode_version}.tar.gz.ab
    wget $args http://sourceforge.net/projects/fusioncatcher/files/data/human_${genome_gencode_version}.tar.gz.ac
    wget $args http://sourceforge.net/projects/fusioncatcher/files/data/human_${genome_gencode_version}.tar.gz.ad
    cat human_${genome_gencode_version}.tar.gz.* | tar xz
    rm human_${genome_gencode_version}.tar*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: \$(echo \$(fusioncatcher --version 2>&1))
    END_VERSIONS
    """

    stub:
    """
    mkdir human_v${genome_gencode_version}
    touch human_v${genome_gencode_version}/ensembl_fully_overlapping_genes.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: \$(echo \$(fusioncatcher --version 2>&1))
    END_VERSIONS
    """
}
