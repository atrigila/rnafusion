process CTATSPLICING_PREPGENOMELIB {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/CTAT-SPLICING/ctat_splicing.v0.0.2.simg' :
        'docker.io/trinityctat/ctat_splicing:0.0.2' }"

    input:
    tuple val(meta), path(genome_lib)
    path(cancer_intron_tsv)

    output:
    tuple val(meta), path(genome_lib, includeInputs:true), emit: reference
    path "versions.yml"                                  , emit: versions

    script:
    def VERSION = '0.0.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    /usr/local/src/CTAT-SPLICING/prep_genome_lib/ctat-splicing-lib-integration.py \\
        --cancer_introns_tsv cancer_introns.*.tsv.gz \\
        --genome_lib_dir $genome_lib

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ctat-splicing: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = '0.0.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch $genome_lib/refGene.bed
    echo | gzip > $genome_lib/refGene.sort.bed.gz
    touch $genome_lib/refGene.sort.bed.gz.tbi
    mkdir $genome_lib/cancer_splicing_lib
    touch $genome_lib/cancer_splicing_lib/cancer_splicing.idx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ctat-splicing: $VERSION
    END_VERSIONS
    """
}
