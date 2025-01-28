process GET_RRNA_TRANSCRIPTS {
    tag 'get_rrna_bed'
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pirate:1.0.5--hdfd78af_0' :
        'biocontainers/pirate:1.0.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path('rrna.gtf')    , emit: rrnagtf
    tuple val(meta), path('rrna.bed')    , emit: bed
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    $baseDir/bin/get_rrna_transcripts.py ${gtf} rrna.gtf

    $baseDir/bin/gtf2bed rrna.gtf > rrna.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_rrna_transcripts: v1.0
    END_VERSIONS
    """

    stub:
    """
    touch rrna.gtf
    touch rrna.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_rrna_transcripts: v1.0
    END_VERSIONS
    """
}
