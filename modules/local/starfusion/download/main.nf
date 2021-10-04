// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STARFUSION_DOWNLOAD {
    tag 'star-fusion'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::star-fusion=1.10.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/star-fusion:1.10.0--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/star-fusion:1.10.0--hdfd78af_1"
    }

    input:
    path fasta
    path gtf

    output:
    path "*"  , emit: reference

    script:
    def software = getSoftwareName(task.process)
    """
    export TMPDIR=/tmp

    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.hmm.gz
    wget https://github.com/FusionAnnotator/CTAT_HumanFusionLib/releases/download/v0.3.0/fusion_lib.Mar2021.dat.gz -O CTAT_HumanFusionLib_Mar2021.dat.gz
    wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/AnnotFilterRule.pm -O AnnotFilterRule.pm
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3f
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3i
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3m
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3p
    
    gunzip Pfam-A.hmm.gz && hmmpress Pfam-A.hmm

    prep_genome_lib.pl \\
        --genome_fa $fasta \\
        --gtf $gtf \\
        --annot_filter_rule AnnotFilterRule.pm \\
        --fusion_annot_lib CTAT_HumanFusionLib_Mar2021.dat.gz \\
        --pfam_db Pfam-A.hmm \\
        --dfam_db homo_sapiens_dfam.hmm \\
        --max_readlength $params.read_length \\
        --CPU $task.cpus
    """

}
