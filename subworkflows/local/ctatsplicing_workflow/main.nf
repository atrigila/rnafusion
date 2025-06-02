include { CTATSPLICING_STARTOCANCERINTRONS } from '../../../modules/local/ctatsplicing/startocancerintrons'

workflow CTATSPLICING_WORKFLOW {
    take:
    splice_junctions // [ val(meta), path(splice_junctions.SJ.out.tab) ]
    junctions       // [ val(meta), path(junctions.Chimeric.out.junction) ]
    aligned_bams    // [ val(meta), path(aligned_bams.Aligned.sortedByCoord.out.bam), path(aligned_bams.Aligned.sortedByCoord.out.bam.bai) ]
    ctat_genome_lib // [ val(meta2), path(path/to/ctat_genome_lib) ]

    main:
    def ch_versions = Channel.empty()

    def ch_ctatsplicing_input = splice_junctions
        .join(junctions, failOnMismatch:true, failOnDuplicate:true)
        .join(aligned_bams, failOnMismatch:true, failOnDuplicate:true)
        .map { meta, split_junction, junction, bam, bai ->
            [ meta, split_junction, junction, bam, bai ]
        }

    CTATSPLICING_STARTOCANCERINTRONS(
        ch_ctatsplicing_input,
        ctat_genome_lib
    )
    ch_versions = ch_versions.mix(CTATSPLICING_STARTOCANCERINTRONS.out.versions.first())

    emit:
    versions              = ch_versions
}
