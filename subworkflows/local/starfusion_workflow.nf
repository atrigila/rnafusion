include { STAR_ALIGN as STAR_FOR_STARFUSION }    from '../../modules/nf-core/modules/star/align/main'
include { STARFUSION }                           from '../../modules/local/starfusion/detect/main'
include { GET_PATH }                             from '../../modules/local/getpath/main'


workflow STARFUSION_WORKFLOW {
    take:
        reads
        ch_chrgtf
        ch_starindex_ref

    main:
        ch_versions = Channel.empty()
        ch_align = Channel.empty()
        ch_dummy_file = file("$baseDir/assets/dummy_file_starfusion.txt", checkIfExists: true)

        if (params.starfusion || params.all){
            if (params.starfusion_fusions){
                ch_starfusion_fusions = params.starfusion_fusions
            } else {
                STAR_FOR_STARFUSION( reads, ch_starindex_ref, ch_chrgtf, params.star_ignore_sjdbgtf, params.seq_platform, params.seq_center )
                ch_versions = ch_versions.mix(STAR_FOR_STARFUSION.out.versions)
                ch_align = STAR_FOR_STARFUSION.out.bam_sorted

                reads_junction = reads.join(STAR_FOR_STARFUSION.out.junction )

                STARFUSION( reads_junction, params.starfusion_ref)
                ch_versions = ch_versions.mix(STARFUSION.out.versions)

                GET_PATH(STARFUSION.out.fusions)
                ch_starfusion_fusions = GET_PATH.out.file
            }
        }
        else {
            ch_starfusion_fusions = ch_dummy_file
        }
    emit:
        fusions         = ch_starfusion_fusions
        bam_sorted      = ch_align
        versions        = ch_versions.ifEmpty(null)


}

