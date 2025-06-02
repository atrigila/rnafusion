include { STARFUSION                                            }   from '../../../modules/local/starfusion/detect/main'

workflow STARFUSION_WORKFLOW {
    take:
        junctions           // channel: [ meta, junctions ]
        ch_starfusion_ref
        starfusion_fusions

    main:
        def ch_versions = Channel.empty()
        def ch_starfusion_fusions = Channel.empty()

        if (starfusion_fusions){
            fusions = file(starfusion_fusions, checkIfExists:true)
            ch_starfusion_fusions = junctions.map { meta, _junc -> [ meta, fusions ] }
        } else {
            STARFUSION(
                junctions.map { meta, junc -> [ meta, [], junc ]},
                ch_starfusion_ref.map { it -> it[1] }
            )
            ch_versions = ch_versions.mix(STARFUSION.out.versions)
            ch_starfusion_fusions = STARFUSION.out.fusions
        }
    emit:
        fusions               = ch_starfusion_fusions
        versions              = ch_versions
    }

