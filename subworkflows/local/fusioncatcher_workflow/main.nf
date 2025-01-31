include { FUSIONCATCHER_DETECT } from '../../../modules/local/fusioncatcher/detect/main'

// TODO: Remove fusioncatcher_fusions as parameter.
// TODO: remove dummy file. Work with Channel.empty()
// TODO: if the files were already produced and the user want to skip the module because of this, they should be taken them from the sample sheet
// TODO: harmonize `run_fusioncatcher` and `fusioncatcher_only` parameters at main workflow level to activate/skip this one.

workflow FUSIONCATCHER_WORKFLOW {
    take:
        reads                   // channel [ meta, [ fastqs ] ]
        fusioncatcher_ref       // channel [ meta, path       ]
        run_fusioncatcher       // boolean
        all                     // boolean
        fusioninspector_only    // boolean
        fusioncatcher_fusions   // path, string

    main:
        ch_versions   = Channel.empty()
        ch_dummy_file = file("$baseDir/assets/dummy_file_fusioncatcher.txt", checkIfExists: true)

        if (( run_fusioncatcher || all) && !fusioninspector_only ) {
            if (fusioncatcher_fusions){

                ch_fusioncatcher_fusions = reads.combine(Channel.value(file(fusioncatcher_fusions, checkIfExists:true)))
                                            .map { meta, reads, fusions -> [ meta, fusions ] }
            } else {

                FUSIONCATCHER_DETECT (
                    reads,
                    fusioncatcher_ref
                )
                ch_fusioncatcher_fusions = FUSIONCATCHER_DETECT.out.fusions
                ch_versions              = ch_versions.mix(FUSIONCATCHER_DETECT.out.versions)
            }
        }
        else {
            ch_fusioncatcher_fusions = reads.combine(Channel.value(file(ch_dummy_file, checkIfExists:true)))
                                        .map { meta, reads, fusions -> [ meta, fusions ] }
        }

    emit:
        fusions  = ch_fusioncatcher_fusions     // channel [ meta, fusions ]
        versions = ch_versions                  // channel [ versions      ]
    }

