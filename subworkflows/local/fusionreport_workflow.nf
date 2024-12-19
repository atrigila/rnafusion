include { FUSIONREPORT      }     from '../../modules/local/fusionreport/detect/main'


workflow FUSIONREPORT_WORKFLOW {
    take:
        reads
        fusionreport_ref
        arriba_fusions
        starfusion_fusions
        fusioncatcher_fusions

    main:
        ch_versions = Channel.empty()
        ch_report = Channel.empty()
        ch_csv = Channel.empty()

        if (!params.fusioninspector_only) {
            arriba_fusions.view { it -> "arriba fusions:${it}"}
            starfusion_fusions.view { it -> "starfusion fusions:${it}"}
            fusioncatcher_fusions.view { it -> "fusioncatcher fusions:${it}"}

            reads_fusions = reads
            .join(arriba_fusions, failOnMismatch:true, failOnDuplicate:true)
            .join(starfusion_fusions, failOnMismatch:true, failOnDuplicate:true)
            .join(fusioncatcher_fusions, failOnMismatch:true, failOnDuplicate:true)

            reads_fusions.view()

            FUSIONREPORT(reads_fusions, fusionreport_ref, params.tools_cutoff)
            ch_fusion_list = FUSIONREPORT.out.fusion_list
            ch_fusion_list_filtered = FUSIONREPORT.out.fusion_list_filtered
            ch_versions = ch_versions.mix(FUSIONREPORT.out.versions)
            ch_report = FUSIONREPORT.out.report
            ch_csv = FUSIONREPORT.out.csv
        } else {
            ch_fusion_list = reads.combine(Channel.value(file(params.fusioninspector_fusions, checkIfExists:true)))
                            .map { meta, reads, fusions -> [ meta, fusions ] }

            ch_fusion_list_filtered  = ch_fusion_list
        }

    emit:
        versions                 = ch_versions
        fusion_list              = ch_fusion_list
        fusion_list_filtered     = ch_fusion_list_filtered
        report                   = ch_report.ifEmpty(null)
        csv                      = ch_csv.ifEmpty(null)

}

