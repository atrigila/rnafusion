include { ARRIBA_ARRIBA                               } from '../../../modules/nf-core/arriba/arriba/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FOR_ARRIBA } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_FOR_ARRIBA   } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FOR_ARRIBA   } from '../../../modules/nf-core/samtools/view/main'
include { STAR_ALIGN as STAR_FOR_ARRIBA               } from '../../../modules/nf-core/star/align/main'

include { CTATSPLICING_WORKFLOW         }   from '../ctatsplicing_workflow'

workflow ARRIBA_WORKFLOW {
    take:
        reads                           // channel [ meta, [ fastqs ]         ]
        ch_gtf                          // channel [ meta, path_gtf           ]
        ch_fasta                        // channel [ meta, path_fasta         ]
        ch_starindex_ref                // channel [ meta, path_index         ]
        ch_arriba_ref_blacklist         // channel [ meta, path_blacklist     ]
        ch_arriba_ref_cytobands         // channel [ meta, path_cytobands     ]
        ch_arriba_ref_known_fusions     // channel [ meta, path_known_fusions ]
        ch_arriba_ref_protein_domains   // channel [ meta, path_proteins      ]
        ch_starfusion_ref               // channel [ meta, path_starfusion_ref ]
        arriba                          // boolean
        all                             // boolean
        fusioninspector_only            // boolean
        star_ignore_sjdbgtf             // boolean
        ctatsplicing                    // boolean
        seq_center                      // string
        arriba_fusions                  // path
        cram                            // array

    main:

        def ch_versions   = Channel.empty()
        def ch_cram_index = Channel.empty()
        def ch_dummy_file = file("$projectDir/assets/dummy_file_arriba.txt", checkIfExists: true)

        if (( arriba || all ) && !fusioninspector_only) {

            STAR_FOR_ARRIBA(
                reads,
                ch_starindex_ref,
                ch_gtf,
                star_ignore_sjdbgtf,
                '',
                seq_center
            )
            ch_versions = ch_versions.mix(STAR_FOR_ARRIBA.out.versions)

            if ( ctatsplicing || all ) {
                CTATSPLICING_WORKFLOW(
                    STAR_FOR_ARRIBA.out.spl_junc_tab,
                    STAR_FOR_ARRIBA.out.junction,
                    STAR_FOR_ARRIBA.out.bam,
                    ch_starfusion_ref
                )
                ch_versions = ch_versions.mix(CTATSPLICING_WORKFLOW.out.versions)
            }

            if ( arriba_fusions ) {

                ch_arriba_fusions = reads.combine( Channel.value( file( arriba_fusions, checkIfExists: true ) ) )
                    .map { it -> [ it[0], it[2] ] }
                ch_arriba_fusion_fail = ch_dummy_file

            } else {

                ARRIBA_ARRIBA (
                    STAR_FOR_ARRIBA.out.bam,
                    ch_fasta,
                    ch_gtf,
                    ch_arriba_ref_blacklist,
                    ch_arriba_ref_known_fusions,
                    ch_arriba_ref_cytobands,
                    ch_arriba_ref_protein_domains
                )

                ch_versions = ch_versions.mix(ARRIBA_ARRIBA.out.versions)

                ch_arriba_fusions     = ARRIBA_ARRIBA.out.fusions
                ch_arriba_fusion_fail = ARRIBA_ARRIBA.out.fusions_fail.map{ it -> return it[1] }
            }

            if ( cram.contains('arriba') ) {

                SAMTOOLS_SORT_FOR_ARRIBA(STAR_FOR_ARRIBA.out.bam, ch_fasta)
                ch_versions = ch_versions.mix(SAMTOOLS_SORT_FOR_ARRIBA.out.versions )

                SAMTOOLS_VIEW_FOR_ARRIBA(SAMTOOLS_SORT_FOR_ARRIBA.out.bam.map { meta, bam -> [ meta, bam, [] ] }, ch_fasta, [])
                ch_versions = ch_versions.mix(SAMTOOLS_VIEW_FOR_ARRIBA.out.versions )

                SAMTOOLS_INDEX_FOR_ARRIBA(SAMTOOLS_VIEW_FOR_ARRIBA.out.cram)
                ch_versions = ch_versions.mix(SAMTOOLS_INDEX_FOR_ARRIBA.out.versions )

                // Join cram and index files
                ch_cram_index = SAMTOOLS_VIEW_FOR_ARRIBA.out.cram.join(SAMTOOLS_INDEX_FOR_ARRIBA.out.crai)
            }

        } else {

            ch_arriba_fusions = reads
                .combine(Channel.value( file(ch_dummy_file, checkIfExists: true ) ) )
                .map { it -> [ it[0], it[2] ] }

            ch_arriba_fusion_fail = ch_dummy_file
        }

    emit:
        fusions      = ch_arriba_fusions        // channel [ meta, path_fusions   ]
        fusions_fail = ch_arriba_fusion_fail    // channel [ path, fusions_failed ]
        cram_index   = ch_cram_index            // channel [ meta, cram, crai     ]
        versions     = ch_versions              // channel [ versions             ]
    }

