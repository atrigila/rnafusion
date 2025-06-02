

include { FASTP                      } from '../../../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_FOR_FASTP } from '../../../modules/nf-core/fastqc/main'

workflow TRIM_WORKFLOW {

    take:
        reads           // channel [ meta, [ fastq files ] ]
        adapter_fasta   // channel [ path ]
        skip_qc         // boolean

    main:
        ch_versions       = Channel.empty()
        ch_fastp_html     = Channel.empty()
        ch_fastp_json     = Channel.empty()
        ch_fastqc_trimmed = Channel.empty()

        FASTP(reads, adapter_fasta, false, false, false)
        ch_versions = ch_versions.mix(FASTP.out.versions)

        if(!skip_qc) {
            FASTQC_FOR_FASTP(FASTP.out.reads)
            ch_versions       = ch_versions.mix(FASTQC_FOR_FASTP.out.versions)
            ch_fastqc_trimmed = FASTQC_FOR_FASTP.out.zip
        }

        ch_reads_all           = FASTP.out.reads
        ch_reads_fusioncatcher = ch_reads_all
        ch_fastp_html          = FASTP.out.html
        ch_fastp_json          = FASTP.out.json

    emit:
        ch_reads_all            // Channel [ meta, [reads]   ]
        ch_reads_fusioncatcher  // Channel [ meta, [reads]   ]
        ch_fastp_html           // Channel [ meta, path_html ]
        ch_fastp_json           // Channel [ meta, path_json ]
        ch_fastqc_trimmed       // Channel [ meta, path_zip  ]
        versions = ch_versions  // Channel [ versions        ]
    }

