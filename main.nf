#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnafusion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rnafusion
    Website: https://nf-co.re/rnafusion
    Slack  : https://nfcore.slack.com/channels/rnafusion
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_rnafusion_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_rnafusion_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_rnafusion_pipeline'
include { RNAFUSION               } from './workflows/rnafusion'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        args,
        params.outdir,
    )

    def tools = params.tools.tokenize(",")
    if (tools.contains("all")) {
        def json = new groovy.json.JsonSlurper().parseText(file("${projectDir}/nextflow_schema.json").text)
        def pattern = json.get('$defs')?.get('input_output_options')?.get('properties')?.get('tools')?.get('pattern')
        if (!pattern) {
            error("Could not fetch the allowed tools from the JSON schema, please check the code. If you see this as a pipeline user, please contact the developers instead.")
        }
        tools = pattern.replace('^((', "").replace(')?,?)*(?<!,)$', "").tokenize("|") - "all"
    }

    def profiles = workflow.profile
    if ((profiles.contains("conda") || profiles.contains("mamba")) && (tools.contains("ctatsplicing"))) {
        error("Conda or Mamba runs are not supported when ctatsplicing is in `--tools`")
    }

    //
    // WORKFLOW: Run main workflow
    //
    RNAFUSION(
        PIPELINE_INITIALISATION.out.samplesheet,
        tools
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        RNAFUSION.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
