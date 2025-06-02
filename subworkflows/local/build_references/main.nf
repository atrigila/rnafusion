/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { GENCODE_DOWNLOAD }                from '../../../modules/local/gencode_download/main'
include { FUSIONCATCHER_DOWNLOAD }          from '../../../modules/local/fusioncatcher/download/main'
include { FUSIONREPORT_DOWNLOAD }           from '../../../modules/local/fusionreport/download/main'
include { HGNC_DOWNLOAD }                   from '../../../modules/local/hgnc/main'
include { STARFUSION_BUILD }                from '../../../modules/local/starfusion/build/main'
include { GTF_TO_REFFLAT }                  from '../../../modules/local/uscs/custom_gtftogenepred/main'
include { GET_RRNA_TRANSCRIPTS }            from '../../../modules/local/get_rrna_transcript/main'
include { CTATSPLICING_PREPGENOMELIB }      from '../../../modules/local/ctatsplicing/prepgenomelib/main.nf'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

include { ARRIBA_DOWNLOAD }                 from '../../../modules/nf-core/arriba/download/main'
include { SAMTOOLS_FAIDX }                  from '../../../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE }             from '../../../modules/nf-core/star/genomegenerate/main'
include { GATK4_CREATESEQUENCEDICTIONARY }  from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_BEDTOINTERVALLIST }         from '../../../modules/nf-core/gatk4/bedtointervallist/main'
include { SALMON_INDEX }                    from '../../../modules/nf-core/salmon/index/main'
include { GFFREAD }                         from '../../../modules/nf-core/gffread/main'
include { getFileSuffix } from '../../../modules/nf-core/cat/cat/main.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow BUILD_REFERENCES {

    take:
    tools // list of all the tools to create references for

    main:
    def ch_versions = Channel.empty()

    def ch_fasta = Channel.empty()
    def ch_gtf   = Channel.empty()
    if (!exists_not_empty(params.fasta) || !exists_not_empty(params.gtf)){
        GENCODE_DOWNLOAD(params.genome_gencode_version, params.genome)
        ch_versions = ch_versions.mix(GENCODE_DOWNLOAD.out.versions)
        ch_fasta = GENCODE_DOWNLOAD.out.fasta.map { that -> [[id:that.Name], that] }
        ch_gtf = GENCODE_DOWNLOAD.out.gtf.map { that -> [[id:that.Name], that] }
    } else {
        ch_fasta = Channel.fromPath(params.fasta).map { that -> [[id:that.Name], that] }
        ch_gtf = Channel.fromPath(params.gtf).map { that -> [[id:that.Name], that] }
    }

    def ch_fai = Channel.empty()
    if (!exists_not_empty(params.fai)){
        SAMTOOLS_FAIDX(ch_fasta, [[],[]], false)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_fai = SAMTOOLS_FAIDX.out.fai
    } else {
        ch_fai = Channel.fromPath(params.fai).map { that -> [[id:that.Name], that] }
    }

    def ch_hgnc_date = Channel.empty()
    def ch_hgnc_ref  = Channel.empty()
    //TODO: unify as if(tools.contains("fusioninspector")) once nextflow bug fixed
    def run_fusioninspector = tools.contains("fusioninspector")
    if(run_fusioninspector && !params.skip_vcf) {
        if ((!exists_not_empty(params.hgnc_ref) || !exists_not_empty(params.hgnc_date)) && !params.skip_vcf){
            HGNC_DOWNLOAD( )
            ch_versions = ch_versions.mix(HGNC_DOWNLOAD.out.versions)
            ch_hgnc_ref = HGNC_DOWNLOAD.out.hgnc_ref.map { that -> [[id:that.Name], that] }
            ch_hgnc_date = HGNC_DOWNLOAD.out.hgnc_date.map { that -> [[id:that.Name], that] }
        } else {
            ch_hgnc_ref = Channel.fromPath(params.hgnc_ref).map { that -> [[id:that.Name], that] }
            ch_hgnc_date = Channel.fromPath(params.hgnc_date).map { that -> [[id:that.Name], that] }
        }
    }

    def ch_rrna_interval = Channel.empty()
    if (!params.skip_qc) {
        if (!exists_not_empty(params.rrna_intervals)){
            GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
            ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
            GET_RRNA_TRANSCRIPTS(ch_gtf)
            ch_versions = ch_versions.mix(GET_RRNA_TRANSCRIPTS.out.versions)
            GATK4_BEDTOINTERVALLIST(GET_RRNA_TRANSCRIPTS.out.bed, GATK4_CREATESEQUENCEDICTIONARY.out.dict )
            ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions)
            ch_rrna_interval = GATK4_BEDTOINTERVALLIST.out.interval_list
        } else {
            ch_rrna_interval = Channel.fromPath(params.rrna_intervals).map { that -> [[id:that.Name], that] }
        }
    }

    def ch_refflat = Channel.empty()
    if (!params.skip_qc) {
        if (!exists_not_empty(params.refflat)){
            GTF_TO_REFFLAT(ch_gtf)
            ch_versions = ch_versions.mix(GTF_TO_REFFLAT.out.versions)
            ch_refflat = GTF_TO_REFFLAT.out.refflat.map { that -> [[id:that.Name], that] }
        } else {
            ch_refflat = Channel.fromPath(params.refflat).map { that -> [[id:that.Name], that] }
        }
    }

    def ch_salmon_index = Channel.empty()
    if (tools.contains("salmon")) {
        if (!params.skip_qc) {
            if (!exists_not_empty(params.salmon_index)){
                GFFREAD(ch_gtf, ch_fasta.map{ it -> it[1] })
                ch_versions = ch_versions.mix(GFFREAD.out.versions)

                SALMON_INDEX(ch_fasta.map{ it -> it[1] }, GFFREAD.out.gffread_fasta.map{ it -> it[1] })
                ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
                ch_salmon_index = SALMON_INDEX.out.index
            } else {
                ch_salmon_index = Channel.fromPath(params.salmon_index)
            }
        }
    }

    def ch_starindex_ref = Channel.empty()
    def star_index_tools = tools.intersect(["starfusion", "arriba", "ctatsplicing", "stringtie"])
    if (star_index_tools) {
        if (!exists_not_empty(params.starindex_ref)) {
            STAR_GENOMEGENERATE(ch_fasta, ch_gtf)
            ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
            ch_starindex_ref = STAR_GENOMEGENERATE.out.index
        } else {
            ch_starindex_ref = Channel.fromPath(params.starindex_ref).map { that -> [[id:that.Name], that] }
        }
    }

    def ch_arriba_ref_blacklist       = Channel.empty()
    def ch_arriba_ref_cytobands       = Channel.empty()
    def ch_arriba_ref_known_fusions   = Channel.empty()
    def ch_arriba_ref_protein_domains = Channel.empty()
    if (tools.contains("arriba")) {
        if (!exists_not_empty(params.arriba_ref_blacklist) || !exists_not_empty(params.arriba_ref_known_fusions) || !exists_not_empty(params.arriba_ref_protein_domains) || !exists_not_empty(params.arriba_ref_cytobands)) {
            ARRIBA_DOWNLOAD(params.genome)
            ch_versions = ch_versions.mix(ARRIBA_DOWNLOAD.out.versions)
            ch_arriba_ref_blacklist       = ARRIBA_DOWNLOAD.out.blacklist
            ch_arriba_ref_cytobands       = ARRIBA_DOWNLOAD.out.cytobands
            ch_arriba_ref_known_fusions   = ARRIBA_DOWNLOAD.out.known_fusions
            ch_arriba_ref_protein_domains = ARRIBA_DOWNLOAD.out.protein_domains
        } else {
            ch_arriba_ref_blacklist       = Channel.fromPath(params.arriba_ref_blacklist)
            ch_arriba_ref_cytobands       = Channel.fromPath(params.arriba_ref_cytobands)
            ch_arriba_ref_known_fusions   = Channel.fromPath(params.arriba_ref_known_fusions)
            ch_arriba_ref_protein_domains = Channel.fromPath(params.arriba_ref_protein_domains)
        }
    }

    def ch_fusioncatcher_ref = Channel.empty()
    if (tools.contains("fusioncatcher")) {
        if (!exists_not_empty(params.fusioncatcher_ref)) {
            FUSIONCATCHER_DOWNLOAD(params.genome_gencode_version)
            ch_versions = ch_versions.mix(FUSIONCATCHER_DOWNLOAD.out.versions)
            ch_fusioncatcher_ref = FUSIONCATCHER_DOWNLOAD.out.reference
        }
        else {
            ch_fusioncatcher_ref = Channel.fromPath(params.fusioncatcher_ref).map { it -> [[id:it.name], it] }
        }
    }

    def ch_starfusion_ref = Channel.empty()
    if (tools.intersect(["starfusion", "ctatsplicing", "fusioninspector"])) {
        if (!exists_not_empty(params.starfusion_ref)) {
            if(!params.fusion_annot_lib) {
                error("Expected --fusion_annot_lib to be specified when using StarFusion or any tools that depend on it")
            }
            STARFUSION_BUILD(ch_fasta, ch_gtf, params.fusion_annot_lib, params.species, params.dfam_version, params.pfam_version)
            ch_versions = ch_versions.mix(STARFUSION_BUILD.out.versions)
            if (tools.contains("ctatsplicing")) {
                CTATSPLICING_PREPGENOMELIB(
                    STARFUSION_BUILD.out.reference,
                    params.ctatsplicing_cancer_introns
                )
                ch_versions = ch_versions.mix(CTATSPLICING_PREPGENOMELIB.out.versions)
                ch_starfusion_ref = CTATSPLICING_PREPGENOMELIB.out.reference
            } else {
                ch_starfusion_ref = STARFUSION_BUILD.out.reference
            }
        }
        else {
            ch_starfusion_ref = Channel.fromPath(params.starfusion_ref).map { it -> [[id:it.name], it] }
        }
    }

    def ch_fusionreport_ref = Channel.empty()
    if (tools.contains("fusionreport")) {
        if (!exists_not_empty(params.fusionreport_ref)) {
            if (!params.no_cosmic && (!params.cosmic_username || !params.cosmic_passwd)) {
                error('COSMIC username and/or password missing, this is needed to download the fusionreport reference')
            }
            FUSIONREPORT_DOWNLOAD()
            ch_versions = ch_versions.mix(FUSIONREPORT_DOWNLOAD.out.versions)
            ch_fusionreport_ref = FUSIONREPORT_DOWNLOAD.out.fusionreport_ref
        } else {
            ch_fusionreport_ref = Channel.fromPath(params.fusionreport_ref).map { that -> [[id:that.Name], that] }
        }
    }

    emit:
    fasta                       = ch_fasta.collect()
    gtf                         = ch_gtf.collect()
    fai                         = ch_fai.collect()
    hgnc_ref                    = ch_hgnc_ref.collect()
    hgnc_date                   = ch_hgnc_date.collect()
    rrna_interval               = ch_rrna_interval.collect()
    refflat                     = ch_refflat.collect()
    salmon_index                = ch_salmon_index.collect()
    starindex_ref               = ch_starindex_ref.collect()
    arriba_ref_blacklist        = ch_arriba_ref_blacklist.collect()
    arriba_ref_cytobands        = ch_arriba_ref_cytobands.collect()
    arriba_ref_known_fusions    = ch_arriba_ref_known_fusions.collect()
    arriba_ref_protein_domains  = ch_arriba_ref_protein_domains.collect()
    fusioncatcher_ref           = ch_fusioncatcher_ref.collect()
    starfusion_ref              = ch_starfusion_ref.collect()
    fusionreport_ref            = ch_fusionreport_ref.collect()
    versions                    = ch_versions
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

def exists_not_empty(path) {
    def path_to_check = file(path as String)
    // Return false if the path does not exist
    if(!path_to_check.exists()) {
        return false
    }

    // Don't check directories if the path is not local
    def is_local = path_to_check.getScheme() == "file"
    if(!is_local || !path_to_check.toFile().isDirectory()) {
        return !path_to_check.isEmpty()
    }

    // Get the first file in a directory and return whether it is empty or not
    def first_file = null
    path_to_check.toFile().eachFileRecurse(groovy.io.FileType.FILES) { file ->
        first_file = file
        return
    }
    return !first_file.toPath().isEmpty()
}
