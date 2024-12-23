/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { GENCODE_DOWNLOAD }                from '../../modules/local/gencode_download/main'
include { FUSIONCATCHER_BUILD }             from '../../modules/local/fusioncatcher/build/main'
include { FUSIONREPORT_DOWNLOAD }           from '../../modules/local/fusionreport/download/main'
include { HGNC_DOWNLOAD }                   from '../../modules/local/hgnc/main'
include { STARFUSION_BUILD }                from '../../modules/local/starfusion/build/main'
include { GTF_TO_REFFLAT }                  from '../../modules/local/uscs/custom_gtftogenepred/main'
include { GET_RRNA_TRANSCRIPTS }            from '../../modules/local/get_rrna_transcript/main'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/
include { ARRIBA_DOWNLOAD }                 from '../../modules/nf-core/arriba/download/main'
include { SAMTOOLS_FAIDX }                  from '../../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE }             from '../../modules/nf-core/star/genomegenerate/main'
include { GATK4_CREATESEQUENCEDICTIONARY }  from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_BEDTOINTERVALLIST }         from '../../modules/nf-core/gatk4/bedtointervallist/main'
include { SALMON_INDEX }                    from '../../modules/nf-core/salmon/index/main'
include { GFFREAD }                         from '../../modules/nf-core/gffread/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow BUILD_REFERENCES {

    main:
    ch_versions = Channel.empty()

    if (!file(params.fasta).exists() || file(params.fasta).isEmpty() ||
            !file(params.gtf).exists() || file(params.gtf).isEmpty()){
        GENCODE_DOWNLOAD(params.genome_gencode_version, params.genome)
        ch_versions = ch_versions.mix(GENCODE_DOWNLOAD.out.versions)
        ch_fasta = GENCODE_DOWNLOAD.out.fasta.map { that -> [[id:that.Name], that] }
        ch_gtf = GENCODE_DOWNLOAD.out.gtf.map { that -> [[id:that.Name], that] }
    } else {
        ch_fasta = Channel.fromPath(params.fasta).map { that -> [[id:that.Name], that] }
        ch_gtf = Channel.fromPath(params.gtf).map { that -> [[id:that.Name], that] }
    }

    if (!file(params.fai).exists() || file(params.fai).isEmpty()){
        SAMTOOLS_FAIDX(ch_fasta, [[],[]])
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_fai = SAMTOOLS_FAIDX.out.fai
    } else {
        ch_fai = Channel.fromPath(params.fai).map { that -> [[id:that.Name], that] }
    }

    if ((!file(params.hgnc_ref).exists() || file(params.hgnc_ref).isEmpty() ||
            !file(params.hgnc_date).exists() || file(params.hgnc_date).isEmpty()) && !params.skip_vcf){
        HGNC_DOWNLOAD( )
        ch_versions = ch_versions.mix(HGNC_DOWNLOAD.out.versions)
        ch_hgnc_ref = HGNC_DOWNLOAD.out.hgnc_ref
        ch_hgnc_date = HGNC_DOWNLOAD.out.hgnc_date
    } else {
        ch_hgnc_ref = Channel.fromPath(params.hgnc_ref).map { that -> [[id:that.Name], that] }
        ch_hgnc_date = Channel.fromPath(params.hgnc_date).map { that -> [[id:that.Name], that] }
    }

    if (!file(params.rrna_intervals).exists() || file(params.rrna_intervals).isEmpty()){
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

    if (!file(params.refflat).exists() || file(params.refflat).isEmpty()){
        GTF_TO_REFFLAT(ch_gtf)
        ch_versions = ch_versions.mix(GTF_TO_REFFLAT.out.versions)
        ch_refflat = GTF_TO_REFFLAT.out.refflat.map { that -> [[id:that.Name], that] }
    } else {
        ch_refflat = Channel.fromPath(params.refflat).map { that -> [[id:that.Name], that] }
    }

    if (!file(params.salmon_index).exists() || file(params.salmon_index).isEmpty() ||
        !file(params.salmon_index_stub_check).exists() || file(params.salmon_index_stub_check).isEmpty()){ // add condition for qc
        GFFREAD(ch_gtf, ch_fasta.map{ it -> it[1] })
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
        SALMON_INDEX(ch_fasta.map{ it -> it[1] }, GFFREAD.out.gffread_fasta.map{ it -> it[1] })
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
        ch_salmon_index = SALMON_INDEX.out.index
    } else {
        ch_salmon_index = Channel.fromPath({params.salmon_index})
    }

    if ((params.starindex || params.all || params.starfusion || params.arriba) &&
            (!file(params.starindex_ref).exists() || file(params.starindex_ref).isEmpty() ||
            !file(params.starindex_ref_stub_check).exists() || file(params.starindex_ref_stub_check).isEmpty() )) {
        STAR_GENOMEGENERATE(ch_fasta, ch_gtf)
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        ch_starindex_ref = STAR_GENOMEGENERATE.out.index
    } else {
        ch_starindex_ref = Channel.fromPath(params.starindex_ref).map { that -> [[id:that.Name], that] }
    }

    if ((params.arriba || params.all) &&
            (!file(params.arriba_ref_blacklist).exists() || file(params.arriba_ref_blacklist).isEmpty() ||
            !file(params.arriba_ref_known_fusions).exists() || file(params.arriba_ref_known_fusions).isEmpty() ||
            !file(params.arriba_ref_protein_domains).exists() || file(params.arriba_ref_protein_domains).isEmpty())) {
        ARRIBA_DOWNLOAD(params.genome)
        ch_versions = ch_versions.mix(ARRIBA_DOWNLOAD.out.versions)
        ch_arriba_ref_blacklist = ARRIBA_DOWNLOAD.out.blacklist
        ch_arriba_ref_cytobands = ARRIBA_DOWNLOAD.out.cytobands
        ch_arriba_ref_known_fusions = ARRIBA_DOWNLOAD.out.known_fusions
        ch_arriba_ref_protein_domains = ARRIBA_DOWNLOAD.out.protein_domains
    } else {
        ch_arriba_ref_blacklist = Channel.fromPath(params.arriba_ref_blacklist)
        ch_arriba_ref_cytobands = Channel.fromPath(params.arriba_ref_cytobands)
        ch_arriba_ref_known_fusions = Channel.fromPath(params.arriba_ref_known_fusions)
        ch_arriba_ref_protein_domains = Channel.fromPath(params.arriba_ref_protein_domains)
    }


    if ((params.fusioncatcher || params.all) &&
            (!file(params.fusioncatcher_ref).exists() || file(params.fusioncatcher_ref).isEmpty() ||
            !file(params.fusioncatcher_ref_stub_check).exists() || file(params.fusioncatcher_ref_stub_check).isEmpty() )) {
            FUSIONCATCHER_BUILD(params.genome_gencode_version)
            ch_versions = ch_versions.mix(FUSIONCATCHER_BUILD.out.versions)
            ch_fusioncatcher_ref = FUSIONCATCHER_BUILD.out.reference
    }
    else {
        ch_fusioncatcher_ref = Channel.fromPath(params.fusioncatcher_ref)
    }


    if ((params.starfusion || params.all) &&
            (!file(params.starfusion_ref).exists() || file(params.starfusion_ref).isEmpty() ||
            !file(params.starfusion_ref_stub_check).exists() || file(params.starfusion_ref_stub_check).isEmpty() )) {
            STARFUSION_BUILD(ch_fasta, ch_gtf, params.fusion_annot_lib, params.species)
            ch_versions = ch_versions.mix(STARFUSION_BUILD.out.versions)
            ch_starfusion_ref = STARFUSION_BUILD.out.reference
    }
    else {
        ch_starfusion_ref = Channel.fromPath(params.starfusion_ref)
    }


    if ((params.fusionreport || params.all) &&
            (!file(params.fusionreport_ref).exists() || file(params.fusionreport_ref).isEmpty() ||
            !file(params.fusionreport_ref_stub_check).exists() || file(params.fusionreport_ref_stub_check).isEmpty())) {
        if (!params.no_cosmic && (!params.cosmic_username || !params.cosmic_passwd)) { exit 1, 'COSMIC username and/or password missing' }
        FUSIONREPORT_DOWNLOAD()
        ch_versions = ch_versions.mix(FUSIONREPORT_DOWNLOAD.out.versions)
        ch_fusionreport_ref = FUSIONREPORT_DOWNLOAD.out.fusionreport_ref
    } else {
        ch_fusionreport_ref = Channel.fromPath(params.fusionreport_ref).map { that -> [[id:that.Name], that] }
    }

    emit:
    ch_fasta
    ch_gtf
    ch_fai
    ch_hgnc_ref
    ch_hgnc_date
    ch_rrna_interval
    ch_refflat
    ch_salmon_index
    ch_starindex_ref
    ch_arriba_ref_blacklist
    ch_arriba_ref_cytobands
    ch_arriba_ref_known_fusions
    ch_arriba_ref_protein_domains
    ch_fusioncatcher_ref
    ch_starfusion_ref
    ch_fusionreport_ref
    versions        = ch_versions
}

/*
========================================================================================
    THE END
========================================================================================
*/