#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
name = 'gatk-mutect2-variant-calling'
short_name = 'gatk-mutect2'
version = '4.1.8.0-2.0'


/*
========================================================================================
                    ICGC ARGO GATK Mutect2 Variant Calling Workflow
========================================================================================
#### Homepage / Documentation
https://github.com/icgc-argo/gatk-mutect2-variant-calling
#### Authors
Junjun Zhang @junjun-zhang <junjun.zhang@oicr.on.ca>
Linda Xiang @lindaxiang <linda.xiang@oicr.on.ca>
----------------------------------------------------------------------------------------

Required Parameters (no default):
--study_id                              SONG study ID
--tumour_aln_analysis_id                Tumour WGS sequencing_alignment SONG analysis ID
--normal_aln_analysis_id                Normal WGS sequencing_alignment SONG analysis ID
--tumour_aln_metadata                   Tumour WGS alignment metadata JSON file
--tumour_extra_info                     Tumour sample extra info TSV including submitter ID to uniform ID mapping
--tumour_aln_cram                       Tumour WGS aligned CRAM file (index file .crai is also required)
--normal_aln_metadata                   Normal WGS alignment metadata JSON file
--normal_extra_info                     Normal sample extra info TSV including submitter ID to uniform ID mapping
--normal_aln_cram                       Normal WGS aligned CRAM file (index file .crai is also required)
--ref_fa                                Reference genome '.fa' file, secondary files ('.fa.fai', '.dict') are expected to be under the same folder
--song_url                              SONG server URL
--score_url                             SCORE server URL
--api_token                             SONG/SCORE API Token

General Parameters (with defaults):
--cpus                                  cpus given to all process containers (default 1)
--mem                                   memory (GB) given to all process containers (default 1)

Download Parameters (object):
--download
{
    song_container_version              song docker container version, defaults set below
    score_container_version             score docker container version, defaults set below
    song_url                            song url for download process (defaults to main song_url param)
    score_url                           score url for download process (defaults to main score_url param)
    api_token                           song/score API token for download process (defaults to main api_token param)
    song_cpus
    song_mem
    score_cpus
    score_mem
    score_transport_mem
}

bqsr Parameters (object):
--bqsr
{
    cpus                                cpus for seqDataToLaneBam container, defaults to cpus parameter
    mem                                 memory (GB) for seqDataToLaneBam container, defaults to mem parameter
    dbsnp_vcf                           vcf for known germline variants from dbsnp
    known_indels_sites_vcfs             vcf for known germline indel sites
    ref_dict                            reference genome .dict file
    ref_fa                              reference genome .fa file
}

calculateContamination Parameters (object):
--calculateContamination
{
    cpus                                cpus for bwaMemAligner container, defaults to cpus parameter
    mem                                 memory (GB) for jvm_mem used by gatk jar
    ref_dict                            reference genome .dict file
    ref_fa                              reference genome .fa file
    variants_for_contamination          reference vcf file from exac_common
}

mutect2 Parameters (object):
--mutect2
{
    container_version                   docker container version, defaults to unset
    cpus                                cpus for bamMergeSortMarkdup container, defaults to cpus parameter
    mem                                 memory (GB) for bamMergeSortMarkdup container, defaults to mem parameter
}

Upload Parameters (object):
--upload
{
    song_container_version              song docker container version, defaults set below
    score_container_version             score docker container version, defaults set below
    song_url                            song url for upload process (defaults to main song_url param)
    score_url                           score url for upload process (defaults to main score_url param)
    api_token                           song/score API token for upload process (defaults to main api_token param)
    song_cpus                           cpus for song container, defaults to cpus parameter
    song_mem                            memory (GB) for song container, defaults to mem parameter
    score_cpus                          cpus for score container, defaults to cpus parameter
    score_mem                           memory (GB) for score container, defaults to mem parameter
    score_transport_mem                 memory (GB) for score_transport, defaults to mem parameter
    extract_cpus                        cpus for extract container, defaults to cpus parameter
    extract_mem                         memory (GB) extract score container, defaults to mem parameter
}

*/

params.study_id = ""

// aligned seq will be downloaded from SONG/SCORE
params.tumour_aln_analysis_id = ""
params.normal_aln_analysis_id = ""

// the following params if provided local files will be used
params.tumour_aln_metadata = "NO_FILE1"
params.tumour_aln_cram = "NO_FILE2"
params.tumour_extra_info = "NO_FILE3"
params.normal_aln_metadata = "NO_FILE4"
params.normal_aln_cram = "NO_FILE5"
params.normal_extra_info = "NO_FILE6"

// dir for outputs, must be set when running in local mode
params.publish_dir = ""

params.perform_bqsr = true  // default to true

params.ref_fa = "tests/reference/tiny-grch38-chr11-530001-537000.fa"

params.mutect2_scatter_interval_files = "assets/mutect2.scatter_by_chr/chr*.interval_list"
params.bqsr_recal_grouping_file = "assets/bqsr.sequence_grouping.grch38_hla_decoy_ebv.csv"
params.bqsr_apply_grouping_file = "assets/bqsr.sequence_grouping_with_unmapped.grch38_hla_decoy_ebv.csv"

// Allele frequency only, pass-only gnomAD vcf file
params.germline_resource_vcfs = []  // "tests/data/HCC1143-mini-Mutect2-calls/HCC1143.mutect2.copy.vcf.gz"

params.panel_of_normals = "NO_FILE"  // optional PoN VCF file, for now we use GATK's public PoN based on 1000 genome project

params.contamination_variants = ""

params.api_token = ""
params.song_url = ""
params.score_url = ""
params.cleanup = true

params.cpus = 2
params.mem = 4

params.download = [:]
params.bqsr = [:]
params.mutect2 = [:]
params.ref_dict = ""
params.gatherPileupSummaries = [:]
params.learnReadOrientationModel = [:]
params.mergeVcfs = [:]
params.mergeMutectStats = [:]
params.filterMutectCalls = [:]

params.mutect2_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.mutect2 ?: [:])
]

params.gatherPileupSummaries_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'ref_dict': params.ref_dict,
    *:(params.gatherPileupSummaries ?: [:])
]

params.calculateContamination = [
    'variants_for_contamination': 'NO_FILE'
]

params.upload = [:]

download_params = [
    'song_cpus': params.cpus,
    'song_mem': params.mem,
    'score_cpus': params.cpus,
    'score_mem': params.mem,
    'song_url': params.song_url,
    'score_url': params.score_url,
    'api_token': params.api_token,
    *:(params.download ?: [:])
]

mutect2_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'ref_fa': params.ref_fa,
    *:(params.mutect2 ?: [:])
]

calculateContamination_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'ref_dict': params.ref_dict,
    'ref_fa': params.ref_fa,
    *:(params.calculateContamination ?: [:])
]

learnReadOrientationModel_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.learnReadOrientationModel ?: [:])
]

mergeVcfs_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.mergeVcfs ?: [:])
]

mergeMutectStats_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.mergeMutectStats ?: [:])
]

filterMutectCalls_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.filterMutectCalls ?: [:])
]

upload_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'song_url': params.song_url,
    'score_url': params.score_url,
    'api_token': params.api_token,
    *:(params.upload ?: [:])
]


// Include all modules and pass params

include { songScoreDownload as dnldT; songScoreDownload as dnldN } from './song-score-utils/song-score-download' params(download_params)
include { bqsr as bqsrT; bqsr as bqsrN } from './bqsr/bqsr'
include { gatkMutect2 as Mutect2 } from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-mutect2.4.1.8.0-2.2/tools/gatk-mutect2/gatk-mutect2' params(mutect2_params)
include { getSecondaryFiles as getSec } from './wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/helper-functions@1.0.0/main'
include { gatkLearnReadOrientationModel as learnROM } from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-learn-read-orientation-model.4.1.8.0-2.0/tools/gatk-learn-read-orientation-model/gatk-learn-read-orientation-model' params(learnReadOrientationModel_params)
include { gatkMergeVcfs as mergeVcfs } from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-merge-vcfs.4.1.8.0-2.0/tools/gatk-merge-vcfs/gatk-merge-vcfs' params(mergeVcfs_params)
include { gatkMergeMutectStats as mergeMS } from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-merge-mutect-stats.4.1.8.0-2.0/tools/gatk-merge-mutect-stats/gatk-merge-mutect-stats' params(mergeMutectStats_params)
include { calculateContamination as calCont } from './calculate-contamination/calculate-contamination' params(calculateContamination_params)
include { gatkFilterMutectCalls as filterMC } from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-filter-mutect-calls.4.1.8.0-2.2/tools/gatk-filter-mutect-calls/gatk-filter-mutect-calls' params(filterMutectCalls_params)
include { gatkSelectVariants as excIndel; gatkSelectVariants as selIndel } from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-select-variants.4.1.8.0-1.0/tools/gatk-select-variants/gatk-select-variants'
include { payloadGenVariantCalling as pGenVarSnv; payloadGenVariantCalling as pGenVarIndel; payloadGenVariantCalling as pGenQc } from "./modules/raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/payload-gen-variant-calling.0.3.6.0/tools/payload-gen-variant-calling/payload-gen-variant-calling"
include { prepMutect2Qc as prepQc } from './modules/raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/prep-mutect2-qc.0.1.2.0/tools/prep-mutect2-qc/prep-mutect2-qc'
include { songScoreUpload } from './song-score-utils/song-score-upload' params(upload_params)
include { songScoreUpload as upSnv; songScoreUpload as upIndel; songScoreUpload as upQc} from './song-score-utils/song-score-upload' params(upload_params)
include { cleanupWorkdir as cleanupM2; cleanupWorkdir as cleanupBqsr } from './wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/cleanup-workdir@1.0.0/main'
include { payloadAddUniformIds as pAddIdT; payloadAddUniformIds as pAddIdN } from './wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/payload-add-uniform-ids@0.1.1/main'


workflow M2 {
    take:
        study_id
        ref_fa
        ref_fa_2nd
        ref_fa_dict
        ref_fa_img
        germline_resource_vcfs
        germline_resource_indices
        panel_of_normals
        panel_of_normals_idx
        contamination_variants
        contamination_variants_indices
        tumour_aln_analysis_id
        normal_aln_analysis_id
        tumour_aln_metadata
        tumour_extra_info
        tumour_aln_cram
        normal_aln_metadata
        normal_extra_info
        normal_aln_cram
        mutect2_scatter_interval_files
        perform_bqsr
        bqsr_recal_grouping_file
        bqsr_apply_grouping_file

    main:
        local_mode = false

        tumour_aln_seq = Channel.from()
        tumour_aln_seq_idx = Channel.from()
        normal_aln_seq = Channel.from()
        normal_aln_seq_idx = Channel.from()

        Channel
            .fromPath(mutect2_scatter_interval_files, checkIfExists: true)
            .set{ mutect2_scatter_interval_files_ch }

        if (tumour_aln_analysis_id && normal_aln_analysis_id) {
            // download tumour aligned seq and metadata from song/score (analysis type: sequencing_alignment)
            dnldT(study_id, tumour_aln_analysis_id)
            tumour_aln_seq = dnldT.out.files.flatten().first()
            tumour_aln_seq_idx = dnldT.out.files.flatten().last()
            tumour_aln_meta = dnldT.out.song_analysis

            // download normal aligned seq and metadata from song/score (analysis type: sequencing_alignment)
            dnldN(study_id, normal_aln_analysis_id)
            normal_aln_seq = dnldN.out.files.flatten().first()
            normal_aln_seq_idx = dnldN.out.files.flatten().last()
            normal_aln_meta = dnldN.out.song_analysis
        } else if (
            !tumour_aln_metadata.startsWith('NO_FILE') && \
            !tumour_extra_info.startsWith('NO_FILE') && \
            !tumour_aln_cram.startsWith('NO_FILE') && \
            !normal_aln_metadata.startsWith('NO_FILE') && \
            !normal_extra_info.startsWith('NO_FILE') && \
            !normal_aln_cram.startsWith('NO_FILE')
        ) {
            if (!params.publish_dir) {
                exit 1, "When use local inputs, params.publish_dir must be specified."
            } else {
                log.info "Use local inputs, outputs will be in: ${params.publish_dir}"
            }

            local_mode = true

            tumour_aln_seq = file(tumour_aln_cram)
            tumour_aln_seq_idx = Channel.fromPath(getSec(tumour_aln_cram, ['crai', 'bai']))
            pAddIdT(file(tumour_aln_metadata), file(tumour_extra_info))
            tumour_aln_meta = pAddIdT.out.payload

            normal_aln_seq = file(normal_aln_cram)
            normal_aln_seq_idx = Channel.fromPath(getSec(normal_aln_cram, ['crai', 'bai']))
            pAddIdN(file(normal_aln_metadata), file(normal_extra_info))
            normal_aln_meta = pAddIdN.out.payload
        } else {
            exit 1, "To download input aligned seq files from SONG/SCORE, please provide `params.tumour_aln_analysis_id` and `params.normal_aln_analysis_id`.\n" +
                "Or please provide `params.tumour_aln_metadata`, `params.tumour_extra_info`, `params.tumour_aln_cram`, `params.normal_aln_metadata`, `params.normal_extra_info` and `params.normal_aln_cram` to use local files as input."
        }


        if (perform_bqsr) {

            Channel
                .fromPath(bqsr_recal_grouping_file)
                .splitCsv()
                .map{ row -> tuple(row[0].toInteger(), row[1].trim()) }
                .set{ bqsr_recal_grouping_ch }

            Channel
                .fromPath(bqsr_apply_grouping_file)
                .splitCsv()
                .map{ row -> tuple(row[0].toInteger(), row[1].trim()) }
                .set{ bqsr_apply_grouping_ch }

            // BQSR Tumour
            bqsrT(
                tumour_aln_seq,  // aln seq
                tumour_aln_seq_idx.collect(),   // aln idx
                ref_fa,
                ref_fa_2nd,
                germline_resource_vcfs,  // use gnomAD as known_sites
                germline_resource_indices,
                bqsr_recal_grouping_ch,
                bqsr_apply_grouping_ch,
                'tumour.recalibrated_bam'
            )
            tumour_aln_seq = bqsrT.out.bqsr_bam
            tumour_aln_seq_idx = bqsrT.out.bqsr_bam_bai

            // BQSR Normal
            bqsrN(
                normal_aln_seq,  // aln seq
                normal_aln_seq_idx.collect(),   // aln idx
                ref_fa,
                ref_fa_2nd,
                germline_resource_vcfs,  // use gnomAD as known_sites
                germline_resource_indices,
                bqsr_recal_grouping_ch,
                bqsr_apply_grouping_ch,
                'normal.recalibrated_bam'
            )
            normal_aln_seq = bqsrN.out.bqsr_bam
            normal_aln_seq_idx = bqsrN.out.bqsr_bam_bai

        }

        // Mutect2
        Mutect2(
            tumour_aln_seq,
            tumour_aln_seq_idx.collect(),
            normal_aln_seq,
            normal_aln_seq_idx.collect(),
            ref_fa,
            ref_fa_2nd.collect(),
            germline_resource_vcfs.collect(),
            germline_resource_indices.collect(),
            panel_of_normals.collect(),
            panel_of_normals_idx.collect(),
            mutect2_scatter_interval_files_ch.flatten()
        )

        // learnROM
        learnROM(
            Mutect2.out.f1r2_counts
        )

        // mergeVcfs
        mergeVcfs(
            Mutect2.out.output_vcf.collect(),
            tumour_aln_seq.name  // use tumour seq basename for output
        )

        // mergeMS
        mergeMS(
            Mutect2.out.mutect_stats.collect(),
            'merged-mutect-stats'
        )

        // calCont
        calCont(
            tumour_aln_seq,
            tumour_aln_seq_idx.collect(),
            normal_aln_seq,
            normal_aln_seq_idx.collect(),
            ref_fa,
            ref_fa_2nd,
            ref_fa_dict,
            contamination_variants,
            contamination_variants_indices,
            mutect2_scatter_interval_files_ch
        )

        // filterMC
        filterMC(
            mergeVcfs.out.output_vcf,
            mergeVcfs.out.output_tbi,
            ref_fa,
            ref_fa_2nd.collect(),
            calCont.out.tumour_contamination_metrics,
            calCont.out.tumour_segmentation_metrics,
            learnROM.out.artifact_prior_table.collect(),
            mergeMS.out.merged_stats,
            ''  // nothing for m2_extra_filtering_args
        )

        excIndel (
            filterMC.out.filtered_vcf,
            filterMC.out.filtered_vcf_tbi,
            '',  // variant type to include
            'INDEL',  // variant type to exclude
            'mutect2-snv'
        )

        selIndel (
            filterMC.out.filtered_vcf,
            filterMC.out.filtered_vcf_tbi,
            'INDEL',  // variant type to include
            '',  // variant type to exclude
            'mutect2-indel'
        )

        // genPayloadSNV
        pGenVarSnv(
            normal_aln_meta, tumour_aln_meta,
            excIndel.out.output.collect(),
            name, short_name, version
        )

        // genPayloadIndel
        pGenVarIndel(
            normal_aln_meta, tumour_aln_meta,
            selIndel.out.output.collect(),
            name, short_name, version
        )

        // prepQc
        prepQc(calCont.out.tumour_contamination_metrics.concat(
                    calCont.out.tumour_segmentation_metrics,
                    calCont.out.normal_contamination_metrics,
                    calCont.out.normal_segmentation_metrics,
                    filterMC.out.filtering_stats,
                    mergeMS.out.merged_stats
               ).collect())

        // genPayloadQc
        pGenQc(
            normal_aln_meta, tumour_aln_meta,
            prepQc.out.qc_metrics_tar.collect(),
            name, short_name, version
        )

        // skip upload if in local_mode
        if (!local_mode) {
            // uploadVariant
            upSnv(study_id, pGenVarSnv.out.payload, pGenVarSnv.out.files_to_upload)
            upIndel(study_id, pGenVarIndel.out.payload, pGenVarIndel.out.files_to_upload)

            // upQc
            upQc(study_id, pGenQc.out.payload, pGenQc.out.files_to_upload)
        }

        // cleanup, skip cleanup when running in local mode
        if (params.cleanup) {
            if (local_mode) {
                cleanupM2(
                    Mutect2.out.output_vcf.concat(
                        learnROM.out, mergeVcfs.out, mergeMS.out, calCont.out,
                        filterMC.out, excIndel.out, selIndel.out
                    ).collect(),
                    true
                )
            } else {
                cleanupM2(
                    dnldT.out.files.concat(
                        dnldN.out, Mutect2.out, learnROM.out, mergeVcfs.out, mergeMS.out, calCont.out,
                        filterMC.out, excIndel.out, selIndel.out
                    ).collect(),
                    upSnv.out.analysis_id.concat(upIndel.out.analysis_id, upQc.out.analysis_id).collect()
                )
            }

            if (params.perform_bqsr) {
                cleanupBqsr(
                    bqsrT.out.bqsr_bam.concat(bqsrN.out).collect(),
                    upSnv.out.analysis_id.concat(upIndel.out.analysis_id, upQc.out.analysis_id).collect()
                )
            }
        }

}


workflow {
    germline_resource_vcfs = Channel.fromPath(params.germline_resource_vcfs)
    germline_resource_indices = germline_resource_vcfs.flatMap { v -> getSec(v.name, ['tbi']) }

    panel_of_normals = Channel.fromPath(params.panel_of_normals)
    panel_of_normals_idx = panel_of_normals.flatMap { v -> getSec(v.name, ['tbi']) }

    contamination_variants = Channel.fromPath(params.contamination_variants)
    contamination_variants_indices = contamination_variants.flatMap { v -> getSec(v.name, ['tbi']) }

    M2(
        params.study_id,
        file(params.ref_fa),
        Channel.fromPath(getSec(params.ref_fa, ['^dict', 'fai']), checkIfExists: true).collect(),
        Channel.fromPath(getSec(params.ref_fa, ['^dict']), checkIfExists: true).collect(),
        Channel.fromPath(getSec(params.ref_fa, ['img'])),
        germline_resource_vcfs,
        germline_resource_indices,
        panel_of_normals,
        panel_of_normals_idx,
        contamination_variants,
        contamination_variants_indices,
        params.tumour_aln_analysis_id,
        params.normal_aln_analysis_id,
        params.tumour_aln_metadata,
        params.tumour_extra_info,
        params.tumour_aln_cram,
        params.normal_aln_metadata,
        params.normal_extra_info,
        params.normal_aln_cram,
        params.mutect2_scatter_interval_files,
        params.perform_bqsr,
        params.bqsr_recal_grouping_file,
        params.bqsr_apply_grouping_file
    )
}
