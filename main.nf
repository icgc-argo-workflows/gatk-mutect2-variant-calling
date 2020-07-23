#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
name = 'gatk-mutect2-variant-calling'
short_name = 'gatk-mutect2'
version = '4.1.8.0-1.0-dev'


/*
========================================================================================
                    ICGC-ARGO GATK Mutect2 Variant Calling Workflow
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
--tumour_aln_cram                       Tumour WGS aligned CRAM file (index file .crai is also required)
--normal_aln_metadata                   Normal WGS alignment metadata JSON file
--normal_aln_cram                       Normal WGS aligned CRAM file (index file .crai is also required)
--ref_genome_fa                         Reference genome '.fa' file, secondary file ('.fa.fai') is expected to be under the same folder
--sanger_ref_genome_tar                 Tarball containing reference genome files from the same genome build
--sanger_vagrent_annot                  Tarball containing VAGrENT annotation reference
--sanger_ref_cnv_sv_tar                 Tarball containing CNV/SV reference
--sanger_ref_snv_indel_tar              Tarball containing SNV/Indel reference
--sanger_qcset_tar                      Tarball containing QC Genotype reference
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

// if provided local files will be used
params.tumour_aln_metadata = "NO_FILE"
params.tumour_aln_cram = "NO_FILE"
params.normal_aln_metadata = "NO_FILE"
params.normal_aln_cram = "NO_FILE"

params.ref_fa = "NO_FILE"
params.api_token = ""
params.song_url = ""
params.score_url = ""
params.cleanup = true

params.cpus = 2
params.mem = 4

params.download = [:]
params.bqsr = [
    'dbsnp_vcf_gz': 'NO_FILE',
    'known_indels_sites_vcf_gzs': 'NO_FILE'
]
params.splitInterval
params.mutect2 = [
    'germline_resource': 'NO_FILE',
    'pon': 'NO_FILE'
]
params.gatherPileupSummaries = [
    'ref_dict': 'NO_FILE'
]
params.calculateContamination = [
    'variants_for_contamination': 'NO_FILE'
]
params.filterAlignmentArtifacts = [
    'bwa_mem_index_image': 'NO_FILE'
]
params.upload = [:]

download_params = [
    'song_url': params.song_url,
    'score_url': params.score_url,
    'api_token': params.api_token,
    *:(params.download ?: [:])
]

bqsr_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'ref_dict': params.ref_dict,
    'ref_fa': params.ref_fa,
    *:(params.bqsr ?: [:])
]

splitIntervals_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'ref_fa': params.ref_fa,
    *:(params.splitIntervals ?: [:])
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

filterMutectCalls_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.filterMutectCalls ?: [:])
]

filterAlignmentArtifacts_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.calculateContamination ?: [:])
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
include { bqsr as bqsrT; bqsr as bqsrN } from './bqsr/bqsr' params(bqsr_params)
include splitIntervals as splitItv from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-split-intervals.4.1.4.1-1.0/tools/gatk-split-intervals/gatk-split-intervals' params(splitIntervals_params)
// include mutect2 from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-mutect2.4.1.8.0-2.0/tools/gatk-mutect2/gatk-mutect2' params(sangerWxsVariantCaller_params)
// include learnReadOrientationModel as learnROM from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-learn-read-orientation-model.4.1.8.0-2.0/tools/gatk-learn-read-orientation-model/gatk-learn-read-orientation-model'
// include mergeVcfs from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-merge-vcfs.4.1.8.0-2.0/tools/gatk-merge-vcfs/gatk-merge-vcfs'
// include mergeMutectStats as mergeMS from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-merge-mutect-stats.4.1.8.0-2.0/tools/gatk-merge-mutect-stats/gatk-merge-mutect-stats'
// include calculateContamination as calCont from './calculate-contamination/calculate-contamination' params(calculateContamination_params)
// include filterMutectCalls as filterMC from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-filter-mutect-calls.4.1.8.0-2.0/tools/gatk-filter-mutect-calls/gatk-filter-mutect-calls' params(filterMutectCalls_params)
// include filterAlignmentArtifacts as filterAA from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-filter-alignment-artifacts.4.1.8.0-2.0/tools/gatk-filter-alignment-artifacts/gatk-filter-alignment-artifacts' params(filterAlignmentArtifacts_params)
include cleanupWorkdir as cleanup from './modules/raw.githubusercontent.com/icgc-argo/nextflow-data-processing-utility-tools/1.1.5/process/cleanup-workdir'


def getSecondaryFiles(main_file, exts){  //this is kind of like CWL's secondary files
    def all_files = []
    for (ext in exts) {
        all_files.add(main_file + ext)
    }
    return all_files
}


workflow M2 {
    take:
        study_id
        ref_fa
        tumour_aln_analysis_id
        normal_aln_analysis_id
        tumour_aln_metadata
        tumour_aln_cram
        normal_aln_metadata
        normal_aln_cram

    main:
        local_mode = false

        if (tumour_aln_analysis_id && normal_aln_analysis_id) {
            // download tumour aligned seq and metadata from song/score (analysis type: sequencing_alignment)
            dnldT(study_id, tumour_aln_analysis_id)

            // download normal aligned seq and metadata from song/score (analysis type: sequencing_alignment)
            dnldN(study_id, normal_aln_analysis_id)
        } else if (
            tumour_aln_metadata != 'NO_FILE' && \
            tumour_aln_cram != 'NO_FILE' && \
            normal_aln_metadata != 'NO_FILE' && \
            normal_aln_cram != 'NO_FILE'
        ) {
            local_mode = true
        } else {
            exit 1, "To download input aligned seq files from SONG/SCORE, please provide `params.tumour_aln_analysis_id` and `params.normal_aln_analysis_id`.\n" +
                "Or please provide `params.tumour_aln_metadata`, `params.tumour_aln_cram`, `params.normal_aln_metadata` and `params.normal_aln_cram` to use local files as input."
        }

        // splitInterval

        // BQSR Tumour

        // BQSR Normal

        // Mutect2

        // learnROM

        // mergeVcfs

        // mergeMS

        // calCont

        // filterMC

        // filterAA

        // genPayloadVariant

        // uploadVariant

        // genPayloadQC

        // upQC

        // genSuppl

        // upSuppl

        // cleanup

}


workflow {
    M2(
        params.study_id,
        params.ref_fa,
        params.tumour_aln_analysis_id,
        params.normal_aln_analysis_id,
        params.tumour_aln_metadata,
        params.tumour_aln_cram,
        params.normal_aln_metadata,
        params.normal_aln_cram
    )
}
