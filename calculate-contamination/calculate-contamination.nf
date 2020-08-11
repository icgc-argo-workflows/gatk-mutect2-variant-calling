#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
name = 'gatk-calculate-contamination'
short_name = 'contamination'
version = '4.1.8.0-1.0-dev'


/*
  ==============================================================================================
                    GATK Calculate Contamination Nextflow Workflow
  ==============================================================================================
  #### Homepage / Documentation
  https://github.com/icgc-argo/gatk-mutect2-variant-calling/tree/master/calculate-contamination
  #### Authors
  Junjun Zhang @junjun-zhang <junjun.zhang@oicr.on.ca>
  Linda Xiang @lindaxiang <linda.xiang@oicr.on.ca>
  ----------------------------------------------------------------------------------------------

Required Parameters (no default):
--tumour_aln_seq                Aligned tumour sequence file
--normal_aln_seq                Aligned normal sequence file
--ref_fa                        Reference genome '.fa' file, secondary file ('.fa.fai') is expected to be under the same folder
--known_variants                Known variants in a VCF file

*/

params.cpus = 2
params.mem = 4


params.aln_seq = "NO_FILE"
params.match_aln_seq = "NO_FILE"
params.ref_genome_fa = "NO_FILE"
params.variants_resources = "NO_FILE"
params.interval_files = []
params.tumour_normal = ""

params.getPileupSummaries = [:]
params.gatherPileupSummaries = [:]
params.calculateContamination = [:]

getPileupSummaries_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.getPileupSummaries ?: [:]) 
]

gatherPileupSummaries_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.gatherPileupSummaries ?: [:]) 
]

calculateContamination_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.calculateContamination ?: [:]) 
]

// Include all modules and pass params
include { gatkGetPileupSummaries as getPS; gatkGetPileupSummaries as getPSM } from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-get-pileup-summaries.4.1.8.0-2.0/tools/gatk-get-pileup-summaries/gatk-get-pileup-summaries' params(getPileupSummaries_params)
include { gatkGatherPileupSummaries as gatherPS; gatkGatherPileupSummaries as gatherPSM } from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-gather-pileup-summaries.4.1.8.0-2.0/tools/gatk-gather-pileup-summaries/gatk-gather-pileup-summaries' params(gatherPileupSummaries_params)
include { gatkCalculateContamination as calCont } from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-calculate-contamination.4.1.8.0-2.0/tools/gatk-calculate-contamination/gatk-calculate-contamination' params(calculateContamination_params)


workflow calculateContamination {
    take:
        aln_seq
        aln_seq_idx
        match_aln_seq
        match_aln_seq_idx
        ref_genome_fa
        ref_genome_fa_2nd
        ref_genome_fa_dict
        variants_resources
        variants_resources_indices
        interval_files
        tumour_normal

    main:
        // getPS
        getPS(
            aln_seq,
            aln_seq_idx.collect(),
            ref_genome_fa,
            ref_genome_fa_2nd.collect(),
            variants_resources.collect(),
            variants_resources_indices.collect(),
            interval_files.flatten()
        )

        // gatherPS
        gatherPS(
            ref_genome_fa_dict.collect(),
            getPS.out.pileups_metrics.collect()
        )

        if (match_aln_seq.name != 'NO_FILE') {
            // getPSM
            getPSM(
                match_aln_seq,
                match_aln_seq_idx.collect(),
                ref_genome_fa,
                ref_genome_fa_2nd.collect(),
                variants_resources.collect(),
                variants_resources_indices.collect(),
                interval_files.flatten()
            )

            // gatherPSM
            gatherPSM(
                ref_genome_fa_dict.collect(),
                getPSM.out.pileups_metrics.collect()
            )
            
            // calCont
            calCont(
                gatherPS.out.merged_pileups_metrics,
                gatherPSM.out.merged_pileups_metrics,
                'tumour'
            )
        } else {
            // calCont
            calCont(
                gatherPS.out.merged_pileups_metrics,
                file("NO_FILE"),
                tumour_normal
            )
        }

    emit:
        segmentation_metrics = calCont.out.segmentation_metrics
        contamination_metrics = calCont.out.contamination_metrics

}

