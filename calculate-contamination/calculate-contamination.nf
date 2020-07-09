#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
name = 'gatk-calculate-contamination'
short_name = 'contamination'
version = '4.1.7.0-1.0-dev'


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

params.ref_fa = "NO_FILE"
params.known_variants = "NO_FILE"
params.cpus = 2
params.mem = 4


// Include all modules and pass params

include { getPileupSummaries as getPST; getPileupSummaries as getPSN } from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-get-pileup-summaries.4.1.7.0-2.0/tools/gatk-get-pileup-summaries/gatk-get-pileup-summaries' params(getPileupSummaries_params)
include gatherPileupSummaries as gatherPS from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-gather-pileup-summaries.4.1.7.0-2.0/tools/gatk-gather-pileup-summaries/gatk-gather-pileup-summaries' params(gatherPileupSummaries_params)
include calculateContamination as calCont from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-calculate-contamination.4.1.7.0-2.0/tools/gatk-calculate-contamination/gatk-calculate-contamination'

def getSecondaryFiles(main_file, exts){  //this is kind of like CWL's secondary files
    def all_files = []
    for (ext in exts) {
        all_files.add(main_file + ext)
    }
    return all_files
}


workflow calculateContaminationWf {
    take:
        intervals
        tumour_aln_seq
        normal_aln_seq


    main:
        // getPST

        // getPSN

        // gatherPS

        // calCont

}


workflow {
    calculateContaminationWf(
        params.study_id,
        params.tumour_aln_analysis_id,
        params.normal_aln_analysis_id
    )
}
