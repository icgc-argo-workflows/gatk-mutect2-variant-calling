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

getPileupSummaries_params = [
    *:(params.getPileupSummaries ?: [:]) 
]

gatherPileupSummaries_params = [
    *:(params.gatherPileupSummaries ?: [:]) 
]

// Include all modules and pass params
include { gatkGetPileupSummaries as getPS; gatkGetPileupSummaries as getPSM } from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-get-pileup-summaries.4.1.8.0-2.0/tools/gatk-get-pileup-summaries/gatk-get-pileup-summaries' params(getPileupSummaries_params)
include { gatkGatherPileupSummaries as gatherPS; gatkGatherPileupSummaries as gatherPSM } from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-gather-pileup-summaries.4.1.8.0-2.0/tools/gatk-gather-pileup-summaries/gatk-gather-pileup-summaries' params(gatherPileupSummaries_params)
include { gatkCalculateContamination as calCont } from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-calculate-contamination.4.1.8.0-2.0/tools/gatk-calculate-contamination/gatk-calculate-contamination'


def getSecondaryFiles(main_file, exts){  //this is kind of like CWL's secondary files
  def secondaryFiles = []
  for (ext in exts) {
    if (ext.startsWith("^")) {
      ext = ext.replace("^", "")
      parts = main_file.split("\\.").toList()
      parts.removeLast()
      secondaryFiles.add((parts + [ext]).join("."))
    } else {
      secondaryFiles.add(main_file + '.' + ext)
    }
  }
  return secondaryFiles
}


workflow calculateContaminationWf {
    take:
        aln_seq
        match_aln_seq
        ref_genome_fa
        variants_resources
        interval_files
        tumour_normal

    main:
        // getPS
        getPS(
            file(aln_seq),
            Channel.fromPath(getSecondaryFiles(aln_seq, ['bai', 'crai'])).collect(),
            file(ref_genome_fa),
            Channel.fromPath(getSecondaryFiles(ref_genome_fa, ['fai', '^dict']), checkIfExists: true).collect(),
            file(variants_resources),
            Channel.fromPath(getSecondaryFiles(variants_resources, ['tbi']), checkIfExists: true).collect(),
            interval_files.flatten()
        )

        // gatherPS
        gatherPS(
            Channel.fromPath(getSecondaryFiles(ref_genome_fa, ['^dict']), checkIfExists: true).collect(),
            getPS.out.pileups_metrics.collect()
        )

        if (match_aln_seq != 'NO_FILE') {
            // getPSM
            getPSM(
                file(match_aln_seq),
                Channel.fromPath(getSecondaryFiles(match_aln_seq, ['bai', 'crai'])).collect(),
                file(ref_genome_fa),
                Channel.fromPath(getSecondaryFiles(ref_genome_fa, ['fai', '^dict']), checkIfExists: true).collect(),
                file(variants_resources),
                Channel.fromPath(getSecondaryFiles(variants_resources, ['tbi']), checkIfExists: true).collect(),
                interval_files.flatten()
            )

            // gatherPSM
            gatherPSM(
                Channel.fromPath(getSecondaryFiles(ref_genome_fa, ['^dict']), checkIfExists: true).collect(),
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


workflow {
    calculateContaminationWf(
        params.aln_seq,
        params.match_aln_seq,
        params.ref_genome_fa,
        params.variants_resources,
        Channel.fromPath(params.interval_files),
        params.tumour_normal
    )
}
