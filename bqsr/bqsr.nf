#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
name = 'gatk-bqsr'
short_name = 'bqsr'
version = '4.1.8.0-1.0-dev'

/*
  ========================================================================================
                GATK Base Quality Score Recalibration Nextflow Workflow
  ========================================================================================
  #### Homepage / Documentation
  https://github.com/icgc-argo/gatk-mutect2-variant-calling/tree/master/bqsr
  #### Authors
  Junjun Zhang @junjun-zhang <junjun.zhang@oicr.on.ca>
  Linda Xiang @lindaxiang <linda.xiang@oicr.on.ca>
  ----------------------------------------------------------------------------------------

Required Parameters (no default):
--aln_seq                       Aligned sequence file
--ref_fa                        Reference genome '.fa' file, secondary file ('.fa.fai') is expected to be under the same folder
--dbsnp_vcf_gz                  dbSNP variants in a VCF file
--known_indels_sites_vcf_gzs    VCFs for known indel sites
--sequence_group_interval       intervals used to perform scattered execution
*/

params.ref_fa = "NO_FILE"
params.ref_dict = "NO_FILE"
params.dbsnp_vcf_gz = "NO_FILE"
params.known_indels_sites_vcf_gzs = []
params.sequence_group_interval = []
params.cpus = 2
params.mem = 4


// Include all modules and pass params
include BaseRecalibrator as recalibrate from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-base-recalibrator.4.1.8.0-1.0/tools/gatk-base-recalibrator/gatk-base-recalibrator'
include GatherBqsrReports as gatherBS from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-gather-bqsr-reports.4.1.8.0-1.0/tools/gatk-gather-bqsr-reports/gatk-gather-bqsr-reports'
include ApplyBQSR as applyBQSR from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-apply-bqsr.4.1.8.0-1.0/tools/gatk-apply-bqsr/gatk-apply-bqsr'
include GatherBamFiles as gatherBAMs from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-gather-bam-files.4.1.8.0-1.0/tools/gatk-gather-bam-files/gatk-gather-bam-files'


def getSecondaryFiles(main_file, exts){  //this is kind of like CWL's secondary files
    def all_files = []
    for (ext in exts) {
        all_files.add(main_file + ext)
    }
    return all_files
}


workflow bqsr {
    take:
        sequence_group_interval
        aln_seq


    main:
        // BaseRecalibrator
        recalibrate()

        // GatherBqsrReports

        // ApplyBQSR

        // GatherBamFiles

}


workflow {
    bqsr(
        params.sequence_group_interval,
        params.aln_seq
    )
}
