#!/usr/bin/env nextflow

/*
 * Copyright (c) 2019-2020, Ontario Institute for Cancer Research (OICR).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/*
 * author Junjun Zhang <junjun.zhang@oicr.on.ca>
 *         Linda Xiang <linda.xiang@oicr.on.ca>
 */

nextflow.preview.dsl = 2
version = '4.1.8.0-2.0'

params.unfiltered_vcf = "NO_FILE"
params.ref_genome_fa = "NO_FILE"
params.contamination_table = "NO_FILE_cont"
params.segmentation_table = "NO_FILE_seg"
params.artifact_priors_tar_gz = "NO_FILE_ob"
params.mutect_stats = "NO_FILE_stats"
params.m2_extra_filtering_args = ""

params.container_version = ""
params.cpus = 1
params.mem = 1  // in GB


def getSecondaryFiles(main_file, exts){
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

process gatkFilterMutectCalls {
  container "quay.io/icgc-argo/gatk-filter-mutect-calls:gatk-filter-mutect-calls.${params.container_version ?: version}"
  cpus params.cpus
  memory "${params.mem} GB"
  publishDir "output/filterMutectCalls"

  input:
    path unfiltered_vcf
    path unfiltered_vcf_tbi
    path ref_genome_fa
    path ref_genome_secondary_file
    path contamination_table
    path segmentation_table
    path artifact_priors_tar_gz
    path mutect_stats
    val m2_extra_filtering_args

  output:
    path "*.filtering-stats", emit: filtering_stats
    path "*.filtered.vcf.gz", emit: filtered_vcf
    path "*.filtered.vcf.gz.tbi", emit: filtered_vcf_tbi

  script:
    arg_contamination_table = contamination_table.name == 'NO_FILE_cont' ? "" : " --contamination-table ${contamination_table}"
    arg_segmentation_table = segmentation_table.name == 'NO_FILE_seg' ? "" : " --tumor-segmentation ${segmentation_table}"
    arg_artifact_priors_tar_gz = artifact_priors_tar_gz.name == 'NO_FILE_ob' ? "" : " --ob-priors ${artifact_priors_tar_gz}"
    arg_m2_extra_filtering_args = m2_extra_filtering_args == '' ? "" : " -e ${m2_extra_filtering_args}"
    
    """
    gatk-filter-mutect-calls.py -V ${unfiltered_vcf} \
                      -R ${ref_genome_fa} \
                      -j ${(int) (params.mem * 1000)} \
                      --stats ${mutect_stats} \
                      ${arg_contamination_table} \
                      ${arg_segmentation_table} \
                      ${arg_artifact_priors_tar_gz} \
                      ${arg_m2_extra_filtering_args}
    """
}
