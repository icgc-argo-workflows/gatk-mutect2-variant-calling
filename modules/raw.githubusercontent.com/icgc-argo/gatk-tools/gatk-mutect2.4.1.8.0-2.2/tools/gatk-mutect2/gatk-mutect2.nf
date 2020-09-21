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
 */

nextflow.enable.dsl = 2
version = '4.1.8.0-2.2'

params.tumour_reads = "NO_FILE1"
params.normal_reads = "NO_FILE2"
params.interval_file = "NO_FILE3"
params.ref_genome_fa = "NO_FILE4"
params.germline_resource = "NO_FILE5"
params.panel_of_normals = "NO_FILES6"

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



process gatkMutect2 {
  container "quay.io/icgc-argo/gatk-mutect2:gatk-mutect2.${params.container_version ?: version}"
  cpus params.cpus
  memory "${params.mem} GB"
  //publishDir "output/mutect2"

  input:
    path tumour_reads
    path tumour_idx
    path normal_reads
    path normal_idx
    path ref_genome_fa
    path ref_genome_secondary_file
    path germline_resource
    path germline_resource_idx
    path panel_of_normals
    path panel_of_normals_idx
    path interval_file

  output:
    path "*${arg_output_prefix}.vcf.gz", emit: output_vcf
    path "*${arg_output_prefix}.vcf.gz.tbi", emit: output_vcf_tbi
    path "*${arg_output_prefix}.vcf.gz.stats", emit: mutect_stats
    path "*${arg_output_prefix}.bamout.bam", emit: bam_output
    path "*${arg_output_prefix}.bamout.bai", emit: bam_output_bai
    path "*${arg_output_prefix}.f1r2.tar.gz", emit: f1r2_counts

  script:
    arg_interval_file = interval_file.name.startsWith('NO_FILE') ? "" : "-L ${interval_file}"
    arg_output_prefix = tumour_reads
    arg_germline_resource = germline_resource.name.startsWith('NO_FILE') ? "" : "-g ${germline_resource}"
    arg_panel_of_normals = panel_of_normals.name.startsWith('NO_FILE') ? "" : "-p ${panel_of_normals}"

    """
    gatk-mutect2.py \
      -j ${(int) (params.mem * 1000)} \
      -R ${ref_genome_fa} \
      -t ${tumour_reads} -n ${normal_reads} \
      -O ${arg_output_prefix}.vcf.gz \
      -b ${arg_output_prefix}.bamout.bam \
      -f ${arg_output_prefix}.f1r2.tar.gz \
      ${arg_interval_file} \
      ${arg_germline_resource} \
      ${arg_panel_of_normals}
    """
}
