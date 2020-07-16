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

nextflow.preview.dsl = 2
version = '4.1.8.0-1.0'

params.seq = "NO_FILE"
params.interval_file = "NO_FILE"
params.ref_genome_fa = "NO_FILE"
params.recalibration_report = "NO_FILE"
params.cpus = 1
params.mem = 1  // in GB
params.container_version = ""


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



process gatkApplyBQSR {
  container "quay.io/icgc-argo/gatk-apply-bqsr:gatk-apply-bqsr.${params.container_version ?: version}"
  cpus params.cpus
  memory "${params.mem} GB"

  input:
    path seq
    path seq_idx
    path ref_genome_fa
    path ref_genome_secondary_file
    path recalibration_report
    path interval_file
    val output_bam_basename

  output:
    path "${output_bam_basename}.bam", emit: recalibrated_bam
    path "${output_bam_basename}.bam.bai", emit: recalibrated_bam_bai
    path "${output_bam_basename}.bam.md5", emit: recalibrated_bam_md5

  script:
    arg_interval_file = interval_file.name == 'NO_FILE' ? "" : "-i ${interval_file}"

    """
    gatk-apply-bqsr.py -s ${seq} \
                      -r ${ref_genome_fa} \
                      -m ${(int) (params.mem * 1000)} \
                      -c ${recalibration_report} \
                      -o ${output_bam_basename} ${arg_interval_file}

    ln -s ${output_bam_basename}.bai ${output_bam_basename}.bam.bai
    """
}
