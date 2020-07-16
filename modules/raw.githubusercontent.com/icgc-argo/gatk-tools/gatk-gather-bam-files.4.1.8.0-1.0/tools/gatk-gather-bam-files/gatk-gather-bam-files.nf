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

params.input_bams = "NO_FILE"
params.output_bam_basename = "final_bam"
params.compression_level = 5

params.cpus = 1
params.mem = 1  // in GB
params.container_version = ""


process gatkGatherBamFiles {
  container "quay.io/icgc-argo/gatk-gather-bam-files:gatk-gather-bam-files.${params.container_version ?: version}"
  cpus params.cpus
  memory "${params.mem} GB"

  input:
    path input_bams
    val output_bam_basename
    val compression_level

  output:
    path "${output_bam_basename}.bam", emit: merged_bam
    path "${output_bam_basename}.bam.bai", emit: merged_bam_bai
    path "${output_bam_basename}.bam.md5", emit: merged_bam_md5

  script:
    arg_compression_level = compression_level != "null" ? compression_level : 5

    """
    gatk-gather-bam-files.py -m ${(int) (params.mem * 1000)} \
      -i ${input_bams} -o ${output_bam_basename} \
      -c ${arg_compression_level}

    ln -s ${output_bam_basename}.bai ${output_bam_basename}.bam.bai
    """
}
