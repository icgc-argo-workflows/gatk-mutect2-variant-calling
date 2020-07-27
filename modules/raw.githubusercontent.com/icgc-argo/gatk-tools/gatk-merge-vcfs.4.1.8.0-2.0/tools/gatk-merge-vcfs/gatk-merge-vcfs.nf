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
version = '4.1.8.0-2.0'

params.input_vcfs = "NO_FILE"
params.output_basename = "NO_FILE"

params.container_version = ""
params.cpus = 1
params.mem = 1  // in GB


process gatkMergeVcfs {
  container "quay.io/icgc-argo/gatk-merge-vcfs:gatk-merge-vcfs.${params.container_version ?: version}"
  cpus params.cpus
  memory "${params.mem} GB"

  input:
    path input_vcfs
    val output_basename

  output:
    path "${output_basename}.vcf.gz", emit: output_vcf
    path "${output_basename}.vcf.gz.tbi", emit: output_tbi

  script:
    """
    gatk-merge-vcfs.py \
                      -j ${(int) (params.mem * 1000)} \
                      -I ${input_vcfs} \
                      -O ${output_basename}.vcf.gz
    """
}
