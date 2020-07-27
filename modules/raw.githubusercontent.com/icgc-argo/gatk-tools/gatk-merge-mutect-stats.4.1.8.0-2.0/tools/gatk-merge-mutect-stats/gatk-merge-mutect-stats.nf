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

params.input_stats = "NO_FILE"
params.output_basename = "NO_FILE"

params.container_version = ""
params.cpus = 1
params.mem = 1  // in GB


process gatkMergeMutectStats {
  container "quay.io/icgc-argo/gatk-merge-mutect-stats:gatk-merge-mutect-stats.${params.container_version ?: version}"
  cpus params.cpus
  memory "${params.mem} GB"

  input:
    path input_stats
    val output_basename

  output:
    path "${output_basename}.stats", emit: merged_stats

  script:
    """
    gatk-merge-mutect-stats.py \
                      -j ${(int) (params.mem * 1000)} \
                      -I ${input_stats} \
                      -O ${output_basename}.stats
    """
}
