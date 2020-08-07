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
params.intervals = []
params.container_version = ""
params.ref_genome_fa = "NO_FILE"
params.known_sites_vcfs = "NO_FILE"
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



process gatkBaseRecalibrator {
  container "quay.io/icgc-argo/gatk-base-recalibrator:gatk-base-recalibrator.${params.container_version ?: version}"
  cpus params.cpus
  memory "${params.mem} GB"

  input:
    path seq
    path seq_idx
    path ref_genome_fa
    path ref_genome_secondary_file
    path known_sites_vcfs
    path known_sites_indices
    val intervals

  output:
    path "${arg_output_prefix}.recal_data.csv", emit: recalibration_report

  script:
    if (!intervals) {
      arg_intervals = ""
      arg_output_prefix = seq
    } else {
      (index, intervals) = intervals
      arg_intervals = "-i ${intervals}"
      arg_output_prefix = "${index.toString().padLeft(6, '0')}.${seq}"
    }
    arg_known_sites = known_sites_vcfs.name == 'NO_FILE' ? "" : "-k ${known_sites_vcfs}"

    """
    gatk-base-recalibrator.py -s ${seq} \
                      -r ${ref_genome_fa} \
                      -m ${(int) (params.mem * 1000)} \
                      -o ${arg_output_prefix}.recal_data.csv ${arg_intervals} ${arg_known_sites}
    """
}
