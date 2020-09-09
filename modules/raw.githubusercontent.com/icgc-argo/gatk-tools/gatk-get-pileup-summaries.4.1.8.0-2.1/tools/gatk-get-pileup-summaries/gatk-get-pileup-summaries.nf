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
 *        Linda Xiang  <linda.xiang@oicr.on.ca>
 */

nextflow.preview.dsl = 2
version = '4.1.8.0-2.1'

params.seq = "NO_FILE1"
params.ref_genome_fa = "NO_FILE2"
params.variants_resources = "NO_FILE3"
params.interval_file = "NO_FILE4"

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

process gatkGetPileupSummaries {
  container "quay.io/icgc-argo/gatk-get-pileup-summaries:gatk-get-pileup-summaries.${params.container_version ?: version}"
  cpus params.cpus
  memory "${params.mem} GB"
  //publishDir "output/getPileupSummaries"

  input:
    path seq
    path seq_idx
    path ref_genome_fa
    path ref_genome_secondary_file
    path variants_resources
    path variants_secondary_file
    path interval_file

  output:
    path "*.pileups_metrics.txt", emit: pileups_metrics

  script:
    arg_ref_genome_fa = ref_genome_fa.name.startsWith('NO_FILE') ? "" : "-R ${ref_genome_fa}"
    arg_interval_file = interval_file.name.startsWith('NO_FILE') ? "" : "-L ${interval_file}"
    """
    gatk-get-pileup-summaries.py -I ${seq} \
                      -V ${variants_resources} \
                      ${arg_interval_file} \
                      -j ${(int) (params.mem * 1000)} \
                      ${arg_ref_genome_fa}
    """
}
