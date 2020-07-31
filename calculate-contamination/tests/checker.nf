#!/usr/bin/env nextflow

/*
 * Copyright (c) 2020, Ontario Institute for Cancer Research (OICR).
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
 *        Linda Xiang <linda.xiang@oicr.on.ca>
 */

nextflow.preview.dsl = 2

params.aln_seq = "NO_FILE"
params.match_aln_seq = "NO_FILE"
params.ref_genome_fa = "NO_FILE"
params.variants_resources = "NO_FILE"
params.interval_files = []
params.tumour_normal = ""

params.cpus = 2
params.mem = 4

include { calculateContaminationWf as calCont; getSecondaryFiles } from '../calculate-contamination'
include { gatkSplitIntervals as splitItvls } from '../modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-split-intervals.4.1.4.1-1.0/tools/gatk-split-intervals/gatk-split-intervals'

Channel
  .fromPath(getSecondaryFiles(params.ref_genome_fa, ['^dict', 'fai']), checkIfExists: true)
  .set { ref_genome_fai_ch }

workflow {
  main:

    if (params.interval_files.size() == 0) {
        splitItvls(params.scatter_count, file(params.ref_genome_fa), ref_genome_fai_ch.collect(), file('NO_FILE'))
        interval_files = splitItvls.out.interval_files
    } else {
        interval_files = Channel.fromPath(params.interval_files)
    }

    calCont(
      params.aln_seq,
      params.match_aln_seq,
      params.ref_genome_fa,
      params.variants_resources,
      interval_files,
      params.tumour_normal
    )
}
