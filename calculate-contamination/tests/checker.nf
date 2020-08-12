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

include { calculateContamination as calCont } from '../calculate-contamination'
include { getSecondaryFiles } from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-base-recalibrator.4.1.8.0-1.0/tools/gatk-base-recalibrator/gatk-base-recalibrator'
include { gatkSplitIntervals as splitItvls } from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-split-intervals.4.1.4.1-1.0/tools/gatk-split-intervals/gatk-split-intervals'

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

    variants_resources = Channel.fromPath(params.variants_resources)

    variants_resources_indices = variants_resources.flatMap { v -> getSecondaryFiles(v, ['tbi']) }

    calCont(
      file(params.aln_seq),
      Channel.fromPath(getSecondaryFiles(params.aln_seq, ['bai', 'crai'])).collect(),
      file(params.match_aln_seq),
      Channel.fromPath(getSecondaryFiles(params.match_aln_seq, ['bai', 'crai'])).collect(),
      file(params.ref_genome_fa),
      Channel.fromPath(getSecondaryFiles(params.ref_genome_fa, ['^dict', 'fai'])).collect(),
      Channel.fromPath(getSecondaryFiles(params.ref_genome_fa, ['^dict'])).collect(),
      variants_resources.collect(),
      variants_resources_indices.collect(),
      interval_files.flatten(),
      params.tumour_normal
    )
}
