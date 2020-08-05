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
 */

nextflow.preview.dsl = 2

params.aln_seq = "NO_FILE"
params.ref_genome_fa = "NO_FILE"
params.known_sites_vcfs = []
params.sequence_group_interval = []

params.cpus = 2
params.mem = 4

include bqsr from '../bqsr'
include getSecondaryFiles from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-base-recalibrator.4.1.8.0-1.0/tools/gatk-base-recalibrator/gatk-base-recalibrator'
include gatkSplitIntervals as splitItvls from './modules/raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-split-intervals.4.1.4.1-1.0/tools/gatk-split-intervals/gatk-split-intervals'


Channel
  .fromPath(getSecondaryFiles(params.aln_seq, ['bai', 'crai']))
  .set { seq_crai_ch }

Channel
  .fromPath(getSecondaryFiles(params.ref_genome_fa, ['^dict', 'fai']), checkIfExists: true)
  .set { ref_genome_fai_ch }

known_sites_vcfs = Channel.fromPath(params.known_sites_vcfs)

known_sites_indices = known_sites_vcfs.flatMap { v -> getSecondaryFiles(v, ['tbi']) }


workflow {
  main:
    if (params.sequence_group_interval.size() == 0) {
        splitItvls(params.scatter_count, file(params.ref_genome_fa), ref_genome_fai_ch.collect(), file('NO_FILE'))
        interval_files = splitItvls.out.interval_files
    } else {
       Channel.fromPath(params.sequence_group_interval).set{ interval_files }
    }

    bqsr(
      file(params.aln_seq),
      seq_crai_ch.collect(),
      file(params.ref_genome_fa),
      ref_genome_fai_ch.collect(),  // secondary files: .fai and .dict
      known_sites_vcfs.collect(),
      known_sites_indices.collect(),
      interval_files.flatten(),
      'recalibrated_bam'
    )

}
