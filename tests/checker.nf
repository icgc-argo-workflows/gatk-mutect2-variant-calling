#!/bin/bash nextflow

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

params.tumour_aln_analysis_id = ""
params.normal_aln_analysis_id = ""
params.ref_fa = "reference/tiny-grch38-chr11-530001-537000.fa"
params.mutect2 = [
  'germline_resource': 'reference/tiny-chr11-exac_common_3.hg38.vcf.gz'
]
params.gatherPileupSummaries = [
  'ref_dict': 'reference/tiny-grch38-chr11-530001-537000.dict'
]
params.calculateContamination = [
  'variants_for_contamination': 'reference/tiny-chr11-exac_common_3.hg38.vcf.gz'
]
params.filterAlignmentArtifacts = [
    'bwa_mem_index_image': 'reference/tiny-grch38-chr11-530001-537000.fa.img'
]
params.api_token = ""
params.song_url = ""
params.score_url = ""
params.cpus = 1
params.mem = 1  // GB
params.cleanup = true

include BroadMutect2 from "../main" params(params)


workflow {
  BroadMutect2(
    params.study_id,
    params.tumour_aln_analysis_id,
    params.normal_aln_analysis_id
  )
}
