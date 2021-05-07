#!/usr/bin/env nextflow

/*
  Copyright (C) 2021,  Ontario Institue for Cancer Research

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Authors:
    Junjun Zhang
*/

/*
 This is an auto-generated checker workflow to test the generated main template workflow, it's
 meant to illustrate how testing works. Please update to suit your own needs.
*/

nextflow.enable.dsl = 2
version = '4.1.8.0-6.0'  // package version

// universal params
params.publish_dir = ""
params.container = ""
params.container_registry = ""
params.container_version = ""

// tool specific parmas go here, add / change as needed
params.study_id = ""

// aligned seq will be downloaded from SONG/SCORE
params.tumour_aln_analysis_id = ""
params.normal_aln_analysis_id = ""

// the following params if provided local files will be used
params.tumour_aln_metadata = "NO_FILE1"
params.tumour_aln_cram = "NO_FILE2"
params.tumour_extra_info = "NO_FILE3"
params.normal_aln_metadata = "NO_FILE4"
params.normal_aln_cram = "NO_FILE5"
params.normal_extra_info = "NO_FILE6"

params.perform_bqsr = true  // default to true

params.ref_fa = "reference/tiny-grch38-chr11-530001-537000.fa"

params.mutect2_scatter_interval_files = "assets/mutect2.scatter_by_chr/chr*.interval_list"
params.bqsr_recal_grouping_file = "assets/bqsr.sequence_grouping.grch38_hla_decoy_ebv.csv"
params.bqsr_apply_grouping_file = "assets/bqsr.sequence_grouping_with_unmapped.grch38_hla_decoy_ebv.csv"

// Allele frequency only, pass-only gnomAD vcf file
params.germline_resource_vcfs = []  // "tests/data/HCC1143-mini-Mutect2-calls/HCC1143.mutect2.copy.vcf.gz"

params.panel_of_normals = "NO_FILE"  // optional PoN VCF file, for now we use GATK's public PoN based on 1000 genome project

params.contamination_variants = ""

params.api_token = ""
params.song_url = ""
params.score_url = ""
params.cleanup = true

params.cpus = 2
params.mem = 4

include { M2 } from '../main'
include { getSecondaryFiles as getSec } from './wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/helper-functions@1.0.1/main'


workflow {
    germline_resource_vcfs = Channel.fromPath(params.germline_resource_vcfs)
    germline_resource_indices = germline_resource_vcfs.flatMap { v -> getSec(v, ['tbi']) }

    panel_of_normals = Channel.fromPath(params.panel_of_normals)
    panel_of_normals_idx = panel_of_normals.flatMap { v -> getSec(v, ['tbi']) }

    contamination_variants = Channel.fromPath(params.contamination_variants)
    contamination_variants_indices = contamination_variants.flatMap { v -> getSec(v, ['tbi']) }

    M2(
        params.study_id,
        file(params.ref_fa),
        Channel.fromPath(getSec(params.ref_fa, ['^dict', 'fai']), checkIfExists: true).collect(),
        Channel.fromPath(getSec(params.ref_fa, ['^dict']), checkIfExists: true).collect(),
        Channel.fromPath(getSec(params.ref_fa, ['img'])),
        germline_resource_vcfs,
        germline_resource_indices,
        panel_of_normals,
        panel_of_normals_idx,
        contamination_variants,
        contamination_variants_indices,
        params.tumour_aln_analysis_id,
        params.normal_aln_analysis_id,
        params.tumour_aln_metadata,
        params.tumour_extra_info,
        params.tumour_aln_cram,
        params.normal_aln_metadata,
        params.normal_extra_info,
        params.normal_aln_cram,
        params.mutect2_scatter_interval_files,
        params.perform_bqsr,
        params.bqsr_recal_grouping_file,
        params.bqsr_apply_grouping_file
    )
}
