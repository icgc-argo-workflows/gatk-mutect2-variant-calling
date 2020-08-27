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

nextflow.enable.dsl = 2
version = '4.1.8.0-1.0'

params.input_vcf = "NO_FILE"
params.select_type_to_include = ""
params.select_type_to_exclude = ""
params.output_basename = "my_variants"

params.container_version = ""
params.cpus = 1
params.mem = 1  // in GB


process gatkSelectVariants {
  container "quay.io/icgc-argo/gatk-select-variants:gatk-select-variants.${params.container_version ?: version}"
  cpus params.cpus
  memory "${params.mem} GB"

  input:
    path input_vcf
    path input_vcf_tbi
    val select_type_to_include
    val select_type_to_exclude
    val output_basename

  output:
    tuple path("${output_basename}.vcf.gz"), path("${output_basename}.vcf.gz.tbi"), emit: output

  script:
    variant_types = ["INDEL", "SNP", "MIXED", "MNP", "SYMBOLIC", "NO_VARIATION"]
    if (!select_type_to_include && !select_type_to_exclude) {
      exit 1, "Please specify either sselect_type_to_include or select_type_to_exclude"
    } else if (select_type_to_include && select_type_to_exclude){
      exit 1, "Please specify either select_type_to_include or select_type_to_exclude"
    } else if (select_type_to_include && !variant_types.contains(select_type_to_include)) {
      exit 1, "select_type_to_include must be one of: ${variant_types}"
    } else if (select_type_to_exclude && !variant_types.contains(select_type_to_exclude)) {
      exit 1, "select_type_to_exclude must be one of: ${variant_types}"
    }

    arg_select = "-i ${select_type_to_include}"
    if (select_type_to_exclude) {
        arg_select = "-x ${select_type_to_exclude}"
    }

    """
    gatk-select-variants.py \
                      -m ${(int) (params.mem * 1000)} \
                      -v ${input_vcf} \
                      -o ${output_basename}.vcf.gz ${arg_select}

    """
}
