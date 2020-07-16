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

params.input_bqsr_reports = "NO_FILE"
params.output_report_filename = "bqsr_report.csv"

params.cpus = 1
params.mem = 1  // in GB
params.container_version = ""


process gatkGatherBQSRReports {
  container "quay.io/icgc-argo/gatk-gather-bqsr-reports:gatk-gather-bqsr-reports.${params.container_version ?: version}"
  cpus params.cpus
  memory "${params.mem} GB"

  input:
    path input_bqsr_reports
    val output_report_filename

  output:
    path "${output_report_filename}", emit: bqsr_report

  script:
    """
    gatk-gather-bqsr-reports.py -m ${(int) (params.mem * 1000)} \
      -i ${input_bqsr_reports} -o ${output_report_filename}
    """
}
