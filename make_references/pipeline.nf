#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process curl_vcf {
  publishDir "references/", mode: 'copy'

  when:
    !file("references/${ params.mouse_vcf.split('/')[-1] }").exists()

  output:
    path "${ params.mouse_vcf.split('/')[-1] }"

  script:
  """
    curl -O ${params.mouse_vcf}
  """
}

workflow {
  curl_vcf()
}

