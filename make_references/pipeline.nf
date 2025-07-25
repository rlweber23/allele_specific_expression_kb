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

process download_igvf_fasta {
  publishDir "references/", mode: 'copy'

  when:
    !file("references/${params.fasta_IGVF_acession}.fasta.gz").exists()

  output:
    path "${params.fasta_IGVF_acession}.fasta.gz"

  script:
  """
  python ${params.download_igvf_portal} --acc ${params.fasta_IGVF_acession} 
  """
}


process remove_chr_fasta {
  publishDir "references/", mode: 'copy'

  input:
    path fasta_file

  output:
    path "${fasta_file.simpleName}.noCHR.fasta"

  when:
      !file("${params.topDir}/references/${params.fasta_IGVF_acession}.noCHR.fasta").exists()

  script:
  """
  sed 's/chr//g' "${fasta_file}" > "${fasta_file.simpleName}.noCHR.fasta"
  """
}



process download_igvf_gtf {
  publishDir "references/", mode: 'copy'

  when:
    !file("references/${params.gtf_IGVF_acession}.gtf.gz").exists()

  output:
    path "${params.gtf_IGVF_acession}.gtf.gz"

  script:
  """
  python ${params.download_igvf_portal} --acc ${params.gtf_IGVF_acession} 
  """
}



process remove_chr_gtf {
  publishDir "references/", mode: 'copy'

  input:
    path gtf_file

  output:
    path "${params.gtf_IGVF_acession}.noCHR.gtf"

  when:
      !file("${params.topDir}/references/${params.gtf_IGVF_acession}.noCHR.gtf").exists()

  script:
  """
  sed 's/chr//g' "${gtf_file}" > "${gtf_file.simpleName}.noCHR.gtf"
  """
}



workflow {
  curl_vcf()

  fasta_file = download_igvf_fasta()
  remove_chr_fasta(fasta_file)
  
  gtf_file = download_igvf_gtf()
  remove_chr_gtf(gtf_file)
}

