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
    path "${params.fasta_IGVF_acession}.fasta.gz"

  when:
    !file("references/${params.fasta_IGVF_acession}.noCHR.fasta").exists()

  output:
    path "${params.fasta_IGVF_acession}.noCHR.fasta"

  script:
  """
  zcat "${params.fasta_IGVF_acession}.fasta.gz" | sed 's/chr//g' > "${params.fasta_IGVF_acession}.noCHR.fasta"
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



workflow {
  curl_vcf()

  fasta_file = download_igvf_fasta()
  remove_chr_fasta(fasta_file)
  
  download_igvf_gtf()
}

