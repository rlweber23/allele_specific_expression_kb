#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process curl_vcf {
  storeDir "references/"

  //when:
  //  !file("references/${ params.mouse_vcf.split('/')[-1] }").exists()

  output:
    path "${ params.mouse_vcf.split('/')[-1] }"

  script:
  """
    curl -O ${params.mouse_vcf}
  """
}

process download_igvf_fasta {
  storeDir "references/"

  when:
    !file("references/${params.fasta_IGVF_acession}.fasta.gz").exists()

  //output:
  //  path "${params.fasta_IGVF_acession}.fasta.gz"

  script:
  """
  python ${params.download_igvf_portal} --acc ${params.fasta_IGVF_acession} 
  """
}


process remove_chr_fasta {
  storeDir "references/"

  input:
    path fasta_file

  output:
    path "${fasta_file.simpleName}.noCHR.fasta"

  //when:
  //    !file("${params.topDir}/references/${params.fasta_IGVF_acession}.noCHR.fasta").exists()

  script:
  """
  sed 's/chr//g' "${fasta_file}" > "${fasta_file.simpleName}.noCHR.fasta"
  """
}



process download_igvf_gtf {
  storeDir "references/"

  //when:
  //  !file("references/${params.gtf_IGVF_acession}.gtf.gz").exists()

  output:
    path "${params.gtf_IGVF_acession}.gtf.gz"

  script:
  """
  python ${params.download_igvf_portal} --acc ${params.gtf_IGVF_acession} 
  """
}



process remove_chr_gtf {
  storeDir "references/"

  input:
    path gtf_file

  output:
    path "${params.gtf_IGVF_acession}.noCHR.gtf"

  //when:
    //  !file("${params.topDir}/references/${params.gtf_IGVF_acession}.noCHR.gtf").exists()

  script:
  """
  sed 's/chr//g' "${gtf_file}" > "${gtf_file.simpleName}.noCHR.gtf"
  """
}


process vcf2vci{
  publishDir "references/", mode: 'copy'

    input:
    each strain                 // fan-out dimension
    path vcf_file               // singleton
    path fasta_file_nchr        // singleton

  
  output:
    path "${strain}.vci"

  script:
  """
    g2gtools vcf2vci -p 1 -i ${vcf_file} -o ${strain}.vci -s ${strain} -f ${fasta_file_nchr}

  """

}



workflow {
  // Fetch inputs
  vcf_ch        = curl_vcf()
  fasta_gz_ch   = download_igvf_fasta()
  fasta_nchr_ch = remove_chr_fasta(fasta_gz_ch)

  // Per-strain fan-out
  strains_ch = Channel.of(params.strains).flatten()

  vcf2vci(
    strains_ch,        // `each strain`
    vcf_ch,            // singleton VCF
    fasta_nchr_ch      // singleton FASTA (no chr)
  )

  // GTF side chain
  gtf_gz_ch   = download_igvf_gtf()
  gtf_nochr_ch= remove_chr_gtf(gtf_gz_ch)
}

