#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// these processes download and prepare reference files

process curl_vcf {
  storeDir "references/"

  output:
    path "${ params.mouse_vcf.split('/')[-1] }"

  script:
  """
    curl -O ${params.mouse_vcf}
  """
}


process download_igvf_fasta {
  storeDir "references/"

  output:
    path "${params.fasta_IGVF_acession}.fasta.gz"

  script:
  """
    python ${params.download_igvf_portal} --acc ${params.fasta_IGVF_acession} 
  """
}


process download_igvf_gtf {
  storeDir "references/"

  output:
    path "${params.gtf_IGVF_acession}.gtf.gz"

  script:
  """
    python ${params.download_igvf_portal} --acc ${params.gtf_IGVF_acession} 
  """
}

// these processes remove 'chr' from the fasta and gtf files

process remove_chr_fasta {
  storeDir "references/"

  input:
    path fasta_file

  output:
    path "${fasta_file.simpleName}.noCHR.fasta"

  script:
  """
    zcat "${fasta_file}" | sed 's/chr//g' > "${fasta_file.simpleName}.noCHR.fasta"
  """
}


process remove_chr_gtf {
  storeDir "references/"

  input:
    path gtf_file

  output:
    path "${gtf_file.simpleName}.noCHR.gtf"

  script:
  """
    zcat "${gtf_file}" | sed 's/chr//g'  > "${gtf_file.simpleName}.noCHR.gtf"
  """
}

// these processes create the strain 'personalized' fasta and gtf files
// using g2gtools

process vcf2vci{
  storeDir "references/${strain}/"

  input:
    each strain                 // fan-out dimension
    path vcf_file               // singleton
    path fasta_file_nchr        // singleton

  output:
    path "${strain}.vci.gz"
    path "${strain}.vci.gz.tbi"
    val strain

  script:
  """
    g2gtools vcf2vci -p 1 -i ${vcf_file} -o ${strain}.vci -s ${strain} -f ${fasta_file_nchr}

  """
}


process patch_fasta{
  storeDir "references/${strain}/"
    
  input:
    path fasta                  
    path vci
    path vci_index
    val strain                   
  
  output:
    path "mm39.${strain}.patch.fa"
    path vci
    path vci_index
    val strain

  script:
  """

    g2gtools patch -p 1 -i ${fasta} -c ${vci} -o mm39.${strain}.patch.fa

  """
}


process transform_fasta{
  storeDir "references/${strain}/"
    
  input:
    path patch_fasta                  
    path vci
    path vci_index
    val strain                   
  
  output:
    path "mm39.${strain}.unnamed.fa"
    val strain

  script:
  """
    g2gtools transform -p 1 -i ${patch_fasta} -c ${vci} -o mm39.${strain}.unnamed.fa
  """
}


process rename_fasta{
  storeDir "references/${strain}/"
    
  input:
    path unnamed_fasta
    val strain
  
  output:
    path "mm39.${strain}.fa"
    val strain

  script:
  """
    sed "s/^>/>${strain}_/" ${unnamed_fasta} > mm39.${strain}.fa
  """
}


process cat_fastas{
  storeDir "references/${strain}/"
    
  input:
    path fasta                  
    path renamed_fasta
    val strain

  output:
    path "mm39_B6J_${strain}_genome.fa"
    
  script:
  """
    cat ${fasta} ${renamed_fasta} > mm39_B6J_${strain}_genome.fa
  """
}


process convert_gtf{
  storeDir "references/${strain}/"
    
  input:
    path gtf                  
    path vci
    path vci_index
    val strain
  
  output:
    path "mm39.${strain}.unnamed.gtf"
    val strain
  
  script:
  """
    g2gtools convert -c ${vci} -i ${gtf} -o mm39.${strain}.unnamed.gtf
  """
}


process rename_gtf{
  storeDir "references/${strain}/"
    
  input:
    path unnamed_gtf
    val strain
  
  output:
    path "mm39.${strain}.gtf"
    val strain
  
  script:
  """
    sed -E \
      -e 's/ENSMUS/'"${strain}"'_ENSMUS/g' \
      -e 's/gene_name "([^"]+)"/gene_name "'"${strain}"'_\\1"/g' \
      -e 's/transcript_name "([^"]+)"/transcript_name "'"${strain}"'_\\1"/g' \
      -e '/^#/! s/^/'"${strain}"'_/' \
      "${unnamed_gtf}" > "mm39.${strain}.gtf"
  """
}


process cat_gtfs{
  storeDir "references/${strain}/"
    
  input: 
    path gtf                  
    path renamed_gtf
    val strain

  output:
    path "mm39_B6J_${strain}_gtf.gtf"
    val strain

  script:
  """
    cat ${gtf} mm39.${strain}.gtf > mm39_B6J_${strain}_gtf.gtf
  """
}

process kb_index{
  storeDir "references/${strain}/kb_index/"
    
  input: 
    path fasta                  
    path gtf
    val strain

  output:
    path index.idx
    path t2g.t2g
    path c1.c1
    path c2.c2
    path fasta.fa
    path output.na
    val strain

  script:
  """
    kb ref \
        --workflow=nac \
        -i index.idx \
        -g t2g.t2g \
        -c1 c1.c1 \
        -c2 c2.c2 \
        -f1 fasta.fa \
        -f2 output.na \
        --tmp ${strain}tmp \
        ${fasta} \
        ${gtf}
  """
}


workflow make_references {
  // Fetch inputs
  vcf_ch        = curl_vcf()
  fasta_gz_ch   = download_igvf_fasta()
  fasta_nchr_ch = remove_chr_fasta(fasta_gz_ch)
  
  // GTF 
  gtf_gz_ch   = download_igvf_gtf()
  gtf_nochr_ch= remove_chr_gtf(gtf_gz_ch)

  // Per-strain fan out
  strains_ch = Channel.of(params.strains).flatten()

  vci = vcf2vci(
    strains_ch,        
    vcf_ch,            
    fasta_nchr_ch      
  )

  fasta_patch = patch_fasta(
    fasta_nchr_ch,
    vci
  )
  fasta_transform = transform_fasta(fasta_patch)
  fasta_rename = rename_fasta(fasta_transform)
  fasta_cat = cat_fastas(fasta_nchr_ch, fasta_rename)

  gtf_convert = convert_gtf(
    gtf_nochr_ch,
    vci
  )
  gtf_rename = rename_gtf(gtf_convert)

  gtf_cat = cat_gtfs(gtf_nochr_ch, gtf_rename)

  emit:
    fasta_cat
    gtf_cat[0]
    gtf_cat[1]



}

workflow make_index {
  take:
    fasta_cat
    gtf_cat
    strain

  main:
    if (params.readType == 'RNA') {
      kb_index(
        fasta_cat,
        gtf_cat
      )
    } else {
      error "ATAC currently unsupported"
    }
}

workflow {
  make_references()
  make_index(make_references.out)
}
