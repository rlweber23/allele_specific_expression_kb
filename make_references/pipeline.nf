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
    tuple path("mm39_B6J_${strain}_genome.fa"), val(strain)
    
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
    tuple path("mm39_B6J_${strain}_gtf.gtf"), val(strain)

  script:
  """
    cat ${gtf} mm39.${strain}.gtf > mm39_B6J_${strain}_gtf.gtf
  """
}

process kb_index{
  storeDir "references/${strain}/kb_index/"
    
  input: 
    tuple val(strain), path(fasta), path(gtf)

  output:
    tuple val(strain), path("index.idx"), path("t2g.t2g"), path("c1.c1"), path("c2.c2"), path("fasta.fa"), path("output.na"), emit: kb_index_files
    // path "index.idx"
    // path "t2g.t2g"
    // path "c1.c1"
    // path "c2.c2"
    // path "fasta.fa"
    // path "output.na"
    // val strain

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
        --tmp "${strain}tmp" \
        ${fasta} \
        ${gtf}
  """
}

process chromap_index {
  storeDir "references/${strain}/chromap_index/"

  input:
    tuple val(strain), path(fasta)

  output:
    tuple val(strain), path("chromap_index"), emit: chromap_index_files

  script:
  """
    chromap -i  -r ${fasta} -o ${strain}_chromap_index
  """
}


samples_csv_ch = Channel.fromPath( params.file_csv )


process kb_count {
  tag        { "${strain}-${subpool}" }
  publishDir { "kb_count/${plate}/${strain}" }, mode: 'copy'

  input:
    // tuple val(strain), val(subpool), val(plate), val(pairs)
    tuple val(strain), val(subpool), val(plate), val(pairs), path("index.idx"), path("t2g.t2g"), path("c1.c1"), path("c2.c2"), path("fasta.fa"), path("output.na")


  output:
    tuple val(strain), val(subpool), val(plate), path("${subpool}"), emit: subpool_dir

  script:
  """
    set -euo pipefail

    # make subdirs (nova1, nova2, …) inside the work dir
    for fq in ${pairs.join(' ')}; do
      subdir=\$(basename \$(dirname \"\$fq\"))
      mkdir -p \"\$subdir\"
      ln -s \"\$fq\" \"\$subdir/\$(basename \"\$fq\")\"
    done

    #  build the exact list of links in the order we got them
    files_to_use=""
    for fq in ${pairs.join(' ')}; do
      subdir=\$(basename \$(dirname \"\$fq\"))
      files_to_use="\$files_to_use \$subdir/\$(basename \"\$fq\")"
    done


    # finally call kb count on those links
    kb count \
      --h5ad \
      --gene-names \
      --sum=total \
      --strand=forward \
      -r ${params.kb_replacement_bcs} \
      -w ${params.kb_onlist} \
      --workflow=nac \
      -g t2g.t2g \
      -x SPLIT-SEQ \
      -i index.idx \
      -t 4 \
      -c1 c1.c1 \
      -c2 c2.c2 \
      -o ${subpool} \
      \$files_to_use
  """
}



process make_adata {

  input:
    tuple val(strain), val(subpool), val(plate), path(subpool_dir)

  output:
    tuple val(strain), val(subpool), val(plate), path("${subpool_dir.getName()}.h5ad"), emit: h5ad

  script:
  """
  python ${params.make_adata_script} \
    --directory ${subpool_dir} \
    --output ${subpool_dir.getName()}.h5ad
  """
}


process cellbender {

  errorStrategy { 
    task.attempt <= 6 ? 'retry'  : 'ignore' 
  }

  publishDir { "cellbender/${plate}/${strain}/" }, mode: 'copy'

  input:
    tuple val(strain), val(subpool), val(plate), path(total_adata)

  output:
    tuple val(strain), val(subpool), val(plate),
      path("${total_adata.getName()}_cellbender_filtered.h5"),
      path("${total_adata.getName()}_cellbender_report.html"),
      emit: cb_outputs


  script:
  """
    cellbender remove-background \
        --input ${total_adata} \
        --output ${total_adata}_cellbender.h5 \
        --total-droplets-included 200000 \
        --learning-rate 0.000025 \
        --expected-cells 67000 \
        --cuda
  """
}


process cb_h5_to_h5ad {
  publishDir { "cellbender/${plate}/${strain}/" }, mode: 'copy'

  input:
    tuple val(strain), val(subpool), val(plate), path(cb_h5), path(_)


  output:
    tuple val(strain), val(subpool), val(plate),
      path("${cb_h5.getName().replaceFirst(/\.h5$/, '.h5ad')}"), emit: cb_h5ad

  script:
  """
  python ${params.cb_h5_to_h5ad_script} \
    --input ${cb_h5} \
    --output ${cb_h5.getName().replaceFirst(/\.h5$/, '.h5ad')}
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

  ref_channel = fasta_cat.join(gtf_cat, by: [1])//.view()

  emit:
    ref_channel
}

workflow make_index {
  take:
    ref_channel

  main:
    if (params.readType == 'RNA') {
      index_ch = kb_index(ref_channel)
    } else if (params.readType == 'ATAC') {
      index_ch = chromap_index(ref_channel)
      }

  emit:
    index_ch
}

workflow RNA_ASE{
  take:
  index_ch

  main:
    config_ch = Channel
      .fromPath(params.file_csv)
      .splitCsv(header:true)
      .collect()
      .flatMap { rowsList ->
        def groups = rowsList.groupBy { [ it.strain, it.subpool, it.plate ] }
        groups.collect { key, groupRows ->
          def (strain, subpool, plate) = key
          groupRows = groupRows.sort { it.lane }   // simple lex‐sort on "L001","L002"
          def pairs = groupRows.collectMany { row ->
            [ file(row.fastq_r1), file(row.fastq_r2) ]
          }
          tuple(strain, subpool, plate, pairs)
        }
      }

    config_ch//.view()    
    //
    // 2) one kb_count run per tuple, with real file paths
    //

    kb_count_channel = config_ch.combine(index_ch, by: [0]).view()


    subpool_dirs = kb_count(kb_count_channel)

    //
    // 3) downstream as before
    //
    total_adata = make_adata(subpool_dirs)
    cellbender_outputs = cellbender(total_adata)
    cb_h5_to_h5ad(cellbender_outputs.cb_outputs)

}


workflow {
  make_references()
  make_index(make_references.out)

  if (params.readType == 'RNA') {
    RNA_ASE(make_index.out)
    } 
  
}
