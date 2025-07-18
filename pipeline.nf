#!/usr/bin/env nextflow
nextflow.enable.dsl=2


samples_csv_ch = Channel.fromPath( params.file_csv )



process kb_count {
  tag        { "${strain}-${subpool}" }
  publishDir { "kb_count/${plate}/${strain}" }, mode: 'copy'

  input:
    tuple val(strain), val(subpool), val(plate), val(pairs)


  output:
    tuple val(strain), val(subpool), val(plate), path("${subpool}"), emit: subpool_dir

  script:
  """
  # make subdirs (nova1, nova2, …) inside the work dir
  for fq in ${pairs.join(' ')}; do
    subdir=\$(basename \$(dirname \"\$fq\"))
    mkdir -p \"\$subdir\"
    ln -s \"\$fq\" \"\$subdir/\$(basename \"\$fq\")\"
  done

  # now build the exact list of links in the order we got them
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
    -g ${params.kb_index}/${strain}/t2g.t2g \
    -x SPLIT-SEQ \
    -i ${params.kb_index}/${strain}/index.idx \
    -t 4 \
    -c1 ${params.kb_index}/${strain}/c1.c1 \
    -c2 ${params.kb_index}/${strain}/c2.c2 \
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
    task.attempt <= 4 ? 'retry'  : 'ignore' 
  }

  tag { total_adata[3].getName() }
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
  tag { cb_h5.getName() }
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


workflow {
  //
  // 1) read & group in pure Groovy
  //
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

  config_ch.view()    // should now print exactly one 4‐element tuple

  //
  // 2) one kb_count run per tuple, with real file paths
  //
  subpool_dirs = kb_count(config_ch)

  //
  // 3) downstream as before
  //
  total_adata = make_adata(subpool_dirs)
  cellbender_outputs = cellbender(total_adata)
  cb_h5_to_h5ad(cellbender_outputs.cb_outputs)

}





