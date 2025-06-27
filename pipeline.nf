#!/usr/bin/env nextflow
nextflow.enable.dsl=2


samples_csv_ch = Channel.fromPath( params.file_csv )



process kb_count {
  tag        { "${strain}-${subpool}" }
  publishDir { "kb_count/${plate}/${strain}" }, mode: 'copy'

  input:
    tuple val(strain), val(subpool), val(plate),
          val(pairs)


  output:
    path "${subpool}", emit: subpool_dir

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
  tag { subpool_dir.getName() }

  input:
    path subpool_dir

  output:
    path "${subpool_dir.getName()}.h5ad", emit: h5ad

  script:
  """
  python ${params.make_adata_script} \
    --directory ${subpool_dir} \
    --output ${subpool_dir.getName()}.h5ad
  """
}


process cellbender {
  tag { total_adata.getName() }

  input:
    path(total_adata)

  output:
    path "${total_adata.getName()}_cellbender.h5", emit: cb_h5

  script:
  """
    cellbender remove-background \
        --input ${total_adata} \
        --output ${total_adata}_cellbender.h5 \
        --total-droplets-included 200000 \
        --learning-rate 0.00005 \
        --expected-cells 67000 \
        --cuda
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
  subpool_dirs = kb_count(config_ch).subpool_dir

  //
  // 3) downstream as before
  //
  total_adata = make_adata(subpool_dirs)
  cellbender(total_adata)
}





