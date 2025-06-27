#!/usr/bin/env nextflow
nextflow.enable.dsl=2


samples_csv_ch = Channel.fromPath( params.file_csv )



process kb_count {
  tag        { "${strain}-${subpool}" }
  publishDir { "kb_count/${plate}/${strain}/" }, mode: 'copy'

  input:
    tuple \
      val(strain), \
      val(subpool), \
      val(plate), \
      path(pairs), \
      stageAs: { f -> "${f.parent.name}/${f.name}" }




  output:
    path "${subpool}", emit: subpool_dir


  script:
  """

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
    ${fastq_r1.join(' ')} ${fastq_r2.join(' ')}

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
  config_ch = Channel
    .fromPath(params.file_csv)
    .splitCsv(header:true)
    .map { row ->
      // first element is the grouping key
      tuple(
        [ row.strain, row.subpool, row.plate ],
        file(row.fastq_r1),
        file(row.fastq_r2)
      )
    }
    .groupTuple()
    .map { key, fastq_r1, fastq_r2 ->
      def (strain, subpool, plate) = key
      tuple(strain, subpool, fastq_r1, fastq_r2, plate)
    }
    config_ch.view()
    

    subpool_dirs = kb_count(config_ch).subpool_dir

    total_adata = make_adata(subpool_dirs)

    cellbender(total_adata)

}












