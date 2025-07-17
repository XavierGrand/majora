version = "0.8.2"
container_url = "xgrand/dorado:${version}"

process pod5convert{
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$fast5_folder"
  if (params.pod5convert_out != "") {
    publishDir "results/${params.pod5convert_out}", mode: 'copy'
  }

  input:
    path(fast5_folder)
  
  output:
    path("converted.pod5"), emit: pod5
  
  script:
"""
pod5 convert fast5 \
     --output converted.pod5 \
     --threads ${task.cpus} \
     --recursive \
     ${fast5_folder}
"""
}


process basecalling {
    container = "${container_url}"
    label "gpus"
    tag "$pod5"
    if (params.basecalling_out != "") {
    publishDir "results/${params.basecalling_out}", mode: 'copy'
    }
    input:
    path pod5

    output:
    path ("basecalled/basecalled.bam"), emit : basecalled
// Files are kept in bam because demultiplexing bam files is more has been more efficient than demultiplexing fastq.
    script:
    """
    dorado basecaller \
           --device ${params.cuda} \
           ${params.model} \
           ${pod5} \
           --kit-name ${params.kit_barcoding} \
           --no-trim \
           --output-dir basecalled/
    cd basecalled
    mv *.bam basecalled.bam
    """
}

process basecalling_nobc {
    container = "${container_url}"
    label "gpus"
    tag "$fast5_folder"
    if (params.basecalling_nobc_out != "") {
    publishDir "results/${params.basecalling_nobc_out}", mode: 'copy'
    }
    input:
    path fast5_folder

    output:
    path ("basecalled_nobc/*.fastq"), emit : basecalled_nobc // a vérifier si un ou plsueiurs fastq

    script:
    """
    dorado basecaller \
           --device ${params.cuda} \
           ${params.model} \
           --emit-fastq \
           ${fast5_folder} \
           --output-dir basecalled_nobc/
    """
}

process demux {
    container = "${container_url}"
    label "big_mem_multi_cpus"
    tag "$basecalled"
    if (params.demux_out != "") {
    publishDir "results/${params.demux_out}", mode: 'copy'
    }
    input:
    path basecalled

    output:
    path("demuxed/*.bam"), emit : demuxed
    path("demuxed/barcoding_summary.txt"), emit: barcoding_summary
    path("basecalled_summary.tsv"), emit: sequence_summary
    
    script:
    """
    #!/usr/bin/bash

  dorado summary ${basecalled} > basecalled_summary.tsv
  dorado demux \
         --emit-summary \
         --threads ${task.cpus} \
         --kit-name ${params.kit_barcoding} \
         --output-dir demuxed/ \
         ${basecalled}

  cd demuxed/
# Loop removes digits and "_" before "barcodenumber.bam"
  for file in *.bam
    do
      mv "\$file" "\${file##*_}"
    done
  
"""
}



