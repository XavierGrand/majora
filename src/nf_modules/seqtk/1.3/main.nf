version = "1.3"
container_url = "xgrand/seqtk:${version}"

//params.sample_fastq_out = ""
params.reads_number = 100
process sample_fastq {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${barcode}"
  if (params.sample_fastq_out != "") {
    publishDir "results/${params.sample_fastq_out}", mode: 'copy'
  }

  input:
    tuple val(barcode), path(fastq)

  output:
    tuple val(barcode), path("${barcode}/${barcode}_sampled.fasta"), optional: true, emit: sampled_fastq
    path("${barcode}/*.txt"), optional: true

  script:
    """
    mkdir ${barcode}
    if ! zgrep -q "^@" ${fastq}
      then
        echo "The file ${barcode}_filtered.fastq.gz does not contain any sequence, it will not be emitted as a fastq or a fasta file" > ${barcode}/summary.txt
      else
        seq_count=\$(zgrep -c "^@" ${fastq})

        if [ \${seq_count} -gt 100 ]
          then
            seqtk sample -s100 ${fastq} ${params.reads_number} > ${barcode}/${barcode}_sampled.fastq
        else
            echo "The file ${barcode}_filtered.fastq.gz 100 or less sequences, all of them will be taken to create the file ${barcode}_sampled.fasta" > ${barcode}/summary.txt
            cp ${fastq} ${barcode}/${barcode}_sampled.fastq
        fi

        seqtk seq -a ${barcode}/${barcode}_sampled.fastq > ${barcode}/${barcode}_sampled.fasta
    fi
    """
}
