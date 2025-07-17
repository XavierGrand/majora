version = "2.25.0"
container_url = "lbmc/bedtools:${version}"

// Demultiplexing in done using bam files to maintain barcoding efficiency, thus bam demultiplexed files are converted into fastq files before filtration.
process convertfastq {
  container = "${container_url}"
  label "small_mem_multi_cpus"
  tag "demuxed.bam"
  if (params.convert_fastq != "") {
    publishDir "results/${params.convert_fastq}", mode: 'copy'
  }
  
  input:
    path(ready_fastq)
  output:
    path("converted/*.fastq"), emit : converted_fastq
  
  """
  mkdir converted
  for file in *.bam
    do
      baseName=\$(basename "\$file" .bam)
      bedtools bamtofastq -i \$file -fq converted/\$baseName.fastq
    done
  """
}