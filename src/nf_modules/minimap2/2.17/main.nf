version = "2.17"
container_url = "lbmc/minimap2:${version}"

params.mapping_hbv_genome = "-ax map-ont"
process mapping_hbv_genome {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${barcode}"
  if (params.mapping_hbv_genome_out != "") {
    publishDir "results/${params.mapping_hbv_genome_out}", mode: 'copy'
  }

  input:
  tuple val(barcode), path(fastq), path(genome)

  output:
  tuple val(barcode), path("${barcode}_res.bam"), emit: bam
  tuple val(barcode), path("${barcode}_unmapped_res.bam")

  script:
  memory = "${task.memory}" - ~/\s*GB/
  memory = memory.toInteger() / (task.cpus + 1.0)
  """
  minimap2 ${params.mapping_hbv_genome} -t ${task.cpus} -K ${memory} ${genome} ${fastq} | tee >(samtools view -Shb -F4 - > ${barcode}_res.bam) >(samtools view -Shb -f4 - > ${barcode}_unmapped_res.bam)
   
  
  """
  // Found in doc: -S
  //Ignored for compatibility with previous samtools versions. Previously this option was required if input was in SAM format, but now the correct format is automatically detected by examining the first few characters of input. 

  // -F4
  // filter the alignments that will be included in the output to only those alignments that match certain criteria
  // Do not output alignments with any bits set in FLAG present in the FLAG field
}

