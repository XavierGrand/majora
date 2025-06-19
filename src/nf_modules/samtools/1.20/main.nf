version = "1.20"
container_url = "xgrand/samtools:${version}"

process index_fasta {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$barcode"
  if (params.index_fasta_out != "") {
    publishDir "results/${params.index_fasta_out}", mode: 'copy'
  }

  input:
    tuple val(barcode), path(fasta)
  output:
    tuple val(barcode), path("*.fai"), emit: index
    tuple val(barcode), path("*.fasta"), emit: toMap

  script:
"""
if gzip -t ${fasta}; then
  zcat ${fasta} > ${fasta.simpleName}.fasta
  samtools faidx ${params.index_fasta} ${fasta.simpleName}.fasta
  cp ${fasta.simpleName}.fasta ${barcode}.fasta
else
  samtools faidx ${params.index_fasta} ${fasta}
  cp ${fasta} ${barcode}.fasta
fi
"""
}



params.filter_bam_mapped = "-F 4"
process filter_bam_mapped {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$barcode"
  if (params.filter_bam_mapped_out != "") {
    publishDir "results/${params.filter_bam_mapped_out}", mode: 'copy'
  }

  input:
    tuple val(barcode), path(bam)

  output:
    tuple val(barcode), path("*_mapped.bam"), emit: bam
  script:
"""
samtools view -@ ${task.cpus} ${params.filter_bam_mapped} -hb ${bam} > \
  ${bam.simpleName}_mapped.bam
"""
}


process sort_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$barcode"
  if (params.sort_bam_out != "") {
    publishDir "results/${params.sort_bam_out}", mode: 'copy'
  }

  input:
    tuple val(barcode), path(bam)

  output:
    tuple val(barcode), path("*.bam*"), emit: bam

  script:
"""
samtools sort -@ ${task.cpus} ${params.sort_bam} -O BAM -o ${bam.simpleName}_sorted.bam ${bam}
"""
}


process index_bam {
  container = "${container_url}"
  label "small_mem_mono_cpus"
  tag "$barcode"
  if (params.index_bam_out != "") {
    publishDir "results/${params.index_bam_out}", mode: 'copy'
  }

  input:
    tuple val(barcode), path(bam)

  output:
    tuple val(barcode), path("${bam}"), path("*.bam.bai"), emit: bam_idx

  script:
"""
samtools index ${params.index_bam} ${bam}
"""
}

process consensus {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$barcode"
  if (params.consensus_out != "") {
    publishDir "results/${params.consensus_out}", mode: 'copy'
  }

  input:
    tuple val(barcode), path(bam)

  output:
    tuple val(barcode), path("*_consensus.fasta"), emit: fasta

  script:
"""
samtools consensus -@ ${task.cpus} -f "fasta" ${bam} -o ${barcode}_consensus.fasta
"""
}

process bam2sam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$barcode"
  if (params.bam2sam_out != "") {
    publishDir "results/${params.bam2sam_out}", mode: 'copy'
  }

  input:
    tuple val(barcode), path(bam)
    tuple val(barcode), path(idx)
  output:
    tuple val(barcode), path("${barcode}_calls_to_ref.sam"), emit: sam

  script:
"""
samtools view -h ${bam} > ${barcode}_calls_to_ref.sam
"""
}

process bam2fasta {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$barcode"
  if (params.bam2fasta_out != "") {
    publishDir "results/${params.bam2fasta_out}", mode: 'copy'
  }

  input:
    tuple val(barcode), path(bam)
  output:
    tuple val(barcode), path("${barcode}_calls_to_ref.fasta"), emit: fasta

  script:
"""
samtools fasta ${bam} > ${barcode}_calls_to_ref.fasta
"""
}