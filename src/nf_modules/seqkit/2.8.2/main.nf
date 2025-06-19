version = "2.8.2"
container_url = "xgrand/seqkit:${version}"

process doublefastaref {
  container = "${container_url}"
  label "small_mem_mono_cpus"
  tag "${fasta.baseName}"
  if (params.doublefastaref_out != "") {
    publishDir "results/${params.doublefastaref_out}", mode: 'copy'
  }
  
  input:
    path(fasta)

  output:
    tuple val(params.hbvdb), path("doubled/${fasta.baseName}_doubled.fasta"), emit: doubledfasta
  // to concatenate sequences with same ID from multiple files     


  script:
    """
    #!/usr/bin/bash 
    mkdir doubled
    seqkit concat ${fasta} ${fasta} -o doubled/${fasta.baseName}_doubled.fasta

    """
}

process extractref {
  container = "${container_url}"
  label "small_mem_mono_cpus"
  tag "${bestref.baseName}"
  if (params.extractref_out != "") {
    publishDir "results/${params.extractref_out}", mode: 'copy'
  }
  
  input:
    tuple val(barcode), path(bestref)
    tuple val(genotype), path(multifasta)

  output:
    tuple val(barcode), path("${barcode}/${barcode}_toMap.fasta"), emit: referenceseq

  script:
    """
    echo "bestref: ${bestref}"
    echo "multifasta: ${multifasta}"
    mkdir ${barcode}
    seqkit grep -r -f ${bestref} ${multifasta} -o ${barcode}/${barcode}_toMap.fasta
    #adding -r to enable partly matching
    """
}

process grep_primer {
  container = "${container_url}"
  label "small_mem_multi_cpus"
  tag "${barcode}"
  if (params.grep_primer_out != "") {
    publishDir "results/${params.grep_primer_out}", mode: 'copy'
  }
  
  input:
    tuple val(barcode), path(ready_fastq)
    val(primer)
    val(primerloc)
  output:
    tuple val(barcode), path("${barcode}/${barcode}_${primerloc}_filtered.fastq.gz"), emit: filtered_fastq
    tuple val(barcode), path("${barcode}/*.tsv"), emit: filtered_tsv
    path("${barcode}/*.txt")
  script:
  length = Math.round(primer.size().div(10)) // allowed mismatch size = 1/10 bp length of primer size
    """
    mkdir ${barcode}
    cd ${barcode}/
    echo "mismatch allowed to ${primerloc} primer search: ${length}" > ${barcode}_${primerloc}_mismatch.txt
    echo ${primer} >> primer.txt
    #Find sequences that respect authorized mismatch length
    seqkit grep -i -f primer.txt -m ${length} ../${ready_fastq} -o ${barcode}_${primerloc}_filtered.fastq -j ${task.cpus}
   

    gzip ${barcode}_${primerloc}_filtered.fastq

    # Gets file extension, data type, number of sequences, number of bases/residues, length of the longest and shortest sequence
    seqkit stats ../${ready_fastq} -T -j ${task.cpus} >> ${barcode}_seq_stats.tsv
    seqkit stats ${barcode}_${primerloc}_filtered.fastq.gz -T -j ${task.cpus} | tail -n1 >> ${barcode}_seq_stats.tsv

    rm primer.txt
    """
}


process filterbylength {
  container = "${container_url}"
  label "small_mem_multi_cpus"
  tag "${barcode}"
  if (params.filterbylength_out != "") {
    publishDir "results/${params.filterbylength_out}", mode: 'copy'
  }
  
  input:
    tuple val(barcode), path(fastq)
    tuple val(barcode), path(tsv)

  output:
    tuple val(barcode), path("${barcode}/${barcode}_length_filtered.fastq"), emit: length_filtered_fastq
    path("${barcode}/*.txt")
    tuple val(barcode), path("${barcode}/${barcode}_seq_stats.tsv"), emit: length_filtered_tsv

  script:
    """
    mkdir ${barcode}
    cd ${barcode}/
    minlength=\$(awk -F '\t' 'NR==3{A=\$7} END {res=int((A*0.9)); print res}' ../${tsv})
    maxlength=\$(awk -F '\t' 'NR==3{B=\$7} END {res=int((B*1.1)); print res}' ../${tsv})
    echo "Calculated minimal length is: \$minlength bp and calculated maximal length is: \$maxlength bp" > length.txt
    
    seqkit seq --min-len \$minlength --max-len \$maxlength --remove-gaps ../${fastq} -j ${task.cpus} > ${barcode}_length_filtered.fastq        
    seqkit stats ${barcode}_length_filtered.fastq -T -j ${task.cpus} >> ${barcode}_seq_stats.tsv

    """
}

process read_stats {
  container = "${container_url}"
  label "small_mem_multi_cpus"
  tag "${barcode}"
  if (params.read_stats_out != "") {
    publishDir "results/${params.read_stats_out}", mode: 'copy'
  }
  
  input:
    tuple val(barcode), path(fastq)

  output:
    tuple val(barcode), path("${barcode}/${barcode}_seq_stats.tsv")

  script:
    """
    mkdir ${barcode}
    cd ${barcode}/  
    seqkit stats ../${fastq} -T -j ${task.cpus} >> ${barcode}_seq_stats.tsv

    """
}