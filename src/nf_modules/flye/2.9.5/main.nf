version = "1.0"
container_url = "lousahra/flye:${version}"

process de_novo_assembly {
  errorStrategy 'ignore'
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$barcode"
  if (params.de_novo_assembly_out != "") {
    publishDir "results/${params.de_novo_assembly_out}", mode: 'copy'
  }

  input:
  tuple val(barcode), path(fastq)
  output:
  tuple val(barcode), path("${barcode}/${barcode}_de_novo_assembly/${barcode}_assembly.fasta"), optional: true, emit: assembly_fasta
  tuple val(barcode), path("${barcode}/${barcode}_de_novo_assembly"), optional: true
  path("${barcode}/${barcode}_summary.txt")
  path("${barcode}/${barcode}.vcf"), optional: true
  script:
"""
mkdir ${barcode}
num_seq=\$(grep -c "@" ${fastq})
echo "\$num_seq"
if (( \$num_seq < 20 )); then
  echo "The fastq file for ${barcode} contains less than 20 reads, de novo assembly cannot be processed" > ${barcode}/${barcode}_summary.txt
else
  flye --nano-hq ${fastq} -g 3000 --asm-coverage 40 --keep-haplotypes -o ${barcode}/${barcode}_de_novo_assembly
  mv ${barcode}/${barcode}_de_novo_assembly/assembly.fasta ${barcode}/${barcode}_de_novo_assembly/${barcode}_assembly.fasta
  echo "The fastq file for ${barcode} contains 20 or more reads, de novo assembly has been processed" > ${barcode}/${barcode}_summary.txt
fi
"""
}