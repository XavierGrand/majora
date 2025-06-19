version = "1.0"
container_url = "lousahra/index_vcf:${version}"

process reindex_vcf {
  container = "${container_url}"
  label "small_mem_multi_cpus"
  tag "$barcode"
  if (params.reindex_vcf_out != "") {
    publishDir "results/${params.reindex_vcf_out}", mode: 'copy'
  }
  
  input:
    tuple val(barcode), path(fasta), path(id), path(vcf)
    val(fwd_primer)
    val(rev_primer)
  output:
    tuple val(barcode), path("${barcode}/*")
    tuple val(barcode), path("${barcode}/final_vcf.tsv"), emit : tsv
    tuple val(barcode), path("${barcode}/positioned_ref.fasta"), emit : ref_fasta

  """
  bash /app/reindex_vcf.sh ${fasta} ${id} ${vcf} ${fwd_primer} ${rev_primer}
  mkdir ${barcode}
  mv *.fasta *.txt *.vcf *.tsv *.bed ${barcode}/
  """
}