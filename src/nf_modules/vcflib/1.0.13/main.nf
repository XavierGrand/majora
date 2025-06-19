version = "1.0"
container_url = "lousahra/vcflib:${version}"

process vcf_filter {
  label "small_mem_multi_cpus"
  container = "${container_url}"
  tag "$barcode"
  if (params.vcf_filter_out!= "") {
    publishDir "results/${params.vcf_filter_out}", mode: 'copy'
  }
  
  input:
    tuple val(barcode), path(vcf)
  output:
    tuple val(barcode), path("${barcode}_filtered.vcf"), optional: true, emit: vcf_filtered
    tuple val(barcode), path("${vcf}")
  script:
    """
    vcffilter  -f "( DP > 19 & QUAL > 19 )" ${vcf} > ${barcode}_filtered.vcf
    """
}


