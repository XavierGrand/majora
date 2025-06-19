version = "1.0"
container_url = "lousahra/rlollipop:${version}"

process lollipop_vcf {
  label "small_mem_multi_cpus"
  container = "${container_url}"
  tag "${barcode}"
  if (params.lollipop_vcf_out != "") {
    publishDir "results/${params.lollipop_vcf_out}", mode: 'copy'
  }
  
  input:
    tuple val(barcode), path(tsv)
  output:
    tuple val(barcode), path("${barcode}/*.pdf"), optional: true
    tuple val(barcode), path("${barcode}/*.html"), optional: true
    tuple val(barcode), path("${barcode}/empty.txt"), optional: true
  script:
    """
  mkdir ${barcode}
 if [ ! -f ${tsv} ] || [ ! -s ${tsv} ]; then
    touch empty.txt
    echo "The vcf file from ${barcode} was empty" > empty.txt
    mv empty.txt ${barcode}
else
    Rscript /circular_lollipop.R "${tsv}"
    mv *.pdf ${barcode}/
    mv *.html ${barcode}/
fi
    """
}
