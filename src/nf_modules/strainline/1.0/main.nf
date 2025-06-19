version = "1.0"
container_url = "lousahra/strainline:${version}"

process strainline {
  errorStrategy 'ignore'
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$barcode"
  if (params.strainline_out != "") {
    publishDir "results/${params.strainline_out}", mode: 'copy'
  }

  input:
  tuple val(barcode), path(ref)
  output:
  tuple val(barcode), path("${barcode}_strainline_out/"), optional: true
  script:
"""
/app/src/strainline.sh -i ${ref} -o ${barcode}_strainline_out -p ont --minSeedLen 2000 -t 32
"""
}
