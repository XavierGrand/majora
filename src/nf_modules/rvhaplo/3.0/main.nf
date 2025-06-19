version = "1.0"
container_url = "lousahra/rvhaplo:${version}"

process rv_haplo {
  //errorStrategy 'ignore'
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$barcode"
  if (params.rv_haplo_out != "") {
    publishDir "results/${params.rv_haplo_out}", mode: 'copy'
  }

  input:
  tuple val(barcode), path(sam), path(ref)
  output:
  tuple val(barcode), path("${barcode}_rvhaplo_out/"), optional: true
  script:
"""
/app/rvhaplo.sh -i ${sam} -r ${ref} -p ${barcode} -o ${barcode}_rvhaplo_out
"""
}