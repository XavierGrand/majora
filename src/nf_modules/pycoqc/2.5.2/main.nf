version = "1.0"
container_url = "lousahra/pycoqc:${version}"
// Version for numpy : 1.24.14 last version using numpy.uint64 before switching to numpy.uint32 accordind to the following stackoverflow post :
// https://bioinformatics.stackexchange.com/questions/23229/pycoqc-n50-is-zero-or-null-overflow-encountered-in-scalar-add
// Newer versions (released afer 1.26.14) cause a "overflow encountered in scalar add" error that stops the process


process qualitycontrol {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$barcode"
  if (params.qualitycontrol_out != "") {
    publishDir "results/${params.qualitycontrol_out}", mode: 'copy'
  }

  input:
  tuple val(barcode), path(bam), path(bai), path(tsv)


  output:
  path("${barcode}/${barcode}_QC.html")

  """
  mkdir ${barcode}
  pycoQC -f ${tsv} -a ${bam} -o ${barcode}_QC.html
  mv ${barcode}_QC.html ${barcode}
  """
}
