version = "1.0"
container_url = "lousahra/ultimate_filter:${version}"

process ultimate_read_filter {
  label "small_mem_multi_cpus"
  container = "${container_url}"
  tag "${barcode}"
  if (params.ultimate_read_filter_out != "") {
    publishDir "results/${params.ultimate_read_filter_out}", mode: 'copy'
  }
  
  input:
    tuple val(barcode), path(tsv)
    tuple val(barcode), path(fastq)
  output:
    tuple val(barcode), path("${barcode}_length_filtered.fastq"), optional: true, emit: ultimately_filtered_fastq
    path("summary.txt")
  script:
    """
    python /app/ultimate_filter.py ${tsv} ${fastq} ${barcode}
    """
}

//         system("gzip $fastq > ${barcode}_length_filtered.fastq.gz");             system('rm ' $fastq); 
/*

    awk 'NR==2 { 
        if (\$4<20) { 
            print "${barcode} contains less than 20 reads, it will be suppressed"; 
            print "${fastq}" > ${barcode}_to_delete.txt
        } else {
            print "${barcode} contains 20 or more reads, process continues"
        }
    }' $tsv >> summary.txt

    if [ -f ${barcode}_to_delete.txt ]; then
        while read fastq_file; do
            echo "Deleting \$fastq_file"
            rm \$fastq_file
        done < ${barcode}_to_delete.txt
    else
        echo "There is no need to delete files for $barcode"
    fi


    /app/ultimate_filter.py $tsv $fastq $barcode
    */