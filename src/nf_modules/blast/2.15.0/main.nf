version = "2.15.0"
container_url = "ncbi/blast:${version}"

//params.hbvdb = ""
process dl_hbvdb {
  container = "${container_url}"
  label "small_mem_mono_cpus"
  tag "${params.hbvdb}_Genomes"
  if (params.dl_hbvdb_out != "") {
    publishDir "results/${params.dl_hbvdb_out}", mode: 'copy'
  }

  input:
    
  output:
    //tuple val(params.hbvdb), path("*.fasta"), emit: reference_db
    path("*.fasta"), emit: reference_db

  script:
  if (params.hbvdb == "all" || params.hbvdb == "A" || params.hbvdb == "B" || params.hbvdb == "C" || 
      params.hbvdb == "D" || params.hbvdb == "E" || params.hbvdb == "F" || params.hbvdb == "G" || 
      params.hbvdb == "H") {
    link = "https://hbvdb.lyon.inserm.fr/data/nucleic/fasta/${params.hbvdb}_Genomes.fas"
    output_name = "${params.hbvdb}_Genomes.fasta"
  }
  /*
  else if(params.hbvdb == "hdvdb"){
    link = 
    output_name = "${params.hbvdb}_Genomes.fasta"
  }
  */
"""
wget --quiet --no-check-certificate -O ${output_name} ${link}
"""
}


//params.makeblastdb_out = ""
process makeblastdb {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${file_id}_Genomes"
  if (params.makeblastdb_out != "") {
    publishDir "results/${params.makeblastdb_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(ref_fasta)

  output:
    tuple val(ref_fasta.baseName), path("*.fasta.n*"), emit: blastdb

  script:

"""
makeblastdb \
    -in ${ref_fasta} \
    -input_type 'fasta' \
    -dbtype 'nucl' \
    -parse_seqids
"""
}

params.blast_them_all_out = ""
process blast_them_all {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${barcode}_Genomes"
  if (params.blast_them_all_out != "") {
    publishDir "results/${params.blast_them_all_out}", mode: 'copy'
  }

  input:
    tuple val(barcode), path(fastq)
    tuple val(genotype), path(blastdb)

  output:
    tuple val(barcode), path("${barcode}/${barcode}_hits.txt"), emit: blasthits
    tuple val(barcode), path("${barcode}/${barcode}_hits_counts.txt"), emit: counthits
    tuple val(barcode), path("${barcode}/${barcode}_best_ref.txt"), emit: bestref

  script:
"""mkdir ${barcode}
blastn -db ${genotype}.fasta -query ${fastq} \
       -task megablast \
       -max_target_seqs 1 \
       -max_hsps 1 \
	   -outfmt "6 qseqid sseqid evalue bitscore slen qlen length pident" \
	   -out ${barcode}/${barcode}_hits.txt -num_threads ${params.blasthreads}

cut -f2 ${barcode}/${barcode}_hits.txt | sort | uniq -c | sort -k 1,1 -r | head -n1 | awk -F'|' '{print \$3}' | grep -o '^[^_ ]*' > ${barcode}/${barcode}_best_ref.txt
# grep -o pour extraire tout ce qu'il y a avant le premier _ ou le premier espace

cut -f2 ${barcode}/${barcode}_hits.txt | sort | uniq -c | sort -k 1,1 -r > ${barcode}/${barcode}_hits_counts.txt

"""
}


process makerefdb {
  label "big_mem_multi_cpus"
  // mettre en tag l'ID de la séquence_Genomes
  if (params.makerefdb_out != "") {
    publishDir "results/${params.makerefdb_out}", mode: 'copy'
  }

  input:
  path(ref_file)

  output: 
  path("combined_ref.fasta"), emit: ref_db
  path("*.txt")
// à tester sans le shebang
  script:
"""
  #!/usr/bin/bash 
  
  sed '1d' ${ref_file}  > temp1.txt
  awk -F'\t' '{print "https://hbvdb.lyon.inserm.fr/tmp/hbvdb_dat/"\$2"/"\$2"_sequence.txt"}'  temp1.txt > temp2.txt
  mkdir refdb
  wget -P refdb  -i temp2.txt
  touch combined_ref.fasta

  for fichier in refdb/*_sequence.txt
  do
    cat "\$fichier" >> combined_ref.fasta
    echo "" >> combined_ref.fasta
  done
  sed -i 's/>/>gnl|hbvnuc|/' combined_ref.fasta
 #rm temp*.txt
"""
}
