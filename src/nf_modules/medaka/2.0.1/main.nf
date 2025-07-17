version = "1.0"
container_url = "lousahra/medaka:${version}"

process call_hap_variant {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$barcode"
  if (params.call_hap_variant_out != "") {
    publishDir "results/${params.call_hap_variant_out}", mode: 'copy'
  }
  input:
    tuple val(barcode), path(ref), path(fastq)
  output:
    tuple val(barcode), path("${barcode}/${barcode}_medaka.annotated.vcf"), emit: vcf
    tuple val(barcode), path("${barcode}/${barcode}_calls_to_ref.bam"), emit: calls
    tuple val(barcode), path("${barcode}/${barcode}_calls_to_ref.bam.bai"), emit: calls_idx
    path("${barcode}/${barcode}_inputs.txt")
    path("${barcode}/${barcode}_medaka_log.txt")
  script:
"""
  echo "The reference file is $ref and the bam file is $fastq" > ${barcode}_inputs.txt
  medaka_variant -i $fastq -r $ref -o ${barcode} > ${barcode}_medaka_log.txt 2>&1
  cd ${barcode}
  mv medaka.annotated.vcf ${barcode}_medaka.annotated.vcf
  mv calls_to_ref.bam ${barcode}_calls_to_ref.bam
  mv calls_to_ref.bam.bai ${barcode}_calls_to_ref.bam.bai
  mv ../${barcode}_inputs.txt .
  mv ../${barcode}_medaka_log.txt .
"""
}

// Process needed to add elements into INFO field of VCF file. These informations are needed to calculate variations frequencies.
process annotate_variant {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$barcode"
  if (params.annotate_variant_out != "") {
    publishDir "results/${params.annotate_variant_out}", mode: 'copy'
  }

  input:
    tuple val(barcode), path(vcf), path(ref)
    tuple val(barcode), path(calls)
    tuple val(barcode), path(calls_idx)
  output:
    tuple val(barcode), path("${barcode}_calling_supp.vcf"), emit: annotated_vcf
    tuple val(barcode), path("${barcode}_annotate.txt")
  script:
"""
echo "VCF file is $vcf" >> ${barcode}_annotate.txt
echo "Reference fasta file is $ref" >> ${barcode}_annotate.txt
echo "Bam file is $calls" >> ${barcode}_annotate.txt
echo "Indexed bam file is $calls_idx" >> ${barcode}_annotate.txt
medaka tools annotate --dpsp $vcf $ref $calls ${barcode}_calling_supp.vcf
"""
}

process annotate_variant_de_novo {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$barcode"
  if (params.annotate_variant_de_novo_out != "") {
    publishDir "results/${params.annotate_variant_de_novo_out}", mode: 'copy'
  }

  input:
    tuple val(barcode), path(vcf), path(ref)
    tuple val(barcode), path(calls)
    tuple val(barcode), path(calls_idx)
  output:
    tuple val(barcode), path("${barcode}_calling_supp.vcf"), emit: annotated_vcf
    tuple val(barcode), path("${barcode}_annotate.txt")
  script:
"""
echo "VCF file is $vcf" >> ${barcode}_annotate.txt
echo "Reference fasta file is $ref" >> ${barcode}_annotate.txt
echo "Bam file is $calls" >> ${barcode}_annotate.txt
echo "Indexed bam file is $calls_idx" >> ${barcode}_annotate.txt
medaka tools annotate --dpsp $vcf $ref $calls ${barcode}_calling_supp.vcf
"""
}

process medaka_consenus {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$barcode"
  if (params.medaka_consenus_out != "") {
    publishDir "results/${params.medaka_consenus_out}", mode: 'copy'
  }

  input:
    tuple val(barcode), path(calls), path(ref)
  output:
    tuple val(barcode), path("${barcode}/${barcode}_consensus.fasta"), emit: consensus_fasta
    tuple val(barcode), path("${barcode}/")
  script:
"""
echo "Reference fasta file is $ref" >> ${barcode}_annotate.txt
echo "Bam file is $calls" >> ${barcode}_annotate.txt
# Nécessite de modifier le header du fasta de référence
sed 's/:/ :/' ${ref} > new_ref.fasta
medaka_consensus -i ${calls} -d new_ref.fasta -o ${barcode} > ${barcode}_medaka_log.txt 2>&1
cd ${barcode}
mv consensus.fasta ${barcode}_consensus.fasta
mv consensus.fasta.gaps_in_draft_coords.bed ${barcode}_consensus.fasta.gaps_in_draft_coords.bed
mv consensus_probs.hdf ${barcode}_consensus_probs.hdf
rm calls_to_draft.bam
rm calls_to_draft.bam.bai
"""
}