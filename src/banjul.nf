#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
/*
========================================================
                        Banjul
========================================================

banjul pipeline :
*ajouter une explication

********************************************************
                    Help message definition
********************************************************
*/

def help() {
  log.info'''
    Usage:
    Command for running the pipeline:

    ICI AJOUTER LA COMMANDE FINALE POUR RUN LE PIPELINE

    Help:
        --help | --h      Displays this help  message

    AJOUTER PAR LA SUITE

    '''.stripIndent() // voir Groovy documentation, mise en forme dans le terminal
}

// Show help message
params.help = '' // initialisation du paramètre help à une chaine 
params.h = ''

if (params.help || params.h) {
  help()
  exit 0
}

/*
********************************************************
                   Default Parameters
********************************************************
*/
/**** Params in ****/
params.blasthreads = 18
params.cuda = "cuda:all"
params.fast5 = ''
params.fastq = ''
params.fwprimer = 'CTACTGTTCAAGCCTCCAAGC' // Denomination P3(1857-1877)
params.hbvdb = ''
params.index_bam = ''
params.index_fasta = ''
params.kit_barcoding = '' // Barcoding kit : "EXP-PBC001"
params.model = "fast" // other params: sup, hac
params.pod5 = ''
params.primerloc = ''
params.ref = ''
params.ref_db = ''
params.rvprimer = 'CGCAGACCAATTTATGCCTAC' // Denomination P4(1783-1803)
params.skipBC = true
params.sort_bam = ""
params.minimap = "skip"
/**** Params out ****/
params.pod5convert_out = '01_pod5convert/'
params.basecalling_out = '02_basecalling/'
params.basecalling_nobc_out = '02_basecalling_nobc/'
params.demux_out = '03_demux/'
params.read_stats_out = '04_reads_stats'
params.convert_fastq = '05_convert_fastq/'
params.grep_primer_out = '06_grep_primer/'
params.filterbylength_out = '07_filterbylength/'
params.ultimate_read_filter_out = '08_ultimate_read_filter'
params.sample_fastq_out = '09_sample_fastq/'
params.dl_hbvdb_out =  '10_dl_hdvdb/'
params.makerefdb_out = '11_makerefdb/'
params.doublefastaref_out = '12_double_fatsaref/'
params.makeblastdb_out = '013_makeblastdb/'
params.blast_them_all_out = '14_blast_them_all/'
params.extractref_out = '15_extractref/'
params.index_fasta_out = '12_index_fasta/'
params.mapping_hbv_genome_out = '13_mapping_hbv_genome/'
params.filter_bam_mapped_out = '14_filter_bam_mapped/'
params.sort_bam_out = '15_sort_bam/'
params.index_bam_out = '16_index_bam_out/'
params.qualitycontrol_out = '17_qualitycontrol/'
params.call_hap_variant_out = '18_call_hap_variant_out/'
params.annotate_variant_out = '19_annotate_variant/'
params.de_novo_assembly_out = '20_de_novo_assembly/'
params.call_hap_variant_out_de_novo = '21_call_de_novo_hap_variant_out/'
params.vcf_filter_out = '22_vcf_filter/'
params.reindex_vcf_out = '23_reindex_vcf/'
params.lollipop_vcf_out = '24_lollipop_vcf/'

params.annotate_variant_de_novo_out = '29_annotate_variant/'


// TEST FOR HAPLOTYPING
params.bam2sam_out = 'a_bam2sam/'
params.rv_haplo_out = 'b_rv_haplo/'
params.strainline_out = 'e_strainline/'
params.medaka_consenus_out = 'd_medaka_consenus/'
params.bam2fasta_out = 'c_bam2fasta'
/*
 ****************************************************************
                        Channel Definitions
 ****************************************************************
  .branch {file->
        fast5: file == 'fast5'
        pod5: file == 'pod5'
  }
*/

// Create input channel depending on input file extension

// FASTQ input
if (params.fastq != '' ) {
  Channel.fromPath(params.fastq )
                 .ifEmpty { error 'No input folder defined.' }
                 .map { file -> tuple('fastq', file) }
                 .set { ready_fastq }
  log.info "Input fastq folder : ${params.fastq}"

}

// POD5 input
else if (params.pod5 != '') {
  Channel.fromPath(params.pod5 )
                .ifEmpty { error 'No input folder defined.' }
                .set { pod5_input}
  log.info "Input pod5 folder : ${params.pod5}"
}

// FAST5 input
else if (params.fast5 != '') {
  Channel.fromPath( params.fast5 )
                 .ifEmpty { error 'No input folder defined.' }
                 .set { fast5_input}
  log.info "Input fast5 folder : ${params.fast5}"
}

/*
To detect input's extension
Channel
    .fromPath(params.input+'**', type: 'file')
    .filter { it.extension in ['fastq', 'fast5', 'pod5'] }
    .collect()
    .map { file -> file.extension }
    .unique()
    .set{file_ext}
*/

// To download hbv or hdv data bases (all or selected genome)
if (params.hbvdb != '') {
  Channel
    .fromPath(params.hbvdb)
    .set { hbvdb }
}

// To make the user use its own multi-fasta file
if (params.ref != '') {
  Channel
    .fromPath(params.ref)
    .set { user_ref }
}

if (params.ref_db != '') {
  Channel
    .fromPath(params.ref_db)
    .set { ref_db }
}

/*
 ****************************************************************
                          Imports
 ****************************************************************
*/
// ********** FILE COLLECTION CONDITIONAL IMPORTS *******************

  // FAST5 + Barcoding = TRUE
  if (params.kit_barcoding &&  params.fast5 != '') {
  include { pod5convert } from './nf_modules/dorado/0.8.2/main.nf'
  include { basecalling } from './nf_modules/dorado/0.8.2/main.nf'
  include { demux } from './nf_modules/dorado/0.8.2/main.nf'
  include { convertfastq } from './nf_modules/bedtools/2.25.0/main.nf'
  }
  // FAST5 + Barcoding = FALSE
  if (!params.kit_barcoding &&  params.fast5 != '') {
  include { pod5convert } from './nf_modules/dorado/0.8.2/main.nf'
  include { basecalling_nobc } from './nf_modules/dorado/0.8.2/main.nf'
  }
  // POD5 + Barcoding = TRUE
  else if (params.kit_barcoding && params.pod5 != '') {
  include { basecalling } from './nf_modules/dorado/0.8.2/main.nf'
  include { demux } from './nf_modules/dorado/0.8.2/main.nf'
  include { convertfastq } from './nf_modules/bedtools/2.25.0/main.nf'
    }
  // POD5 + Barcoding = FALSE
  else if (!params.kit_barcoding && params.pod5 != '') {
  include { basecalling_nobc } from './nf_modules/dorado/0.8.2/main.nf'
  }

  // FASTQ or FQ + Barcoding = TRUE
  // VERSION UTILISATEUR PEU ATTENTIF
  //else if ( params.kit_barcoding && params.file_type == "fastq" || params.file_type == "fastq.gz" || params.file_type == "fq") || params.file_type == "fq.gz"{
  //}
  else if (params.kit_barcoding && params.pod5 == '' && params.fast5 == '') {
  include { demux } from './nf_modules/dorado/0.8.2/main.nf'
  }

   // FASTQ or FQ + Barcoding = FALSE
// ********** PRE PROCESSING IMPORTS *******************

include { read_stats } from './nf_modules/seqkit/2.8.2/main.nf'
include { grep_primer as pick_fw_primer } from './nf_modules/seqkit/2.8.2/main.nf' // include as pour appliquer le même process à deux fichers/variables distincts, ici primers fwd et rev
include { grep_primer as pick_rv_primer } from './nf_modules/seqkit/2.8.2/main.nf'
include { filterbylength } from './nf_modules/seqkit/2.8.2/main.nf'

include { sample_fastq } from './nf_modules/seqtk/1.3/main.nf'
include { dl_hbvdb } from './nf_modules/blast/2.15.0/main.nf'
include { doublefastaref  } from './nf_modules/seqkit/2.8.2/main.nf'
include { makeblastdb  } from './nf_modules/blast/2.15.0/main.nf'
include { makerefdb  } from './nf_modules/blast/2.15.0/main.nf'
include { blast_them_all  } from './nf_modules/blast/2.15.0/main.nf'
include { extractref } from './nf_modules/seqkit/2.8.2/main.nf'
include { index_fasta  } from './nf_modules/samtools/1.20/main.nf'

if (params.minimap != 'skip') {
  include { mapping_hbv_genome } from './nf_modules/minimap2/2.17/main.nf'
  include { filter_bam_mapped } from './nf_modules/samtools/1.20/main.nf'
  include { sort_bam } from './nf_modules/samtools/1.20/main.nf'
  include { index_bam } from './nf_modules/samtools/1.20/main.nf'
}


include { qualitycontrol } from './nf_modules/pycoqc/2.5.2/main.nf'
include { call_hap_variant as vc1 } from './nf_modules/medaka/2.0.1/main.nf' params(call_hap_variant_out: '18_call_hap_variant_out/')
include { call_hap_variant as vc2 } from './nf_modules/medaka/2.0.1/main.nf' params(call_hap_variant_out: '21_call_de_novo_hap_variant_out/')

include { annotate_variant } from './nf_modules/medaka/2.0.1/main.nf'
include { de_novo_assembly } from './nf_modules/flye/2.9.5/main.nf'
include { annotate_variant_de_novo } from './nf_modules/medaka/2.0.1/main.nf'


include { reindex_vcf as reindex_vcf1 } from './nf_modules/bash/main.nf' params(reindex_vcf_out: '24_reindexed_vcf/')
include { reindex_vcf as reindex_vcf2 } from './nf_modules/bash/main.nf' params(reindex_vcf_out: '25_reindexed_de_novo_vcf/')

include { lollipop_vcf as lollipop_vcf1} from './nf_modules/R/main.nf' params(lollipop_vcf_out: '26_lollipop_vcf/')
include { lollipop_vcf as lollipop_vcf2} from './nf_modules/R/main.nf' params(lollipop_vcf_out: '27_lollipop_de_novo_vcf/')

include { vcf_filter as filtervcf1 } from './nf_modules/vcflib/1.0.13/main.nf' params(vcf_filter_out: '22_vcf_filter_out/')
include { vcf_filter as filtervcf2 } from './nf_modules/vcflib/1.0.13/main.nf' params(vcf_filter_out: '23_de_novo_vcf_filter_out/')


// TEST FOR HAPLOTYPING
include { medaka_consenus } from './nf_modules/medaka/2.0.1/main.nf'
include { rv_haplo } from './nf_modules/rvhaplo/3.0/main.nf'
include { bam2sam } from './nf_modules/samtools/1.20/main.nf'
include { strainline } from './nf_modules/strainline/1.0/main.nf'
include { bam2fasta } from './nf_modules/samtools/1.20/main.nf'
/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/
workflow {

// *****File Collection*****

  //Multiplexed samples fast5
  if ( params.kit_barcoding != '' && params.fast5 != '' ) {
    pod5convert(fast5_input)
    basecalling(pod5convert.out.pod5)
    demux(basecalling.out.basecalled)
    demux.out.demuxed.collect()
                     .set {ready_bam}      
    convertfastq(ready_bam)
    convertfastq.out.converted_fastq.flatten()
                                    .map { it -> [it.getSimpleName(), it] }
                                    .set { ready_fastq }
  }

// Multiplexed samples pod5
  else if ( params.kit_barcoding != '' && params.pod5 != '' ) {
    basecalling(pod5_input)
    demux(basecalling.out.basecalled)
    demux.out.demuxed.collect()
                     .set {ready_bam}
    convertfastq(ready_bam)
    convertfastq.out.converted_fastq.flatten()
                                    .map { it -> [it.getSimpleName(), it] }
                                    .set { ready_fastq }
  }

// Multiplexed samples fastq
  else if ( params.kit_barcoding != '' && params.pod5 == '' && params.fast5 == '' ) {
    demux(fastq_input)
    demux.out.demuxed.flatten()
                           .map { it -> [it.getSimpleName(), it] }
                           .set { ready_fastq }
  }

// Single sample fast5
  else if ( params.kit_barcoding == '' && params.fast5 != '' ) {
    pod5convert(fast5_input)
    basecalling_nobc(pod5convert.out.pod5)
    basecalling_nobc.out.basecalled_nobc.map { it -> [it.getSimpleName(), it] }
                                               .set { ready_fastq }
  }

// Single sample pod5
  else if ( params.kit_barcoding == '' && params.pod5 != '' ) {
    basecalling_nobc(pod5_input)
    basecalling_nobc.out.basecalled_nobc.map { it -> [it.getSimpleName(), it] }
                                               .set { ready_fastq }
  }
  

// Single sample fastq  --> Pre-processing


//***********************Pre-processing**********************

// If multiple single sample fastq files provided need to concatenate in one single file
// Add process concatenate from banjul_xgrand

// Filter primers

 read_stats(ready_fastq)
 pick_fw_primer(ready_fastq, params.fwprimer, 'fw')

 pick_rv_primer(pick_fw_primer.out.filtered_fastq, params.rvprimer, 'rv')



 filterbylength(pick_rv_primer.out.filtered_fastq, pick_rv_primer.out.filtered_tsv)


  // take a sample of 100 sequences from reads
 sample_fastq(filterbylength.out.length_filtered_fastq)

//***********************Download / upload reference **********************

  if ( params.ref_db !== '' ) {
      makerefdb(ref_db)
      doublefastaref(makerefdb.out.ref_db)
  }

  else if ( params.hbvdb !== '' ) {
    dl_hbvdb()
    dl_hbvdb .out.reference_db.set { genotype }
    doublefastaref(genotype)
  }

  else if ( params.ref !== '' ) {
   doublefastaref(user_ref)
   makeblastdb(doublefastaref.out.doubledfasta)
  }

//*********************** Blast **********************
 makeblastdb(doublefastaref.out.doubledfasta)
 blast_them_all(sample_fastq.out.sampled_fastq, makeblastdb.out.blastdb.collect())
/*
//*********************** Extract best results and corresponding reference sequence **********************
 extractref(blast_them_all.out.bestref, doublefastaref.out.doubledfasta.collect())
 //index_fasta(extractref.out.referenceseq)




//*********************** Filter mapping results *****************
if ( params.minimap != 'skip' ) {
  //********************** Align reads on reference sequence ************
 mapping_hbv_genome (filterbylength.out.length_filtered_fastq.combine(extractref.out.referenceseq, by: 0))
 filter_bam_mapped(mapping_hbv_genome.out.bam)
 sort_bam(filter_bam_mapped.out.bam)
 index_bam(sort_bam.out.bam)
 qualitycontrol(demux.out.sequence_summary, index_bam.out.bam_idx)
  }


 //********************** Variant calling *************

 vc1(extractref.out.referenceseq.combine(filterbylength.out.length_filtered_fastq, by: 0))
 annotate_variant(vc1.out.vcf.combine(extractref.out.referenceseq, by: 0), vc1.out.calls, vc1.out.calls_idx)
 filtervcf1(annotate_variant.out.annotated_vcf)

 reindex_vcf1(extractref.out.referenceseq.combine(blast_them_all.out.counthits, by: 0).combine(filtervcf1.out.vcf_filtered, by: 0), params.fwprimer, params.rvprimer)
 lollipop_vcf1(reindex_vcf1.out.tsv)
 */
 //********************* de novo assembly *************************
 /*
 de_novo_assembly(filterbylength.out.length_filtered_fastq)
 vc2(de_novo_assembly.out.assembly_fasta.combine(filterbylength.out.length_filtered_fastq, by: 0))
 annotate_variant_de_novo(vc2.out.vcf.combine(extractref.out.referenceseq, by: 0), vc2.out.calls, vc2.out.calls_idx)

 filtervcf2(vc2.out.vcf)
 reindex_vcf2(extractref.out.referenceseq.combine(blast_them_all.out.bestref, by: 0).combine(filtervcf2.out.vcf_filtered, by: 0), params.fwprimer, params.rvprimer)
 lollipop_vcf2(reindex_vcf2.out.tsv)
 */
/* 
if (params.pod5 != '' || params.fast5 != '' ) {
qualitycontrol(vc1.out.calls.combine(vc1.out.calls_idx, by: 0).combine(demux.out.sequence_summary))
}
*/

// TEST HAPLOTYPONG
/*
bam2sam(vc1.out.calls, vc1.out.calls_idx)
medaka_consenus(vc1.out.calls.combine(reindex_vcf1.out.ref_fasta, by: 0))
rv_haplo(bam2sam.out.sam.combine(medaka_consenus.out.consensus_fasta, by: 0))
bam2fasta(vc1.out.calls)
strainline(bam2fasta.out.fasta)
*/
}