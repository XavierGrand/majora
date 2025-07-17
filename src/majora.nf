#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
/*
========================================================
                        Majora
========================================================

Majora pipeline :
  - Analyzes raw DNA nanopore sequencing data from hepatits B virus
  - Designed for HBV genotyping and variant calling
  - Generates visual outputs and summary table in csv (tab separated) format

********************************************************
                    Help message definition
********************************************************
*/

def help() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:

      nextflow ./src/majora.nf -c ./src/nextflow.config -profile singularity --fast5 or --pod5 or --fastq <PATH_TO_INPUT_FOLDER>

    Nextflow parameters:
      -profile [str]                  Configuration profile to use.
                                      Available: docker, singularity, podman, psmn, ccin2p3, etc.
                                      User can use his own, associated with a nextflow configuration file.
                                      Refer to https://www.nextflow.io/docs/latest/config.html#configuration-file. 

    Input file:                       Only one of the three input types can be specified.
      --fast5 [path]                  Path to the folder containing fast5 files.
      --pod5  [path]                  Path to the folder containing pod5 files.
      --fastq [path]                  Path to the fastq FILE. Only ONE SINGLE fastq file can be used at a time.      
                     

    Nanopore basecalling:
      --skipBC [boolean]              Use if input data are not multiplxed. Thus, dorado basecalling will perform trimming. 
                                      Default is "false"

      --kit_barcoding [str]           To provide Nanopore barcoding kit. 
                                      No default.
            
      --model [str]                   Basecalling model : sup, hac, fast. 
                                      Default: "fast"
    
    Primers:
      --fwprimer [str]                Sequence of forward primer.
                                      Default: "CTACTGTTCAAGCCTCCAAGC"
      --rvprimer [str]                Sequence of reverse primer.
                                      Default: "CGCAGACCAATTTATGCCTAC"
    
    References:                       These sequences will be used during BLAST alignment.
    --ref_db [path]                   To provide csv (tab separated) file containg hbvdb sequence IDs in the second column.
                                      
    --hbvdb  [str]                    Should be followed by either "all" or the uppercase letter corresponding to the desired genotype.
                                      
    --ref_user [path]                 To provide fasta/multifasta file containg hbvdb sequences.

    Advanced parameters:
    --blasthreads [int]               To specify the number of CPU threads to use for faster parallel processing during BLAST alignment.
                                      Default: 18.
    --cuda [str]                      To specify how many GPUs should be used during basecalling. Should follow this synthax: "cuda:<all/integer>".
                                      Default: "cuda:all"
    --minimap [str]                   During variant calling, medaka already performs reads alignment using minimap2. To perform additionnal aligment with minimap2 use "align".
                                      Default: "skip".
    --reads_number [int]              To specify the number of reads that will be sampled for the BLAST alignment.
                                      Default: 100.

    Help:
      --help | --h                    Display this help message.
    
    """.stripIndent()
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
params.ref_user = ''
params.ref_db = ''
params.rvprimer = 'CGCAGACCAATTTATGCCTAC' // Denomination P4(1783-1803)
params.skipBC = false
params.sort_bam = ""
params.minimap = "skip"
params.reads_number = 100

/**** Params out ****/
params.pod5convert_out = '01_pod5convert/'
params.basecalling_out = '02_basecalling/'
params.basecalling_nobc_out = '02_basecalling_nobc/'
params.demux_out = '03_demux/'
params.read_stats_out = '04_reads_stats'
params.convert_fastq = '05_convert_fastq/'
params.grep_primer_out = '06_grep_primer/'
params.filterbylength_out = '07_filterbylength/'
params.sample_fastq_out = '08_sample_fastq/'
params.dl_hbvdb_out =  '09_dl_hdvdb/'
params.makerefdb_out = '10_makerefdb/'
params.doublefastaref_out = '11_double_fatsaref/'
params.makeblastdb_out = '12_makeblastdb/'
params.blast_them_all_out = '13_blast_them_all/'
params.extractref_out = '14_extractref/'
params.index_fasta_out = '15_index_fasta/'
params.mapping_hbv_genome_out = '16_mapping_hbv_genome/'
params.filter_bam_mapped_out = '17_filter_bam_mapped/'
params.sort_bam_out = '16_sort_bam/'
params.index_bam_out = '17_index_bam_out/'
params.qualitycontrol_out = '18_qualitycontrol/'
params.call_hap_variant_out = '19_call_hap_variant  /'
params.annotate_variant_out = '20_annotate_variant/'
params.de_novo_assembly_out = '21_de_novo_assembly/'
params.call_hap_variant_out_de_novo = '22_call_de_novo_hap_variant_out/'
params.vcf_filter_out = '23_vcf_filter/'
params.reindex_vcf_out = '24_reindex_vcf/'
params.lollipop_vcf_out = '25lollipop_vcf/'

params.annotate_variant_de_novo_out = '29_annotate_variant/'


/*
 ****************************************************************
                        Channel Definitions
 ****************************************************************
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


// To download hbv or hdv data bases (all or selected genome)
if (params.hbvdb != '') {
  Channel
    .fromPath(params.hbvdb)
    .set { hbvdb }
}

// To make the user use its own multi-fasta file
if (params.ref_user != '') {
  Channel
    .fromPath(params.ref_user)
    .set { user_ref }
}

// To use a csv file containing reference sequences ID from hbvdb
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

  // Multiplexed FAST5 
  if (params.skipBC  == false &&  params.fast5 != '') {
  include { pod5convert } from './nf_modules/dorado/0.8.2/main.nf'
  include { basecalling } from './nf_modules/dorado/0.8.2/main.nf'
  include { demux } from './nf_modules/dorado/0.8.2/main.nf'
  include { convertfastq } from './nf_modules/bedtools/2.25.0/main.nf'
  }
  // Non-mutiplexed FAST5
  if (params.skipBC  == true &&  params.fast5 != '') {
  include { pod5convert } from './nf_modules/dorado/0.8.2/main.nf'
  include { basecalling_nobc } from './nf_modules/dorado/0.8.2/main.nf'
  }
  // Multiplexed POD5
  else if (params.skipBC  == false && params.pod5 != '') {
  include { basecalling } from './nf_modules/dorado/0.8.2/main.nf'
  include { demux } from './nf_modules/dorado/0.8.2/main.nf'
  include { convertfastq } from './nf_modules/bedtools/2.25.0/main.nf'
    }
  // Non-mutiplexed POD5
  else if (params.skipBC  == true && params.pod5 != '') {
  include { basecalling_nobc } from './nf_modules/dorado/0.8.2/main.nf'
  }

  // Multiplexed FASTQ or FQ  
  else if (params.skipBC  == false && params.pod5 == '' && params.fast5 == '') {
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


if (params.minimap == 'align') {
  include { mapping_hbv_genome } from './nf_modules/minimap2/2.28/main.nf'
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

include { lollipop_vcf as lollipop_vcf1} from './nf_modules/r-majora/main.nf' params(lollipop_vcf_out: '26_lollipop_vcf/')
include { lollipop_vcf as lollipop_vcf2} from './nf_modules/r-majora/main.nf' params(lollipop_vcf_out: '27_lollipop_de_novo_vcf/')

include { vcf_filter as filtervcf1 } from './nf_modules/vcflib/1.0.13/main.nf' params(vcf_filter_out: '22_vcf_filter_out/')
include { vcf_filter as filtervcf2 } from './nf_modules/vcflib/1.0.13/main.nf' params(vcf_filter_out: '23_de_novo_vcf_filter_out/')


/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/
workflow {

// *****File Collection*****

  // Multiplexed FAST5 
  if ( params.skipBC  == false &&  params.fast5 != '' ) {
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

// Multiplexed POD5
  else if ( params.skipBC  == false && params.pod5 != '' ) {
    basecalling(pod5_input)
    demux(basecalling.out.basecalled)
    demux.out.demuxed.collect()
                     .set {ready_bam}
    convertfastq(ready_bam)
    convertfastq.out.converted_fastq.flatten()
                                    .map { it -> [it.getSimpleName(), it] }
                                    .set { ready_fastq }
  }

// Multiplexed  FASTQ
  else if ( params.skipBC  == false && params.pod5 == '' && params.fast5 == '' ) {
    demux(fastq_input)
    demux.out.demuxed.flatten()
                           .map { it -> [it.getSimpleName(), it] }
                           .set { ready_fastq }
  }

// Non-multiplexed FAST5
  else if ( params.skipBC  == true && params.fast5 != '' ) {
    pod5convert(fast5_input)
    basecalling_nobc(pod5convert.out.pod5)
    basecalling_nobc.out.basecalled_nobc.map { it -> [it.getSimpleName(), it] }
                                               .set { ready_fastq }
  }

// Non-mutliplexed POD5
  else if ( params.skipBC  == true && params.pod5 != '' ) {
    basecalling_nobc(pod5_input)
    basecalling_nobc.out.basecalled_nobc.map { it -> [it.getSimpleName(), it] }
                                               .set { ready_fastq }
  }
  

// Non-mutliplexed FASTQ  --> Pre-processing


//***********************Pre-processing**********************

// If multiple single sample fastq files provided need to concatenate in one single file
// Add process concatenate from banjul_xgrand

// Reads statistics before filtering
 read_stats(ready_fastq)

// Reads filtering on primer + on length
 pick_fw_primer(ready_fastq, params.fwprimer, 'fw')
 pick_rv_primer(pick_fw_primer.out.filtered_fastq, params.rvprimer, 'rv')
 filterbylength(pick_rv_primer.out.filtered_fastq, pick_rv_primer.out.filtered_tsv)

// Randomly draw filtered reads for BLAST alignment (Default = 100 reads)
 sample_fastq(filterbylength.out.length_filtered_fastq)


//***********************Download / upload reference depending on used option **********************

  if ( params.ref_db !== '' ) {
      makerefdb(ref_db)
      doublefastaref(makerefdb.out.ref_db)
  }

  else if ( params.hbvdb !== '' ) {
    dl_hbvdb()
    dl_hbvdb .out.reference_db.set { genotype }
    doublefastaref(genotype)
  }

  else if ( params.ref_user !== '' ) {
   doublefastaref(user_ref)
   makeblastdb(doublefastaref.out.doubledfasta)
  }

//*********************** Blast **********************
 makeblastdb(doublefastaref.out.doubledfasta)
 blast_them_all(sample_fastq.out.sampled_fastq, makeblastdb.out.blastdb.collect())

//*********************** Extract best results and corresponding reference sequence **********************
 extractref(blast_them_all.out.bestref, doublefastaref.out.doubledfasta.collect())
 //index_fasta(extractref.out.referenceseq)

//*********************** Optionnal alignment with minimap2 *****************
if ( params.minimap == 'align' ) {
 mapping_hbv_genome (filterbylength.out.length_filtered_fastq.combine(extractref.out.referenceseq, by: 0))
 filter_bam_mapped(mapping_hbv_genome.out.bam)
 sort_bam(filter_bam_mapped.out.bam)
 index_bam(sort_bam.out.bam)
 qualitycontrol(demux.out.sequence_summary, index_bam.out.bam_idx)
  }


 //********************** Variant calling *************
 vc1(extractref.out.referenceseq.combine(filterbylength.out.length_filtered_fastq, by: 0))

 // To get additionnal informations about variations. These will be used to calculate variation frenquencies
 annotate_variant(vc1.out.vcf.combine(extractref.out.referenceseq, by: 0), vc1.out.calls, vc1.out.calls_idx)

 // Filtering identified variation on quality (must be at least 20 to be kept) 
 filtervcf1(annotate_variant.out.annotated_vcf)

 // Reindexing variations positions according to actual EcoRI cutting site
 reindex_vcf1(extractref.out.referenceseq.combine(blast_them_all.out.counthits, by: 0).combine(filtervcf1.out.vcf_filtered, by: 0), params.fwprimer, params.rvprimer)

 // Visualization
 lollipop_vcf1(reindex_vcf1.out.tsv)
 

 //********************* de novo assembly *************************
 /*
 de_novo_assembly(filterbylength.out.length_filtered_fastq)
 vc2(de_novo_assembly.out.assembly_fasta.combine(filterbylength.out.length_filtered_fastq, by: 0))
 annotate_variant_de_novo(vc2.out.vcf.combine(extractref.out.referenceseq, by: 0), vc2.out.calls, vc2.out.calls_idx)

 filtervcf2(vc2.out.vcf)
 reindex_vcf2(extractref.out.referenceseq.combine(blast_them_all.out.bestref, by: 0).combine(filtervcf2.out.vcf_filtered, by: 0), params.fwprimer, params.rvprimer)
 lollipop_vcf2(reindex_vcf2.out.tsv)
 */
 
if (params.pod5 != '' || params.fast5 != '' ) {
qualitycontrol(vc1.out.calls.combine(vc1.out.calls_idx, by: 0).combine(demux.out.sequence_summary))
}


}