#!/usr/bin/env nextflow
import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.sample_sheet = false
params.fastq_folder = false
params.fastq_suffix = ".fastq.gz"
params.output = false
params.genome = false
params.genome_json = false
params.genome_key = false
params.barcodes = false

// Quality trimming
params.min_qvalue = 20
params.min_align_score = 40
params.RD1_ADAPTER_3P = "GATCGGAAGAGCACACGTCTGAACTCCAGTC"
params.RD2_ADAPTER_3P = "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

// Unique molecular tags
params.barcode_length = 6
params.barcode_max_homopolymer = 6
params.max_barcode_mismatch = 2

// Trim a fixed amount from the 5' of both reads
params.trim_length = 5

// Split up each specimen into shards for parallel processing
params.n_shards = 100

// Minimum proportion of bases needed to call an SSC base
params.min_base_prop = 0.7

// Maximum distance that an SSC may change position after realignment
params.max_realign_offset = 5

// Minimum number of reads needed for EACH SSC to keep a DSC
// Note that we will keep a DSC from either end of a molecule,
// even if the other end doesn't have enough data to use
params.min_reads = 3

// If specified, mask everything outside of a user-provided target region (BED file)
params.target_regions_bed = false

// If specified, a CSV file containing the set of coordinates which
// should be ignored during mutation/variant calling
params.ignore_coordinates = false

// Set the containers to use for each component
params.container__cutadapt = "quay.io/biocontainers/cutadapt:3.5--py36hc5360cc_0"
params.container__fastqc = "quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1"
params.container__multiqc = "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
params.container__pandas = "quay.io/fhcrc-microbiome/python-pandas:0fd1e29"
params.container__python_plotting = "quay.io/hdc-workflows/python-plotting:b50a842"
params.container__bwa = "quay.io/hdc-workflows/bwa-samtools:93deeda"
params.container__bcftools = "quay.io/biocontainers/bcftools:1.13--h3a49de5_0"

// Import sub-workflows
include { manifest_wf } from './modules/manifest'
include { quality_wf } from './modules/quality'
include { barcodes_wf } from './modules/barcodes'
include { trim_wf } from './modules/trim'
include { align_wf } from './modules/align'
include { family_wf } from './modules/family'
include { variants_wf } from './modules/variants'


// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run FredHutch/maddd-seq <ARGUMENTS>

Input Data:
  --sample_sheet        CSV file listing samples with headers: specimen, R1, and R2
  -or-
  --fastq_folder        Folder containing paired-end FASTQ data
                        Note that all FASTQ files must contain either 'R1' or 'R2'
                        in the filename, with the suffix noted by --fastq_suffix,
                        described below.

Required Arguments:
  --output              Folder to write output files to
  --barcodes            Path to text file containing the barcode whitelist

  --genome              Reference genome indexed for alignment by BWA.
  -or-
  --genome_json         JSON encoded dict mapping keys to `genome` (required) and `target_regions_bed` (optional).
  --genome_key          Key corresponding to entry in `genome_json`

Optional Arguments:
  --min_qvalue          Minimum quality score used to trim data (default: ${params.min_qvalue})
  --min_align_score     Minimum alignment score (default: ${params.min_align_score})
  --barcode_length      Length of the barcodes ligated to the 5' end of each read
                        This should correspond to the length of every entry in
                        the barcode text file specified with --barcodes.
                        (default: ${params.barcode_length})
  --barcode_max_homopolymer
                        Maximum number of repeated bases in barcode sequences
                        (default: ${params.barcode_max_homopolymer})
  --max_barcode_mismatch
                        Maximum number of nucleotide differences used for barcode
                        error correction
  --trim_length         Number of bases to trim from the 5' of each read after
                        ligated barcode sequences are removed
  --n_shards            Number of parallel processes to use for creating families
                        (default: ${params.n_shards})
  --min_base_prop       Minimum proportion of bases needed to call one base of a SSC
                        (default: ${params.min_base_prop})
  --min_reads           Minimum number of reads needed for EACH SSC to keep a DSC
                        Note that we will keep a DSC from either end of a molecule,
                        even if the other end doesn't have enough data to use.
                        (default: ${params.min_reads})
  --max_realign_offset  Maximum distance that an SSC may change position after realignment
                        (default: ${params.max_realign_offset})
  --RD1_ADAPTER_3P      Sequence of the universal Illumina adapter found at the 3'
                        of the R1 (default: ${params.RD1_ADAPTER_3P})
  --RD2_ADAPTER_3P      Sequence of the universal Illumina adapter found at the 3'
                        of the R2 (default: ${params.RD2_ADAPTER_3P})
  --target_regions_bed  Optional BED file describing a set of target genomic regions
                        to use for analysis, e.g. exome target capture.
  --ignore_coordinates  Optional CSV listing locations which should be ignored from
                        SNP/adduct calling, with columns 'chr' and 'pos'.
                        Note that positional coordinates in this CSV should be 1-based,
                        as is the convention with genomic coordinates, rather than the
                        0-based coordinates used more commonly in programming languages.
  --fastq_suffix        When input data is specified by --fastq_folder (instead of --sample_sheet),
                        only files with this suffix will be used as inputs (default: ${params.fastq_suffix})


Sample Sheet:
  The sample sheet is a CSV listing all of the duplex sequencing data to be analyzed.
  The sample sheet must contain the column headers: specimen,R1,R2
  The R1 and R2 file paths may start with ftp://, s3://, or even just be a path to a local file.
  The specimen must be unique, and only contain a-z, A-Z, 0-9, or _.
    """.stripIndent()
}


// Main workflow
workflow {

    // Show help message if the user specifies the --help flag at runtime
    // or if --genome and --output are not provided
    if ( params.help || params.output == false ){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    }

    // If neither --sample_sheet or --fastq_folder are provided
    if ( !params.sample_sheet && !params.fastq_folder ){
        log.info"""
        ERROR: Must provide either --sample_sheet or --fastq_folder.
        View help text with --help for more details.
        """.stripIndent()
        exit 1
    }

    // If both --sample_sheet and --fastq_folder are provided
    if ( params.sample_sheet && params.fastq_folder ){
        log.info"""
        ERROR: Must provide either --sample_sheet or --fastq_folder, but not both.
        View help text with --help for more details.
        """.stripIndent()
        exit 1
    }

    // Check on the combinations of genome parameter inputs
    // If both --genome and --genome_json are provided
    if ( params.genome && params.genome_json ){
        log.info"""
        ERROR: Must provide only one of --genome or --genome_json, not both.
        View help text with --help for more details.
        """.stripIndent()
        exit 1
    }

    // If neither --genome and --genome_json are provided
    if ( !params.genome && !params.genome_json ){
        log.info"""
        ERROR: Must provide either --genome or --genome_json
        View help text with --help for more details.
        """.stripIndent()
        exit 1
    }

    // If --genome_json is provided but not --genome_key
    if ( params.genome_json && !params.genome_key ){
        log.info"""
        ERROR: Must provide --genome_key when using --genome_json
        View help text with --help for more details.
        """.stripIndent()
        exit 1
    }

    if ( params.genome_json && !file(params.genome_json).exists() ) {
        log.info"""
        ERROR: --genome_json file ${params.genome_json} does not exist
        """.stripIndent()
        exit 1
    }

    // If the user provided a genome json mapping
    if ( params.genome_json ) {
        genome_map = new JsonSlurper().parseText( file(params.genome_json).text )
        genome_map_partial = genome_map["${params.genome_key}"]
        if ( ! genome_map["${params.genome_key}"] ) {
            log.info"""
            ERROR: --genome_key ${params.genome_key} does not exist in --genome_json ${params.genome_json}
            """.stripIndent()
            exit 1
        }
        params.genome_path = genome_map["${params.genome_key}"]['genome']

        //String configJSON = new File("${params.wfconfig}").text
        //def wfi = jsonSlurper.parseText(configJSON)
        
        // Use the genome json map's target region bed only
        // if the parameter was not passed in
        if ( !params.target_regions_bed ) {
            params.target_regions_bed_path = genome_map["${params.genome_key}"]['target_regions_bed']
        } else {
            params.target_regions_bed_path = params.target_regions_bed
        }
        log.info"""
        Using --genome_json argument
        """.stripIndent()
    } else {
        // If the user provided a genome via --genome
        params.genome_path = params.genome
        params.target_regions_bed_path = params.target_regions_bed
    }

    // If the user provided a sample sheet
    if ( params.sample_sheet ){

        // Parse the manifest (sample sheet)
        manifest_wf(
            params.sample_sheet
        )
        // output:
        //   reads:
        //     tuple val(specimen), path(read_1), path(read_2)

        // Set fastq_ch with the output of that subworkflow
        fastq_ch = manifest_wf.out.reads

    // Otherwise, the user must have provided a FASTQ folder
    } else {

        // Set fastq_ch with the pairs of files from the FASTQ folder
        // which end with the suffix, and which vary only by R1/R2
        fastq_ch = Channel
            .fromFilePairs("${params.fastq_folder}/*R{1,2}*${params.fastq_suffix}")
            .map{
                [it[0], it[1][0], it[1][1]]
            }
        // Cardinality: tuple val(specimen), path(read_1), path(read_2)

    }


    // Perform quality trimming on the ends of all reads
    quality_wf(
        fastq_ch
    )
    // output:
    //   reads:
    //     tuple val(specimen), path(read_1), path(read_2)
    // publish:
    //   1_input_data/input/multiqc_report.html
    //   1_input_data/quality_trimmed/multiqc_report.html

    // Clip the barcodes from the ends of the FASTQ files
    // This part of the workflow will also perform error
    // correction on the barcode sequences.
    barcodes_wf(
        quality_wf.out.reads,
        file(params.barcodes)
    )
    // output:
    //   bam:
    //     tuple val(specimen), path(bam)
    // publish:
    //   2_barcode_trimmed/{specimen}/barcode_counts.csv.gz
    //   2_barcode_trimmed/{specimen}/barcode_corrections.csv.gz
    //   2_barcode_trimmed/multiqc_report.html

    // Trim a fixed number of bases from the beginning of each read
    trim_wf(
        barcodes_wf.out.reads
    )
    // output:
    //   bam:
    //     tuple val(specimen), path(bam)
    // publish:
    //   3_end_trimmed/multiqc_report.html

    // The pre-compiled genome reference is provided as an input
    // It must contain a wildcard capturing all required index files
    genome_ref = Channel.fromPath("${params.genome_path}").toSortedList()

    // Align the barcode-clipped reads to the reference genome
    align_wf(
        // WGS data in FASTQ format
        trim_wf.out.reads,
        // Reference genome, indexed for alignment with BWA
        genome_ref,
        // Table linking each uncorrected barcode to its corrected sequence
        barcodes_wf.out.csv
    )
    // output:
    //   bam:
    //     tuple val(specimen), path(bam)
    // publish:
    //   4_aligned/{specimen}/{specimen}.flagstats
    //   4_aligned/multiqc_report.html

    // Group reads into families based on barcodes and alignment position
    // This sub-workflow will also collapse and summarize SSCs and DSCs
    family_wf(
        align_wf.out.bam,
        genome_ref
    )
    // input:
    //   tuple val(specimen), path(bam), path(barcodes_csv_gz)
    // output:
    //   bam: 
    //     tuple val(specimen), path(ssc_bam), path(dsc_bam)
    // publish:
    //   5_families/family_summary.csv
    //   5_families/<specimen>/aligned.bam[.bai]
    //   6_all_SSC/<specimen>/POS.SSC.bam[.bai]
    //   6_all_SSC/<specimen>/NEG.SSC.bam[.bai]
    //   6_all_SSC/<specimen>/SSC.details.csv.gz

    // Call variants and adducts
    variants_wf(
        // Channel with aligned SSCs and DSCs
        family_wf.out.bam,
        // Genome sequence reference
        genome_ref,
        // CSV with the number of reads per barcode, for each specimen
        barcodes_wf.out.counts
    )
    // publish:
    //   7_filtered_SSC/<specimen>/POS.SSC.bam[.bai]
    //   7_filtered_SSC/<specimen>/NEG.SSC.bam[.bai]
    //   7_filtered_SSC/<specimen>/SSC.details.csv.gz
    //   8_variants/<specimen>/variants.vcf
    //   8_variants/<specimen>/adducts.vcf
    //   8_variants/variant_summary.csv

}