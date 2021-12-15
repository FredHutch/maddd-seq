#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.genome_fasta = false
params.output = false

// If specified, do not run RepeatMasker on the genome before indexing
params.skip_repeatmasker = false
params.repeatmasker_species = 'human'

// If specified, publish the intermediate result files
params.save_intermediates = false


// Set the containers to use for each component
params.container__bwa = "quay.io/hdc-workflows/bwa-samtools:93deeda"
params.container__repeatmasker = "quay.io/biocontainers/repeatmasker:4.1.2-p1"

// Import sub-workflows
include { repeatmasker; bwa_index } from './modules/index'


// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run FredHutch/maddd-seq/index.nf <ARGUMENTS>

Required Arguments:
  --genome_fasta        Reference genome in FASTA format
  --output              Folder to write output files to

Optional Arguments:
  --skip_repeatmasker   If specified, do not mask repetitive elements with the RepeatMasker tool
                        (default: ${params.skip_repeatmasker})
  --save_intermediates  If specified, publish all 'intermediate' files.
                        These are files created by various steps but not usually published
                        (default: ${params.save_intermediates})

  --repeatmasker_species
                        Species name provided to RepeatMasker (default: ${params.repeatmasker_species})
                        Allowed: human mouse rattus "ciona savignyi" arabidopsis,
                        mammal, carnivore, rodentia, rat, cow, pig, cat, dog, chicken, fugu,
                        danio, "ciona intestinalis" drosophila, anopheles, worm, diatoaea,
                        artiodactyl, arabidopsis, rice, wheat, maize
    """.stripIndent()
}


// Main workflow
workflow {

    // Show help message if the user specifies the --help flag at runtime
    // or if --genome_fasta and --output are not provided
    if ( params.help || params.output == false || params.genome_fasta == false ){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    }

    // If the user requested --skip_repeatmasker
    if ( params.skip_repeatmasker ){

        // Index the masked genome
        bwa_index(file(params.genome_fasta))

    // Otherwise, run RepeatMasker
    } else {

        // Run RepeatMasker on the genome FASTA
        repeatmasker(file(params.genome_fasta))

        // Index the input
        bwa_index(repeatmasker.out)

    }

}