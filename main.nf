#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters
params {
    input = "samples.txt"           // Input file with sample_id, R1, R2 paths
    outdir = "./results"            // Output directory
    publish_dir_mode = "copy"       // Copy outputs to outdir
    kmer_size = 25                  // For normalization
    target_depth = 100              // Normalization target
    assembler = "rnaspades"         // Options: "rnaspades" or "trinity"
}

// Include modules
include { FASTP } from './modules/local/fastp'
include { REPAIR } from './modules/local/repair'
include { MERGE_READS } from './modules/local/merge_reads'
include { BBNORM } from './modules/local/bbnorm'
include { RNASPADES } from './modules/local/rnaspades'
include { TRINITY } from './modules/local/trinity'
include { BUSCO } from './modules/local/busco'

// Input channel from samples.txt
Channel
    .fromPath(params.input, checkIfExists: true)
    .splitText()
    .map { line ->
        def (sample_id, r1, r2) = line.trim().split(/\s+/)
        [sample_id, [file(r1, checkIfExists: true), file(r2, checkIfExists: true)]]
    }
    .set { read_pairs }

// Workflow
workflow {
    // Trimming
    FASTP(read_pairs)

    // Parallel repair on each trimmed pair
    REPAIR(FASTP.out.trimmed_reads)

    // Merge repaired reads
    MERGE_READS(REPAIR.out.repaired_reads.collect())

    // Normalize merged reads
    BBNORM(MERGE_READS.out.merged_reads)

    // Assembly (conditional based on params.assembler)
    if (params.assembler == "rnaspades") {
        RNASPADES(BBNORM.out.normalized_reads)
        assembly_ch = RNASPADES.out.assembly
    } else if (params.assembler == "trinity") {
        TRINITY(BBNORM.out.normalized_reads)
        assembly_ch = TRINITY.out.assembly
    } else {
        error "Invalid assembler specified: ${params.assembler}. Use 'rnaspades' or 'trinity'."
    }

    // Quality assessment with BUSCO
    BUSCO(assembly_ch)
} 