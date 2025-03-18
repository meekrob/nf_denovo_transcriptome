#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters are defined in nextflow.config
// This section should be removed from main.nf
// params {
//     input = "samples.txt"
//     outdir = "./results"
//     publish_dir_mode = "copy"
//     kmer_size = 25
//     target_depth = 100
//     assembler = "rnaspades"
// }

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
        if (line.trim().startsWith('#')) return  // Skip comment lines
        def fields = line.trim().split(/\s+/)
        if (fields.size() >= 3) {
            def sample_id = fields[0]
            def r1 = file(fields[1], checkIfExists: true)
            def r2 = file(fields[2], checkIfExists: true)
            return [sample_id, [r1, r2]]
        }
    }
    .filter { it != null }  // Filter out skipped lines
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