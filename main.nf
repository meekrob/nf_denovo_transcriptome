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
//     skip_to_merge = false
//     repaired_reads_dir = "${params.outdir}/repair"
// }

// Include modules
include { FASTP } from './modules/local/fastp'
include { REPAIR } from './modules/local/repair'
include { MERGE_READS } from './modules/local/merge_reads'
include { BBNORM } from './modules/local/bbnorm'
include { RNASPADES } from './modules/local/rnaspades'
include { TRINITY } from './modules/local/trinity'
include { BUSCO } from './modules/local/busco'

// Create log file for missing files
def logMissingFile(file_path) {
    def logFile = file("${params.outdir}/missing_files.log")
    if(!logFile.parent.exists()) {
        logFile.parent.mkdirs()
    }
    logFile.append("Missing file: ${file_path}\n")
    log.warn "Missing file: ${file_path} - skipping this sample"
}

// Input channel from samples.txt
Channel
    .fromPath(params.input, checkIfExists: true)
    .splitText()
    .map { line ->
        if (line.trim().startsWith('#')) return  // Skip comment lines
        def fields = line.trim().split(/\s+/)
        if (fields.size() >= 3) {
            def sample_id = fields[0]
            def r1_path = fields[1]
            def r2_path = fields[2]
            
            // Check if files exist without causing pipeline failure
            def r1 = file(r1_path)
            def r2 = file(r2_path)
            
            if (!r1.exists()) {
                logMissingFile(r1_path)
                return null
            }
            
            if (!r2.exists()) {
                logMissingFile(r2_path)
                return null
            }
            
            return [sample_id, [r1, r2]]
        }
        return null
    }
    .filter { it != null }  // Filter out skipped lines
    .set { read_pairs }

// Workflow
workflow {
    // Log the number of samples that will be processed
    if (!params.skip_to_merge) {
        read_pairs.count().subscribe { count ->
            log.info "Processing ${count} samples"
        }
        
        // Trimming
        FASTP(read_pairs)
    
        // Parallel repair on each trimmed pair
        REPAIR(FASTP.out.trimmed_reads)
        
        // Extract just file paths
        REPAIR.out.repaired_reads
            .map { sample_id, reads -> reads }
            .collect()
            .set { all_repaired_reads }
    }
    else {
        // If skipping to merge, collect already repaired files
        Channel
            .fromPath("${params.repaired_reads_dir}/*_repaired_R{1,2}.fastq.gz")
            .map { file -> file }
            .collect()
            .set { all_repaired_reads }
        
        log.info "Skipping trimming and repair. Using existing files in ${params.repaired_reads_dir}"
    }
    
    // Merge repaired reads
    MERGE_READS(all_repaired_reads)
    
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