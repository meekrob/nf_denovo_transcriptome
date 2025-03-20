process MERGE_R2 {
    publishDir "${params.outdir}/merging", mode: params.publish_dir_mode
    cpus 8
    memory '64 GB'
    time '24h'

    input:
    path trimmed_reads_r2

    output:
    path "merged_R2.fastq.gz", emit: merged_r2

    script:
    """
    # Merge R2 files
    cat ${trimmed_reads_r2} > temp_R2.fastq.gz
    pigz -p ${task.cpus} -dc temp_R2.fastq.gz | pigz -p ${task.cpus} > merged_R2.fastq.gz
    """
} 