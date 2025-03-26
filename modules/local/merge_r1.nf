process MERGE_R1 {
    publishDir "${params.outdir}/merging", mode: params.publish_dir_mode
    cpus 8
    memory '64 GB'
    time '24h'

    input:
    path cleaned_reads_r1

    output:
    path "merged_R1.fastq.gz", emit: merged_r1

    script:
    """
    # Merge R1 files
    cat ${cleaned_reads_r1} > temp_R1.fastq.gz
    pigz -p ${task.cpus} -dc temp_R1.fastq.gz | pigz -p ${task.cpus} > merged_R1.fastq.gz
    """
} 