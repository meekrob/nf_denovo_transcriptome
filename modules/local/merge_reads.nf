process MERGE_READS {
    publishDir "${params.outdir}/merging", mode: params.publish_dir_mode
    cpus 8
    memory '64 GB'
    time '24h'

    input:
    path repaired_reads

    output:
    tuple path("merged_R1.fastq.gz"), path("merged_R2.fastq.gz"), emit: merged_reads

    script:
    """
    # Merge R1 files
    cat *_repaired_R1.fastq.gz > temp_R1.fastq.gz
    pigz -p ${task.cpus} -dc temp_R1.fastq.gz | pigz -p ${task.cpus} > merged_R1.fastq.gz

    # Merge R2 files
    cat *_repaired_R2.fastq.gz > temp_R2.fastq.gz
    pigz -p ${task.cpus} -dc temp_R2.fastq.gz | pigz -p ${task.cpus} > merged_R2.fastq.gz
    """
} 