process SEQKIT_QC {
    tag "$sample_id"
    publishDir "${params.outdir}/seqkit_qc", mode: params.publish_dir_mode
    cpus 64
    memory '128 GB'
    time '24h'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_clean_R{1,2}.fastq.gz"), emit: cleaned_reads

    script:
    """
    # Run seqkit sana on R1
    seqkit sana -j ${task.cpus} ${reads[0]} -o ${sample_id}_sana_R1.fastq.gz

    # Run seqkit sana on R2
    seqkit sana -j ${task.cpus} ${reads[1]} -o ${sample_id}_sana_R2.fastq.gz

    # Run seqkit pair on sana outputs
    seqkit pair -j ${task.cpus} ${sample_id}_sana_R1.fastq.gz ${sample_id}_sana_R2.fastq.gz -1 ${sample_id}_clean_R1.fastq.gz -2 ${sample_id}_clean_R2.fastq.gz
    """
} 