process SEQKIT_QC {
    tag "$sample_id"
    publishDir "${params.outdir}/seqkit_qc", mode: params.publish_dir_mode
    cpus 4
    memory '16 GB'
    time '24h'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_clean_R{1,2}.fastq.gz"), emit: cleaned_reads

    script:
    """
    # Create temporary directory for processing
    mkdir -p temp_dir

    # Run seqkit sana on R1
    seqkit sana -j ${task.cpus} ${reads[0]} -o temp_dir/${sample_id}_sana_R1.fastq.gz

    # Run seqkit sana on R2
    seqkit sana -j ${task.cpus} ${reads[1]} -o temp_dir/${sample_id}_sana_R2.fastq.gz

    # Create output directory for paired reads
    mkdir -p temp_dir/paired_output

    # Run seqkit pair on sana outputs
    seqkit pair -j ${task.cpus} \
        -1 temp_dir/${sample_id}_sana_R1.fastq.gz \
        -2 temp_dir/${sample_id}_sana_R2.fastq.gz \
        -O temp_dir/paired_output -u

    # Move the outputs to the work directory
    mv temp_dir/paired_output/${sample_id}_sana_R1.fastq.gz ${sample_id}_clean_R1.fastq.gz
    mv temp_dir/paired_output/${sample_id}_sana_R2.fastq.gz ${sample_id}_clean_R2.fastq.gz

    # Clean up
    rm -rf temp_dir
    """
} 