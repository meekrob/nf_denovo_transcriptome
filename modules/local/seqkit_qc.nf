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
    # Create temporary directory
    mkdir -p temp

    # Debug - list input files
    ls -la
    echo "Input files: ${reads[0]}, ${reads[1]}"

    # Run seqkit sana on R1 - note order of parameters
    seqkit sana -j ${task.cpus} -o temp/${sample_id}_sana_R1.fastq.gz ${reads[0]}

    # Run seqkit sana on R2 - note order of parameters
    seqkit sana -j ${task.cpus} -o temp/${sample_id}_sana_R2.fastq.gz ${reads[1]}

    # Create output directory for paired reads
    mkdir -p temp/paired

    # Run seqkit pair on sana outputs
    seqkit pair -j ${task.cpus} \
        -1 temp/${sample_id}_sana_R1.fastq.gz \
        -2 temp/${sample_id}_sana_R2.fastq.gz \
        -O temp/paired -u

    # Move the outputs to final filenames
    mv temp/paired/${sample_id}_sana_R1.fastq.gz ${sample_id}_clean_R1.fastq.gz
    mv temp/paired/${sample_id}_sana_R2.fastq.gz ${sample_id}_clean_R2.fastq.gz

    # Clean up
    rm -rf temp
    """
} 