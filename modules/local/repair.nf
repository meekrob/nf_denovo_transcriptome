process REPAIR {
    tag "$sample_id"
    publishDir "${params.outdir}/repair", mode: params.publish_dir_mode
    cpus 4
    memory '32 GB'
    time '24h'

    input:
    tuple val(sample_id), path(trimmed_reads)

    output:
    tuple val(sample_id), path("${sample_id}_repaired_R{1,2}.fastq.gz"), emit: repaired_reads
    path "${sample_id}_singletons.fastq.gz", emit: singletons

    script:
    """
    repair.sh \\
        in1=${trimmed_reads[0]} \\
        in2=${trimmed_reads[1]} \\
        out1=${sample_id}_repaired_R1.fastq.gz \\
        out2=${sample_id}_repaired_R2.fastq.gz \\
        outs=${sample_id}_singletons.fastq.gz \\
        threads=${task.cpus} \\
        -Xmx${task.memory.toGiga()-5}g \\
        verbose=t \\
        overwrite=t

    # Verify read counts match
    r1_count=\$(zcat ${sample_id}_repaired_R1.fastq.gz | grep -c '^@')
    r2_count=\$(zcat ${sample_id}_repaired_R2.fastq.gz | grep -c '^@')
    if [ "\$r1_count" -ne "\$r2_count" ]; then
        echo "Error: Read counts mismatch for ${sample_id}: R1=\$r1_count, R2=\$r2_count" >&2
        exit 1
    fi
    """
} 