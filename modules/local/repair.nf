process REPAIR {
    publishDir "${params.outdir}/repair", mode: params.publish_dir_mode
    cpus 16
    memory '128 GB'
    time '24h'

    input:
    path merged_r1
    path merged_r2

    output:
    tuple path("repaired_R1.fastq.gz"), path("repaired_R2.fastq.gz"), emit: repaired_reads
    path "singletons.fastq.gz", emit: singletons

    script:
    """
    repair.sh \\
        in1=${merged_r1} \\
        in2=${merged_r2} \\
        out1=repaired_R1.fastq.gz \\
        out2=repaired_R2.fastq.gz \\
        outs=singletons.fastq.gz \\
        threads=${task.cpus} \\
        -Xmx${task.memory.toGiga()-5}g \\
        verbose=t \\
        overwrite=t

    # Verify read counts match
    r1_count=\$(zcat repaired_R1.fastq.gz | grep -c '^@')
    r2_count=\$(zcat repaired_R2.fastq.gz | grep -c '^@')
    if [ "\$r1_count" -ne "\$r2_count" ]; then
        echo "Error: Read counts mismatch: R1=\$r1_count, R2=\$r2_count" >&2
        exit 1
    fi
    """
} 