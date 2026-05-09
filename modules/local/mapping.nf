process MAP_READS {
    tag { sample }

    label 'process_high'

    publishDir "${params.outdir}/mapped_reads", mode: 'copy', pattern: '*.sam'

    input:
    tuple val(sample), path(reads)
    path ref
    path bwt
    path pac
    path ann
    path amb
    path sa

    output:
    tuple val(sample), path("${sample}.sam"), emit: sam

    script:
    def read_args = reads.size() == 2
        ? "${reads[0]} ${reads[1]}"
        : "${reads[0]}"
    """
    bwa mem -t ${task.cpus} "${ref}" ${read_args} > "${sample}.sam"
    """
}

process FILTER_SORT_INDEX {
    tag { sample }

    label 'process_high'

    publishDir "${params.outdir}/mapped_reads", mode: 'copy'

    input:
    tuple val(sample), path(sam)

    output:
    tuple val(sample), path("${sample}_sorted.bam"),     emit: bam
    tuple val(sample), path("${sample}_sorted.bam.csi"), emit: bai

    script:
    """
    samtools view "${sam}" -b -F256 -F2048 -F4 -q60 \\
        | samtools sort -o "${sample}_sorted.bam" --write-index -@ ${task.cpus}
    """
}

process COUNT_READS {
    tag { sample }

    label 'process_high'

    publishDir "${params.outdir}/read_counts", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path ref
    path fai
    path sites

    output:
    path "${sample}_read_counts.tsv", emit: read_counts

    script:
    """
    count_reads.py \\
        --bam       "${bam}" \\
        --ref       "${ref}" \\
        --positions "${sites}" \\
        --output    "\$PWD" \\
        --threads   ${task.cpus}
    """
}
