process FAIDX_REF {
    label 'process_single'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path ref

    output:
    path "${ref}.fai", emit: fai
    path ref,          emit: ref  // pass through so callers can access both

    script:
    """
    samtools faidx "${ref}"
    """
}

process BWA_INDEX {
    label 'process_medium'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path ref

    output:
    path "${ref}.bwt", emit: bwt
    path "${ref}.pac", emit: pac
    path "${ref}.ann", emit: ann
    path "${ref}.amb", emit: amb
    path "${ref}.sa",  emit: sa
    path ref,          emit: ref  // pass-through for map_reads

    script:
    """
    bwa index "${ref}"
    """
}
