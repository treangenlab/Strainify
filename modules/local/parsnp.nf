process RUN_PARSNP {
    label 'process_high'

    publishDir "${params.outdir}/parsnp_results", mode: 'copy'

    input:
    path renamed_genomes  // staged .fna files in task working dir

    output:
    path "parsnp.maf", emit: maf
    path "*.ref",      emit: ref_file
    path "parsnp.log", emit: log
    path "parsnp_out", emit: parsnp_dir

    script:
    def flags = params.parsnp_flags ?: '-c'
    """
    set -o pipefail
    parsnp -r ! -o parsnp_out ${flags} -p ${task.cpus} --fo -d ./*.fna 2>&1 | tee parsnp.log
    cp parsnp_out/parsnp.maf ./parsnp.maf
    find parsnp_out -name '*.ref' -exec cp {} . ';'
    """
}

process MAF2VCF {
    label 'process_medium'

    publishDir "${params.outdir}/parsnp_results", mode: 'copy'

    input:
    path maf

    output:
    path "merged.vcf", emit: vcf

    script:
    """
    maf2vcf.sh "${maf}"
    """
}

process FILTER_VARIANTS {
    label 'process_single'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path maf
    path vcf

    output:
    path "significantly_enriched_windows.tsv", emit: windows
    path "filtered_variant_matrix.csv",        emit: variant_matrix
    path "sites.txt",                          emit: sites

    script:
    def window_args = "--window_size ${params.window_size} --window_overlap ${params.window_overlap}"
    def filter_arg  = params.filter_off ? '--filter_off' : ''
    """
    filter_variants.py \
        --maf  "${maf}" \
        --vcf  "${vcf}" \
        --output_dir "\$PWD" \
        ${window_args} \
        ${filter_arg}
    """
}

process GET_REF {
    label 'process_single'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path ref_file   // the *.ref file produced by parsnp

    output:
    path "reference.fna", emit: ref

    script:
    """
    cp "${ref_file}" reference.fna
    """
}