process COMPUTE_ABUNDANCES {
    label 'process_high'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path read_counts_files  // collected list of *_read_counts.tsv
    path variant_matrix

    output:
    path "abundance_estimates_combined.csv",      emit: abundances
    path "abundance_estimates_bootstrap_CIs.csv", emit: bootstrap_cis, optional: true

    script:
    def weight_flag     = params.weight_by_entropy  ? '--weight_by_entropy'  : ''
    def bootstrap_flag  = params.bootstrap          ? '--bootstrap'          : ''
    def bootstrap_iters = params.bootstrap ? "--bootstrap_iterations ${params.bootstrap_iterations}" : ''
    // Gather all read count files into a dedicated subdirectory so the
    // script's glob of *_read_counts.tsv works correctly.
    """
    mkdir -p read_counts
    for f in ${read_counts_files}; do
        ln -sf "\$PWD/\$f" read_counts/
    done
    compute_abundances.py \\
        --read_counts_dir read_counts \\
        --filtered_variants "${variant_matrix}" \\
        --threads ${task.cpus} \\
        ${weight_flag} \\
        ${bootstrap_flag} \\
        ${bootstrap_iters}
    mv read_counts/../abundance_estimates_combined.csv . 2>/dev/null || true
    """
}
