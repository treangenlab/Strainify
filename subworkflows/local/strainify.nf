include { RENAME_FASTA }                                  from '../../modules/local/rename_fasta.nf'
include { RUN_PARSNP; MAF2VCF; FILTER_VARIANTS; GET_REF } from '../../modules/local/parsnp.nf'
include { FAIDX_REF; BWA_INDEX }                          from '../../modules/local/reference.nf'
include { MAP_READS; FILTER_SORT_INDEX; COUNT_READS }     from '../../modules/local/mapping.nf'
include { COMPUTE_ABUNDANCES }                            from '../../modules/local/abundance.nf'

workflow STRAINIFY {
    take:
    genomes_ch  // channel of individual genome Path objects
    reads_ch    // channel of tuple(val(sample), list(path(reads)))

    main:

    // ----------------------------------------------------------------
    // Branch: full run vs. precomputed variants
    // ----------------------------------------------------------------
    def variant_matrix_ch
    def sites_ch
    def ref_ch
    def fai_ch
    def bwa_index_ch  // tuple: (ref, bwt, pac, ann, amb, sa)

    if (!params.use_precomputed_variants) {
        // --- Step 1: rename genome headers ---
        def renamed_ch = RENAME_FASTA(genomes_ch)

        // --- Step 2: parsnp whole-genome alignment ---
        // Sort the collected genomes by file name (not by Path, whose string
        // includes the random work-dir hash) so parsnp receives a deterministic,
        // reproducible input set/order across runs.
        def parsnp_out = RUN_PARSNP(renamed_ch.toSortedList { a, b -> a.name <=> b.name })

        // --- Step 3: MAF → VCF ---
        def maf2vcf_out = MAF2VCF(parsnp_out.maf)

        // --- Step 4: filter recombinant windows ---
        def fv_out = FILTER_VARIANTS(parsnp_out.maf, maf2vcf_out.vcf)
        variant_matrix_ch = fv_out.variant_matrix
        sites_ch          = fv_out.sites

        // --- Step 5: extract reference FASTA ---
        def ref_out = GET_REF(parsnp_out.ref_file)
        ref_ch = ref_out.ref
    } else {
        // Use precomputed outputs from a previous run
        if (!params.precomputed_dir) {
            error "--precomputed_dir must be set when --use_precomputed_variants is true"
        }
        def pdir = params.precomputed_dir
        // Wrap as VALUE channels (channel.value), not channel.fromPath (which is a
        // single-item QUEUE that gets consumed after one sample). Value channels
        // are reused for every sample, matching the non-precomputed branch -- so
        // mapping/counting fan out across all samples instead of just the first.
        variant_matrix_ch = channel.value(file("${pdir}/filtered_variant_matrix.csv", checkIfExists: true))
        sites_ch          = channel.value(file("${pdir}/sites.txt",                   checkIfExists: true))
        ref_ch            = channel.value(file("${pdir}/reference.fna",               checkIfExists: true))
    }

    // --- Step 6: index the reference ---
    def faidx_out  = FAIDX_REF(ref_ch)
    def bwaidx_out = BWA_INDEX(faidx_out.ref)

    // --- Step 7: per-sample mapping ---
    def sam_ch = MAP_READS(
        reads_ch,
        bwaidx_out.ref,
        bwaidx_out.bwt,
        bwaidx_out.pac,
        bwaidx_out.ann,
        bwaidx_out.amb,
        bwaidx_out.sa
    )

    // --- Step 8: filter, sort, index BAM ---
    def bam_ch = FILTER_SORT_INDEX(sam_ch.sam)

    // tuple(sample, bam, bai) for count_reads
    def bam_with_bai = bam_ch.bam.join(bam_ch.bai.map { sample, bai -> tuple(sample, bai) })

    // --- Step 9: per-site read counting ---
    def counts_ch = COUNT_READS(
        bam_with_bai,
        faidx_out.ref,
        faidx_out.fai,
        sites_ch
    )

    // --- Step 10: compute strain abundances ---
    COMPUTE_ABUNDANCES(
        counts_ch.read_counts.collect(),
        variant_matrix_ch
    )

    emit:
    abundances = COMPUTE_ABUNDANCES.out.abundances
}