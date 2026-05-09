include { PREFILTER_GENOMES } from './subworkflows/local/prefilter_genomes.nf'
include { STRAINIFY }         from './subworkflows/local/strainify.nf'

workflow {
    main:

    // Input validation
    if (!params.genome_folder) {
        error "Please provide --genome_folder pointing to a directory of reference genome FASTA files (.fna)"
    }
    if (!params.fastq_folder) {
        error "Please provide --fastq_folder pointing to a directory of FASTQ files"
    }

    // ----------------------------------------------------------------
    // Build genome channel
    // ----------------------------------------------------------------
    def genomes_ch = channel.fromPath("${params.genome_folder}/*.{fna,fa,fasta}", checkIfExists: true)

    // ----------------------------------------------------------------
    // Build reads channel:
    //   paired: tuple(sample, [r1, r2])  from *_r{1,2}.fq or *.fq.gz variants
    //   single: tuple(sample, [r])
    // ----------------------------------------------------------------
    def reads_ch
    if (params.read_type == 'paired') {
        // Build a combined channel from lowercase and uppercase R1/R2 patterns,
        // handling both uncompressed and gzip-compressed FASTQ files.
        def r_lower = channel.fromFilePairs("${params.fastq_folder}/*_r{1,2}.fq{,.gz}", flat: false)
        def r_upper = channel.fromFilePairs("${params.fastq_folder}/*_R{1,2}.fq{,.gz}", flat: false)
        reads_ch = r_lower.mix(r_upper)
    } else {
        reads_ch = channel.fromPath(
            "${params.fastq_folder}/*.f{ast,}q{,.gz}",
            checkIfExists: true
        ).filter { p ->
            !(p.name =~ /_(r|R)[12]\./)
        }.map { p -> tuple(p.simpleName, [p]) }
    }

    // ----------------------------------------------------------------
    // Optional genome prefiltering step
    // ----------------------------------------------------------------
    def input_genomes_ch
    if (params.prefilter) {
        input_genomes_ch = PREFILTER_GENOMES(genomes_ch, reads_ch).filtered
    } else {
        input_genomes_ch = genomes_ch
    }

    // ----------------------------------------------------------------
    // Core Strainify workflow
    // ----------------------------------------------------------------
    STRAINIFY(input_genomes_ch, reads_ch)

    onComplete:
    def status = workflow.success ? 'completed successfully' : 'failed'
    log.info "Strainify ${status}. Results in: ${params.outdir}"
}
