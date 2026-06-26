include { PREFILTER_GENOMES } from './subworkflows/local/prefilter_genomes.nf'
include { STRAINIFY }         from './subworkflows/local/strainify.nf'

workflow {
    main:

    // ----------------------------------------------------------------
    // Input validation
    // ----------------------------------------------------------------
    if (!params.genome_folder) {
        error "Please provide --genome_folder pointing to a directory of reference genome FASTA files (.fna)"
    }

    // ----------------------------------------------------------------
    // Build genome channel (the full candidate panel)
    // ----------------------------------------------------------------
    def genomes_ch = channel.fromPath("${params.genome_folder}/*.{fna,fa,fasta,fna.gz,fa.gz,fasta.gz}", checkIfExists: true)

    if (params.prefilter) {
        // ================================================================
        // FILTER-RUN MODE (MAGNET prefilter -> per-set STRAINIFY)
        //
        // One set of fastqs per pipeline run. MAGNET filters this set against a
        // shared, build-once bwa index of the reference panel (params.magnet_ref_dir),
        // and STRAINIFY runs ONCE on this set using only its own filtered genomes.
        //
        // To process many sets, loop the pipeline once per set against the SAME
        // --magnet_ref_dir (the index is built on the first run and reused by the
        // rest). See bin/run_filter_batch.sh.
        // ================================================================
        // magnet_ref_dir is optional: PREFILTER_GENOMES defaults it to a
        // launch-directory-stable path when unset (and warns). Set it explicitly
        // to control the build-once index location / for distributed filesystems.
        if (!params.fastq1) {
            error """\
                |--prefilter expects ONE set of fastqs per run, given explicitly:
                |    --fastq1 <R1> [--fastq2 <R2>] [--sample_id <name>]
                |This is required because strainify.nf builds one reference per run, so
                |each set of fastqs gets its own filtered genomes and its own strainify run.
                |To process a whole folder, use bin/run_filter_batch.sh, which loops this
                |pipeline once per set sharing the same --magnet_ref_dir.
                |""".stripMargin()
        }

        // Single-set reads channel. Single-end when --fastq2 is omitted: only
        // --fastq1 is provided and MAGNET maps that one fastq to the shared index.
        def sid   = params.sample_id ?: file(params.fastq1).simpleName
        def reads = params.fastq2 ? [ file(params.fastq1), file(params.fastq2) ]
                                  : [ file(params.fastq1) ]
        def reads_ch = channel.of( tuple(sid, reads) )

        log.info "MAGNET prefilter: ${params.fastq2 ? 'paired-end' : 'single-end'} set '${sid}'"

        // MAGNET prefilter -> per-sample filtered genome set.
        def pf = PREFILTER_GENOMES(genomes_ch, reads_ch)

        // Re-attach this set's reads and split into the inputs STRAINIFY expects:
        //   genomes_ch : individual filtered genome Paths (this set's references)
        //   reads_ch   : tuple(sample, [reads]) for this same set
        def joined = pf.per_sample.join(reads_ch)   // (sample, [filtered_files], [reads])

        def filtered_genomes_ch = joined.flatMap { sample, files, rds ->
            (files instanceof List) ? files : [files]
        }
        def filtered_reads_ch = joined.map { sample, files, rds ->
            tuple(sample, rds)
        }

        // STRAINIFY (unchanged) on this set with its own filtered genomes.
        STRAINIFY(filtered_genomes_ch, filtered_reads_ch)
    }
    else {
        // ================================================================
        // ORIGINAL MODE (no prefilter): one reference from all genomes,
        // every sample in --fastq_folder mapped against it.
        // ================================================================
        if (!params.fastq_folder) {
            error "Please provide --fastq_folder pointing to a directory of FASTQ files (or use --prefilter with --fastq1)."
        }

        // Eagerly inspect the folder so a read_type/layout mismatch fails LOUDLY
        // here, instead of silently yielding an empty reads channel and a
        // near-empty "successful" run. file(glob) is evaluated now; with no
        // checkIfExists a non-match is an empty list rather than an error.
        // (Normalize inline -- the strict parser forbids calling a local closure.)
        def _fq       = file("${params.fastq_folder}/*.f{ast,}q{,.gz}")
        def all_fq    = (_fq instanceof List) ? _fq : (_fq != null ? [_fq] : [])
        def paired_fq = all_fq.findAll { p -> (p.name =~ /_(r|R)[12]\./) }
        def single_fq = all_fq.findAll { p -> !(p.name =~ /_(r|R)[12]\./) }

        def reads_ch
        if (params.read_type == 'paired') {
            if (paired_fq.isEmpty()) {
                def hint = single_fq ? " Found ${single_fq.size()} file(s) that look single-end (no _R1/_R2) -- did you mean --read_type single?" : ""
                error "No paired-end FASTQs found in '${params.fastq_folder}' (expected *_R1/_R2 or *_r1/_r2 with .fq/.fastq[.gz]).${hint}"
            }
            def r_lower = channel.fromFilePairs("${params.fastq_folder}/*_r{1,2}.f{ast,}q{,.gz}", flat: false)
            def r_upper = channel.fromFilePairs("${params.fastq_folder}/*_R{1,2}.f{ast,}q{,.gz}", flat: false)
            reads_ch = r_lower.mix(r_upper)
        } else {
            if (single_fq.isEmpty()) {
                def hint = paired_fq ? " Found ${paired_fq.size()} file(s) that look paired-end (_R1/_R2) -- did you mean --read_type paired (the default)?" : ""
                error "No single-end FASTQs found in '${params.fastq_folder}' (expected .fq/.fastq[.gz] without _R1/_R2).${hint}"
            }
            reads_ch = channel.fromPath(
                "${params.fastq_folder}/*.f{ast,}q{,.gz}",
                checkIfExists: true
            ).filter { p ->
                !(p.name =~ /_(r|R)[12]\./)
            }.map { p -> tuple(p.simpleName, [p]) }
        }

        STRAINIFY(genomes_ch, reads_ch)
    }
}