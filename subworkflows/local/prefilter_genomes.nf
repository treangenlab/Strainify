// Prefilter reference genomes with MAGNET.
//
// Replaces the old bwa-mem + `samtools coverage > 0` prefilter. Instead of
// keeping any genome touched by a single read, we run MAGNET against the panel
// and keep only the genomes MAGNET calls Present.
//
// MAGNET is invoked with the SAME settings as the standalone driver script:
//   -m illumina  --no-cluster  --shared-reference <dir>
//   --min-covscore 0.5  --min-breadth 0.005  --min-reads 0
//   --rescue-min-reads 100  --rescue-min-ani 0.97
// (--rescue-min-breadth and --min-mapq are intentionally NOT passed, matching
// the script; every value above is overridable via params -- see below.)
//
// genome_id handling also mirrors the script: each genome_id is the FASTA
// basename with its extension stripped (.gz, then .fna/.fa/.fasta) and any
// character outside [A-Za-z0-9._-] replaced by '_', de-duplicated with a _N
// suffix. That is byte-for-byte what MAGNET emits in 'Assembly Accession ID',
// so the Present-call -> local-file mapping is exact.
//
// Design (matches `magnet --shared-reference`):
//   1. MAGNET_BUILD_INDEX builds ONE shared reference panel + bwa index ONCE
//      into the persistent `params.magnet_ref_dir` -- the single reference
//      index every sample reuses (marker-guarded, so re-running is a no-op).
//   2. MAGNET_FILTER runs magnet for EACH set of fastqs against that SAME
//      shared index, producing a per-sample Present/Absent table.
//   3. MAGNET_SELECT resolves each sample's Present genome_ids back to the
//      actual genome FASTAs (exact genome_id -> path lookup in the genome list)
//      to hand to STRAINIFY as that set's reference genomes.
//
// Tunable params (defaults shown match the driver script):
//   magnet_min_covscore     = 0.5
//   magnet_min_breadth      = 0.005
//   magnet_min_reads        = 0
//   magnet_rescue_min_reads = 100
//   magnet_rescue_min_ani   = 0.97
//   magnet_mode             = 'illumina'
//   magnet_cmd              = "python3 ${projectDir}/bin/magnet/magnet.py"
//                             (the in-repo MAGNET script; override with e.g.
//                              --magnet_cmd magnet if it is installed on PATH)
//   magnet_threads          = task.cpus   // set e.g. --magnet_threads 36
// Optional extras (OFF unless set, the script does not use them):
//   magnet_rescue_min_breadth, magnet_min_mapq, magnet_extra
//
// Emits a PER-SAMPLE channel: tuple(sample, [filtered_genome_files]).

process MAGNET_BUILD_INDEX {
    label 'process_high'

    input:
    path genomes        // collected list of genome FASTAs (.fna/.fa/.fasta[.gz])
    val  magnet_ref_dir // absolute path to the shared reference directory

    output:
    path "index_ready.txt", emit: ready

    script:
    def magnet_cmd  = params.containsKey('magnet_cmd')     && params.magnet_cmd     ? params.magnet_cmd     : "python3 ${projectDir}/bin/magnet/magnet.py"
    def magnet_mode = params.containsKey('magnet_mode')    && params.magnet_mode    ? params.magnet_mode    : 'illumina'
    def threads     = params.containsKey('magnet_threads') && params.magnet_threads ? params.magnet_threads : task.cpus
    """
    set -euo pipefail
    # Build the panel dir + genome_list.tsv with a portable python3 step instead
    # of bash associative arrays, so this runs identically on Linux and on macOS
    # (whose /bin/bash is 3.2 and lacks `declare -A`). genome_id logic is
    # byte-for-byte unchanged: basename, strip (.gz then .fna/.fa/.fasta), sanitize
    # to [A-Za-z0-9._-], de-duplicate with a _N suffix -- exactly what MAGNET emits
    # as 'Assembly Accession ID'. Genomes are copied to stable panel paths.
    python3 - "${magnet_ref_dir}" ${genomes} <<'PY'
import os, re, sys, shutil
ref_dir = sys.argv[1]
genomes = sys.argv[2:]
gdir    = os.path.join(ref_dir, "genomes")
os.makedirs(gdir, exist_ok=True)
seen = {}
with open(os.path.join(ref_dir, "genome_list.tsv"), "w") as gl:
    gl.write("genome_path\\tgenome_id\\n")
    for genome in genomes:
        fname = os.path.basename(genome)
        stem  = fname
        if stem.endswith(".gz"):     stem = stem[:-3]
        if stem.endswith(".fna"):    stem = stem[:-4]
        if stem.endswith(".fa"):     stem = stem[:-3]
        if stem.endswith(".fasta"):  stem = stem[:-6]
        suffix = fname[len(stem):]
        gid = re.sub(r"[^A-Za-z0-9._-]", "_", stem) or "genome"
        if gid in seen:
            seen[gid] += 1
            gid = "%s_%d" % (gid, seen[gid])
        else:
            seen[gid] = 0
        dest = os.path.join(gdir, gid + suffix)
        if not (os.path.exists(dest) and os.path.getsize(dest) > 0):
            shutil.copy(genome, dest)
        gl.write("%s\\t%s\\n" % (dest, gid))
PY

    # Build the shared concatenated reference + bwa index ONCE (idempotent:
    # marker-guarded inside magnet, a no-op if already present).
    ${magnet_cmd} --build-index-only \\
        -g "${magnet_ref_dir}/genome_list.tsv" \\
        --shared-reference "${magnet_ref_dir}" \\
        --no-cluster \\
        -m ${magnet_mode} \\
        --threads ${threads}

    # Readiness sentinel so downstream filtering waits for the index.
    cp "${magnet_ref_dir}/genome_list.tsv" index_ready.txt
    """
}

process MAGNET_FILTER {
    tag { sample }

    label 'process_high'

    publishDir { "${params.outdir}/magnet/${sample}" }, mode: 'copy',
        pattern: '*.present_genomes.txt'
    publishDir { "${params.outdir}/magnet/${sample}" }, mode: 'copy',
        pattern: 'cluster_representative.csv'

    input:
    tuple val(sample), path(reads)
    val  magnet_ref_dir
    path ready          // ordering dependency on MAGNET_BUILD_INDEX

    output:
    tuple val(sample), path("${sample}.present_genomes.txt"), emit: present
    path "cluster_representative.csv"

    script:
    def magnet_cmd  = params.containsKey('magnet_cmd')     && params.magnet_cmd     ? params.magnet_cmd     : "python3 ${projectDir}/bin/magnet/magnet.py"
    def magnet_mode = params.containsKey('magnet_mode')    && params.magnet_mode    ? params.magnet_mode    : 'illumina'
    def threads     = params.containsKey('magnet_threads') && params.magnet_threads ? params.magnet_threads : task.cpus

    // Single- vs paired-end: normalize to a list so detection is robust whether
    // Nextflow hands us a bare path (1 file) or a list. Single-end -> only -i,
    // i.e. the one fastq is mapped to the shared reference index; paired -> -i/-I.
    def read_list = (reads instanceof List) ? reads : [reads]
    def r1 = read_list[0]
    def r2 = (read_list.size() > 1) ? read_list[1] : null
    def read_args = r2 ? "-i ${r1} -I ${r2}" : "-i ${r1}"

    // Present/Absent thresholds. Defaults match the driver script; override
    // any of them via the matching param. (Use != null so an explicit 0 sticks.)
    def min_covscore     = (params.containsKey('magnet_min_covscore')     && params.magnet_min_covscore     != null) ? params.magnet_min_covscore     : 0.5
    def min_breadth      = (params.containsKey('magnet_min_breadth')      && params.magnet_min_breadth      != null) ? params.magnet_min_breadth      : 0.005
    def min_reads        = (params.containsKey('magnet_min_reads')        && params.magnet_min_reads        != null) ? params.magnet_min_reads        : 0
    def rescue_min_reads = (params.containsKey('magnet_rescue_min_reads') && params.magnet_rescue_min_reads != null) ? params.magnet_rescue_min_reads : 100
    def rescue_min_ani   = (params.containsKey('magnet_rescue_min_ani')   && params.magnet_rescue_min_ani   != null) ? params.magnet_rescue_min_ani   : 0.97

    def thr = "--min-covscore ${min_covscore}" +
              " --min-breadth ${min_breadth}" +
              " --min-reads ${min_reads}" +
              " --rescue-min-reads ${rescue_min_reads}" +
              " --rescue-min-ani ${rescue_min_ani}"

    // Extras the driver script does NOT pass; included only if explicitly set.
    if (params.containsKey('magnet_rescue_min_breadth') && params.magnet_rescue_min_breadth != null)
        thr += " --rescue-min-breadth ${params.magnet_rescue_min_breadth}"
    if (params.containsKey('magnet_min_mapq') && params.magnet_min_mapq != null)
        thr += " --min-mapq ${params.magnet_min_mapq}"
    if (params.containsKey('magnet_extra') && params.magnet_extra)
        thr += " ${params.magnet_extra}"
    """
    set -euo pipefail

    ${magnet_cmd} \\
        -g "${magnet_ref_dir}/genome_list.tsv" \\
        --shared-reference "${magnet_ref_dir}" \\
        --no-cluster \\
        -m ${magnet_mode} \\
        ${read_args} \\
        -o magnet_out \\
        ${thr} \\
        --threads ${threads}

    cp magnet_out/cluster_representative.csv cluster_representative.csv

    # Extract genomes MAGNET called Present (by 'Assembly Accession ID', which
    # equals our genome_id). csv module handles the quoted, comma-bearing fields;
    # header lookup is case-insensitive, matching the driver script's parser.
    python3 - <<'PY'
import csv
present = []
with open("cluster_representative.csv", newline="") as fh:
    reader = csv.DictReader(fh)
    fields = {(c or "").strip().lower(): c for c in (reader.fieldnames or [])}
    id_col = fields.get("assembly accession id")
    pa_col = fields.get("presence/absence")
    if id_col is None:
        raise SystemExit("No 'Assembly Accession ID' column in cluster_representative.csv")
    for row in reader:
        if pa_col is not None and str(row.get(pa_col, "")).strip().lower() != "present":
            continue
        gid = str(row.get(id_col, "")).strip()
        if gid:
            present.append(gid)
with open("${sample}.present_genomes.txt", "w") as out:
    out.write("\\n".join(present))
    if present:
        out.write("\\n")
PY
    """
}

process MAGNET_SELECT {
    tag { sample }

    label 'process_single'

    publishDir { "${params.outdir}/magnet/${sample}" }, mode: 'copy',
        pattern: 'filtered_*'

    input:
    tuple val(sample), path(present_file)
    val  magnet_ref_dir

    output:
    tuple val(sample), path("filtered_${sample}/*"), emit: filtered

    script:
    // Resolve each Present genome_id to its stable panel path via exact lookup
    // in the genome list (column 2 == genome_id), mirroring the driver script.
    // We COPY (not symlink) so Nextflow can stage real files into STRAINIFY.
    """
    set -euo pipefail
    mkdir -p "filtered_${sample}"

    while IFS= read -r gid; do
        [[ -n "\$gid" ]] || continue
        genome="\$(awk -F'\\t' -v id="\$gid" 'NR>1 && \$2==id {print \$1; exit}' "${magnet_ref_dir}/genome_list.tsv")"
        if [[ -z "\$genome" || ! -f "\$genome" ]]; then
            echo "WARN: Present genome_id has no local file: \$gid" >&2
            continue
        fi
        cp -f "\$genome" "filtered_${sample}/\$(basename "\$genome")"
    done < "${present_file}"

    n=\$(find "filtered_${sample}" -maxdepth 1 -type f | wc -l)
    if [[ "\$n" -lt 2 ]]; then
        echo "ERROR: not enough genomes to run strainify for sample '${sample}': MAGNET called \$n genome(s) Present, but strainify needs at least 2 reference genomes. Loosen the MAGNET presence thresholds (e.g. --magnet_min_covscore, --magnet_rescue_min_reads, --magnet_rescue_min_ani) or provide more candidate genomes." >&2
        exit 1
    fi
    """
}

workflow PREFILTER_GENOMES {
    take:
    genomes_ch  // channel of individual genome Path objects
    reads_ch    // channel of tuple(val(sample), list(path(reads)))

    main:
    // Persistent shared panel + build-once index, reused across every set. When
    // not set, default to a launch-directory-stable path: the batch driver
    // launches every per-set run from the same cwd, so all sets resolve to the
    // SAME default dir and reuse one index (rather than each rebuilding its own).
    def magnet_ref_dir = file(params.magnet_ref_dir ?: "${launchDir}/magnet_shared_reference").toAbsolutePath().toString()
    if (!params.magnet_ref_dir) {
        log.warn "PREFILTER: --magnet_ref_dir not set; using default '${magnet_ref_dir}'. " +
                 "Set it explicitly to control where the build-once index lives and to " +
                 "guarantee a shared path on distributed/cluster filesystems."
    }

    def collected_genomes = genomes_ch.collect()

    // Build the shared reference panel + bwa index ONCE.
    def built = MAGNET_BUILD_INDEX(collected_genomes, magnet_ref_dir)

    // Filter EACH set of fastqs against that SAME prebuilt index.
    def filtered = MAGNET_FILTER(reads_ch, magnet_ref_dir, built.ready)

    // Resolve each set's Present calls back to genome FASTAs.
    def selected = MAGNET_SELECT(filtered.present, magnet_ref_dir)

    emit:
    // Per-sample reference set: tuple(sample, [filtered_genome_files]).
    per_sample = selected.filtered
}