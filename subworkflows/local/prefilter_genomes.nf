// Prefilter genomes by mapping reads and removing genomes with zero coverage.
// Replaces filter_genomes.sh.

process PREFILTER_CONCAT_INDEX {
    label 'process_medium'

    publishDir "${params.outdir}/prefilter", mode: 'copy', pattern: 'contig_to_genome.tsv'

    input:
    path genomes  // collected list of genome FASTA files (.fa / .fna / .fasta)

    output:
    path "concat_genome.fna",    emit: concat_fasta
    path "concat_genome.fna.*",  emit: index_files
    path "contig_to_genome.tsv", emit: contig_map
    path "all_input_genomes.txt",emit: genome_list

    script:
    """
    mkdir -p renamed
    : > contig_to_genome.tsv
    : > all_input_genomes.txt

    for genome in ${genomes}; do
        genome_name="\$(basename "\$genome" | sed 's/\\.[^.]*\$//')"
        echo "\$genome_name" >> all_input_genomes.txt
        awk -v g="\$genome_name" -v map="contig_to_genome.tsv" '
            /^>/{
                split(substr(\$0,2), a, " ")
                contig = a[1]
                new = g "|" contig
                print new "\\t" g >> map
                print ">" new
                next
            }
            { print }
        ' "\$genome" > "renamed/\${genome_name}.renamed.fna"
    done

    sort -u all_input_genomes.txt -o all_input_genomes.txt
    cat renamed/*.fna > concat_genome.fna
    bwa index concat_genome.fna
    samtools faidx concat_genome.fna
    """
}

process PREFILTER_MAP_SAMPLE {
    tag { sample }

    label 'process_high'

    input:
    tuple val(sample), path(reads)
    path concat_fasta
    path index_files

    output:
    path "${sample}.covered_genomes.txt", emit: covered

    script:
    def read_args = reads.size() == 2
        ? "${reads[0]} ${reads[1]}"
        : "${reads[0]}"
    """
    bwa mem -t ${task.cpus} "${concat_fasta}" ${read_args} \\
        | samtools view -b -F256 -F2048 -F4 -q60 \\
        | samtools sort -o "${sample}.bam" --write-index -@ ${task.cpus}

    samtools coverage "${sample}.bam" \\
        | awk 'NR>1 && \$5>0 { split(\$1,a,"|"); print a[1] }' \\
        > "${sample}.covered_genomes.txt" || true
    """
}

process PREFILTER_SELECT_GENOMES {
    label 'process_single'

    publishDir "${params.outdir}/prefilter", mode: 'copy'

    input:
    path genome_list
    path covered_files  // collected list of *.covered_genomes.txt
    path original_genomes  // collected list of original genome files

    output:
    path "filtered_genomes/*.{fna,fa,fasta}", emit: filtered, optional: true
    path "genomes_with_zero_coverage_across_all_samples.txt", emit: zero_coverage

    script:
    """
    mkdir -p filtered_genomes

    sort -u ${covered_files} > any_covered.txt

    comm -23 <(sort "${genome_list}") <(sort any_covered.txt) \\
        > genomes_with_zero_coverage_across_all_samples.txt || true

    while read -r g; do
        [[ -n "\$g" ]] || continue
        for genome in ${original_genomes}; do
            base="\$(basename "\$genome" | sed 's/\\.[^.]*\$//')"
            if [[ "\$base" == "\$g" ]]; then
                cp "\$genome" filtered_genomes/
                break
            fi
        done
    done < any_covered.txt

    n=\$(ls filtered_genomes/*.fna filtered_genomes/*.fa filtered_genomes/*.fasta 2>/dev/null | wc -l || echo 0)
    if [[ "\$n" -eq 0 ]]; then
        echo "ERROR: Genome prefiltering removed all genomes. None of the input genomes appear to be present in the metagenomic data." >&2
        exit 1
    fi
    """
}

workflow PREFILTER_GENOMES {
    take:
    genomes_ch  // channel of individual genome Path objects
    reads_ch    // channel of tuple(val(sample), list(path(reads)))

    main:
    def collected_genomes = genomes_ch.collect()

    def ci = PREFILTER_CONCAT_INDEX(collected_genomes)

    def mapped = PREFILTER_MAP_SAMPLE(reads_ch, ci.concat_fasta, ci.index_files.collect())

    def selected = PREFILTER_SELECT_GENOMES(
        ci.genome_list,
        mapped.covered.collect(),
        collected_genomes
    )

    emit:
    filtered = selected.filtered.flatten()
}
