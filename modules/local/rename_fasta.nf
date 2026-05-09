process RENAME_FASTA {
    tag { genome_fasta.baseName }

    label 'process_single'

    publishDir "${params.outdir}/renamed_genomes", mode: 'copy'

    input:
    path genome_fasta

    output:
    path "${genome_fasta.baseName}.fna"

    script:
    def out = "${genome_fasta.baseName}.fna"
    """
    awk '/^>/ {
        split(\$0, a, " "); name = substr(a[1], 2);
        gsub(/[^A-Za-z0-9_]/, "_", name);
        print ">" name; next
    } { print }' "${genome_fasta}" > "${out}.tmp"
    mv "${out}.tmp" "${out}"
    """
}
