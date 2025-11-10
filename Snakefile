import os
import glob
import shutil
import re


def rename_header(header):
    return re.sub(r'[^\w]', '_', header.split()[0])


def create_renamed_fasta(in_path, out_path):
    with open(in_path, 'r') as infile, open(out_path, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                sanitized = rename_header(line[1:].strip())
                outfile.write(f">{sanitized}\n")
            else:
                outfile.write(line)


def find_ref_file(output_dir):
    search_dir = os.path.join(output_dir, "parsnp_results")
    for f in os.listdir(search_dir):
        if f.endswith(".ref"):
            return os.path.join(search_dir, f)
    raise ValueError(f"No .ref file found in {search_dir}.")


# Load config values
genome_folder = config["genome_folder"]
fastq_folder = config["fastq_folder"]
output_dir = config["output_dir"]

# Automatically detect genome and sample names
raw_fasta_files = glob.glob(os.path.join(genome_folder, "*.fna"))
GENOME_IDS = [os.path.splitext(os.path.basename(f))[0] for f in raw_fasta_files]

read_type = config.get("read_type", "paired")
use_precomputed_variants = config.get("use_precomputed_variants", False)

# Allow user to define a separate directory when using precomputed variants
precomputed_output_dir = config.get("precomputed_output_dir", f"{output_dir}/precomputed_results")
original_dir = output_dir
output_dir = precomputed_output_dir if use_precomputed_variants else original_dir

# Define benchmark directory under output_dir
#BENCH_DIR = os.path.join(output_dir, "benchmarks")


if read_type == "paired":
    fastq_files = glob.glob(os.path.join(fastq_folder, "*_r1.fq"))
    SAMPLES = [os.path.basename(f).replace("_r1.fq", "") for f in fastq_files]
elif read_type == "single":
    fastq_files = glob.glob(os.path.join(fastq_folder, "*.fq"))
    SAMPLES = [os.path.basename(f).replace(".fq", "") for f in fastq_files]


rule all:
    input:
        expand(os.path.join(original_dir, "renamed_genomes", "{sample1}.fna"), sample1=GENOME_IDS),
        expand(os.path.join(output_dir, "mapped_reads/{sample}.sam"), sample=SAMPLES),
        expand(os.path.join(output_dir, "mapped_reads/{sample}_sorted.bam"), sample=SAMPLES),
        expand(os.path.join(output_dir, "read_counts/{sample}_read_counts.tsv"), sample=SAMPLES),
        os.path.join(output_dir, "abundance_estimates_combined.csv"),
        *([] if use_precomputed_variants else [
            os.path.join(output_dir,"abundance_estimates_combined.csv"),
            os.path.join(output_dir, "significantly_enriched_windows.tsv"),
            os.path.join(output_dir, "filtered_variant_matrix.csv"),
            os.path.join(output_dir, "sites.txt"),
        ])

if not use_precomputed_variants:
    rule make_output_dir:
        output:
            directory(output_dir)
        # benchmark:
        #     os.path.join(BENCH_DIR, "make_output_dir.tsv")
        shell:
            "mkdir -p {output}"

    rule rename_fasta_headers:
        input:
            fasta=lambda wildcards: os.path.join(genome_folder, f"{wildcards.sample1}.fna")
        output:
            renamed = os.path.join(output_dir, "renamed_genomes", "{sample1}.fna")
        # benchmark:
        #     os.path.join(BENCH_DIR, "rename_fasta_headers/{sample1}.tsv")
        run:
            os.makedirs(os.path.dirname(output.renamed), exist_ok=True)
            create_renamed_fasta(input.fasta, output.renamed)

    rule run_parsnp:
        input:
            genomes = expand(os.path.join(output_dir, "renamed_genomes", "{sample1}.fna"), sample1=GENOME_IDS)
        output:
            maf = os.path.join(output_dir, "parsnp_results", "parsnp.maf")
        params:
            outdir = os.path.join(output_dir, "parsnp_results"),
            genome_dir = os.path.join(output_dir, "renamed_genomes"),
            parsnp_flags = lambda wildcards: config.get("parsnp_flags", "-c")
        threads: 12
        # benchmark:
        #     os.path.join(BENCH_DIR, "run_parsnp.tsv")
        shell:
            """
            mkdir -p {params.outdir}
            parsnp -r ! -o {params.outdir} {params.parsnp_flags} -p {threads} --fo -d {params.genome_dir}
            """

    rule maf2vcf:
        input:
            maf = rules.run_parsnp.output.maf
        output:
            vcf = os.path.join(output_dir, "parsnp_results", "merged.vcf")
        # benchmark:
        #     os.path.join(BENCH_DIR, "maf2vcf.tsv")
        shell:
            "bash maf2vcf.sh {input.maf} && mv merged.vcf {output.vcf}"

    rule filter_variants:
        input:
            maf = rules.run_parsnp.output.maf,
            vcf = rules.maf2vcf.output.vcf
        params:
            modify_windows = lambda wildcards: config.get("modify_windows", "")
        output:
            recomb_windows = os.path.join(output_dir, "significantly_enriched_windows.tsv"),
            filtered_variant_matrix = os.path.join(output_dir, "filtered_variant_matrix.csv"),
            sites = os.path.join(output_dir, "sites.txt")
        # benchmark:
        #     os.path.join(BENCH_DIR, "filter_variants.tsv")
        shell:
            "python filter_variants_v2.py --maf {input.maf} --vcf {input.vcf} --output_dir {output_dir} {params.modify_windows}"

    rule get_ref:
        input:
            maf = rules.run_parsnp.output.maf
        output:
            ref = os.path.join(output_dir, "reference.fna")
        # benchmark:
        #     os.path.join(BENCH_DIR, "get_ref.tsv")
        run:
            ref_path = find_ref_file(output_dir)
            shutil.copy(ref_path, output.ref)

    rule faidx_ref:
        input:
            ref = rules.get_ref.output.ref
        output:
            fai = rules.get_ref.output.ref + ".fai"
        # benchmark:
        #     os.path.join(BENCH_DIR, "faidx_ref.tsv")
        shell:
            "samtools faidx {input.ref}"

else:
    rule get_ref:
        input:
            filtered_variant_matrix = os.path.join(original_dir, "filtered_variant_matrix.csv"),
            original_ref = os.path.join(original_dir, "reference.fna")
        output:
            ref=os.path.join(output_dir, "reference.fna")
        # benchmark:
        #     os.path.join(BENCH_DIR, "get_ref_precomputed.tsv")
        run:
            shutil.copy(input.original_ref,output.ref)

    rule faidx_ref:
        input:
            ref = rules.get_ref.output.ref
        output:
            fai = rules.get_ref.output.ref + ".fai"
        # benchmark:
        #     os.path.join(BENCH_DIR, "faidx_ref_precomputed.tsv")
        shell:
            "samtools faidx {input.ref}"

rule bwa_index:
    input:
        ref = rules.get_ref.output.ref
    output:
        bwt = rules.get_ref.output.ref + ".bwt",
        pac = rules.get_ref.output.ref + ".pac",
        ann = rules.get_ref.output.ref + ".ann",
        amb = rules.get_ref.output.ref + ".amb",
        sa  = rules.get_ref.output.ref + ".sa"
    # benchmark:
    #     os.path.join(BENCH_DIR, "bwa_index.tsv")
    shell:
        "bwa index {input.ref}"


rule map_reads:
    input:
        ref = rules.get_ref.output.ref,
        bwt = rules.bwa_index.output.bwt
    output:
        sam = os.path.join(output_dir, "mapped_reads", "{sample}.sam")
    threads: 12
    params:
        fastq = lambda wildcards: (
            f"{os.path.join(fastq_folder, wildcards.sample)}_r1.fq {os.path.join(fastq_folder, wildcards.sample)}_r2.fq"
            if read_type == "paired"
            else f"{os.path.join(fastq_folder, wildcards.sample)}.fq"
        )
    # benchmark:
    #     os.path.join(BENCH_DIR, "map_reads/{sample}.tsv")
    shell:
        """
        mkdir -p $(dirname {output.sam})
        bwa mem -t {threads} {input.ref} {params.fastq} > {output.sam}
        """


rule filter_and_index_sam:
    input:
        sam = rules.map_reads.output.sam
    output:
        sorted_bam = os.path.join(output_dir, "mapped_reads", "{sample}_sorted.bam")
    threads: 8
    # benchmark:
    #     os.path.join(BENCH_DIR, "filter_and_index_sam/{sample}.tsv")
    shell:
        "samtools view {input.sam} -b -F256 -F2048 -F4 -q60 | samtools sort -o {output.sorted_bam} --write-index -@ {threads}"


rule count_reads:
    input:
        bam = rules.filter_and_index_sam.output.sorted_bam,
        ref = rules.get_ref.output.ref,
        fai = rules.faidx_ref.output.fai,
        positions = (
            os.path.join(original_dir, "sites.txt")
            if use_precomputed_variants
            else rules.filter_variants.output.sites
        )
    output:
        read_counts = os.path.join(output_dir, "read_counts", "{sample}_read_counts.tsv")
    threads: 12
    # benchmark:
    #     os.path.join(BENCH_DIR, "count_reads/{sample}.tsv")
    shell:
        """
        mkdir -p $(dirname {output.read_counts})
        python {workflow.basedir}/count_reads_parallelized_v2.py \
            --bam $(realpath {input.bam}) \
            --ref $(realpath {input.ref}) \
            --positions $(realpath {input.positions}) \
            --output $(realpath $(dirname {output.read_counts})) \
            --threads {threads}
        """


rule finalize_read_counts:
    input:
        expand(os.path.join(output_dir,"read_counts","{sample}_read_counts.tsv"),sample=SAMPLES)
    output:
        marker=os.path.join(output_dir,"read_counts","done.txt")
    # benchmark:
    #     os.path.join(BENCH_DIR, "finalize_read_counts.tsv")
    shell:
        "touch {output.marker}"


rule compute_abundances:
    input:
        marker = rules.finalize_read_counts.output.marker,
        filtered_variant_matrix = (
            os.path.join(original_dir, "filtered_variant_matrix.csv")
            if use_precomputed_variants
            else rules.filter_variants.output.filtered_variant_matrix
        )
    output:
        abundances = os.path.join(output_dir, "abundance_estimates_combined.csv")
    params:
        read_count_dir = os.path.join(output_dir, "read_counts"),
        weight_flag = "--weight_by_entropy" if config.get("weight_by_entropy", False) else ""
    # benchmark:
    #     os.path.join(BENCH_DIR, "compute_abundances.tsv")
    shell:
        "python compute_abundances_all.py --read_counts_dir {params.read_count_dir} --filtered_variants {input.filtered_variant_matrix} {params.weight_flag}"


