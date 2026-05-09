# Strainify
<p align="center">
  <img src="images/Strainify_logo.png" alt="Strainify logo" width="50%">
</p>

Strainify is an accurate strain-level abundance analysis tool for short-read metagenomics.

## Strainify Workflow
<p align="center">
  <img src="images/Strainify_workflow.png" alt="Strainify workflow diagram" width="100%">
</p>

## Installation

```bash
git clone https://github.com/treangenlab/Strainify.git
cd Strainify

# Create the conda environment (Linux / macOS)
conda env create -f environment.yml
conda activate strainify

# Run directly from the repository root
./strainify --help
```

---

## Usage

All parameters are passed as Nextflow `--param` flags. No YAML config file is required, though
you can use Nextflow's built-in `-params-file params.yml` for reproducibility.

### Basic run (paired-end reads)

```bash
strainify \
  --genome_folder path/to/genomes \
  --fastq_folder  path/to/fastqs \
  --outdir        results \
  -profile conda
```

### Single-end reads

```bash
strainify \
  --genome_folder path/to/genomes \
  --fastq_folder  path/to/fastqs \
  --read_type     single \
  --outdir        results
```

### Pre-filter genomes before running

If you suspect some reference genomes are absent from the metagenome, use `prefilter-run` to
remove zero-coverage genomes first:

```bash
strainify prefilter-run \
  --genome_folder path/to/genomes \
  --fastq_folder  path/to/fastqs \
  --outdir        results
```

Filtered genomes are written to `results/prefilter/filtered_genomes/`.

### Use precomputed variants (re-run with new samples)

If you already have `filtered_variant_matrix.csv`, `reference.fna`, and `sites.txt` from a
previous run, skip the parsnp step:

```bash
strainify \
  --genome_folder          path/to/genomes \
  --fastq_folder           path/to/new_fastqs \
  --use_precomputed_variants \
  --precomputed_dir        path/to/previous_results \
  --outdir                 new_results
```

---

## Parameters

| Parameter | Description | Default |
|---|---|---|
| `--genome_folder` | **(required)** Directory of reference genome FASTA files (`.fna`, `.fa`, `.fasta`) | — |
| `--fastq_folder` | **(required)** Directory of FASTQ files. Paired-end: `*_r1.fq[.gz]` / `*_r2.fq[.gz]`. Single-end: `*.fq[.gz]` | — |
| `--outdir` | Output directory | `results` |
| `--read_type` | `paired` or `single` | `paired` |
| `--parsnp_flags` | Extra flags passed to parsnp | `-c` |
| `--window_size` | Window size for variant filtering (positive integer or `average_LCB_length`) | `500` |
| `--window_overlap` | Overlap fraction between windows (0–1) | `0` |
| `--filter_off` | Skip the recombination filter (recommended when genomes differ by <500 variants) | `false` |
| `--weight_by_entropy` | Weight variants by Shannon entropy when estimating abundances | `false` |
| `--bootstrap` | Compute 95% bootstrap confidence intervals | `false` |
| `--bootstrap_iterations` | Number of bootstrap iterations | `100` |
| `--use_precomputed_variants` | Skip parsnp/variant-filtering; use existing outputs | `false` |
| `--precomputed_dir` | Directory containing `filtered_variant_matrix.csv`, `reference.fna`, `sites.txt` | — |
| `--prefilter` | Pre-filter genomes by coverage before the main run | `false` |

### Profiles (`-profile`)

| Profile | Description |
|---|---|
| `standard` | Local execution (default) |
| `conda` | Activates conda environment defined in `environment.yml` |
| `test` | Runs against the bundled `example/` data |

### Controlling resources

```bash
# Cap to 8 CPUs and 64 GB RAM
strainify ... --max_cpus 8 --max_memory 64.GB
```

---

## Output

| File / Directory | Description |
|---|---|
| `renamed_genomes/` | Reference genomes with sanitised FASTA headers |
| `parsnp_results/parsnp.maf` | Whole-genome alignment (MAF format) |
| `parsnp_results/merged.vcf` | Multi-sample VCF |
| `filtered_variant_matrix.csv` | Filtered bi-allelic variant matrix |
| `significantly_enriched_windows.tsv` | Putative recombinant windows |
| `sites.txt` | Variant positions used for read counting |
| `reference.fna` | Reference genome used for mapping |
| `mapped_reads/*.sam` | Per-sample SAM alignments |
| `mapped_reads/*_sorted.bam` | Filtered, sorted BAM files |
| `read_counts/*_read_counts.tsv` | Per-site read counts per sample |
| `abundance_estimates_combined.csv` | **Main output** – relative strain abundances |
| `abundance_estimates_bootstrap_CIs.csv` | Bootstrap 95% CIs (only with `--bootstrap`) |
| `pipeline_info/` | Nextflow execution reports, timeline, and trace |

---

## Tutorial / Examples

For step-by-step instructions using the example data, see:

[Strainify Tutorial](documentation/tutorial.md)

Run the test profile to verify your installation:

```bash
strainify -profile test,conda
```

---

## Reproducing paper analyses

See the scripts and documentation at:
<https://github.com/treangenlab/Strainify_paper>

## Preprint

Strainify: Strain-Level Microbiome Profiling for Low-Coverage Short-Read Metagenomic Datasets  
<https://www.biorxiv.org/content/10.1101/2025.10.10.681738v2>

## Questions / Contact

For questions or suggestions, open an issue or contact Rossie Luo at [rl152@rice.edu](mailto:rl152@rice.edu).
