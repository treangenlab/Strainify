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
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/strainify/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/strainify/badges/version.svg)](https://anaconda.org/bioconda/strainify)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/strainify/badges/downloads.svg)](https://anaconda.org/bioconda/strainify)

### Installing via conda
```bash
conda install -c conda-forge -c bioconda strainify
conda activate strainify
strainify --help
```

### Installing via git

```bash
git clone https://github.com/treangenlab/Strainify.git
cd Strainify

# Create the conda environment (Linux / macOS)
conda env create -f environment.yml
conda activate strainify

# Run directly from the repository root
./strainify --help
```

> **macOS (Apple Silicon) note:** `parsnp` and `harvesttools` do not yet have
> native `osx-arm64` builds. Conda must resolve them through Rosetta 2 using
> the `osx-64` sub-architecture. Set this before creating the environment and
> keep it in your shell profile:
>
> ```bash
> export CONDA_SUBDIR=osx-64
> conda env create -f environment.yml
> ```
>
> To persist across sessions, add `export CONDA_SUBDIR=osx-64` to your
> `~/.zshrc` (or `~/.bash_profile`). Rosetta 2 must be installed
> (`softwareupdate --install-rosetta`).

---

## Usage

Note: if installed via conda, replace `./strainify` with `strainify` in the following commands.

### Basic run (paired-end reads)

```bash
./strainify \
  --genome_folder path/to/genomes \
  --fastq_folder  path/to/fastqs \
  --outdir        results
```

### Single-end reads

```bash
./strainify \
  --genome_folder path/to/genomes \
  --fastq_folder  path/to/fastqs \
  --read_type     single \
  --outdir        results
```

### Pre-filter genomes with MAGNET (`filter-run`)

If you suspect some reference genomes are absent from the metagenome, use `filter-run`. It runs
MAGNET to keep only the genomes called **present** in each sample, then runs Strainify on those
genomes — in a single command. Point `--fastq_folder` at a folder of sets (they share one
build-once index) or give a single set with `--fastq1`/`--fastq2`:

```bash
./strainify filter-run \
  --genome_folder  path/to/genomes \
  --fastq_folder   path/to/fastqs \
  --magnet_ref_dir shared_panel \
  --jobs           1 \
  --outdir         results
```

Each sample's results go under `results/<sample>/`, with MAGNET's present/absent calls and the
genomes it kept under `results/<sample>/magnet/<sample>/`. Note: `--magnet_ref_dir` is optional — it
defaults to a path under the launch directory — but setting it is recommended so every set
reuses the same build-once index. `--jobs` is also optional and sets the number of samples to process concurrently (defaults to 1 if not set, i.e. samples are run serially). Each sample's heavy steps cap at ~min(12,max_cpus) cores, so set N near (cores / 12) to fill a big node without oversubscribing.

#### About the bundled MAGNET
 
Strainify bundles a **modified version of [MAGNET](https://github.com/treangenlab/magnet)** (the
Treangen Lab's whole-genome read-mapping presence/absence caller), adapted for strain-level
filtering. Beyond a "bring-your-own-genomes" input mode and a build-once shared reference reused
across samples, it adds two presence-call refinements that are **enabled by default** in
`filter-run`:
 
- an **absolute evidence gate** that removes false positives a genome can otherwise earn from a
  saturated breadth ÷ expected-breadth ratio off only a handful of mapped reads (it downgrades
  *present → absent* unless there is enough absolute primary evidence), and
- a **primary-evidence rescue** that recovers genuinely-present **low-abundance** strains the
  coverage-score threshold would otherwise drop, using uniquely-mapped (primary) read support
  plus consensus ANI (it upgrades *absent → present*).
See [`bin/magnet/README.md`](bin/magnet/README.md) for the full description, attribution, and the
`--magnet_*` parameters that tune presence sensitivity (and their defaults). If too few genomes
pass (Strainify needs at least 2 references), loosen those thresholds — e.g. lower
`--magnet_min_covscore` or `--magnet_rescue_min_reads`.

### Use precomputed variants (re-run with new samples)

If you already have `filtered_variant_matrix.csv`, `reference.fna`, and `sites.txt` from a
previous run, skip the parsnp step:

```bash
./strainify \
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
| `--prefilter` | Enable the MAGNET present/absent genome prefilter (set automatically by the `filter-run` command) | `false` |

### Profiles (`-profile`)

Strainify runs the tools from your **active conda environment** (built once from
`environment.yml` — see [Installation](#installation)), so a normal run needs no
`-profile` at all.

| Profile | Description |
|---|---|
| `standard` | Local execution from your active environment (default; applied automatically) |
| `test` | Runs against the bundled `example/` data |

### Controlling resources

```bash
# Cap to 8 CPUs and 64 GB RAM
./strainify ... --max_cpus 8 --max_memory 64.GB
```

---

## Tutorial / Examples

For step-by-step instructions using the example data, see:

[Strainify Tutorial](documentation/tutorial.md)

Run the test profile to verify your installation (from your activated environment):

```bash
./strainify -profile test
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
