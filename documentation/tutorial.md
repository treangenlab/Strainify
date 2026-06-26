# Quick Start and Tutorial
## Table of content
- [Running Strainify without a precomputed variant matrix](#running-strainify-without-a-precomputed-variant-matrix)
- [Running Strainify with a precomputed variant matrix](#running-strainify-with-a-precomputed-variant-matrix)
- [Example Output](#example-output)
- [Clustering input reference genomes](#clustering-input-reference-genomes)
- [Unknown strains in metagenomic samples](#unknown-strains-in-metagenomic-samples)


## Running Strainify without a precomputed variant matrix
Strainify includes example input data to help you get started quickly. 

Build and activate the environment once (see the README), then run Strainify on the paired-end
example data. 

```bash
conda activate strainify
./strainify \
  --genome_folder example/genomes \
  --fastq_folder  example/fastqs/paired \
  --outdir        example/results
```
`paired` is the default read type, so `--read_type` is omitted here; for single-end data pass
`--read_type single` (as in the next section). The output is written to `example/results`.

## Running Strainify with a precomputed variant matrix
If you are running Strainify again on the same set of genomes, you can reuse the precomputed
variant matrix by passing `--use_precomputed_variants` along with `--precomputed_dir` (the
directory that **contains** the previous run's `filtered_variant_matrix.csv`, `reference.fna`,
and `sites.txt`) and a fresh `--outdir` for the new results. Reusing the matrix from the
previous example:

```bash
./strainify \
  --genome_folder example/genomes \
  --fastq_folder  example/fastqs/single \
  --read_type     single \
  --use_precomputed_variants \
  --precomputed_dir example/results \
  --outdir          example/new_output
```
The output will be written to the `example/new_output` directory.

## Example Output

The abundance estimates are stored in a CSV file named `abundance_estimates_combined.csv`.

- Each row corresponds to a strain.
- Each column (after the first) corresponds to a sample.
- Values represent the estimated relative abundances.

Example CSV output:
```bash
strain name,10x_ratio_2
E24377A,23.7987
H10407,24.8376
Sakai,24.9406
UTI89,26.423
```
>Note: these numbers are percentages. 

Other important output files:
- `sites.txt` contains a list of variant positions that passed the filter. Read counts supporting the allele and reference base at these positions are then obtained and used as input to the MLE model. 
- `filtered_variant_matrix.csv` contains the filtered variant matrix. Confounding variants (potential recombination sites) have been removed. For metagenomic samples that share the same set of strains (i.e. query genomes), this file can be reused to avoid rerunning the genome alignment and variant filtering steps. For more details, see instructions above for running Strainify with a precomputed variant matrix. 
- `significantly_enriched_windows.tsv` contains the start and end coordinates of windows that are flagged as potential recombination sites. The z-score and p-values of each window are also shown. Variants in these windows are removed from downstream analysis (i.e. excluded from the filtered variant matrix).
- `abundance_estimates_bootstrap_CIs.csv` (if `--bootstrap` is applied) contains the 95% confidence intervals for each estimated relative abundance value. 

## Clustering input reference genomes
By default, Strainify defines each input reference genome as a strain, and it is able to operate with strains that are highly similar provided a reference genome is available for each one. However, strains can also be defined at a less granular level by clustering similar reference genomes and selecting a representative sequence from each cluster as input for Strainify. 

A simple tool that can be used for clustering is [TreeCluster](https://github.com/niemasd/TreeCluster). You can first run [Parsnp](https://github.com/marbl/parsnp) (included in Strainify's dependencies, no need to install separately) with your reference genomes to obtain a tree (one of Parsnp's default outputs). For example:
```bash
parsnp -r reference_genome.fna -d path/to/genome/folder/*.fna
```
>Tip: if you would like Parsnp to randomly select a genome from the input as the reference for alignment, you can simply replace `reference_genome.fna` with `!`. 

Look for  `parsnp.tree` in Parsnp's output, and provide it to TreeCluster. Here is an example command:
```bash
TreeCluster.py -i parsnp.tree -m max_clade -t 0.001 -o clusters.tsv
```
You can then select a representative genome from each cluster and provide them to Strainify. Strainify will then estimate relative abundances at the cluster level.

## Unknown strains in metagenomic samples
When strains in a metagenomic sample to be analyzed are unknown (i.e. reference genomes are not available), Strainify can be paired with an upstream strain identification tool such as [StrainGST](https://strainge.readthedocs.io/en/latest/straingst.html). Just run the upstream tool and collect the genomes that it identified in your metagenomic samples, and then run Strainify with these genomes as references. For longitudinal/time-series studies, we recommend running the upstream tool on metagenomic samples from <i>all</i> timepoints, and provide <i>all</i> of the identified genomes to Strainify.

### Using the bundled StrainGST helper scripts

Two helper scripts in `helpers/` automate the StrainGST hand-off.

> **You must install [StrainGE](https://github.com/broadinstitute/StrainGE) yourself.** Install it in a **separate conda environment** and
> **activate that environment before running these scripts** (the NCBI download route also needs
> `ncbi-genome-download` in that environment). If you'd rather not activate it, the scripts accept
> the optional `--strainge-env <name>` (and `--download-env <name>` for the NCBI step) to run the
> tools in the named env via `conda run -n` instead.

**1. Build a StrainGST database (one-time).** Either download reference genomes from NCBI by
genus/species, or point at a folder of genomes you already have:

```bash
# from NCBI:
helpers/build_straingst_db.sh \
  --taxa "Escherichia,Shigella" \
  --threads 16 --outdir straingst_db

# or from a custom genome folder:
helpers/build_straingst_db.sh \
  --genome-folder my_reference_genomes/ \
  --outdir straingst_db
```

This produces `straingst_db/pan-genome-db.hdf5` (the database) and `straingst_db/db/` (the
reference genomes). Tune with `--kmer-size`, `--cluster-jaccard`, and `--cluster-subset`
(run with `--help` for all options).

**2. Discover strains per sample and collect their genomes.** Run StrainGST on your FASTQs; for
each sample it lists the discovered strains and copies their genome FASTAs into a folder ready
for Strainify:

```bash
helpers/run_straingst.sh \
  --db      straingst_db/pan-genome-db.hdf5 \
  --ref-dir straingst_db \
  --fastq-folder path/to/fastqs --read-type single \
  --outdir straingst_results
```

For each `<sample>` this writes `straingst_results/<sample>/<sample>.strains.txt` (the discovered
strains), a combined `straingst_results/discovered_strains.tsv` summary, and
`straingst_results/<sample>/genomes/` (the genome FASTAs).

StrainGST parameters are fully tunable: pass any `straingst run` flags through `--run-args "..."`
(for example `-i` to cap how many strains are reported, or `--minscore` to raise the minimum
score), any `straingst kmerize` flags through `--kmerize-args "..."`, and set the k-mer size with
`--kmer-size` (it must match the database). Run the script with `--help` for the full list.

**3. Run Strainify on the discovered genomes.** Point `--genome_folder` at a sample's collected
genomes:

```bash
conda activate strainify
./strainify run \
  --genome_folder straingst_results/<sample>/genomes \
  --fastq_folder  path/to/fastqs \
  --read_type     single \
  --outdir        results/<sample>
```

For longitudinal/time-series studies, pool the genomes discovered across *all* timepoints into a
single folder and pass that to Strainify, so every timepoint is profiled against the same strain set.
