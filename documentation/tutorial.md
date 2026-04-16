# Quick Start and Tutorial
## Table of content
- [Running Strainify without a precomputed variant matrix](#running-strainify-without-a-precomputed-variant-matrix)
- [Running Strainify with a precomputed variant matrix](#running-strainify-with-a-precomputed-variant-matrix)
- [Example Output](#example-output)
- [Clustering input reference genomes](#clustering-input-reference-genomes)
- [Unknown strains in metagenomic samples](#unknown-strains-in-metagenomic-samples)


## Running Strainify without a precomputed variant matrix
Strainify includes example input data to help you get started quickly. To run the pipeline on the example data:

1. Unzip the compressed FASTQ files. The following command line can be used for Linux systems. 
```bash
gunzip example/fastq/paired/*.gz
gunzip example/fastq/single/*.gz
```

2. Make sure your `config.yaml` is set like this:

```yaml
genome_folder: example/genomes
fastq_folder: example/fastqs/paired
output_dir: example/results
read_type: paired
modify_windows: --window_size 500 --window_overlap 0
weight_by_entropy: false
use_precomputed_variants: false
```

3. Run Strainify:
```bash
./strainify run --cores 12 --configfile config.yaml
```
The output will be written to the `example/results` directory.

## Running Strainify with a precomputed variant matrix
If you are running Strainify again on the same set of genomes, you can use the precomputed variant matrix by doing the following:

1. In the `config.yaml` file, set `use_precomputed_variants` to `true` and `precomputed_output_dir` to your desired directory for storing the new output. If you do not provide a path for `precomputed_output_dir`, your output directory will be automatically set to `output_dir/precomputed_results`. Make sure you set `output_dir` to the directory that contains the precomputed variant matrix. For example, let's use the variant matrix from the previous example, and the `config.yaml` file would look like this:

```yaml
genome_folder: example/genomes
fastq_folder: example/fastqs/single
output_dir: example/results
read_type: single
modify_windows: --window_size 500 --window_overlap 0
weight_by_entropy: false
use_precomputed_variants: true
precomputed_output_dir: example/new_output
```
2. Run Strainify:
```bash
./strainify run --cores 12 --configfile config.yaml
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

