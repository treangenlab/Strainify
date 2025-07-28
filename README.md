# Strainify

Strainify is an accurate strain-level abundance analysis tool for short-read metagenomics.


## Installation


### Clone the repository

```bash
git clone https://github.com/treangenlab/Strainify.git
cd Strainify
```

### Set up conda environment
```bash
# Install conda and Python dependencies
conda env create -f environment.yaml

# Activate conda environment
conda activate strainify
```

### Install other dependencies
Strainify requires the following external dependencies:

- [wgatools 0.1.0](https://github.com/wjwei-handsome/wgatools.git)  
- [bcftools 1.21](https://www.htslib.org/download/)
- [htslib 1.21](https://www.htslib.org/download/)
- [bwa 0.7.18](https://github.com/lh3/bwa.git)
- [samtools 1.13](https://www.htslib.org/download/)

## Usage
Strainify uses a `config.yaml` file to manage input files and parameters.

### Parameters
These fields can be set in `config.yaml`:

| Parameter                  | Description                                                                                                         | Default                         | Accepted Values                              |
|----------------------------|---------------------------------------------------------------------------------------------------------------------|----------------------------------|-----------------------------------------------|
| `genome_folder`            | **(Required)** Path to the folder containing reference genome files (FASTA format).                                                | —                                | A valid directory path                        |
| `fastq_folder`             | **(Required)** Path to the folder containing input FASTQ files (must be unzipped and named using `*_r1.fq`, `*_r2.fq`).             | —                                | A valid directory path                        |
| `output_dir`               | **(Required)** Directory where all output files will be saved.                                                                     | —                       | A valid directory path                        |
| `read_type`                | Type of sequencing reads.                                                                                           | `paired`                         | `paired` or `single`                          |
| `window_size`              | Size of the genomic window for variant grouping. Can also be set to `average_LCB_length` to use the average length of local colinear blocks from Parsnp. | `500`                            | Any positive integer or `average_LCB_length` |
| `window_overlap`           | Proportion of overlap between consecutive windows.                                                                  | `0`                              | Float between `0` and `1` (e.g., `0.5`)       |
| `weight_by_entropy`        | Whether to weight variants by their Shannon entropy when estimating strain abundances.                              | `false`                          | `true` or `false`                             |
| `use_precomputed_variants`| Use existing filtered variant matrix instead of recomputing from scratch.                                           | `false`                          | `true` or `false`                             |
| `precomputed_output_dir`  | Path to directory where new output will be saved.                                                        | `output_dir/precomputed_results` | A valid directory path                        |



### Edit the `config.yaml`
Open `config.yaml` in a text editor. Modify the fields to match your input and desired options. Example:

```yaml
genome_folder: path/to/genomes
output_dir: path/to/output
fastq_folder: path/to/fastqs
read_type: paired
modify_windows: --window_size 500 --window_overlap 0
weight_by_entropy: false
use_precomputed_variants: false
precomputed_output_dir: path/to/new_output
```

### Running Strainify
Use Snakemake to run the pipeline:
```bash
snakemake --cores 12 --configfile config.yaml
```
>Tip: Replace `12` with the number of CPU cores you want to allocate.
>You can set `--cores` to any positive integer, depending on your system’s available resources.


## Examples

Strainify includes example input data to help you get started quickly.

To run the pipeline on the example data:

1. Make sure your `config.yaml` is set like this:

```yaml
genome_folder: example/genomes
fastq_folder: example/paired
output_dir: example/results
read_type: paired
modify_windows: --window_size 500 --window_overlap 0
weight_by_entropy: false
use_precomputed_variants: false
```
2. Run Strainify with Snakemake:
```bash
snakemake --cores 12 --configfile config.yaml
```
The output will be written to the `example/results` directory.

## Example Output

The abundance estimates are stored in a CSV file named `abundance_estimates_combined.csv`.

- Each row corresponds to a strain.
- Each column (after the first) corresponds to a sample.
- Values represent the estimated relative abundances.

## Questions / Contact

For questions or suggestions, open an issue or contact Rossie Luo at [rl152@rice.edu](mailto:rl152@rice.edu).

