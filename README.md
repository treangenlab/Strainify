# Strainify
<img src="images/strainify.png" alt="Strainify diagram" style="width:30%;"/>

Strainify is an accurate strain-level abundance analysis tool for short-read metagenomics.


## Installation


### Clone the repository

Strainify uses [Git LFS](https://github.com/git-lfs/git-lfs.git) to manage large files (e.g., genomes, test data). Please ensure that it is set up **before cloning the repository**.

```bash
git clone https://github.com/treangenlab/Strainify.git
cd Strainify
```

Use the following commands to initialize `git lfs` and download the example files: 
```bash
git lfs install
git lfs checkout
```

### Set up conda environment
```bash
# Install conda and Python dependencies
conda env create -f environment.yaml

# Activate conda environment
conda activate strainify
```

## Usage
Strainify uses a `config.yaml` file to manage input files and parameters.

### Parameters
These fields can be set in `config.yaml`:

| Parameter                  | Description                                                                                                         | Default                         | Accepted Values                              |
|----------------------------|---------------------------------------------------------------------------------------------------------------------|----------------------------------|-----------------------------------------------|
| `genome_folder`            | **(Required)** Path to the folder containing reference genome files (FASTA format).                                                | —                                | A valid directory path                        |
| `fastq_folder`             | **(Required)** Path to the folder containing input FASTQ files (must be unzipped and named using `*_r1.fq`, `*_r2.fq`). Multiple samples can be added to this folder (for paired-end reads, make sure the two FASTQ files for each sample have matching names).             | —                                | A valid directory path                        |
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

### Running Strainify without a precomputed variant matrix:
To run the pipeline on the example data:

1. Unzip the compressed FASTQ files. The following command line can be used for Linux systems. 
```bash
gunzip example/fastq/paired/*.gz
gunzip example/fastq/single/*.gz
```

2. Make sure your `config.yaml` is set like this:

```yaml
genome_folder: example/genomes
fastq_folder: example/fastq/paired
output_dir: example/results
read_type: paired
modify_windows: --window_size 500 --window_overlap 0
weight_by_entropy: false
use_precomputed_variants: false
```

3. Run Strainify with Snakemake:
```bash
snakemake --cores 12 --configfile config.yaml
```
The output will be written to the `example/results` directory.

### Running Strainify with a precomputed variant matrix:
If you are running Strainify again on the same set of genomes, you can use the precomputed variant matrix by doing the following:

1. In the `config.yaml` file, set `use_precomputed_variants` to `true` and `precomputed_output_dir` to your desired directory for storing the new output. If you do not provide a path for `precomputed_output_dir`, your output directory will be automatically set to `output_dir/precomputed_results`. Make sure you set `output_dir` to the directory that contains the precomputed variant matrix. For example, let's use the variant matrix from the previous example, and the `config.yaml` file would look like this:

```yaml
genome_folder: example/genomes
fastq_folder: example/fastq/single
output_dir: example/results
read_type: single
modify_windows: --window_size 500 --window_overlap 0
weight_by_entropy: false
use_precomputed_variants: true
precomputed_output_dir: example/new_output
```
2. Run Strainify with Snakemake:
```bash
snakemake --cores 12 --configfile config.yaml
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

## Questions / Contact

For questions or suggestions, open an issue or contact Rossie Luo at [rl152@rice.edu](mailto:rl152@rice.edu).

