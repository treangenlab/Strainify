# MAGnet

MAGnet is a computational tool that uses reference based method to refine taxonomic assignment and provide accurate abundance estimation. This pipeline handles reference genome selection, download using NCBI Datasets, read alignment using Minimap2, and breadth and depth of genome coverage calculations.

## Installation

Please make sure to have a working conda installation first. If you don't, please select the appropriate version found [here](https://github.com/conda-forge/miniforge).

```
# Clone the respository
git clone https://github.com/treangenlab/magnet.git

# Go to repo
cd magnet

# Create conda environment
conda env create -n magnet -f magnet_env.yml

# Activate conda environment
conda activate magnet

# Test install with 
magnet -h
```

## Getting started

```
usage: magnet [-h] -c CLASSIFICATION -i FASTQ [-I FASTQ2] [-m {ont,illumina}] -o OUTPUT [-t TAXID_IDX]
              [-a ABUNDANCE_IDX] [--min-abundance MIN_ABUNDANCE] [--min-mapq MIN_MAPQ]
              [--min-covscore MIN_COVSCORE] [--threads THREADS] [--include-mag] [--subspecies]

Universal Taxonomic Classification Verifier.

optional arguments:
  -h, --help            show this help message and exit
  -c CLASSIFICATION, --classification CLASSIFICATION
                        Path to the Taxonomic Classification Report. Accepting csv/tsv file format, other text formats
                        are treated as tsv.
  -i FASTQ, --fastq FASTQ
                        Path to the first fastq file.
  -I FASTQ2, --fastq2 FASTQ2
                        Path to the second fastq file for paired-end reads.
  -m {ont,illumina}, --mode {ont,illumina}
                        Modes for different sequencing platforms [ont, illumina]. Default:[ont]
  -o OUTPUT, --output OUTPUT
                        Path to the output directory.
  -t TAXID_IDX, --taxid-idx TAXID_IDX
                        The column index (0-based) of the taxids. Default:[0]
  -a ABUNDANCE_IDX, --abundance-idx ABUNDANCE_IDX
                        The column index (0-based) of the abundance. Default:[None]
  --min-abundance MIN_ABUNDANCE
                        Minimum abundance (0-1) for pre-filtering, exclude taxa below the threshold.
  --min-mapq MIN_MAPQ   Minimum MAPQ for primary alignments. Default:[20]
  --min-covscore MIN_COVSCORE
                        Minimum Coverage Score for supplementary alignments. Default:[0.7]
  --threads THREADS     Number of threads for Multi-threading. Default:[1]
  --include-mag         Include metagenomic assemble genomes. Default:[False]
  --subspecies          Verify taxonomic classification at subspecies rank. Default:[False]
```

## Example of Using Lemur Report as Input for ONT datasets
Since Lemur report uses a format that taxid is the first column, we set '-t' to '0', which is already the default setting, and can be omit. Using the following command to run Magnet. 

```
magnet -c {lemur_report_file} -i {input_fastq_file} -o {output_path} -m ont
```

## Example of Using Sylph Merged Report as Input
Sylph uses GTDB databases for it's default so we need to convert these GTDB taoxnomies to NCBI taxonomies to use Magnet. We wrote a helper script that will take in a merged sylph-tax report in the metaphlann style format and will output a .tsv for use within Magnet. An full example workflow can be seen below:

We already created the sylph_gtdb226_2_ncbi.tsv.gz file (unzip before using) but for future use follow the bash script in the magnet/tax dir similarly.

```
# Activate your conda env
conda activate magnet

# Unzip the acc-taxmap file
gzip -d /path/to/magnet/tax/sylph_gtdb226_2_ncbi.tsv.gz

# Now use our helper script (column name for sample from merged sylph output)
python /path/to/magnet/helper/sylph2magnet.py \
  -s /path/to/sylph_mpa_merged_abundance_table.tsv \
  --acc-taxmap /path/to/magnet/tax/sylph_gtdb226_2_ncbi.tsv \
  --sample "sample1.fastq.gz" \
  -o sample1_magnet_species.tsv \
  --min-abundance 0.0 \
  --scale-to-fraction

# Run magnet for sample1
magnet -c sample1_magnet_species.tsv \
  -i /path/to/sample1.fastq.gz \
  -o /path/to/magnet_sample1 \
  -m ont \
  -t 0 \
  -a 1 \
  --min-abundance 0.001
```

With this format, taxid is the first column, we set '-t' to '0', and since the abundance estimation is located in the second column, we set '-a' to '1'. To filter out potiential noise at low abundance, we  set '--min-abundance' to 0.001.

## Example Output

The following table is an example of the output csv file 'cluster_representative.csv'.

| Taxonomy ID | Assembly Accession ID | Source Database | Is Representative | Assembly Level | Organism of Assembly | Strain | Total Length | Downloaded | Species | Cluster Representative | Cluster Members | Secondary Breadth | Secondary Expected | Secondary Score | Secondary Depth | Primary Breadth | Primary Expected | Primary Score | Primary Depth | Consensus ANI | Combined PS and ANI (Sqrt(ANI)xPSx100) | Presence/Absence |
|-------------|-----------------------|-----------------|-------------------|----------------|----------------------|--------|--------------|------------|---------|------------------------|-----------------|-------------------|--------------------|-----------------|-----------------|-----------------|------------------|---------------|---------------|---------------|----------------------------------------|-----------------|
| 2678528 | GCF_028867355.1 | SOURCE_DATABASE_REFSEQ | False | Complete Genome | Listeria sp. LM90SB2 | LM90SB2 | 2915834.0 | True | Listeria sp. LM90SB2 | True | GCF_028867355.1,GCF_000196035.1 | 1.0 | 1.0 | 1.0 | 3195.6499999657044 | 1.0 | 1.0 | 1.0 | 3188.0799997530726 | 0.9957277140090639 | 99.79 | Present |
| 96241   | GCF_006094475.1 | SOURCE_DATABASE_REFSEQ | True  | Complete Genome | Bacillus spizizenii ATCC 6633 = JCM 2499 | ATCC 6633 | 4045538.0 | True | Bacillus spizizenii | True | GCF_006094475.1 | 0.9999987640704401 | 0.9999999999999607 | 0.9999987640704794 | 30.85103792256793 | 0.9999987640704401 | 0.999999999999958 | 0.999998764070482 | 30.786837976602836 | 0.9952862581301721 | 99.76 | Present |
| 1613    | GCF_029961225.1 | SOURCE_DATABASE_REFSEQ | True  | Complete Genome | Limosilactobacillus fermentum | EFEL6800 | 2103331.0 | True | Limosilactobacillus fermentum | True | GCF_029961225.1 | 0.17646569702852932 | 0.1786636068369798 | 0.9876980553154514 | 1.1147205952452701 | 0.17275099611416372 | 0.17304955697814606 | 0.9982747088799537 | 1.0992270709852185 | 0.853534466277766 | 92.23 | Present |
| 287     | GCF_000006765.1 | SOURCE_DATABASE_REFSEQ | True  | Complete Genome | Pseudomonas aeruginosa PAO1 | PAO1 | 6264404.0 | True | Pseudomonas aeruginosa | True | GCF_024714315.1,GCF_000006765.1 | 0.9736496879830867 | 1.0 | 0.9736496879830867 | 133.57679009269043 | 0.9736029157761855 | 1.0 | 0.9736029157761855 | 133.2863715317914 | 0.9790992928162352 | 96.34 | Present |
| 4932    | GCF_000146045.2 | SOURCE_DATABASE_REFSEQ | True  | Complete Genome | Saccharomyces cerevisiae S288C | S288C | 12071326.0 | True | Saccharomyces cerevisiae | True | GCF_000146045.2 | 0.9702460413067091 | 0.9988842650498445 | 0.9713297878991953 | 7.005664504492986 | 0.9551612822296097 | 0.9987129198063666 | 0.9563922357335672 | 6.96673509016021 | 0.8611351208262882 | 88.75 | Present |
| 28901   | GCF_000006945.2 | SOURCE_DATABASE_REFSEQ | True  | Complete Genome | Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 | LT2 | 4951383.0 | True | Salmonella enterica | True | GCF_000006945.2,GCF_014492125.1,GCF_022556855.1,GCF_014491785.1 | 0.7829200506438564 | 0.8320860674724839 | 0.9409123421835766 | 2.2776524817939543 | 0.7817768582280826 | 0.8295483061281343 | 0.9424126991193291 | 2.2617706866885654 | 0.9299024245256468 | 90.88 | Present |
| 562     | GCF_000008865.2 | SOURCE_DATABASE_REFSEQ | True  | Complete Genome | Escherichia coli O157:H7 str. Sakai | Sakai substr. RIMD 0509952 | 5594605.0 | True | Escherichia coli | True | GCF_000008865.2 | 0.6591498747494352 | 0.7795045795270928 | 0.8456010292451726 | 2.2926196898951465 | 0.6569027483105632 | 0.7772529211094938 | 0.8451595747917728 | 2.284965162259128 | 0.8996365155332352 | 80.16 | Present |
| 1913579 | GCF_001875785.1 | SOURCE_DATABASE_REFSEQ | False | Scaffold        | Bacillus sp. FMQ74 | FMQ74 | 4161595.0 | True | Bacillus sp. FMQ74 | True | GCF_001875785.1 | 0.003557530225790833 | 0.03775324221474807 | 0.0942311180998676 | 10.81661600810537 | 0.002592275317516481 | 0.02381212312913672 | 0.10886367853291296 | 9.295791620318873 | 0.905627438899158 | 10.36 | Absent |

The table contains the following information on the reference genome that are being selected for alignment:

- Taxonomy ID: The taxonomy ID of the species that users queries
- Assembly Accession ID: The assembly accession ID of the genome that is being selected as reference genome for the alignment process
- Source Database: source of the assembly, part of the metadata from NCBI Database
- Is Representative: Whether or not the assembly is labeled as representative genome in NCBI database
- Assembly Level: the assembly level of the assembly, can be one of the following: complete genome, chromosome, contig, scaffold
- Organism of Assembly: The description of the assembly, may contain strain level information
- Strain: The specific strain used for alignment
- Species: taxonomy name of the species
- Cluster Representative: Whether the genome is a cluster representative after ANI clustering
- Cluster Members: Members of the same cluster

After the reference genomes are being selected, the pipeline create a reference in multi-fasta format by concatenate the reference genomes, and aligns the reads (in fasta format) to this concatenated reference, and calculated the following:

- Secondary Breadth Coverage: The proportion of the genome has reads aligned, including secondary alignment and alignments with low MAPQ. 
- Secondary Expected Coverage: An estimation of the expected breadth of coverage based on the assumption that the alignment of the read is randomly distributed along the entire reference genome, excluding secondary alignment and alignments with low MAPQ. 
- Secondary Coverage Score: Secondary Breadth Coverage / Secondary Expected Coverage 
- Secondary Depth Coverage: The mean depth of coverage, calculated with aligned regions only, excluding secondary alignment and alignments with low MAPQ. 
- Primary Breadth Coverage: The proportion of the genome has reads aligned, excluding secondary alignment and alignments with low MAPQ. 
- Primary Expected Coverage: An estimation of the expected breadth of coverage based on the assumption that the alignment of the read is randomly distributed along the entire reference genome, excluding secondary alignment and alignments with low MAPQ. 
- Primary Coverage Score: Primary Breadth Coverage / Primary Expected Coverage
- Primary Depth Coverage: The mean depth of coverage, calculated with aligned regions only, excluding secondary alignment and alignments with low MAPQ. 
- Consensus ANI: The ANI is calculated by comparing the consensus genome created by the primary alignment with the original reference genome. Unaligned regions are excluded. 
- Combined PS and ANI (Sqrt(ANI)xPSx100)
- Presence/Absence Calls

