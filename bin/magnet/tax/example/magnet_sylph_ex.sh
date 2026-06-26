#! /bin/bash

# Make sure to run in magnet conda env created from taxdev branch install instructions

# Convert the sylph GTDB taxids to the NCBI accessions that Magnet will use
python /dodo/am503/magnet/helper/sylph2magnet.py \
  -s /dodo/am503/rice_mgx_fqs/25_11_10_sylph_mpa_merged_abundance_table.tsv \
  --acc-taxmap /dodo/am503/magnet/tax/sylph_gtdb226_2_ncbi.tsv \
  --sample "/home/tmhagm8/scratch/rice_concuss/rice_prom_data/prep/25_11_10_final_clean_fqs/127_1.chop.fastq.gz" \
  -o /dodo/am503/magnet/tax/example/127_1_magnet_species.tsv \
  --min-abundance 0.0 \
  --scale-to-fraction

# Now run Magnet on the output
python /dodo/am503/magnet/magnet.py \
  -c /dodo/am503/magnet/tax/example/127_1_magnet_species.tsv \
  -i /dodo/am503/rice_mgx_fqs/all_fastqs/127_1.chop.fastq.gz \
  -o /dodo/am503/magnet/tax/example/magnet_127_1 \
  -m ont \
  -t 0 \
  -a 1 \
  --threads 120 \
  --min-abundance 0.001

