#! /bin/bash

# Now use our helper script (column name for sample from merged sylph output)
python /path/to/magnet/helper/sylph2magnet.py \
  -s /path/to/magnet/test/sylph_mpa_merged_abundance_table.tsv \
  --acc-taxmap /path/to/magnet/tax/sylph_gtdb226_2_ncbi.tsv \
  --sample "name_of_sample_in_sylph_abund_header" \
  -o taxtest_magnet_species.tsv \
  --min-abundance 0.0 \
  --scale-to-fraction

# Run magnet
magnet -c taxtest_magnet_species.tsv \
  -i /path/to/sample.fastq.gz \
  -o taxtest_magnet_species_out \
  -m ont \
  -t 0 \
  -a 1 \
  --threads 64 \
  --min-abundance 0.001
