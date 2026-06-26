#! /bin/bash

# accession \t species_name
awk -F'\t' '{
  lineage = $2
  n = split(lineage, a, ";")
  species = a[n]
  sub(/^s__/, "", species)
  print $1 "\t" species
}' gtdb_r226_metadata.tsv > accession_species.tsv

# Adjust this to wherever your NCBI taxdump is
NCBI_DB=/dodo/am503/ncbi_taxdump

#mkdir -p "$NCBI_DB"
#cd "$NCBI_DB"
#wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
#tar xzf new_taxdump.tar.gz

# Now run name2taxid; -i 2 means "use 2nd column as name field"
taxonkit name2taxid \
  --data-dir "$NCBI_DB" \
  -i 2 \
  /dodo/am503/sylph_db/accession_species.tsv \
  > /dodo/am503/sylph_db/accession_species_taxid_raw.tsv

# filter out those without a taxid
awk -F'\t' 'NF>=3 && $3 != ""' /dodo/am503/sylph_db/accession_species_taxid_raw.tsv \
  > /dodo/am503/sylph_db/accession_species_taxid.tsv

