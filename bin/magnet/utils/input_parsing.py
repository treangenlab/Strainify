import os
import argparse
import pathlib
import sys
import subprocess
import math
import warnings
from io import StringIO
from collections import defaultdict
from multiprocessing import Pool, Manager
from itertools import repeat
warnings.filterwarnings("ignore")

import pandas as pd
from ete3 import NCBITaxa
from Bio import SeqIO

def get_species_taxid(taxid, ncbi_taxa_db, valid_kingdom, ret_subspecies):
    try:
        lineage = ncbi_taxa_db.get_lineage(taxid)
    except ValueError:
        return None 
    
    if bool(set(lineage) & valid_kingdom):
        taxid2rank_dict = ncbi_taxa_db.get_rank(lineage)
        species_rank_visited = False
        for lineage_taxid in lineage:
            if species_rank_visited and ret_subspecies:
                return lineage_taxid
            if taxid2rank_dict[lineage_taxid] == 'species':
                species_rank_visited = True
                if not ret_subspecies or lineage_taxid == lineage[-1]:
                    return lineage_taxid
    return None

def filter_input_df(input_df, min_abundance, ncbi_taxa_db, valid_kingdom, ret_subspecies):
    valid_taxids = []
    input_df = input_df.fillna(0)
    for idx, row in input_df.iterrows():
        if min_abundance > 0 and row['abundance'] < min_abundance:
            continue
        if row['tax_id'].isnumeric():
            species_taxid = get_species_taxid(row['tax_id'], ncbi_taxa_db, valid_kingdom, ret_subspecies)
            if species_taxid is not None:
                valid_taxids.append(species_taxid)
        else:
            print(row['tax_id'], 'is not a valid taxid.')
            continue
    
    valid_taxids = list(set(valid_taxids))
    print(f"num_species: {len(valid_taxids)}")
    return valid_taxids

def parsing_input_f(input_file, sep, taxid_col_idx, abundance_col_idx, min_abundance):
    input_df = pd.DataFrame(columns = ['tax_id', 'abundance'])
    taxid_col_values_df = pd.read_csv(input_file, sep=sep, usecols=[taxid_col_idx], dtype=str, index_col=False)
    taxid_col_values = taxid_col_values_df[taxid_col_values_df.columns[0]].values
    input_df['tax_id'] = list(taxid_col_values)
    if abundance_col_idx != None:
        abundance_col_values_df = pd.read_csv(input_file, sep=sep, usecols=[abundance_col_idx])
        abundance_col_values = abundance_col_values_df[abundance_col_values_df.columns[0]].values
        input_df['abundance'] = list(abundance_col_values)
    else:
        min_abundance = 0
    
    print(f'min_abundance: {min_abundance}')
    return input_df, min_abundance

def get_seq2assembly_dict(working_directory, downloaded_assemblies):
    seq2assembly_dict = dict()
    for assembly_id in downloaded_assemblies['Assembly Accession ID']:
        reference_genome_fasta = os.path.join(working_directory, 'reference_genomes', f'{assembly_id}.fasta')
        seq_dict = SeqIO.to_dict(SeqIO.parse(reference_genome_fasta, "fasta"))
        for seq in seq_dict:
            seq2assembly_dict[seq] = assembly_id
            
    return seq2assembly_dict
