import warnings
warnings.filterwarnings("ignore")

import os
import argparse
import sys
import subprocess
import math
import pickle
from collections import defaultdict
import json
import pandas as pd
import math
from ete3 import NCBITaxa
from Bio import SeqIO

from utils.reference_finder import prepare_reference_genomes
from utils.alignment import run_minimap2, sort_samfile, run_bwa, samtools_calculate_depth
from utils.summary import alignment_1_summary, alignment_2_summary, merge_reference_fasta, call_present_absent
from utils.ani import samtools_merged_consensus, ani_summary

def get_species_taxid(taxid, ncbi_taxa_db, valid_kingdom):
    lineage = ncbi_taxa_db.get_lineage(taxid)
    if bool(set(lineage) & valid_kingdom):
        taxid2rank_dict = ncbi_taxa_db.get_rank(lineage)
        for lineage_taxid in taxid2rank_dict:
            if taxid2rank_dict[lineage_taxid] == 'species':
                return lineage_taxid
    return None

def filter_kraken2_taxonomy(kraken2_report, min_frac, ncbi_taxa_db, valid_kingdom):
    classification_result_df = pd.read_csv(kraken2_report, sep='\t', names=['Abundance', 'Cumulative Count', 'Count', 'Rank', 'Taxid', 'Taxname']) 
    taxids = classification_result_df[(classification_result_df['Rank'] == 'S') & (classification_result_df['Abundance'] >= min_frac*100)]['Taxid'].tolist()
    valid_taxids = []
    for taxid in taxids:
        species_taxid = get_species_taxid(taxid, ncbi_taxa_db, valid_kingdom)
        if species_taxid is not None:
            valid_taxids.append(species_taxid)
    return valid_taxids

ncbi_taxa_db = NCBITaxa()

input_fastq_1 = '/data/ww_data/final_ww_data/fastp/wastewater metagenome/SRR11088368_1.fastq.gz'
input_fastq_2 = '/data/ww_data/final_ww_data/fastp/wastewater metagenome/SRR11088368_2.fastq.gz'
kraken2_report = '/data/ww_data/final_ww_data/kraken/wastewater metagenome/SRR11088368_report.txt'

threads = 40
min_mapq = 20
min_coverage_score = 0.7
min_frac = 0.0005

working_directory = '/scratch0/yl181/ww_pathogen/test_run'
if not os.path.exists(working_directory):
    os.mkdir(working_directory)