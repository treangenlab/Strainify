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

def cal_combined_cs2_ani(cs2, ani):
    return round(math.sqrt(ani)*cs2*100,2)

def _ani_summary(downloaded_assemblies, consensus_record_dict, output_directory):
    ani_list = []
    combined_cs2_ani_list = []
    for idx, row in downloaded_assemblies.iterrows():
        if row['CS2'] != 0:
            assembly_id = row['Assembly Accession ID']
            ani = cal_ani(assembly_id, output_directory, consensus_record_dict)
            combined_cs2_ani = cal_combined_cs2_ani(row["CS2"], ani)
            ani_list.append(cal_ani(assembly_id, output_directory, consensus_record_dict))
            combined_cs2_ani_list.append(combined_cs2_ani)
        else:
            ani_list.append(0)
            combined_cs2_ani_list.append(0)
            
    downloaded_assemblies['Consensus ANI'] = ani_list
    downloaded_assemblies["Combined CS2 and ANI (Sqrt(ANI)xCS2x100)"] = combined_cs2_ani_list
    downloaded_assemblies.sort_values(['CS2'], ascending=False).to_csv(os.path.join(output_directory, 'alignment.csv'), index=False)
    
    return downloaded_assemblies

def ani_summary(downloaded_assemblies, consensus_record_dict, output_directory, threads):  
    assembly_ids = []
    for idx, row in downloaded_assemblies.iterrows():
        assembly_ids.append(row['Assembly Accession ID'])
        
    pool = Pool(processes=threads)
    ani_list = pool.starmap(cal_ani, zip(assembly_ids, repeat(output_directory), repeat(consensus_record_dict)))
    pool.close()  
    pool.join()
    
    downloaded_assemblies['Consensus ANI'] = ani_list
    
    combined_cs2_ani_list = []
    for idx, row in downloaded_assemblies.iterrows():
        combined_cs2_ani = cal_combined_cs2_ani(row["Primary Score"], row['Consensus ANI'])
        combined_cs2_ani_list.append(combined_cs2_ani)
            
    
    downloaded_assemblies["Combined PS and ANI (Sqrt(ANI)xPSx100)"] = combined_cs2_ani_list
    downloaded_assemblies.sort_values(['Primary Score'], ascending=False).to_csv(os.path.join(output_directory,
                                                                                              'alignment.csv'),
                                                                                 index=False)
    
    return downloaded_assemblies

def samtools_merged_consensus(output_directory, threads):
    merged_bam = os.path.join(output_directory, 'bam_files', 'merged.sorted.bam')
    subprocess.run(['samtools', 'consensus', 
                    '--show-ins', 'no', 
                    '--show-del', 'yes',
                    '--min-MQ', str(20),
                    '-a',
                    '--mode', "simple",
                    '--threads', str(threads),
                    merged_bam, 
                    '-o', os.path.join(output_directory, 'merged_consensus.fasta')],
                  check=True)
    
    consensus_record_dict = SeqIO.to_dict(SeqIO.parse(os.path.join(output_directory, 'merged_consensus.fasta'), "fasta"))
    return consensus_record_dict

def cal_ani(assembly_id, output_directory, consensus_record_dict, ignore_del=False, del_count_as_match=False):
    reference_fasta = os.path.join(output_directory, 'reference_genomes', f'{assembly_id}.fasta')

    with open(reference_fasta, "r") as reference_handle:
        total_count = 0
        matched_count = 0

        for record in SeqIO.parse(reference_handle, "fasta"):
            if record.id in consensus_record_dict:
                consensus_seq = consensus_record_dict[record.id].seq
                reference_seq = record.seq
                
                for x, y in zip(consensus_seq, reference_seq):
                    if x == y:
                        matched_count += 1
                        
                total_count += len(reference_seq) - consensus_seq.count('N')
                
            else:
                continue
                #print("No alignment found:", record.id, record.description)
                    
    if total_count != 0:
        return matched_count/total_count
    else:
        return 0
    
def _cal_ani(assembly_id, output_directory, consensus_record_dict, ignore_del=False, del_count_as_match=False):
    reference_fasta = os.path.join(output_directory, 'reference_genomes', f'{assembly_id}.fasta')

    with open(reference_fasta, "r") as handle:
        total_count = 0
        matched_count = 0

        for record in SeqIO.parse(handle, "fasta"):
            if record.id in consensus_record_dict:
                for idx, base in enumerate(record.seq):
                    if consensus_record_dict[record.id][idx] != 'N':
                        if consensus_record_dict[record.id][idx] == '*' and ignore_del:
                            continue
                        elif consensus_record_dict[record.id][idx] == '*':
                            total_count += 1
                            if del_count_as_match:
                                matched_count += 1
                        else:
                            total_count += 1
                            if consensus_record_dict[record.id][idx] == base:
                                matched_count += 1
            else:
                continue
                #print("No alignment found:", record.id, record.description)
                    
    if total_count != 0:
        return matched_count/total_count
    else:
        return 0
