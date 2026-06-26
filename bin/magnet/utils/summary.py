import os
import argparse
import pathlib
import sys
import subprocess
import math
import json
from io import StringIO
import warnings
from collections import defaultdict
from multiprocessing import Pool, Manager
from itertools import repeat
warnings.filterwarnings("ignore")

import pandas as pd
from ete3 import NCBITaxa
from Bio import SeqIO

def alignment_1_summary(downloaded_assemblies, output_directory):
    breadth_coverage_list = []
    depth_coverage_list = []
    expected_breadth_coverage_list = []
    coverage_score = []
    for assembly_id in downloaded_assemblies['Assembly Accession ID']:
        breadth_coverage, depth_coverage, expected_breadth_coverage = calculate_depth(assembly_id, output_directory, min_depth=1)
        breadth_coverage_list.append(breadth_coverage)
        depth_coverage_list.append(depth_coverage)
        expected_breadth_coverage_list.append(expected_breadth_coverage)
        if expected_breadth_coverage != 0:
            coverage_score.append(breadth_coverage/expected_breadth_coverage)
        else:
            coverage_score.append(0)
        
    downloaded_assemblies['Breadth Coverage'] = breadth_coverage_list
    downloaded_assemblies['Expected Coverage'] = expected_breadth_coverage_list
    downloaded_assemblies['Coverage Score'] = coverage_score
    downloaded_assemblies['Depth Coverage'] = depth_coverage_list
    
    downloaded_assemblies.to_csv(os.path.join(output_directory, 'alignment.csv'), index=False)
    
    return downloaded_assemblies

def alignment_2_summary(downloaded_assemblies, output_directory):
    breadth_coverage_dict, depth_coverage_dict, expected_breadth_coverage_dict = \
    calculate_depth_merged(list(downloaded_assemblies['Assembly Accession ID']), output_directory, min_depth=1)
    
    breadth_coverage_list = []
    depth_coverage_list = []
    expected_breadth_coverage_list = []
    coverage_score = []
    for assembly_id in downloaded_assemblies['Assembly Accession ID']:
        breadth_coverage_list.append(breadth_coverage_dict[assembly_id])
        depth_coverage_list.append(depth_coverage_dict[assembly_id])
        expected_breadth_coverage_list.append(expected_breadth_coverage_dict[assembly_id])

        if expected_breadth_coverage_dict[assembly_id] != 0:
            coverage_score.append(breadth_coverage_dict[assembly_id]/expected_breadth_coverage_dict[assembly_id])
        else:
            coverage_score.append(0)
            
    downloaded_assemblies['BC2'] = breadth_coverage_list
    downloaded_assemblies['EC2'] = expected_breadth_coverage_list
    downloaded_assemblies['CS2'] = coverage_score
    downloaded_assemblies['DC2'] = depth_coverage_list
    
    downloaded_assemblies.sort_values(['CS2'], ascending=False).to_csv(os.path.join(output_directory, 'alignment.csv'), index=False)
    
    return downloaded_assemblies

def alignment_summary(downloaded_assemblies, output_directory, seq2assembly_dict, include_supp=True):
    if include_supp:
        coverage_file_name = 'secondary_coverage.tsv'
        column_prefix = 'Secondary'
    else:
        coverage_file_name = 'primary_coverage.tsv'
        column_prefix = 'Primary'
    
    columns = [f'{column_prefix} Breadth',
               f'{column_prefix} Expected',
               f'{column_prefix} Score',
               f'{column_prefix} Depth']
    
    coverage_df = pd.read_csv(os.path.join(output_directory,
                                           'coverage_files',
                                           coverage_file_name),
                              sep='\t')
    
    taxa_records = defaultdict(lambda: defaultdict(int))
    for idx, row in coverage_df.iterrows():
        taxa_reference = seq2assembly_dict[row['#rname']]
        taxa_records[taxa_reference]['genome_length'] += row['endpos']
        taxa_records[taxa_reference]['reads_mapped'] += row['numreads']
        taxa_records[taxa_reference]['genome_totol_count'] += int(row['meandepth'] * row['endpos'])
        taxa_records[taxa_reference]['covbases'] += row['covbases']
    
    breadth_coverage_list = []
    depth_coverage_list = []
    expected_breadth_coverage_list = []
    coverage_score = []
    for assembly_id in downloaded_assemblies['Assembly Accession ID']:
        breadth_coverage, depth_coverage, expected_breadth_coverage = calculate_depth(assembly_id, taxa_records)
        breadth_coverage_list.append(breadth_coverage)
        depth_coverage_list.append(depth_coverage)
        expected_breadth_coverage_list.append(expected_breadth_coverage)
        if expected_breadth_coverage != 0:
            coverage_score.append(min(breadth_coverage/expected_breadth_coverage, 1))
        else:
            coverage_score.append(0)
        
    downloaded_assemblies[columns[0]] = breadth_coverage_list
    downloaded_assemblies[columns[1]] = expected_breadth_coverage_list
    downloaded_assemblies[columns[2]] = coverage_score
    downloaded_assemblies[columns[3]] = depth_coverage_list
    
    downloaded_assemblies.to_csv(os.path.join(output_directory, 'alignment.csv'), index=False)
    
    return downloaded_assemblies

def calculate_depth(assembly_id, taxa_records):
    genome_length = taxa_records[assembly_id]['genome_length']
    covbases = taxa_records[assembly_id]['covbases']
    genome_totol_count = taxa_records[assembly_id]['genome_totol_count']
    reads_mapped = taxa_records[assembly_id]['reads_mapped']
    
    if genome_totol_count == 0 or reads_mapped == 0:
        breadth_coverage = 0
        depth_coverage = 0
        expected_breadth_coverage = 0
    else:
        breadth_coverage = covbases/genome_length
        depth_coverage = genome_totol_count/covbases
        expected_breadth_coverage, std = get_expected_coverage(genome_length, reads_mapped, genome_totol_count)
    
    return breadth_coverage, depth_coverage, expected_breadth_coverage

def _calculate_depth(assembly_id, output_directory, min_depth=1):
    depth_file = os.path.join(output_directory, 'depth_files', f"{assembly_id}.depth")
    
    genome_length, genome_ids = parse_reference_fasta(assembly_id, output_directory)
    mapping_stats = get_mapping_stats(assembly_id, output_directory)
    reads_mapped = mapping_stats['reads mapped']

    genome_pos_count = 0
    genome_totol_count = 0
    with open(depth_file, "r") as depth:
        for line in depth.readlines():
            
            genome_id = line.split("\t")[0]
            pos = int(line.split("\t")[1])
            depth = int(line.strip().split("\t")[2])

            if depth >= min_depth and genome_id in genome_ids:
                genome_pos_count += 1
                genome_totol_count += depth
    
    if genome_totol_count == 0 or reads_mapped == 0:
        breadth_coverage = 0
        depth_coverage = 0
        expected_breadth_coverage = 0
    else:
        breadth_coverage = genome_pos_count/genome_length
        depth_coverage = genome_totol_count/genome_pos_count
        expected_breadth_coverage, std = get_expected_coverage(genome_length, reads_mapped, genome_totol_count)
    
    return breadth_coverage, depth_coverage, expected_breadth_coverage

def _calculate_depth_merged(assembly_ids, output_directory, min_depth=1):
    depth_file = os.path.join(output_directory, 'depth_files', f"merged.depth")
    
    genome_id_pos_count = defaultdict(int)
    genome_id_totol_count = defaultdict(int)
    
    with open(depth_file, "r") as depth:
        for line in depth.readlines():            
            genome_id = line.split("\t")[0]
            pos = int(line.split("\t")[1])
            depth = int(line.strip().split("\t")[2])

            if depth >= min_depth:
                genome_id_pos_count[genome_id] += 1
                genome_id_totol_count[genome_id] += depth
    
    breadth_coverage_dict = defaultdict(float)
    depth_coverage_dict = defaultdict(float)
    reads_mapped_dict = defaultdict(int)
    expected_breadth_coverage_dict = defaultdict(float)
    
    for assembly_id in assembly_ids:
        genome_length, genome_ids = parse_reference_fasta(assembly_id, output_directory)
        
        genome_pos_count = 0
        genome_totol_count = 0
        for genome_id in genome_ids:
            genome_pos_count += genome_id_pos_count[genome_id]
            genome_totol_count += genome_id_totol_count[genome_id]
            if genome_id_pos_count[genome_id] > 0:
                mapping_stats = get_mapping_stats('merged', output_directory, genome_id)
                reads_mapped_dict[assembly_id] += mapping_stats['reads mapped']
            
        if genome_totol_count > 0 and reads_mapped_dict[assembly_id] > 0:
            breadth_coverage_dict[assembly_id] = genome_pos_count/genome_length
            depth_coverage_dict[assembly_id] = genome_totol_count/genome_pos_count
            expected_breadth_coverage, std = get_expected_coverage(genome_length, reads_mapped_dict[assembly_id], genome_totol_count)
            expected_breadth_coverage_dict[assembly_id] = expected_breadth_coverage
        elif genome_totol_count > 0 and reads_mapped_dict[assembly_id] == 0:
            print("WARNING: Inconsistency between samtools depth and samtools stats for:", assembly_id)
            print("Genome ID", "\t", "Stats Mapped Read Count", "\t", "# of Covered Positions")
            for genome_id in genome_ids:
                print(genome_id, "\t", get_mapping_stats('merged', output_directory, genome_id)['reads mapped'], "\t", genome_id_totol_count[genome_id])
            
    return breadth_coverage_dict, depth_coverage_dict, expected_breadth_coverage_dict
    
def parse_reference_fasta(assembly_id, output_directory):
    reference_fasta = os.path.join(output_directory, 'reference_genomes', f'{assembly_id}.fasta')
    
    with open(reference_fasta, "r") as handle:
        genome_length = 0
        genome_ids = []
        
        for record in SeqIO.parse(handle, "fasta"):      
            genome_length += len(record.seq)
            genome_ids.append(record.id)
    
    return genome_length, genome_ids

def merge_reference_fasta(assembly_ids, output_directory):
    merged_fasta = os.path.join(output_directory, 'reference_genomes', f'merged.fasta')
    
    seq_records = []
    for assembly_id in assembly_ids:
        reference_fasta = os.path.join(output_directory, 'reference_genomes', f'{assembly_id}.fasta')

        with open(reference_fasta, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):   
                seq_records.append(record)
                    
    with open(merged_fasta, "w") as output_handle:
        SeqIO.write(seq_records, output_handle, "fasta")
        
    return merged_fasta

def get_mapping_stats(assembly_id, output_directory, genome_id=None):
    bam_file = os.path.join(output_directory, 'bam_files', f"{assembly_id}.sorted.bam")
    
    if genome_id is not None:
        stats_res = subprocess.Popen(['samtools', 'stats', bam_file, genome_id],
                               stdout=subprocess.PIPE)
    else:
        stats_res = subprocess.Popen(['samtools', 'stats', bam_file],
                               stdout=subprocess.PIPE)

    grep_res = subprocess.Popen(['grep', '^SN'],
                              stdin=stats_res.stdout,
                              stdout=subprocess.PIPE)

    mapping_res = subprocess.run(['cut', '-f', '2-'],
                                 check=True,
                                 universal_newlines=True,
                                 stdin = grep_res.stdout,
                                 stdout=subprocess.PIPE)

    mapping_stats = dict()
    for line in mapping_res.stdout.strip().split('\n'):
        attribute, value = line.split(':\t')
        mapping_stats[attribute] = float(value.split('\t')[0])
        
    return mapping_stats

def get_expected_coverage(genome_length, reads_mapped, genome_totol_count):
    mean_mapping_length = genome_totol_count/reads_mapped
    
    N = genome_length/mean_mapping_length
    x = reads_mapped
    
    expected_M = N*(1-((1-1/N)**x))
    variance = N*((1-1/N)**x) + (N**2)*(1-1/N)*((1-2/N)**x)-(N**2)*((1-1/N)**(2*x))
    
    expected_coverage = expected_M/N
    try:
        std = math.sqrt(variance)
    except ValueError:
        std = 0
    return expected_coverage, std

def _call_present_absent(downloaded_assemblies):
    status = []
    cs = downloaded_assemblies["Coverage Score"].tolist()
    cs2 = downloaded_assemblies["CS2"].tolist()
    consensus_ani = downloaded_assemblies["Consensus ANI"].tolist()
    bc = downloaded_assemblies["Breadth Coverage"].tolist()
    bc2 = downloaded_assemblies['BC2'].tolist()
    
    for i in range(len(cs)):
        if consensus_ani[i] > 0.97 and bc2[i] > 0.2:
            status.append("Present")
        else:
            if cs[i] - cs2[i] < 0.1 and bc[i] - bc2[i] < 0.1 and consensus_ani[i] >= 0.70:
                if consensus_ani[i] > 0.9:
                    status.append("Present")
                else:
                    status.append("Genus Present")
            else:
                status.append("Absent")
            
    downloaded_assemblies["Presence/Absence"] = status
    return downloaded_assemblies

def call_present_absent(downloaded_assemblies, min_coverage_score):
    status = []
    cs = downloaded_assemblies["Secondary Score"].tolist()
    cs2 = downloaded_assemblies["Primary Score"].tolist()
    consensus_ani = downloaded_assemblies["Consensus ANI"].tolist()
    bc = downloaded_assemblies["Secondary Breadth"].tolist()
    bc2 = downloaded_assemblies['Primary Breadth'].tolist()
    
    for i in range(len(cs)):
        if cs[i] < min_coverage_score:
            status.append("Absent")
        elif consensus_ani[i] > 0.97 and bc2[i] > 0.2:
            status.append("Present")
        elif cs[i] - cs2[i] < 0.1 and bc[i] - bc2[i] < 0.1:
            if consensus_ani[i] >= 0.70:
                status.append("Present")
            else:
                status.append("Absent")
        else:
            status.append("Absent")
    downloaded_assemblies["Presence/Absence"] = status
    return downloaded_assemblies

# def call_present_absent(downloaded_assemblies, min_coverage_score):
#     status = []
#     cs = downloaded_assemblies["Secondary Score"].tolist()
#     cs2 = downloaded_assemblies["Primary Score"].tolist()
#     consensus_ani = downloaded_assemblies["Consensus ANI"].tolist()
#     bc = downloaded_assemblies["Secondary Breadth"].tolist()
#     bc2 = downloaded_assemblies['Primary Breadth'].tolist()
    
#     for i in range(len(cs)):
#         if cs[i] < min_coverage_score:
#             status.append("Absent")
#         elif consensus_ani[i] > 0.97 and bc2[i] > 0.2:
#             status.append("Present")
#         else:
#             if cs[i] - cs2[i] < 0.1 and bc[i] - bc2[i] < 0.1 and consensus_ani[i] >= 0.70:
#                 if consensus_ani[i] > 0.9:
#                     status.append("Present")
#                 else:
#                     status.append("Genus Present")
#             else:
#                 status.append("Absent")
#     downloaded_assemblies["Presence/Absence"] = status
#     return downloaded_assemblies

# def call_present_absent(downloaded_assemblies, min_coverage_score):
#     status = []
#     cs = downloaded_assemblies["Secondary Score"].tolist()
#     cs2 = downloaded_assemblies["Primary Score"].tolist()
#     consensus_ani = downloaded_assemblies["Consensus ANI"].tolist()
#     bc = downloaded_assemblies["Secondary Breadth"].tolist()
#     bc2 = downloaded_assemblies['Primary Breadth'].tolist()
    
#     for i in range(len(cs)):
#         if cs[i] < min_coverage_score:
#             status.append("Absent")
#         elif consensus_ani[i] > 0.97:
#             status.append("Present")
#         else:
#             if cs2[i] - cs[i] > -0.1:
#                 status.append("Present")
#             else:
#                 status.append("Absent")
            
#     downloaded_assemblies["Presence/Absence"] = status
#     return downloaded_assemblies
