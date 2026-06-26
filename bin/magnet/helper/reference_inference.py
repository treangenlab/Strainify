import os
import argparse
import sys
import subprocess
import math
import pickle
import warnings
from collections import defaultdict
import json
import pandas as pd
from ete3 import NCBITaxa
from Bio import SeqIO

def run_datasets_summary(taxid, flags, assembly_level='complete_genome'):
    args = ['datasets', 'summary', 
            'genome',
            'taxon', str(taxid),
            '--assembly-level', assembly_level,
            '--limit', str(1)]
    
    args = args + flags
    
    grepOut  = subprocess.run(args,
                              universal_newlines=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL)
    
    return json.loads(grepOut.stdout.strip())

def run_datasets_download(taxid, assembly_accession, working_dir):
    grepOut  = subprocess.run(['datasets', 'download', 
                               'genome',
                               'accession', str(assembly_accession),
                               '--exclude-genomic-cds',
                               '--exclude-gff3',
                               '--exclude-protein',
                               '--exclude-rna',
                               '--filename', f"{working_dir}/{taxid}.zip",
                               '--no-progressbar'],
                              check=True, stderr=subprocess.DEVNULL)
    
    return grepOut.returncode

def cut_text(text, text_length):
    try:
        if len(text) > text_length-2:
            return text[0:text_length-2] + ".."
        else:
            return text.ljust(text_length)
    except TypeError:
        return "".ljust(text_length)
		
def download_reference_genome(taxid, working_dir):
    # check reference and representative genomes first
    #source = 'RefSeq'
    assembly_accession = None
    res =None
    representative = True # assembly_category
    assembly_level = 'complete'
    
    res = run_datasets_summary(taxid, ['--reference', '--assembly-source', 'refseq'])
    
    # check complete assembly in RefSeq
    if res['total_count'] == 0:
        representative = False
        res = run_datasets_summary(taxid, ['--assembly-source', 'refseq'])

    # check complete assembly
    if res['total_count'] == 0:
        res = run_datasets_summary(taxid, [])
        #source = 'Assembly'
    
    # check chromosome assembly
    if res['total_count'] == 0:
        assembly_level = 'chromosome'
        res = run_datasets_summary(taxid, [], assembly_level)
        
    # check contig assembly
    if res['total_count'] == 0:
        assembly_level='contig'
        res = run_datasets_summary(taxid, [], assembly_level)

    # check scaffold assembly
    if res['total_count'] == 0:
        assembly_level='scaffold'
        res = run_datasets_summary(taxid, [], assembly_level)
        
    if res['total_count'] == 1:
        assembly_accession = res['assemblies'][0]['assembly']['assembly_accession']
        try:
            source = res['assemblies'][0]['assembly']['annotation_metadata']['source']
        except KeyError:
            source = None
        assembly_level_ret = res['assemblies'][0]['assembly']['assembly_level']
        organism = res['assemblies'][0]['assembly']['org']['sci_name']
        try:
            strain = res['assemblies'][0]['assembly']['org']['strain']
        except KeyError:
            strain = None
		
        print(str(taxid).ljust(10), '\t', assembly_accession, '\t', cut_text(organism, 30), '\t', cut_text(strain, 10), '\t', assembly_level_ret)
        
        return_code = run_datasets_download(taxid, assembly_accession, working_dir)
        if not return_code:
            download_completed = True
        else:
            download_completed = False
            
    else:
        print(str(taxid).ljust(10), '\t', 'Genome Not Found.')
        assembly_accession = None
        source = None
        representative = None
        assembly_level_ret = None
        organism = None
        download_completed = False
        
    return taxid, assembly_accession, source, representative, assembly_level_ret, organism, strain, download_completed

def unpack(working_dir, output_dir):
    subprocess.run(['unzip', '-o', '-q', 
					f'{working_dir}/*.zip', 
					'-d', output_dir],
				   check=True)

def cat_reference_genome(reference_metadata, output_dir, reference_genome_path='reference_genomes'):
    if not os.path.exists(reference_genome_path):
        os.mkdir(reference_genome_path)
        
    downloaded_assemblies = reference_metadata[reference_metadata['Downloaded']]
    
    for assembly_id in downloaded_assemblies['Assembly Accession ID']:
        genome_path = os.path.join(output_dir, 'ncbi_dataset', 'data', assembly_id)
        output_fasta = os.path.join(reference_genome_path, f'{assembly_id}.fasta')

        subprocess.Popen(f'cat {genome_path}/*.fna', shell=True, stdout=open(os.path.join(output_fasta), "w"))

def get_species_taxid(taxid, ncbi_taxa_db, valid_kingdom):
    lineage = ncbi_taxa_db.get_lineage(taxid)
    if bool(set(lineage) & valid_kingdom):
        taxid2rank_dict = ncbi_taxa_db.get_rank(lineage)
        for lineage_taxid in taxid2rank_dict:
            if taxid2rank_dict[lineage_taxid] == 'species':
                return lineage_taxid
    return None

def filter_seqscreen_taxonomy(seqscreen_output, min_frac, ncbi_taxa_db, valid_kingdom):
    classification_result_df = pd.read_csv(os.path.join(seqscreen_output, 
                                                        'taxonomic_identification', 
                                                        'taxonomic_assignment', 
                                                        'taxonomic_results.txt'), 
                                           sep='\t')
    total_read_count, _ = classification_result_df.shape
    
    taxid_count_dict = defaultdict(int)
    taxid_species_lookup = dict()
    error_count = 0
    for taxid in classification_result_df['taxid']:
        try:
            taxid = int(taxid)
            try:
                species_taxid = taxid_species_lookup[taxid]
            except KeyError:
                species_taxid = get_species_taxid(taxid, ncbi_taxa_db, valid_kingdom)
                taxid_species_lookup[taxid] = species_taxid

            if species_taxid is not None:
                taxid_count_dict[species_taxid] += 1
        except ValueError:
            error_count += 1

    taxid_queries = []
    for taxid in taxid_count_dict:
        if taxid_count_dict[taxid] >= min_frac * total_read_count:
            taxid_queries.append(taxid)
            
    return taxid_queries

def prepare_reference_genomes(taxid_queries, output_directory, ncbi_taxa_db):
    working_dir = os.path.join(output_directory, 'ncbi_downloads')
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
        
    download_result = []
    for taxid in taxid_queries:
        download_result.append(download_reference_genome(taxid, working_dir))
        
    unpack(working_dir, output_directory)
    reference_metadata = pd.DataFrame(download_result,
                                      columns=['Taxonomy ID', 
                                               'Assembly Accession ID', 
                                               'Source Database', 
                                               'Is Representative', 
                                               'Assembly Level', 
                                               'Organism of Assembly',
                                               'Strain',
                                               'Downloaded'])
    
    taxonomy_name = []
    for taxid in reference_metadata['Taxonomy ID']:
        taxonomy_name.append(ncbi_taxa_db.get_taxid_translator([taxid])[taxid])
    reference_metadata['Species'] = taxonomy_name
    
    reference_metadata.to_csv(os.path.join(output_directory, 'reference_metadata.csv'), index=False)
    cat_reference_genome(reference_metadata, output_directory, reference_genome_path=os.path.join(output_directory, 'reference_genomes'))
    
    return reference_metadata

def sort_samfile(assembly_id, output_dir, min_mapq, num_cores):
    '''converting and sorting alignment files'''
    bam_files = os.path.join(output_dir, "bam_files")
    sam_files = os.path.join(output_dir, "sam_files")
    
    if not os.path.exists(bam_files):
        os.mkdir(bam_files)

    # covert sam file to binary bam file
    samtools_view_res = subprocess.Popen([
        "samtools",
        "view",
        "-@", str(num_cores),
        "--min-MQ", str(min_mapq),
        "-bS", os.path.join(sam_files, f"{assembly_id}.sam")],
        stdout=subprocess.PIPE)

    # sort the bam file 
    subprocess.run([
        "samtools",
        "sort",
        "-@", str(num_cores),
        "-o", os.path.join(bam_files, f"{assembly_id}.sorted.bam"),
        "-O", "BAM"],
        stdin=samtools_view_res.stdout,
        stderr=subprocess.DEVNULL,
        check=True)

    # indexing the sorted bam file 
    subprocess.run([
        "samtools",
        "index",
        os.path.join(bam_files, f"{assembly_id}.sorted.bam")],
                    check=True)

def run_minimap2(input_fasta, reference_file, assembly_id, output_dir, threads=20):
    sam_files = os.path.join(output_dir, "sam_files")

    if not os.path.exists(sam_files):
        os.mkdir(sam_files)

    subprocess.run(["minimap2", 
                    "-ax", "map-ont", 
                    reference_file, 
                    input_fasta, 
                    "--sam-hit-only",
                    "-o", os.path.join(sam_files, f"{assembly_id}.sam"),
                    "-t", str(threads)],
                    check=True, stderr=subprocess.DEVNULL)

def samtools_calculate_depth(assembly_id, output_dir):
    depth_files = os.path.join(output_dir, "depth_files")
    bam_files = os.path.join(output_dir, "bam_files")

    if not os.path.exists(depth_files):
        os.mkdir(depth_files)

    depth_file = os.path.join(depth_files, f"{assembly_id}.depth")
                
    subprocess.run([
        "samtools",
        "depth",
        os.path.join(bam_files, f"{assembly_id}.sorted.bam"),
        "-G", "SUPP"],
                    check=True,
                    stdout=open(os.path.join(depth_file), "w"))

def parse_reference_fasta(assembly_id, output_directory):
    reference_fasta = os.path.join(output_directory, 'reference_genomes', f'{assembly_id}.fasta')
    
    with open(reference_fasta, "r") as handle:
        genome_length = 0
        genome_ids = []
        
        for record in SeqIO.parse(handle, "fasta"):      
            genome_length += len(record.seq)
            genome_ids.append(record.id)
    
    return genome_length, genome_ids

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

def calculate_depth(assembly_id, output_directory, min_depth=1):
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

def calculate_depth_merged(assembly_ids, output_directory, min_depth=1):
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

def cal_ani(assembly_id, output_directory, consensus_record_dict, ignore_del=False, del_count_as_match=False):
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

def samtools_merged_consensus(output_directory, threads):
    merged_bam = os.path.join(output_directory, 'bam_files', 'merged.sorted.bam')
    subprocess.run(['samtools', 'consensus', 
                    '--show-ins', 'no', 
                    '--show-del', 'yes', 
                    '-a',
                    '--mode', "simple",
                    '--threads', str(threads),
                    merged_bam, 
                    '-o', os.path.join(output_directory, 'merged_consensus.fasta')],
                  check=True)
    
    consensus_record_dict = SeqIO.to_dict(SeqIO.parse(os.path.join(output_directory, 'merged_consensus.fasta'), "fasta"))
    return consensus_record_dict

def ani_summary(downloaded_assemblies, consensus_record_dict, output_directory):
    ani_list = []
    for idx, row in downloaded_assemblies.iterrows():
        if row['CS2'] != 0:
            assembly_id = row['Assembly Accession ID']
            ani_list.append(cal_ani(assembly_id, output_directory, consensus_record_dict, ignore_del=True))
        else:
            ani_list.append(0)
            
    downloaded_assemblies['Consensus ANI'] = ani_list
    downloaded_assemblies.sort_values(['CS2'], ascending=False).to_csv(os.path.join(output_directory, 'alignment.csv'), index=False)
    
    return downloaded_assemblies

def build_record_dict(sequences_db_f):
    sequence_ids = set()
    record_dict = dict()
    count = 0
    for record in SeqIO.parse(sequences_db_f, "fasta"):
        count += 1
        if record.id not in sequence_ids:
            sequence_ids.add(record.id)
            record_dict[record.id] = record
    print(f"WARNING: {count-len(record_dict)} duplicated record found among {count} sequences.")
    
    return record_dict

def local_taxid_fetch(taxid, output_directory, taxid2seqid_dict, record_dict):
    reference_genome_path=os.path.join(output_directory, 'reference_genomes')
    if not os.path.exists(reference_genome_path):
        os.mkdir(reference_genome_path)
    
    records = []
    for sequence_id in taxid2seqid_dict[str(taxid)]:
        records.append(record_dict[sequence_id])
        
    if len(records) > 0:
        with open(os.path.join(reference_genome_path, f"taxid_{taxid}.fasta"), "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")
        return taxid, f"taxid_{taxid}", "Local Database", "N/A", "N/A", "N/A", "N/A", True
    else:
        return taxid, "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", False
	
def prepare_reference_genomes_offline(taxid_queries, output_directory, sequences_db_f, mapping_f, ncbi_taxa_db):
    with open(mapping_f, 'rb') as handle:
        taxid2seqid_dict = pickle.load(handle)
        
    record_dict = build_record_dict(sequences_db_f)
    
    download_result = []
    for taxid in taxid_queries:
        download_result.append(local_taxid_fetch(taxid, output_directory, taxid2seqid_dict, record_dict))
        
    reference_metadata = pd.DataFrame(download_result,
                                      columns=['Taxonomy ID', 
                                               'Assembly Accession ID', 
                                               'Source Database', 
                                               'Is Representative', 
                                               'Assembly Level', 
                                               'Organism of Assembly',
                                               'Strain',
                                               'Downloaded'])
    
    taxonomy_name = []
    for taxid in reference_metadata['Taxonomy ID']:
        taxonomy_name.append(ncbi_taxa_db.get_taxid_translator([taxid])[taxid])
    reference_metadata['Species'] = taxonomy_name
    reference_metadata.to_csv(os.path.join(output_directory, 'reference_metadata.csv'), index=False)
    
    return reference_metadata

def main(argv):
	'''main function'''
	warnings.filterwarnings("ignore")

	parser = argparse.ArgumentParser(description="Repeated Read Alignment Module")
	
	parser.add_argument("-i", "--fasta", type=str, help="Sequences in fasta format.", required=True)
	parser.add_argument("-o", "--output", type=str, help="Output directory.", required=True)
	parser.add_argument("-w", "--working", type=str, help="Working directory.", required=True)
	parser.add_argument("-d", "--databases", type=str, help="SeqScreenDB", required=True)
	parser.add_argument("-m", "--min-frac", type=float, default=0.002,
						help="minimum fraction of assigned reads for a species to be included \
						in the first alignment process. [0.002]", )
	parser.add_argument("-c", "--min-coverage-score", type=float, default=0.7,
						help="minimum coverage score for a species to be included \
						in the second alignment process. [0.7]", )
	parser.add_argument("--min-mapq", type=int, default=20,
						help="minimum mapping quality for the alignment step. [20]", )
	parser.add_argument("--online", action='store_true', help="Use online mode for searching reference genome. Requires internet access.")
	parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads. [1]")
	parser.set_defaults(online=False)
	
	args = parser.parse_args()
	# seqscreen_output = "/home/Users/yl181/seqscreen_nano/ZymoBIOMICS.STD.Even.ont.seqscreen"
	seqscreen_output = args.output
	# input_fasta = "/home/Users/yl181/seqscreen_nano/ZymoBIOMICS.STD.Even.ont.raw_sequences/ERR3152364.downsampled.fasta"
	input_fasta = args.fasta
	threads = args.threads
	min_mapq = args.min_mapq
	
	# working_directory = 'working'
	working_directory = args.working
	if not os.path.exists(working_directory):
		os.mkdir(working_directory)
		
	# TODO, point this path to DB file
	ete3db = os.path.join(args.databases, "reference_inference", "taxa.sqlite")
	sequences_db_f = os.path.join(args.databases, "bowtie2", "blacklist.seqs.nt.fna")
	mapping_f = os.path.join(args.databases, "reference_inference", "taxid2seqid.pickle")

	ncbi_taxa_db = NCBITaxa(dbfile=ete3db)
	
	# valid_kingdom = set(bacteria, archaea, viruses, fungi)
	taxid_queries = filter_seqscreen_taxonomy(seqscreen_output, 
											  min_frac=args.min_frac, 
											  ncbi_taxa_db=ncbi_taxa_db, 
											  valid_kingdom={2, 4751, 2157, 10239})
	
	if args.online:
		reference_metadata = prepare_reference_genomes(taxid_queries, working_directory, ncbi_taxa_db)
	else:
		reference_metadata = prepare_reference_genomes_offline(taxid_queries, working_directory, sequences_db_f, mapping_f, ncbi_taxa_db)
	
	downloaded_assemblies = reference_metadata[reference_metadata['Downloaded']]
	for assembly_id in downloaded_assemblies['Assembly Accession ID']:
		reference_fasta = os.path.join(working_directory, 'reference_genomes', f'{assembly_id}.fasta')
		run_minimap2(input_fasta, reference_fasta, assembly_id, working_directory, threads=threads)
		sort_samfile(assembly_id, working_directory, min_mapq, threads)
		samtools_calculate_depth(assembly_id, working_directory)
	
	downloaded_assemblies = alignment_1_summary(downloaded_assemblies, working_directory)
	
	# Filtered the assemblies by coverage score
	filtered_assemblies = list(downloaded_assemblies[downloaded_assemblies['Coverage Score'] >= args.min_coverage_score]['Assembly Accession ID'])
	
	reference_fasta = merge_reference_fasta(filtered_assemblies, working_directory)
	run_minimap2(input_fasta, reference_fasta, 'merged', working_directory, threads=threads)
	sort_samfile('merged', working_directory, min_mapq, threads)
	samtools_calculate_depth('merged', working_directory)

	downloaded_assemblies = alignment_2_summary(downloaded_assemblies, working_directory)
	
	consensus_record_dict = samtools_merged_consensus(working_directory, threads)
	downloaded_assemblies = ani_summary(downloaded_assemblies, consensus_record_dict, working_directory)
	
	downloaded_assemblies.sort_values(['CS2'], ascending=False).to_csv(os.path.join(seqscreen_output, "taxonomic_identification", "taxonomic_assignment",\
                                                              'alignment.csv'), index=False)
	
if __name__ == "__main__":
	main(sys.argv[1:])
	
