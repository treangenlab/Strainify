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

def sort_samfile(assembly_id, aligner_output, output_dir, min_mapq, threads=20):
    '''converting and sorting alignment files'''
    bam_files = os.path.join(output_dir, "bam_files")
    # sam_files = os.path.join(output_dir, "sam_files")
    
    if not os.path.exists(bam_files):
        os.mkdir(bam_files)

    # covert sam file to binary bam file
    samtools_view_res = subprocess.Popen([
        "samtools",
        "view",
        "-@", str(threads),
        "--min-MQ", str(min_mapq),
        "-bS"], #os.path.join(sam_files, f"{assembly_id}.sam")
        stdin=aligner_output.stdout,
        stdout=subprocess.PIPE)

    # sort the bam file 
    subprocess.run([
        "samtools",
        "sort",
        "-@", str(threads),
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

def run_minimap2(input_fastq, reference_file, assembly_id, output_dir, threads=20):
    minimap2_output = subprocess.Popen(["minimap2", 
                                        "-ax", "map-ont", 
                                        reference_file, 
                                        input_fastq,
                                        "-N", str(50),
                                        "--sam-hit-only",
                                        "-t", str(threads)],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.DEVNULL)
    
    return minimap2_output

def run_bwa(input_fastq_1, input_fastq_2, reference_file, assembly_id, output_dir, threads=20):
    '''map the reads to the reference'''
    subprocess.run(["bwa", "index", reference_file], 
                    stdout=open(os.path.join(output_dir, "bwa_mem.log"), "a"),
                    stderr=open(os.path.join(output_dir, "bwa_mem.err"), "a"),
                    check=True)

    if input_fastq_2:
        bwa_mem_output = subprocess.Popen([
            "bwa",
            "mem",
            "-a",
            "-t", str(threads),
            reference_file,
            input_fastq_1,
            input_fastq_2],
                stdout=subprocess.PIPE,
                stderr=open(os.path.join(output_dir, "bwa_mem.err"), "a"))
    else:
        bwa_mem_output = subprocess.Popen([
            "bwa",
            "mem",
            "-a",
            "-t", str(threads),
            reference_file,
            input_fastq_1],
                stdout=subprocess.PIPE,
                stderr=open(os.path.join(output_dir, "bwa_mem.err"), "a"))
    
    return bwa_mem_output

def run_bowtie2(input_fastq_1, input_fastq_2, reference_file, assembly_id, output_dir, threads=20):
    '''map the reads to the reference'''
    subprocess.run(["bowtie2-build", reference_file, os.path.join(output_dir, assembly_id)], 
                    stdout=open(os.path.join(output_dir, "bowtie2.log"), "a"),
                    stderr=open(os.path.join(output_dir, "bowtie2.err"), "a"),
                    check=True)

    if input_fastq_2:
        bowtie2_output = subprocess.Popen([
            "bowtie2",
            "--local",
            "-a", #str(50),
            "--threads", str(threads),
            "-x", os.path.join(output_dir, assembly_id),
            "-1", input_fastq_1,
            "-2", input_fastq_2],
                stdout=subprocess.PIPE,
                stderr=open(os.path.join(output_dir, "bowtie2.err"), "a"))
    else:
        bowtie2_output = subprocess.Popen([
            "bowtie2",
            "--local",
            "-a", #str(50),
            "--threads", str(threads),
            "-x", os.path.join(output_dir, assembly_id),
            "-1", input_fastq_1],
                stdout=subprocess.PIPE,
                stderr=open(os.path.join(output_dir, "bowtie2.err"), "a"))
    
    return bowtie2_output
        
def _samtools_calculate_depth(assembly_id, output_dir, exclude_supp=True):
    depth_files = os.path.join(output_dir, "depth_files")
    bam_files = os.path.join(output_dir, "bam_files")

    if not os.path.exists(depth_files):
        os.mkdir(depth_files)

    depth_file = os.path.join(depth_files, f"{assembly_id}.depth")

    command = ["samtools",
               "depth",
               os.path.join(bam_files, f"{assembly_id}.sorted.bam")]
	
    if exclude_supp:
        command += ["-G", "0x800"]

    subprocess.run(command,
                    check=True,
                    stdout=open(os.path.join(depth_file), "w"))

def samtools_calculate_coverage(output_dir, include_supp=False):
    coverage_files = os.path.join(output_dir, "coverage_files")
    bam_files = os.path.join(output_dir, "bam_files")

    command = ["samtools",
               "coverage",
               os.path.join(bam_files, f"merged.sorted.bam")]
    
    if include_supp:
        # samtools coverage --ff UNMAP,QCFAIL,DUP -q 0 merged.sorted.bam > secondary_coverage.tsv
        coverage_file = os.path.join(coverage_files, f"secondary_coverage.tsv")
        command += ['--ff', 'UNMAP,QCFAIL,DUP', '-q', str(0)]
    else:
        # samtools coverage merged.sorted.bam -q 20 > primary_coverage.tsv
        coverage_file = os.path.join(coverage_files, f"primary_coverage.tsv")
        command += ['-q', str(20)]
        # command += ['--ff', 'UNMAP,QCFAIL,DUP', '-q', str(20)]

    subprocess.run(command,
                    check=True,
                    stdout=open(coverage_file, "w"))
	
def main(argv):
	'''main function'''
	parser = argparse.ArgumentParser(description="Read Alignment and Coverage Calculation.")
	parser.add_argument("-1", "--fastq1", type=str, help="Sequences in fastq format", required=True)
	parser.add_argument("-2", "--fastq2", type=str, help="Sequences in fastq format, for paired-end reads")
	parser.add_argument("-a", "--aligner", type=str, help="Aligner [minimap2, bwa_mem]", required=True)
	parser.add_argument("-o", "--output", type=str, help="Output directory", required=True)
	parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads. [1]")
	
	args = parser.parse_args()
	
	output_dir = args.output
	input_fastq_1 = args.fastq1
	input_fastq_2 = args.fastq2
	threads = args.threads
	
	reference_metadata = pd.read_csv(os.path.join(output_dir, 'reference_metadata.csv'))
	reference_genome_path = os.path.join(output_dir, 'reference_genomes')

	downloaded_assemblies = reference_metadata[reference_metadata['Downloaded']]

	for assembly_id in downloaded_assemblies['Assembly Accession ID']:
		reference_fasta = os.path.join(reference_genome_path, f'{assembly_id}.fasta')
		if args.aligner == 'minimap2':
			run_minimap2(input_fastq_1, reference_fasta, assembly_id, output_dir, threads=threads)
		elif args.aligner == 'bwa_mem':
			run_bwa(input_fastq_1, input_fastq_2, reference_fasta, assembly_id, output_dir, threads=threads)
		else:
			print('Supported Aligners are minimap2 or bwa_mem.')
			return 1
		sort_samfile(assembly_id, output_dir, threads)
		samtools_calculate_depth(assembly_id, output_dir)

if __name__ == "__main__":
	main(sys.argv[1:])
