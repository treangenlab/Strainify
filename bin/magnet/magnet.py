import os
import argparse
import pathlib
import sys
import subprocess
import math
import warnings
from collections import defaultdict
from multiprocessing import Pool, Manager
from itertools import repeat
warnings.filterwarnings("ignore")
import pandas as pd
from ete3 import NCBITaxa
from Bio import SeqIO
from utils.reference_finder import prepare_reference_genomes
from utils.alignment import run_minimap2, run_bwa, sort_samfile
from utils.summary import merge_reference_fasta, call_present_absent
from utils.ani import ani_summary
from utils.input_parsing import parsing_input_f, filter_input_df, get_seq2assembly_dict
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from ast import literal_eval

def samtools_calculate_coverage(output_dir, include_supp=False):
    coverage_files = os.path.join(output_dir, "coverage_files")
    bam_files = os.path.join(output_dir, "bam_files")
    command = ["samtools",
               "coverage", "--no-header",
               os.path.join(bam_files, f"merged.sorted.bam")]
    if include_supp:
        # samtools coverage --ff UNMAP,QCFAIL,DUP -q 0 merged.sorted.bam > secondary_coverage.tsv
        coverage_file = os.path.join(coverage_files, f"secondary_coverage.tsv")
        command += ['--ff', 'UNMAP,QCFAIL,DUP', '-q', str(1)]
    else:
        # samtools coverage merged.sorted.bam -q 20 > primary_coverage.tsv
        coverage_file = os.path.join(coverage_files, f"primary_coverage.tsv")
        command += ['-q', str(20)]
        # command += ['--ff', 'UNMAP,QCFAIL,DUP', '-q', str(20)]

    subprocess.run(command,
                   check=True,
                   stdout=open(coverage_file, "w"))

def parse_fastani_line(line):
    ret_list = []
    for item in line.strip().split('\t')[1:]:
        if item == 'NA':
            ret_list.append(float(0))
        else:
            ret_list.append(float(item))
    return ret_list

def find_representative_genome(fastani_path, fastani_assemblies, downloaded_assemblies):
    fastani_result = os.path.join(fastani_path, f'pairwise_ani.matrix')
    #fastani_assemblies = downloaded_assemblies[downloaded_assemblies['Genus Taxid'] == genus_taxid]['Assembly Accession ID'].values
    ani_matrix = []
    with open(fastani_result, 'r') as fastani_out_f:
        lines = fastani_out_f.readlines()
        num_seqs = int(lines[0].strip())
        for idx, line in enumerate(lines[1:]):
            ani_matrix.append(parse_fastani_line(line)+list(np.ones(num_seqs-idx)*100))

    dist_nparray = 1-np.array(ani_matrix)/100
    dist_df = pd.DataFrame(dist_nparray.T + dist_nparray,
                           columns=fastani_assemblies,
                           index=fastani_assemblies)

    model = AgglomerativeClustering(metric='precomputed', n_clusters=None, compute_full_tree=True,
                                    linkage='complete',
                                    distance_threshold=0.05).fit(dist_df.values)
    cluster_df = dist_df.copy()
    #print(cluster_df)
    cluster_df['Cluster Label'] = model.labels_

    representative_genomes = defaultdict(list)
    member2representative = dict()
    for cluster_idx in cluster_df['Cluster Label'].unique():
        cluster_members = cluster_df[cluster_df['Cluster Label'] == cluster_idx].index.to_list()
        selected_df = downloaded_assemblies[downloaded_assemblies['Assembly Accession ID'].isin(cluster_members)].copy()
        if selected_df[selected_df['Assembly Level'] == 'Complete Genome'].shape[0] == 1:
            representative_genome = selected_df[selected_df['Assembly Level'] == 'Complete Genome']['Assembly Accession ID'].values[0]
        elif selected_df[selected_df['Assembly Level'] == 'Complete Genome'].shape[0] > 1:
            complete_genomes = list(selected_df[selected_df['Assembly Level'] == 'Complete Genome']['Assembly Accession ID'].values)
            representative_genome = dist_df.loc[complete_genomes].sum().idxmin()
        else:
            representative_genome = dist_df.loc[cluster_members].sum().idxmin()
        representative_genomes[representative_genome] = cluster_members
        for member in cluster_members:
            member2representative[member] = representative_genome

    return representative_genomes, member2representative

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
                              sep='\t',
                              header=None,
                              names=['rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq'])

    taxa_records = defaultdict(lambda: defaultdict(int))
    for idx, row in coverage_df.iterrows():
        taxa_reference = seq2assembly_dict[row['rname']]
        taxa_records[taxa_reference]['genome_length'] += row['endpos']
        taxa_records[taxa_reference]['reads_mapped'] += row['numreads']
        taxa_records[taxa_reference]['genome_totol_count'] += int(row['meandepth'] * row['endpos'])
        taxa_records[taxa_reference]['covbases'] += row['covbases']

    breadth_coverage_list = []
    depth_coverage_list = []
    expected_breadth_coverage_list = []
    coverage_score = []
    reads_mapped_list = []
    for assembly_id in downloaded_assemblies['Assembly Accession ID']:
        breadth_coverage, depth_coverage, expected_breadth_coverage = calculate_depth(assembly_id, taxa_records)
        breadth_coverage_list.append(breadth_coverage)
        depth_coverage_list.append(depth_coverage)
        expected_breadth_coverage_list.append(expected_breadth_coverage)
        reads_mapped_list.append(taxa_records[assembly_id]['reads_mapped'])
        if expected_breadth_coverage != 0:
            coverage_score.append(min(breadth_coverage/expected_breadth_coverage, 1))
        else:
            coverage_score.append(0)

    downloaded_assemblies[columns[0]] = breadth_coverage_list
    downloaded_assemblies[columns[1]] = expected_breadth_coverage_list
    downloaded_assemblies[columns[2]] = coverage_score
    downloaded_assemblies[columns[3]] = depth_coverage_list
    downloaded_assemblies[f'{column_prefix} Reads'] = reads_mapped_list
    downloaded_assemblies.to_csv(os.path.join(output_directory, 'alignment.csv'), index=False)

    return downloaded_assemblies


def apply_absolute_gate(representative_df, min_breadth=0.0, min_reads=0,
                        presence_col='Presence/Absence'):
    """
    Absolute evidence gate applied AFTER call_present_absent.

    The Present/Absent call is a ratio (breadth / expected breadth), which
    saturates to ~1.0 at tiny read counts -- so a genome covered by a handful of
    reads (breadth ~1e-4) can score 1.0 and be called Present. This gate requires
    a minimum absolute amount of evidence (primary breadth and/or primary mapped
    reads) before a Present call is trusted. It ONLY downgrades Present -> Absent,
    never the reverse, so lowering --min-covscore to recover borderline true
    positives no longer admits these near-empty false positives.

    Disabled when both thresholds are 0 (default), preserving original behavior.
    """
    if min_breadth <= 0 and min_reads <= 0:
        return representative_df

    representative_df = representative_df.copy()
    passes = pd.Series(True, index=representative_df.index)
    if min_breadth > 0 and 'Primary Breadth' in representative_df.columns:
        passes &= representative_df['Primary Breadth'].astype(float) >= float(min_breadth)
    if min_reads > 0 and 'Primary Reads' in representative_df.columns:
        passes &= representative_df['Primary Reads'].astype(float) >= float(min_reads)

    representative_df['Absolute Gate Pass'] = passes
    if presence_col in representative_df.columns:
        present = representative_df[presence_col].astype(str).str.strip().str.lower() == 'present'
        representative_df.loc[present & (~passes), presence_col] = 'Absent'
    return representative_df


def apply_primary_evidence_rescue(representative_df, min_reads=0, min_breadth=0.0,
                                  min_ani=0.95, presence_col='Presence/Absence'):
    """
    Primary-evidence rescue applied AFTER call_present_absent (and the absolute
    gate). It is the mirror image of the gate: it ONLY upgrades Absent -> Present.

    Rationale: the coverage-score call thresholds on the breadth/expected ratio,
    which normalizes away the signal from a genuinely-present but LOW-ABUNDANCE
    strain. On a near-identical panel such a strain's reads cluster (observed
    breadth << expected), so its score falls below --min-covscore and it is
    called Absent -- even though it has thousands of UNIQUELY-MAPPED (primary)
    reads, i.e. reads that map to it and to no other genome in the panel. That
    unique read count is direct, competition-proof evidence of presence, and it
    cleanly separates true low-abundance strains (thousands of primary reads)
    from near-identical false positives (a handful). Using per-genome alignment
    to "fix" this instead would remove read competition and let every genome's
    conserved core light up, inflating false positives -- which is why this
    operates on the existing competitive (merged) alignment's primary signal.

    A genome called Absent is upgraded to Present when it has at least
    `min_reads` primary mapped reads OR at least `min_breadth` primary breadth,
    AND a Consensus ANI of at least `min_ani`. Set thresholds well above the
    saturation floor used by the absolute gate so the two never conflict.

    Disabled when both read/breadth thresholds are 0 (default).
    """
    if min_reads <= 0 and min_breadth <= 0:
        return representative_df

    representative_df = representative_df.copy()
    evidence = pd.Series(False, index=representative_df.index)
    if min_reads > 0 and 'Primary Reads' in representative_df.columns:
        evidence |= representative_df['Primary Reads'].astype(float) >= float(min_reads)
    if min_breadth > 0 and 'Primary Breadth' in representative_df.columns:
        evidence |= representative_df['Primary Breadth'].astype(float) >= float(min_breadth)
    if 'Consensus ANI' in representative_df.columns:
        evidence &= representative_df['Consensus ANI'].astype(float) >= float(min_ani)

    if presence_col in representative_df.columns:
        absent = representative_df[presence_col].astype(str).str.strip().str.lower() == 'absent'
        rescued = absent & evidence
        representative_df['Primary Evidence Rescue'] = rescued
        representative_df.loc[rescued, presence_col] = 'Present'
    else:
        representative_df['Primary Evidence Rescue'] = False
    return representative_df


_BWA_INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")


def _bwa_index_sentinel(reference_file):
    return str(reference_file) + ".bwa_index_complete"


def bwa_index_is_complete(reference_file):
    '''
    True only if a COMPLETE, UP-TO-DATE bwa index exists for reference_file:
      - the completion sentinel is present (so the build wasn't interrupted),
      - all five index files exist and are non-empty,
      - every index file is at least as new as the reference FASTA.
    This catches partial indexes (killed mid-build) and stale indexes (the
    reference was regenerated after indexing) -- both of which make bwa mem emit
    `bns_fetch_seq` errors and thrash because the BWT/SA no longer matches the
    reference annotation.
    '''
    ref = str(reference_file)
    if not os.path.exists(ref):
        return False
    if not os.path.exists(_bwa_index_sentinel(ref)):
        return False
    try:
        ref_mtime = os.path.getmtime(ref)
        for suf in _BWA_INDEX_SUFFIXES:
            f = ref + suf
            if not os.path.exists(f) or os.path.getsize(f) == 0:
                return False
            if os.path.getmtime(f) < ref_mtime:
                return False
    except OSError:
        return False
    return True


def build_bwa_index(reference_file, output_dir):
    '''Build a bwa index and drop a completion sentinel only on success.'''
    sentinel = _bwa_index_sentinel(reference_file)
    if os.path.exists(sentinel):
        os.remove(sentinel)  # invalidate before (re)building
    subprocess.run(["bwa", "index", reference_file],
                   stdout=open(os.path.join(output_dir, "bwa_mem.log"), "a"),
                   stderr=open(os.path.join(output_dir, "bwa_mem.err"), "a"),
                   check=True)
    with open(sentinel, "w") as fh:
        fh.write("ok\n")


def run_bwa_index_once(input_fastq_1, input_fastq_2, reference_file, assembly_id, output_dir, threads=20):
    '''
    Paired-end bwa-mem alignment, identical to utils.alignment.run_bwa EXCEPT the
    `bwa index` step runs only when a complete, up-to-date index is NOT already
    present (see bwa_index_is_complete). With --shared-reference the merged
    reference and its index live in the shared directory and are built once, so
    every sample skips indexing and maps against the same prebuilt index. A
    stale or partial index is detected here and rebuilt rather than trusted.
    '''
    if not bwa_index_is_complete(reference_file):
        build_bwa_index(reference_file, output_dir)

    bwa_cmd = ["bwa", "mem", "-a", "-t", str(threads), reference_file, input_fastq_1]
    if input_fastq_2:
        bwa_cmd.append(input_fastq_2)

    bwa_mem_output = subprocess.Popen(bwa_cmd,
                                      stdout=subprocess.PIPE,
                                      stderr=open(os.path.join(output_dir, "bwa_mem.err"), "a"))
    return bwa_mem_output


def run_magnet(args):
    input_tsv = args.classification
    genome_list = getattr(args, 'genome_list', None)
    no_cluster = getattr(args, 'no_cluster', False)
    prefix_contigs = not getattr(args, 'no_prefix_contigs', False)
    shared_reference = getattr(args, 'shared_reference', None)
    prebuilt_merged_fasta = None
    input_fastq = args.fastq
    input_fastq2 = args.fastq2
    mode = args.mode
    working_directory = args.output
    taxid_col_idx = args.taxid_idx
    abundance_col_idx = args.abundance_idx
    min_abundance = args.min_abundance
    min_mapq = args.min_mapq
    min_coverage_score = args.min_covscore
    min_breadth = getattr(args, 'min_breadth', 0.0)
    min_reads = getattr(args, 'min_reads', 0)
    rescue_min_reads = getattr(args, 'rescue_min_reads', 0)
    rescue_min_breadth = getattr(args, 'rescue_min_breadth', 0.0)
    rescue_min_ani = getattr(args, 'rescue_min_ani', 0.95)
    build_index_only = getattr(args, 'build_index_only', False)
    threads = args.threads
    valid_kingdom_str = args.kingdom

    valid_kingdom = set()
    for i in valid_kingdom_str.split(','):
        valid_kingdom.add(int(i))

    if args.include_mag:
        mag_flag = 'all'
    else:
        mag_flag = 'exclude'

    if build_index_only:
        # ---- Build step only ----
        # Construct the shared concatenated reference + its bwa index ONCE and
        # exit, without requiring or aligning any reads. The build is marker-
        # guarded, so if a previous run already produced the index this is a
        # no-op. Every later `--shared-reference` sample run reuses this index
        # (run_bwa_index_once sees a complete, up-to-date index and skips
        # `bwa index`), turning the merge + index into a one-time batch cost.
        from utils.shared_reference import build_shared_reference
        build_shared_reference(genome_list, shared_reference,
                               prefix_contigs=prefix_contigs,
                               index_for_bwa=(mode == 'illumina'),
                               threads=threads)
        merged = os.path.join(str(shared_reference), 'merged_reference.fasta')
        print(f"[magnet] shared reference ready: {merged}", file=sys.stderr)
        if mode == 'illumina':
            print(f"[magnet] bwa index present: {merged}.{{amb,ann,bwt,pac,sa}} "
                  f"(sentinel {merged}.bwa_index_complete)", file=sys.stderr)
        else:
            print("[magnet] ont mode: minimap2 indexes at map time; "
                  "merged reference is prebuilt and reused.", file=sys.stderr)
        return

    accession_flag = args.accession
    call_subspecies = args.subspecies

    sep = '\t'
    if input_tsv is not None and str(input_tsv)[-3:] == 'csv':
        sep = ','

    # NCBITaxa() triggers a (large) taxonomy DB download/load on first use.
    # In genome-list mode we never need it, so defer construction.
    ncbi_taxa_db = None

    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    if genome_list is not None:
        # --- Bring-your-own-genomes path ---
        # Build the same `downloaded_assemblies` schema and reference_genomes/
        # layout that the taxid->NCBI path produces, then fall through unchanged.
        if shared_reference is not None:
            # Shared-reference mode: lay out reference_genomes/, build + index the
            # merged reference ONCE in `shared_reference`, and have every sample
            # reuse it (reference_genomes symlinked in; prebuilt merged fasta +
            # bwa index reused). The build is idempotent (marker-guarded) and is
            # safe to call from every per-sample invocation.
            from utils.shared_reference import build_shared_reference, load_shared_reference
            build_shared_reference(genome_list, shared_reference,
                                   prefix_contigs=prefix_contigs,
                                   index_for_bwa=(mode == 'illumina'),
                                   threads=threads)
            downloaded_assemblies, prebuilt_merged_fasta = load_shared_reference(
                shared_reference, working_directory)
        else:
            from utils.genome_list_input import prepare_reference_genomes_from_list
            downloaded_assemblies = prepare_reference_genomes_from_list(
                genome_list, working_directory, prefix_contigs=prefix_contigs)
    else:
        # --- Original classifier-report path ---
        ncbi_taxa_db = NCBITaxa()
        if accession_flag:
            abundance_col_idx = None
            min_abundance = 0
            input_df, min_abundance = parsing_input_f(input_tsv, sep, taxid_col_idx, abundance_col_idx, min_abundance)
            valid_taxids = list(input_df['tax_id'].values)
        else:
            input_df, min_abundance = parsing_input_f(input_tsv, sep, taxid_col_idx, abundance_col_idx, min_abundance)
            # make valid_kingdom a variable?
            valid_taxids = filter_input_df(input_df, min_abundance, ncbi_taxa_db, valid_kingdom=valid_kingdom, ret_subspecies=call_subspecies)

        reference_metadata = prepare_reference_genomes(valid_taxids, working_directory, ncbi_taxa_db, accession_flag=accession_flag, mag_flag=mag_flag)
        downloaded_assemblies = reference_metadata[reference_metadata['Downloaded']]

    downloaded_assemblies = downloaded_assemblies.copy()
    reference_genome_path = os.path.join(working_directory, 'reference_genomes')

    if 'Cluster Representative' in downloaded_assemblies.columns:
        # Cluster assignment already present (e.g. loaded from a shared reference);
        # nothing to do.
        pass
    elif no_cluster or len(downloaded_assemblies) < 2:
        # Skip ANI clustering: treat every reference genome as its own cluster
        # representative. Also avoids running fastANI on a single genome (which
        # the matrix parser cannot handle).
        downloaded_assemblies['Cluster Representative'] = True
        downloaded_assemblies['Cluster Members'] = downloaded_assemblies['Assembly Accession ID'].astype(str)
    else:
        fastani_path = os.path.join(working_directory, 'fastANI')
        if not os.path.exists(fastani_path):
            os.mkdir(fastani_path)

        fastani_assemblies = downloaded_assemblies['Assembly Accession ID'].values
        assemblie_list_f = os.path.join(fastani_path, f"assemblie_list.txt")
        with open(assemblie_list_f, "w") as rl_f:
            for accession in fastani_assemblies:
                reference_genome = os.path.join(reference_genome_path, f'{accession}.fasta')
                rl_f.write(f"{reference_genome}\n")

        subprocess.run(['fastANI',
                        '--rl', assemblie_list_f,
                        '--ql', assemblie_list_f,
                        "--matrix",
                        '--threads', str(threads),
                        '-o', os.path.join(fastani_path, f'pairwise_ani')],
                       check=True,
                       stdout=open(os.path.join(fastani_path, "fastani.log"), "a"),
                       stderr=open(os.path.join(fastani_path, "fastani.err"), "a"))

        representative_genomes, member2representative = find_representative_genome(fastani_path, fastani_assemblies, downloaded_assemblies)

        representative_labels = []
        cluster_members = []
        for idx, row in downloaded_assemblies.iterrows():
            accession = row['Assembly Accession ID']
            if accession in representative_genomes.keys():
                representative_labels.append(True)
                cluster_members.append(','.join(representative_genomes[accession]))
            else:
                representative_labels.append(False)
                cluster_members.append(','.join(representative_genomes[member2representative[accession]]))

        downloaded_assemblies['Cluster Representative'] = representative_labels
        downloaded_assemblies['Cluster Members'] = cluster_members

    representative_df = downloaded_assemblies[downloaded_assemblies['Cluster Representative']]

    seq2assembly_dict = get_seq2assembly_dict(working_directory, representative_df)
    if prebuilt_merged_fasta is not None:
        # Shared mode: reuse the merged reference (and its bwa index) built once.
        reference_fasta = prebuilt_merged_fasta
    else:
        reference_fasta = merge_reference_fasta(list(representative_df['Assembly Accession ID']), working_directory)

    if mode == 'ont':
        aligner_output = run_minimap2(input_fastq, reference_fasta, 'merged', working_directory, threads=threads)
    if mode == 'illumina':
        # Use the guarded wrapper (defined below) which builds the bwa index only
        # if it is missing. With a shared reference the index is built once and
        # every sample reuses it. (The imported run_bwa always re-indexes.)
        aligner_output = run_bwa_index_once(input_fastq, input_fastq2, reference_fasta, 'merged', working_directory, threads=threads)

    sort_samfile('merged', aligner_output, working_directory, min_mapq=0, threads=threads)

    coverage_files = os.path.join(working_directory, "coverage_files")
    if not os.path.exists(coverage_files):
        os.mkdir(coverage_files)

    pool = Pool(processes=threads)
    pool.starmap(samtools_calculate_coverage, zip(repeat(working_directory), [True, False]))
    pool.close()
    pool.join()

    representative_df = alignment_summary(representative_df,
                                          working_directory,
                                          seq2assembly_dict,
                                          include_supp=True)

    representative_df = alignment_summary(representative_df,
                                          working_directory,
                                          seq2assembly_dict,
                                          include_supp=False)

    consensus_record_dict = samtools_merged_consensus(working_directory, threads)
    representative_df = ani_summary(representative_df, consensus_record_dict, working_directory, threads)
    representative_df = call_present_absent(representative_df, min_coverage_score)
    representative_df = apply_absolute_gate(representative_df,
                                            min_breadth=min_breadth,
                                            min_reads=min_reads)
    representative_df = apply_primary_evidence_rescue(representative_df,
                                                      min_reads=rescue_min_reads,
                                                      min_breadth=rescue_min_breadth,
                                                      min_ani=rescue_min_ani)

    representative_df.sort_values(['Primary Score'], ascending=False).to_csv(os.path.join(working_directory,
                                                                                          'cluster_representative.csv'),
                                                                             index=False)

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

def main():
    parser = argparse.ArgumentParser(description="Universal Taxonomic Classification Verifier.")
    parser.add_argument(
        "-c",
        "--classification",
        type=pathlib.Path,
        required=False,
        help="Path to the Taxonomic Classification Report. Accepting csv/tsv file format, other text formats are treated as tsv. Mutually exclusive with -g/--genome-list.",
    )
    parser.add_argument(
        "-i",
        "--fastq",
        type=pathlib.Path,
        required=False,
        help="Path to the first fastq file. Required unless --build-index-only.",
    )
    parser.add_argument(
        "-I",
        "--fastq2",
        type=pathlib.Path,
        required=False,
        help="Path to the second fastq file for paired-end reads.",
    )
    parser.add_argument(
        "-m",
        "--mode",
        type=str,
        required=False,
        choices=["ont", "illumina"],
        help="Modes for different sequencing platforms [ont, illumina]. Default:[ont]",
        default="ont",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=pathlib.Path,
        required=False,
        help="Path to the output directory. Required unless --build-index-only.",
    )
    parser.add_argument(
        "-t",
        "--taxid-idx",
        type=int,
        required=False,
        help="The column index (0-based) of the taxids. Default:[0]",
        default=0,
    )
    parser.add_argument(
        "-a",
        "--abundance-idx",
        type=int,
        required=False,
        help="The column index (0-based) of the abundance. Default:[None]",
    )
    parser.add_argument(
        "--min-abundance",
        type=float,
        required=False,
        help="Minimum abundance (0-1) for pre-filtering, exclude taxa below the threshold.",
        default=0,
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        required=False,
        help="Minimum MAPQ for primary alignments. Default:[20]",
        default=20,
    )
    parser.add_argument(
        "--min-covscore",
        type=float,
        required=False,
        help="Minimum Coverage Score for supplementary alignments. Default:[0.7]",
        default=0.7,
    )
    parser.add_argument(
        "--min-breadth",
        type=float,
        required=False,
        help="Absolute gate: minimum PRIMARY breadth of coverage (0-1) for a "
             "Present call. Applied after the coverage-score call and only "
             "downgrades Present->Absent. Drops genomes covered by a handful of "
             "reads whose ratio score saturates to ~1.0. Default:[0] (off).",
        default=0.0,
    )
    parser.add_argument(
        "--min-reads",
        type=int,
        required=False,
        help="Absolute gate: minimum number of PRIMARY mapped reads for a Present "
             "call. Applied after the coverage-score call and only downgrades "
             "Present->Absent. Default:[0] (off).",
        default=0,
    )
    parser.add_argument(
        "--rescue-min-reads",
        type=int,
        required=False,
        help="Primary-evidence rescue: upgrade Absent->Present when a genome has at "
             "least this many PRIMARY (uniquely-mapped) reads and ANI >= "
             "--rescue-min-ani. Recovers genuine low-abundance strains whose "
             "breadth/expected score fell below --min-covscore. Default:[0] (off).",
        default=0,
    )
    parser.add_argument(
        "--rescue-min-breadth",
        type=float,
        required=False,
        help="Primary-evidence rescue: alternative to --rescue-min-reads, using "
             "PRIMARY breadth (0-1). A genome is rescued if it meets either "
             "threshold and ANI >= --rescue-min-ani. Default:[0] (off).",
        default=0.0,
    )
    parser.add_argument(
        "--rescue-min-ani",
        type=float,
        required=False,
        help="Consensus ANI floor required for a primary-evidence rescue. "
             "Default:[0.95]",
        default=0.95,
    )
    parser.add_argument(
        "--threads",
        type=int,
        required=False,
        help="Number of threads for Multi-threading. Default:[1]",
        default=1,
    )
    parser.add_argument(
        "--kingdom",
        type=str,
        help="A comma separated list of taxids of valid kingdoms. Default:[2,4751,2157,10239]",
        default="2,4751,2157,10239",
    )
    parser.add_argument(
        "--include-mag",
        action="store_true",
        required=False,
        help="Include metagenomic assemble genomes. Default:[off]",
    )
    parser.set_defaults(include_mag=False)
    parser.add_argument(
        "--subspecies",
        action="store_true",
        required=False,
        help="Verify taxonomic classification at subspecies rank. Default:[off]",
    )
    parser.set_defaults(subspecies=False)
    parser.add_argument(
        "--accession",
        action="store_true",
        required=False,
        help="Take accession ids as taxids. Does not work with min-abundance. Default:[off]",
    )
    parser.set_defaults(accession=False)
    parser.add_argument(
        "-g",
        "--genome-list",
        type=pathlib.Path,
        required=False,
        help="Path to a TSV/CSV listing local strain genomes to use directly as "
             "references (column 'genome_path' required; optional genome_id, species, "
             "strain, taxid, assembly_level, abundance). Mutually exclusive with -c.",
    )
    parser.add_argument(
        "--no-cluster",
        action="store_true",
        required=False,
        help="Skip fastANI/ANI clustering; use every reference genome as-is. "
             "Default:[off] (always effectively on for a single genome).",
    )
    parser.set_defaults(no_cluster=False)
    parser.add_argument(
        "--no-prefix-contigs",
        action="store_true",
        required=False,
        help="With --genome-list, do NOT prefix contig headers with the genome id. "
             "Only safe if every contig name is already globally unique. Default:[off]",
    )
    parser.set_defaults(no_prefix_contigs=False)
    parser.add_argument(
        "--shared-reference",
        type=pathlib.Path,
        required=False,
        help="With --genome-list: directory holding ONE shared reference layout "
             "(reference_genomes/ + merged_reference.fasta + bwa index + metadata). "
             "Built once on first use and reused by every sample, so the merge and "
             "bwa index happen a single time across a batch. Requires --no-cluster.",
    )
    parser.add_argument(
        "--build-index-only",
        action="store_true",
        required=False,
        help="Build step only: construct the shared concatenated reference and its "
             "bwa index (for -m illumina) in --shared-reference, then exit WITHOUT "
             "mapping any reads. Run this once up front; subsequent per-sample runs "
             "with the same --shared-reference reuse the prebuilt index. If the index "
             "already exists from a previous run this is a no-op. Requires "
             "-g/--genome-list, --shared-reference and --no-cluster; -i/-o not needed.",
    )
    parser.set_defaults(build_index_only=False)

    args = parser.parse_args()

    if args.build_index_only:
        if args.genome_list is None or args.shared_reference is None:
            parser.error("--build-index-only requires -g/--genome-list and --shared-reference")
    else:
        if args.output is None:
            parser.error("-o/--output is required")
        if args.fastq is None:
            parser.error("-i/--fastq is required")

    if args.classification is None and args.genome_list is None:
        parser.error("one of -c/--classification or -g/--genome-list is required")
    if args.classification is not None and args.genome_list is not None:
        parser.error("-c/--classification and -g/--genome-list are mutually exclusive")
    if args.shared_reference is not None:
        if args.genome_list is None:
            parser.error("--shared-reference is only supported with -g/--genome-list")
        if not args.no_cluster:
            parser.error("--shared-reference currently requires --no-cluster "
                         "(the shared reference is the candidate set used as-is)")

    run_magnet(args)

if __name__ == "__main__":
    # Make run as a console entry point: `magnet ...`
    main()