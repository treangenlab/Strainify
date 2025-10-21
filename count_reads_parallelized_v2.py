import pysam
import pandas as pd
from collections import Counter
import argparse
import os
import multiprocessing as mp
from math import ceil


def read_positions(file_path):
    positions = []
    with open(file_path) as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            chrom, pos = line.strip().split()
            positions.append((chrom.split('.')[0], int(pos)))
    return positions


def process_chunk(chunk, bam_path, reference_path):
    """Process a subset of positions and return a DataFrame."""
    samfile = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(reference_path)

    records = []

    for chrom, pos in chunk:
        ref_base = fasta.fetch(chrom, pos - 1, pos).upper()
        base_counts = Counter()
        insertion_counts = Counter()
        deletion_counts = Counter()

        for pileupcolumn in samfile.pileup(
            chrom, pos - 1, pos,
            truncate=True,
            ignore_overlaps=True,
            stepper="all",
            min_base_quality=20,
            min_mapping_quality=0,
            flag_filter=0,
            max_depth=10000000,
            ignore_orphans=True
        ):
            for pileupread in pileupcolumn.pileups:
                aln = pileupread.alignment
                if pileupread.query_position is not None:
                    base = aln.query_sequence[pileupread.query_position].upper()
                    if base == 'N':
                        continue
                    base_counts[base] += 1

                if pileupread.indel > 0:
                    seq = '+' + aln.query_sequence[
                        pileupread.query_position + 1 : pileupread.query_position + pileupread.indel + 1
                    ]
                    insertion_counts[seq] += 1

                if pileupread.is_del:
                    ref_pointer = aln.reference_start
                    for op, length in aln.cigartuples:
                        if op in (0, 7, 8):
                            ref_pointer += length
                        if ref_pointer == pos - 1 and op == 2:
                            del_seq = '-' + fasta.fetch(chrom, pos - 1, pos - 1 + length).upper()
                            deletion_counts[del_seq] += 1

        # Add all base counts
        records.append((chrom, pos, ref_base, ref_base, base_counts.get(ref_base, 0)))
        for base, count in base_counts.items():
            if base != ref_base:
                records.append((chrom, pos, ref_base, base, count))
        for seq, count in insertion_counts.items():
            records.append((chrom, pos, ref_base, seq, count))
        for del_seq, count in deletion_counts.items():
            records.append((chrom, pos, ref_base, del_seq, count))

    samfile.close()
    fasta.close()

    df = pd.DataFrame(records, columns=["chrom", "position", "ref", "base", "count"])
    return df


def count_bases_and_write_tsv(bam_path, reference_path, positions_file, output_dir, threads=4):
    """Parallelized counting of bases and indels"""
    os.makedirs(output_dir, exist_ok=True)
    positions = read_positions(positions_file)
    bam_basename = os.path.basename(bam_path)
    output_name = os.path.splitext(bam_basename)[0].split('_sorted')[0] + '_read_counts.tsv'
    output_path = os.path.join(output_dir, output_name)

    # Split into chunks
    chunk_size = ceil(len(positions) / threads)
    chunks = [positions[i:i + chunk_size] for i in range(0, len(positions), chunk_size)]

    # Process in parallel, each returns a DataFrame
    with mp.Pool(threads) as pool:
        dfs = pool.starmap(process_chunk, [(chunk, bam_path, reference_path) for chunk in chunks])

    # Combine all DataFrames
    combined_df = pd.concat(dfs, ignore_index=True)

    # Sort and write to disk
    combined_df.sort_values(["chrom", "position"], inplace=True)
    combined_df.to_csv(output_path, sep="\t", index=False)

    print(f"Finished writing read count output: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Parallel base and indel counting")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--ref", required=True, help="Reference FASTA file")
    parser.add_argument("--positions", required=True, help="Positions file (tab-delimited: chrom<tab>pos)")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--threads", type=int, default=8, help="Number of parallel processes (default: 8)")
    args = parser.parse_args()

    count_bases_and_write_tsv(args.bam, args.ref, args.positions, args.output, args.threads)


if __name__ == "__main__":
    main()
