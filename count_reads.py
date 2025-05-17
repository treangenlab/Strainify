import pysam
from collections import Counter
import time
import argparse
import os

def read_positions(file_path):
    """Reads positions from a tab-delimited file: chromosome<tab>position"""
    positions = []
    with open(file_path) as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            chrom, pos = line.strip().split()
            positions.append((chrom.split('.')[0], int(pos)))
    return positions

def update_fasta_header(fasta_path):
    base_name = os.path.basename(fasta_path)
    base_name_no_ext = os.path.splitext(base_name)[0].split('.')[0]
    new_header = f">{base_name_no_ext}\n"
    new_fasta_path = base_name_no_ext + "_renamed.fna"
    if os.path.exists(new_fasta_path):
        print(f"Renamed file already exists: {new_fasta_path}")
        return new_fasta_path

    with open(fasta_path, 'r') as f:
        lines = f.readlines()

    if not lines:
        raise ValueError("FASTA file is empty.")

    lines[0] = new_header  # Replace the first line

    with open(f"{base_name_no_ext}_renamed.fna", 'w') as f:
        f.writelines(lines)
    return new_fasta_path

def count_bases_and_write_tsv(bam_path, reference_path, positions_file, output_dir):
    positions = read_positions(positions_file)
    samfile = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(reference_path)
    bam_basename = os.path.basename(bam_path)
    output_name = os.path.splitext(bam_basename)[0].split('_sorted')[0] + '_read_counts.tsv'
    output_path = os.path.join(output_dir, output_name)
    #output_name = '_read_counts.tsv'
    #output_file = output_file + '_read_counts.tsv'

    with open(output_path, 'w') as out:
        out.write("chrom\tposition\tref\tbase\tcount\n")

        for chrom, pos in positions:
            ref_base = fasta.fetch(chrom, pos - 1, pos).upper()
            base_counts = Counter()
            insertion_counts = Counter()
            deletion_counts = Counter()

            for pileupcolumn in samfile.pileup(
                chrom, pos - 1, pos,
                truncate=True,
                ignore_overlaps=False,
                stepper="all",
                min_base_quality=0,
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

            out.write(f"{chrom}\t{pos}\t{ref_base}\t{ref_base}\t{base_counts.get(ref_base, 0)}\n")
            for base, count in base_counts.items():
                if base != ref_base:
                    out.write(f"{chrom}\t{pos}\t{ref_base}\t{base}\t{count}\n")
            for seq, count in insertion_counts.items():
                out.write(f"{chrom}\t{pos}\t{ref_base}\t{seq}\t{count}\n")
            for del_seq, count in deletion_counts.items():
                out.write(f"{chrom}\t{pos}\t{ref_base}\t{del_seq}\t{count}\n")

    samfile.close()
    fasta.close()

def main():
    parser = argparse.ArgumentParser(description="Count base and indel frequencies at specified positions in a BAM file.")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--ref", required=True, help="Reference FASTA file")
    parser.add_argument("--positions", required=True, help="Positions file (tab-delimited: chrom<tab>pos)")
    parser.add_argument("--output", required=True, help="Output directory")
    args = parser.parse_args()

    #new_path = update_fasta_header(args.ref)
    count_bases_and_write_tsv(args.bam, args.ref, args.positions, args.output)



if __name__ == "__main__":
    main()
