import pysam
from collections import Counter
import time

start = time.time()
def read_positions(file_path):
    """Reads positions from a tab-delimited file: chromosome<tab>position"""
    positions = []
    with open(file_path) as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            chrom, pos = line.strip().split()
            positions.append((chrom, int(pos)))
    return positions
def count_bases_and_write_tsv(bam_path, ref_fasta_path, positions_file, output_file):
    """
    Count bases and indels at positions and write results to a TSV file.

    Args:
        bam_path (str): Path to the BAM file.
        ref_fasta_path (str): Path to the reference FASTA file.
        positions_file (str): Path to the file listing positions (chrom\tpos).
        output_file (str): Path to write the TSV output.
    """
    positions = read_positions(positions_file)
    samfile = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(ref_fasta_path)

    with open(output_file, 'w') as out:
        #out.write("chrom\tposition\tref\tbase\tvaf\tdepth\tcount\n")
        out.write("chrom\tposition\tref\tbase\tcount\n")


        for chrom, pos in positions:
            ref_base = fasta.fetch(chrom, pos - 1, pos).upper()
            base_counts = Counter()
            insertion_counts = Counter()
            deletion_counts = Counter()

            for pileupcolumn in samfile.pileup(chrom, pos - 1, pos, truncate=True, ignore_overlaps=False,stepper="all", min_base_quality=0, min_mapping_quality=0,flag_filter=0,max_depth=10000000,ignore_orphans=True):
                #print(pos, pileupcolumn.n)
                #if pileupcolumn.reference_pos != pos - 1:
                #    continue

                for pileupread in pileupcolumn.pileups:
                    aln = pileupread.alignment
                      # skip 'N' bases
                    if pileupread.query_position is not None:
                        base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                        if base == 'N':
                            continue
                        base_counts[base] += 1

                    if pileupread.indel > 0:
                        # insertion
                        seq = '+' + pileupread.alignment.query_sequence[
                            pileupread.query_position+1 : pileupread.query_position + pileupread.indel+1]
                        insertion_counts[seq] += 1
                        #ref_base = fasta.fetch(chrom, pos - 1, pos).upper()

                    if pileupread.is_del:
                        # Deletions — NEW: Parse from CIGAR directly
                        #if not aln.cigartuples:
                         #   continue

                        ref_pointer = aln.reference_start

                        for op, length in aln.cigartuples:
                            if op in (0, 7, 8):  # match, equal, mismatch
                                ref_pointer += length
                            if ref_pointer == pos-1 and op == 2:
                                del_len = length
                                #del_ref = fasta.fetch(chrom, pos - 2, pos - 1).upper()
                                del_seq = '-' + fasta.fetch(chrom, pos - 1, pos -1 + del_len).upper()
                                deletion_counts[del_seq] += 1
                                #del_tuple = (del_seq, del_ref)
                            #if ref_pointer < pos - 1 < ref_pointer + length and op == 2:
                            #    del_len = ref_pointer + length - (pos - 1)
                            #    #del_start='.'*(pos-1-ref_pointer)
                            #    #del_ref = fasta.fetch(chrom, pos - 1, pos).upper()
                            #    #del_ref = '-'
                            #    del_seq = '-' + fasta.fetch(chrom, pos - 1, pos -1 + del_len).upper()
                                #del_tuple = (del_seq, del_ref)


                        #del_seq = fasta.fetch(chrom, pos - 2, pos + del_len - 1).upper()

                        #deletion_counts[del_seq] += 1
                        #ref_base = fasta.fetch(chrom, pos - 2, pos - 1).upper()


            # Reference base
            ref_count = base_counts.get(ref_base, 0)
            #vaf = ref_count / depth
            #out.write(f"{chrom}\t{pos}\t{ref_base}\t{ref_base}\t{vaf:.4f}\t{depth}\t{ref_count}\n")
            out.write(f"{chrom}\t{pos}\t{ref_base}\t{ref_base}\t{ref_count}\n")

            # Alternate bases
            for base, count in base_counts.items():
                if base != ref_base:
                    #vaf = count / depth
                    #out.write(f"{chrom}\t{pos}\t{ref_base}\t{base}\t{vaf:.4f}\t{depth}\t{count}\n")
                    out.write(f"{chrom}\t{pos}\t{ref_base}\t{base}\t{count}\n")

            # Insertions
            for seq, count in insertion_counts.items():
                #vaf = count / depth
                #out.write(f"{chrom}\t{pos}\t{ref_base}\t{seq}\t{vaf:.4f}\t{depth}\t{count}\n")
                out.write(f"{chrom}\t{pos}\t{ref_base}\t{seq}\t{count}\n")

            # Deletions (with sequences)
            for del_seq, count in deletion_counts.items():
                #vaf = count / depth
                #out.write(f"{chrom}\t{pos}\t{ref_base}\t{del_seq}\t{vaf:.4f}\t{depth}\t{count}\n")
                out.write(f"{chrom}\t{pos}\t{ref_base}\t{del_seq}\t{count}\n")

    samfile.close()
    fasta.close()


count_bases_and_write_tsv(
    bam_path="/Users/Rossie/Strain_phaser/W0008_q0_sorted.bam",
    ref_fasta_path="/Users/Rossie/Strain_phaser/W0008.gb.fna",
    positions_file="./W0008_sites_merged.txt",
    output_file="./W0008_q0.tsv"
)
end = time.time()

print ("Execution time: {:.2f} minutes".format((end - start)/60))