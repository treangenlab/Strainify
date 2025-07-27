import argparse
import os.path
from Bio import AlignIO
import pandas as pd
from scipy.stats import chi2
import numpy as np
from collections import defaultdict


def create_overlapping_windows(intervals, window_size, step_size):
    windows = []
    buffer = []
    current_start = None

    for start, end in intervals:
        pos = start
        while pos <= end:
            if current_start is None:
                current_start = pos

            chunk_size = min(end - pos + 1, window_size - len(buffer))
            buffer.extend(range(pos, pos + chunk_size))
            pos += chunk_size

            if len(buffer) == window_size:
                windows.append(buffer.copy())
                buffer = buffer[step_size:]
                current_start = buffer[0] if buffer else None

    if buffer:
        windows.append(buffer.copy())

    return windows

def vcf_to_dataframe(vcf_file):
    with open(vcf_file, 'r') as vcf:
        data = []
        header = None
        for line in vcf:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.strip().lstrip('#').split("\t")
                continue
            row = line.strip().split("\t")
            data.append(row)

        df = pd.DataFrame(data, columns=header)
        df['POS'] = df.apply(lambda x: int(x['POS']) + 1 if 'SVTYPE=DEL' in str(x['INFO']) else x['POS'], axis=1)
        df['ALT'] = df.apply(lambda x: f"-{x['REF'][1:]}" if 'SVTYPE=DEL' in str(x['INFO']) else x['ALT'], axis=1)
        df['ALT'] = df.apply(lambda x: f"+{x['ALT'][1:]}" if 'SVTYPE=INS' in str(x['INFO']) else x['ALT'], axis=1)
        df['REF'] = df.apply(lambda x: x['REF'][1] if 'SVTYPE=DEL' in str(x['INFO']) else x['REF'], axis=1)
        df = df[df['ALT'] != '<INV>']
        df['CHROM'] = df['CHROM'].apply(lambda x: x.split("#")[-1])
        df.reset_index(drop=True, inplace=True)
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Detect and remove variant-dense regions from a MAF and VCF input.")
    parser.add_argument('--maf', required=True, help="MAF file from multiple alignment")
    parser.add_argument('--vcf', required=True, help="VCF file with variants")
    parser.add_argument('--output_dir', required=True, help="output directory name")
    parser.add_argument('--window_size', default="500", help="Window size (positive integer) or 'average_LCB_length'")
    parser.add_argument('--window_overlap', type=float, default=0.0, help="Window overlap as a float in [0, 1)")

    args = parser.parse_args()

    maf_file = args.maf
    vcf_file = args.vcf
    output_dir = args.output_dir

    if args.window_size == "average_LCB_length":
        window_size = None
    elif args.window_size.isdigit() and int(args.window_size) > 0:
        window_size = int(args.window_size)
    else:
        raise ValueError("Invalid window_size. Use a positive integer or 'average_LCB_length'.")

    if not (0 <= args.window_overlap < 1):
        raise ValueError("window_overlap must be between 0 (inclusive) and 1 (exclusive).")

    window_overlap = args.window_overlap

    # Step 1: Extract core genome intervals per CHROM
    core_intervals_per_chrom = defaultdict(list)

    for multiple_alignment in AlignIO.parse(maf_file, "maf"):
        # Use the first record to identify the reference strain
        reference_strain = multiple_alignment[0].id.split('#')[0]

        for record in multiple_alignment:
            strain, chrom = record.id.split('#')
            if strain == reference_strain:
                start = record.annotations["start"]
                size = record.annotations["size"]
                end = start + size
                core_intervals_per_chrom[chrom].append((start, end))

    # Optional: Sort intervals for each CHROM
    for chrom in core_intervals_per_chrom:
        core_intervals_per_chrom[chrom].sort()# Flatten all intervals across all contigs
    # Flatten all intervals across all contigs
    all_intervals = [end - start for intervals in core_intervals_per_chrom.values() for start, end in intervals]

    if all_intervals:
        avg_interval_size = int(np.mean(all_intervals))
        print(f"Average core genome interval size across all contigs: {avg_interval_size}")
    else:
        print("No core genome intervals found across any contig.")

    if window_size is None:
        window_size = avg_interval_size

    step_size = int(window_size * (1 - window_overlap))
    print(f"Using window size: {window_size}, step size: {step_size}")

    # Create overlapping windows per CHROM
    windows_per_chrom = {
        chrom: create_overlapping_windows(intervals, window_size, step_size)
        for chrom, intervals in core_intervals_per_chrom.items()
    }

    variant_df = vcf_to_dataframe(vcf_file)
    pd.set_option('display.max_columns', None)
    variant_df['POS'] = variant_df['POS'].astype(int)

    variant_df_no_dup = variant_df[~variant_df.duplicated(subset=['CHROM', 'POS'], keep=False)]

    variant_df = variant_df[~variant_df.duplicated(subset=['CHROM', 'POS'], keep="first")]
    variant_df.reset_index(drop=True, inplace=True)
    variant_df = variant_df.iloc[:, list(range(0, 5)) + list(range(9, variant_df.shape[1]))]
    variant_df.iloc[:, 5:] = variant_df.iloc[:, 5:].apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

    total_variants = len(variant_df)
    genome_length = sum(all_intervals)
    expected_density = total_variants / genome_length

    significant_windows = []

    for chrom, windows in windows_per_chrom.items():
        # Subset variant_df for this chromosome
        chrom_df = variant_df[variant_df['CHROM'] == chrom]

        for positions in windows:
            window_len = len(positions)
            start, end = positions[0], positions[-1]
            in_window = chrom_df[(chrom_df['POS'] >= start) & (chrom_df['POS'] <= end)]
            observed = len(in_window)
            expected = expected_density * window_len

            if expected == 0:
                continue

            chi2_stat = ((observed - expected) ** 2) / expected
            p_value = chi2.sf(chi2_stat, df=1)
            z = (observed - expected) / np.sqrt(expected)

            if p_value < 0.05 and observed > expected:
                significant_windows.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'z': z,
                    'pval': round(p_value, 6)
                })

    # Write results
    filtered_df = pd.DataFrame(significant_windows)
    filtered_df.to_csv(os.path.join(output_dir, "significantly_enriched_windows.tsv"), sep='\t', index=False)

    positions_to_exclude = set()

    for win in significant_windows:
        positions_to_exclude.update(
            variant_df_no_dup[
                (variant_df_no_dup['CHROM'] == win['chrom']) &
                (variant_df_no_dup['POS'] >= win['start']) &
                (variant_df_no_dup['POS'] <= win['end'])
                ].index
        )

    # Drop variants in significantly enriched windows
    df_AF = variant_df_no_dup.drop(index=positions_to_exclude)

    # Keep only biallelic variants (ALT with a single allele)
    df_AF = df_AF[df_AF['ALT'].apply(lambda x: len(set(x.split(','))) == 1)]

    # Reset index for cleanliness
    df_AF.reset_index(drop=True, inplace=True)

    # Save filtered variant matrix
    df_AF.to_csv(os.path.join(output_dir, "filtered_variant_matrix.csv"), index=False)

    # Save site positions (CHROM and POS)
    df_sites = df_AF[['CHROM', 'POS']]
    df_sites.to_csv(os.path.join(output_dir, 'sites.txt'), sep='\t', index=False, header=None)

    # Output log
    print("Finished filtering variants. Output files:")
    print("- significantly_enriched_windows.tsv")
    print("- filtered_variant_matrix.csv")
    print("- sites.txt")

