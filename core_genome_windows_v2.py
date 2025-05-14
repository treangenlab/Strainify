from Bio import AlignIO
import pandas as pd
from scipy.stats import chi2
import numpy as np
from collections import defaultdict


def create_overlapping_windows(intervals, window_size, step_size):
    windows = []
    buffer = []  # To hold contiguous bases across intervals
    current_start = None

    for start, end in intervals:
        pos = start
        while pos <= end:
            if current_start is None:
                current_start = pos

            # Fill buffer with current chunk
            chunk_size = min(end - pos + 1, window_size - len(buffer))
            buffer.extend(range(pos, pos + chunk_size))
            pos += chunk_size

            # When we reach full window size, emit a window
            if len(buffer) == window_size:
                windows.append(buffer.copy())
                # Slide forward by step_size
                buffer = buffer[step_size:]
                current_start = buffer[0] if buffer else None

    # Handle remaining buffer
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
        df['CHROM'] = df['CHROM'].apply(lambda x: x.split("#")[0])
        df.reset_index(drop=True, inplace=True)
    return df


# Input values
maf_file = "W0008_as_ref.maf"
#reference_strain = "W0008"

# Step 1: Extract core genome intervals from reference in MAF
core_intervals = []
for multiple_alignment in AlignIO.parse(maf_file, "maf"):
    reference_strain = multiple_alignment[0].id.split('#')[0]
    for record in multiple_alignment:
        if reference_strain in record.id:
            start = record.annotations["start"]
            size = record.annotations["size"]
            core_intervals.append((start, start + size))

core_intervals.sort()
print(core_intervals)

# Compute average interval size
interval_lengths = [end - start for start, end in core_intervals]
avg_interval_size = int(np.mean(interval_lengths))
print(f"Average core genome interval size: {avg_interval_size}")

# Define overlapping window size and step size
#window_size = avg_interval_size
window_size = 500
step_size = int(window_size*1.00)  # percentage overlap = 1 - step_size

# Create overlapping windows
windows = create_overlapping_windows(core_intervals, window_size, step_size)

# Load variant matrix
variant_df = vcf_to_dataframe("./W0008_as_ref.vcf")
pd.set_option('display.max_columns', None)
print(variant_df.head(5))
variant_df['POS'] = variant_df['POS'].astype(int)
variant_df_no_dup = variant_df[~variant_df['POS'].duplicated(keep=False)]
variant_df = variant_df[~variant_df['POS'].duplicated(keep="first")]
variant_df.reset_index(drop=True, inplace=True)
variant_df = variant_df.iloc[:, list(range(0, 5)) + list(range(9, variant_df.shape[1]))]
variant_df.iloc[:, 5:] = variant_df.iloc[:, 5:].apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)
print(variant_df.head(5))

# Compute expected variant density
total_variants = len(variant_df)
genome_length = sum(end - start + 1 for start, end in core_intervals)
expected_density = total_variants / genome_length

# Identify enriched windows
significant_windows = []
for positions in windows:
    window_len = len(positions)
    start, end = positions[0], positions[-1]
    in_window = variant_df[(variant_df['POS'] >= start) & (variant_df['POS'] <= end)]
    observed = len(in_window)
    expected = expected_density * window_len
    if expected == 0:
        continue
    chi2_stat = ((observed - expected) ** 2) / expected
    p_value = chi2.sf(chi2_stat, df=1)
    z = (observed - expected) / np.sqrt(expected)

    if p_value < 0.05 and observed > expected:
        significant_windows.append({
            'start': start,
            'end': end,
            'z': z,
            'pval': round(p_value, 6)
        })

# Save enriched windows
filtered_df = pd.DataFrame(significant_windows)
filtered_df.to_csv("significantly_enriched_windows.tsv", sep='\t', index=False)

# Filter out variants in enriched windows
positions_to_exclude = set()
for win in significant_windows:
    positions_to_exclude.update(
        variant_df_no_dup[(variant_df_no_dup['POS'] >= win['start']) & (variant_df_no_dup['POS'] <= win['end'])].index
    )

df_AF = variant_df_no_dup.drop(index=positions_to_exclude)
df_AF = df_AF[df_AF['ALT'].apply(lambda x: len(set(x.split(','))) == 1)]
df_AF.reset_index(drop=True, inplace=True)
print(df_AF)

df_AF.to_csv("filtered_variant_matrix.csv", index=False)
df_allele_freq = df_AF.drop(columns=['CHROM','ID','REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
strain_count = len(df_allele_freq.columns[1:])
strains = df_allele_freq.columns

# Create sites to count allele frequencies in aligned reads (BAM file)
df_sites = df_AF.iloc[:, list(range(0, 2))]

df_sites.to_csv('sites.txt', sep='\t', index=False, header=None)
'''
import os
import cvxpy
import numpy as np

# Folder containing the files
folder_path = './no_IS'

# Initialize an empty DataFrame to hold the combined results
combined_results = pd.DataFrame()

# Loop through each file in the folder
for file_name in os.listdir(folder_path):
    if file_name.endswith('.tsv'):  # Assuming the files are .tsv
        file_path = os.path.join(folder_path, file_name)

        # Read the data (Replace with your actual df_read_counts and df_AF loading logic)
        print(f"Processing file: {file_name}")

        df_read_counts = pd.read_csv(file_path, sep='\t')
        print(df_read_counts)

        # Reset the indices
        df_read_counts = df_read_counts.drop_duplicates(keep=False)
        df_read_counts = df_read_counts.reset_index(drop=True)
        df_AF = df_AF.reset_index(drop=True)
        #df_AF.to_csv('df_AF.csv')

        print('Filtering read counts')
        df_read_counts['position'] = df_read_counts['position'].astype(int)
        df_AF['POS'] = df_AF['POS'].astype(int)

        # Merging read counts
        var_reads = pd.merge(df_read_counts, df_AF,
                             left_on=['position', 'ref', 'base'],
                             right_on=['POS', 'REF', 'ALT'],
                             how='inner')

        print('Reads filtered')
        print(var_reads)
        var_reads.to_csv('filtered_reads_w_ref.csv')

        # Reference reads merge
        ref_reads = pd.merge(df_read_counts, df_AF,
                             left_on=['position', 'ref', 'base'],
                             right_on=['POS', 'REF', 'REF'],
                             how='inner')
        print(ref_reads)

        # Merge ref and var
        merged_ref_var = pd.merge(ref_reads.iloc[:, :5], var_reads.iloc[:, :5], on='position', how='inner')
        print(merged_ref_var)
        # merged_ref_var.to_csv('merged_ref_var.csv')

        merged_ref_var['depth'] = merged_ref_var['count_x'] + merged_ref_var['count_y']
        merged_ref_var['vaf_x'] = merged_ref_var['count_x'] / merged_ref_var['depth']
        merged_ref_var['vaf_y'] = merged_ref_var['count_y'] / merged_ref_var['depth']
        merged_ref_var['position'] = merged_ref_var['position'].astype(int)
        print(merged_ref_var[merged_ref_var['position'].duplicated(keep=False)])

        # Extract vectors
        observed_ref_vector = merged_ref_var['vaf_x']
        observed_var_vector = merged_ref_var['vaf_y']
        print('This is observed_var_vector')
        print(observed_var_vector)
        print(len(observed_var_vector))
        print('This is observed_ref_vector')
        print(observed_ref_vector)
        print(len(observed_ref_vector))
        # ref_reads.to_csv('ref_reads.csv')

        # Filter frequencies
        df_allele_freq['POS'] = df_allele_freq['POS'].astype(int)
        merged_ref_var['position'] = merged_ref_var['position'].astype(int)
        filtered_freqs = df_allele_freq[df_allele_freq['POS'].isin(merged_ref_var['position'])]

        filtered_freqs.reset_index(drop=True, inplace=True)

        allele_freq_matrix = filtered_freqs.iloc[:, 1:].to_numpy()

        # CVXPY optimization
        A = allele_freq_matrix
        b = observed_var_vector
        c = observed_ref_vector

        x = cvxpy.Variable(strain_count, nonneg=True)
        ps = cvxpy.matmul(A, x)
        print(ps)

        eps = np.nextafter(0, 1)  # smallest positive float
        ps += eps  # avoid Log(0)
        nllh = (-1) * cvxpy.sum((cvxpy.multiply(b, cvxpy.log(ps))) + cvxpy.multiply(c, cvxpy.log(1 - ps)))

        objective = cvxpy.Minimize(nllh)
        constraints = [0 <= x, x <= 1, sum(x) <= 1]
        prob = cvxpy.Problem(objective, constraints)
        try:
            result = prob.solve(verbose=True, solver='ECOS', max_iters=50000, warm_start=True)
            if result is None or result == 'failed':  # Check if the solver failed
                print(f"Solver failed for {file_name}. Skipping to next file.")
                continue  # Skip to the next file
        except cvxpy.error.SolverError:
            print(f"Solver encountered an error for {file_name}. Skipping to next file.")
            continue  # Skip to the next file

        try:
            result = prob.solve(verbose=True, warm_start=True)
        except cvxpy.error.SolverError:
            result = 'failed'

        print(prob)
        print(result)
        x.value = x.value.astype(float)
        rounded_arr = np.array([np.round(a * 100, 4) for a in x.value])

        # Assuming `strains` is available as a global list or define it as needed
        strains = strains[1:]
        strain_names = pd.DataFrame(strains, columns=['strain name'])
        df_results = pd.DataFrame(rounded_arr)
        df_results = pd.concat([strain_names, df_results], axis=1)
        df_results.columns = ['strain name', f'{file_name.replace("_observed_AF.tsv", "")}']
        df_results['strain name'] = df_results['strain name'].str.replace('.gb.fna', '', regex=False)

        # Combine this file's results with the combined results DataFrame
        if combined_results.empty:
            combined_results = df_results
        else:
            # Merge on 'strain name'
            combined_results = pd.merge(combined_results, df_results, on='strain name', how='outer')

# Save the final combined results to a CSV
combined_results.to_csv('W0008_w_recomb_filter_500bp.csv', index=False)
print(combined_results)

'''