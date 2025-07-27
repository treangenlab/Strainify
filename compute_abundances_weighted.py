import sys
import os
import glob
import cvxpy
import numpy as np
import pandas as pd


def compute_entropy_weights(A):
    # Compute alt allele frequency f for each site (row)
    f = np.mean(A, axis=1)  # shape: (n_variants,)

    # Avoid log(0) by clipping f away from 0 and 1
    epsilon = 1e-10
    f_clipped = np.clip(f, epsilon, 1 - epsilon)

    # Binary entropy formula
    entropy = -f_clipped * np.log2(f_clipped) - (1 - f_clipped) * np.log2(1 - f_clipped)

    # Optional: normalize to [0, 1] for numerical stability
    # entropy /= np.max(entropy)

    return entropy  # shape: (n_variants,)

# === Command-line arguments ===
if len(sys.argv) != 3:
    print("Usage: python compute_abundances_weighted.py <read_counts_dir> <filtered_variant_matrix.csv>")
    sys.exit(1)

read_counts_dir = sys.argv[1]
variant_matrix_file = sys.argv[2]

# === Load variant matrix ===
df_AF = pd.read_csv(variant_matrix_file)
#df_AF = df_AF[~df_AF['INFO'].str.contains('SVTYPE=DEL|SVTYPE=INS', regex=True)]
df_allele_freq = df_AF.drop(columns=['ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
print(df_allele_freq)
strain_names = df_allele_freq.columns[2:]
strain_count = len(strain_names)

# === Store all results ===
abundance_df = pd.DataFrame({'strain name': strain_names})
read_count_files = sorted(glob.glob(os.path.join(read_counts_dir, "*_read_counts.tsv")))

for read_counts_file in read_count_files:
    print(f"Processing {read_counts_file}...")
    df_read_counts = pd.read_csv(read_counts_file, sep='\t')
    df_read_counts = df_read_counts.drop_duplicates(keep=False).reset_index(drop=True)
    df_AF = df_AF.reset_index(drop=True)

    df_read_counts['position'] = df_read_counts['position'].astype(int)
    df_AF['POS'] = df_AF['POS'].astype(int)

    var_reads = pd.merge(df_read_counts, df_AF,
                         left_on=['position', 'ref', 'base', 'chrom'],
                         right_on=['POS', 'REF', 'ALT', 'CHROM'],
                         how='inner')

    ref_reads = pd.merge(df_read_counts, df_AF,
                         left_on=['position', 'ref', 'base', 'chrom'],
                         right_on=['POS', 'REF', 'REF', 'CHROM'],
                         how='inner')

    merged_ref_var = pd.merge(ref_reads.iloc[:, :5], var_reads.iloc[:, :5], on=['chrom','position'], how='inner')
    print(merged_ref_var)
    if merged_ref_var.empty:
        print(f"Warning: No overlap found in {read_counts_file}")
        continue

    merged_ref_var['depth'] = merged_ref_var['count_x'] + merged_ref_var['count_y']
    merged_ref_var['vaf_x'] = merged_ref_var['count_x'] / merged_ref_var['depth']
    merged_ref_var['vaf_y'] = merged_ref_var['count_y'] / merged_ref_var['depth']

    # Compute total read depth per position
    merged_ref_var['depth'] = merged_ref_var['count_x'] + merged_ref_var['count_y']

    # Calculate mean and std of depth
    mean_depth = merged_ref_var['depth'].mean()
    std_depth = merged_ref_var['depth'].std()
    #threshold = mean_depth - 2 * std_depth
    threshold = 0

    print(f"Mean depth: {mean_depth:.2f}, Std depth: {std_depth:.2f}")
    print(f"Filtering positions with depth < {threshold:.2f}")

    # Filter out positions with significantly lower than average depth
    merged_ref_var = merged_ref_var[merged_ref_var['depth'] >= threshold]

    # Filter out plasmids in reference
    #merged_ref_var = merged_ref_var[merged_ref_var['chrom']=='NC_007946_1']
    #merged_ref_var = merged_ref_var[merged_ref_var['chrom'] == 'NC_009801_1']

    if merged_ref_var.empty:
        print(f"All positions filtered out for {read_counts_file} due to low depth")
        continue


    observed_ref_vector = merged_ref_var['vaf_x']
    observed_var_vector = merged_ref_var['vaf_y']

    df_allele_freq['POS'] = df_allele_freq['POS'].astype(int)
    print(df_allele_freq)
    #filtered_freqs = df_allele_freq[df_allele_freq['POS'].isin(merged_ref_var['position'])]
    # mask = df_allele_freq.set_index(['CHROM', 'POS']).index.isin(
    #     merged_ref_var.set_index(['chrom', 'position']).index
    # )
    # filtered_freqs = df_allele_freq[mask]

    filtered_freqs = pd.merge(
        df_allele_freq,
        merged_ref_var[['chrom', 'position']],
        left_on=['CHROM', 'POS'],
        right_on=['chrom', 'position'],
        how='inner'
    )
    print(filtered_freqs)
    filtered_freqs.reset_index(drop=True, inplace=True)
    allele_freq_matrix = filtered_freqs.drop(columns=['CHROM', 'POS', 'chrom', 'position'])
    #print(allele_freq_matrix)
    allele_freq_matrix = allele_freq_matrix.to_numpy()


    # === CVXPY optimization ===
    A = allele_freq_matrix
    b = observed_var_vector
    c = observed_ref_vector

    x = cvxpy.Variable(strain_count, nonneg=True)
    ps = cvxpy.matmul(A, x)
    ps += np.nextafter(0, 1)  # avoid log(0)

    # Apply entropy weights to log-likelihood
    w = compute_entropy_weights(A)
    print(w)
    weighted_loglik = cvxpy.multiply(w, cvxpy.multiply(b, cvxpy.log(ps)) + cvxpy.multiply(c, cvxpy.log(1 - ps)))
    nllh = -cvxpy.sum(weighted_loglik)


    #nllh = -1 * cvxpy.sum(cvxpy.multiply(b, cvxpy.log(ps)) + cvxpy.multiply(c, cvxpy.log(1 - ps)))
    objective = cvxpy.Minimize(nllh)
    constraints = [0 <= x, x <= 1, sum(x) == 1]
    prob = cvxpy.Problem(objective, constraints)

    try:
        result = prob.solve(verbose=True, solver='ECOS', max_iters=50000, warm_start=True)
    except cvxpy.error.SolverError:
        print(f"Solver failed for {read_counts_file}")
        continue

    if result is not None:
        x.value = x.value.astype(float)
        rounded_arr = np.round(np.array(x.value) * 100, 4)

        sample_id = os.path.basename(read_counts_file).replace("_read_counts.tsv", "")
        abundance_df[sample_id] = rounded_arr
    else:
        print(f"Optimization failed for {read_counts_file}")

# === Save final merged output ===
abundance_df['strain name'] = abundance_df['strain name'].str.replace(r'\.fna$|\.fasta$', '', regex=True)
abundance_path = os.path.join(os.path.dirname(variant_matrix_file), "abundance_estimates_combined.csv")
abundance_df.to_csv(abundance_path, index=False)

print("\nAll samples processed. Combined output saved to 'abundance_estimates_combined.csv'.")
