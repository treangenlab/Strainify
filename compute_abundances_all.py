import sys
import os
import glob
import cvxpy
import numpy as np
import pandas as pd

# === Command-line arguments ===
if len(sys.argv) != 3:
    print("Usage: python compute_abundances_multi.py <read_counts_dir> <filtered_variant_matrix.csv>")
    sys.exit(1)

read_counts_dir = sys.argv[1]
variant_matrix_file = sys.argv[2]

# === Load variant matrix ===
df_AF = pd.read_csv(variant_matrix_file)
df_allele_freq = df_AF.drop(columns=['CHROM', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
strain_names = df_allele_freq.columns[1:]
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
                         left_on=['position', 'ref', 'base'],
                         right_on=['POS', 'REF', 'ALT'],
                         how='inner')

    ref_reads = pd.merge(df_read_counts, df_AF,
                         left_on=['position', 'ref', 'base'],
                         right_on=['POS', 'REF', 'REF'],
                         how='inner')

    merged_ref_var = pd.merge(ref_reads.iloc[:, :5], var_reads.iloc[:, :5], on='position', how='inner')
    if merged_ref_var.empty:
        print(f"Warning: No overlap found in {read_counts_file}")
        continue

    merged_ref_var['depth'] = merged_ref_var['count_x'] + merged_ref_var['count_y']
    merged_ref_var['vaf_x'] = merged_ref_var['count_x'] / merged_ref_var['depth']
    merged_ref_var['vaf_y'] = merged_ref_var['count_y'] / merged_ref_var['depth']

    observed_ref_vector = merged_ref_var['vaf_x']
    observed_var_vector = merged_ref_var['vaf_y']

    df_allele_freq['POS'] = df_allele_freq['POS'].astype(int)
    filtered_freqs = df_allele_freq[df_allele_freq['POS'].isin(merged_ref_var['position'])]
    filtered_freqs.reset_index(drop=True, inplace=True)
    allele_freq_matrix = filtered_freqs.iloc[:, 1:].to_numpy()

    # === CVXPY optimization ===
    A = allele_freq_matrix
    b = observed_var_vector
    c = observed_ref_vector

    x = cvxpy.Variable(strain_count, nonneg=True)
    ps = cvxpy.matmul(A, x)
    ps += np.nextafter(0, 1)  # avoid log(0)

    nllh = -1 * cvxpy.sum(cvxpy.multiply(b, cvxpy.log(ps)) + cvxpy.multiply(c, cvxpy.log(1 - ps)))
    objective = cvxpy.Minimize(nllh)
    constraints = [0 <= x, x <= 1, sum(x) <= 1]
    prob = cvxpy.Problem(objective, constraints)

    try:
        result = prob.solve(verbose=False, solver='ECOS', max_iters=50000, warm_start=True)
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
abundance_df.to_csv("abundance_estimates_combined.csv", index=False)
print("\nAll samples processed. Combined output saved to 'abundance_estimates_combined.csv'.")
