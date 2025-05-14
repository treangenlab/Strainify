import sys
import os
import cvxpy
import numpy as np
import pandas as pd

# === Command-line arguments ===
if len(sys.argv) != 3:
    print("Usage: python compute_abundances.py <read_counts.tsv> <filtered_variant_matrix.csv>")
    sys.exit(1)

read_counts_file = sys.argv[1]
variant_matrix_file = sys.argv[2]

# === Load data ===
df_AF = pd.read_csv(variant_matrix_file)
df_allele_freq = df_AF.drop(columns=['CHROM', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
strain_count = len(df_allele_freq.columns[1:])
strains = df_allele_freq.columns[1:]

df_read_counts = pd.read_csv(read_counts_file, sep='\t')
df_read_counts = df_read_counts.drop_duplicates(keep=False).reset_index(drop=True)
df_AF = df_AF.reset_index(drop=True)

# === Filter and merge reads ===
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
    result = prob.solve(verbose=True, solver='ECOS', max_iters=50000, warm_start=True)
except cvxpy.error.SolverError:
    print("Solver failed.")
    sys.exit(1)

if result is not None:
    x.value = x.value.astype(float)
    rounded_arr = np.round(np.array(x.value) * 100, 4)

    #strains = strains[1:]
    strain_names = pd.DataFrame(strains, columns=['strain name'])
    df_results = pd.DataFrame(rounded_arr)
    df_results = pd.concat([strain_names, df_results], axis=1)

    sample_id = os.path.basename(read_counts_file).replace(".tsv", "")
    df_results.columns = ['strain name', sample_id]
    #df_results['strain name'] = df_results['strain name'].str.replace('.fna', '', regex=False)
    df_results['strain name'] = df_results['strain name'].str.replace(r'\.fna$|\.fasta$', '', regex=True)

    df_results.to_csv("abundance_estimates.csv", index=False)
    print(df_results)
else:
    print("Optimization failed; no output saved.")

