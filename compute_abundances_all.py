import sys
import os
import glob
import cvxpy
import numpy as np
import pandas as pd
import argparse


def compute_entropy_weights(A):
    # Compute alt allele frequency f for each site (row)
    f = np.mean(A, axis=1)  # shape: (n_variants,)

    # Avoid log(0) by clipping f away from 0 and 1
    epsilon = 1e-10
    f_clipped = np.clip(f, epsilon, 1 - epsilon)

    # Binary entropy formula
    entropy = -f_clipped * np.log2(f_clipped) - (1 - f_clipped) * np.log2(1 - f_clipped)

    return entropy  # shape: (n_variants,)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute abundances")
    parser.add_argument('--read_counts_dir', required=True, help="read counts directory")
    parser.add_argument('--filtered_variants', required=True, help="filtered_variant_matrix.csv")
    parser.add_argument('--weight_by_entropy', action='store_true', help="use weighted likelihood model")

    args = parser.parse_args()
    read_counts_dir = args.read_counts_dir
    variant_matrix_file = args.filtered_variants


    # === Load variant matrix ===
    df_AF = pd.read_csv(variant_matrix_file)
    # Sort by chromosome and position
    if 'CHROM' in df_AF.columns and 'POS' in df_AF.columns:
        df_AF = df_AF.sort_values(['CHROM', 'POS'], ascending=[True, True]).reset_index(drop=True)

    df_allele_freq = df_AF.drop(columns=['ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
    strain_names = df_allele_freq.columns[2:]
    strain_count = len(strain_names)

    # === Store all results ===
    abundance_df = pd.DataFrame({'strain name': strain_names})
    read_count_files = sorted(glob.glob(os.path.join(read_counts_dir, "*_read_counts.tsv")))

    for read_counts_file in read_count_files:
        print(f"Processing {read_counts_file}...")
        df_read_counts = pd.read_csv(read_counts_file, sep='\t')
        # Sort by chrom and position
        if 'chrom' in df_read_counts.columns and 'position' in df_read_counts.columns:
            df_read_counts = df_read_counts.sort_values(['chrom', 'position'], ascending=[True, True]).reset_index(drop=True)
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

        #print(f"Mean depth: {mean_depth:.2f}, Std depth: {std_depth:.2f}")
        #print(f"Filtering positions with depth < {threshold:.2f}")

        # Filter out positions with significantly lower than average depth
        merged_ref_var = merged_ref_var[merged_ref_var['depth'] >= threshold]

        if merged_ref_var.empty:
            print(f"All positions filtered out for {read_counts_file} due to low depth")
            continue

        observed_ref_vector = merged_ref_var['vaf_x']
        observed_var_vector = merged_ref_var['vaf_y']

        df_allele_freq['POS'] = df_allele_freq['POS'].astype(int)


        filtered_freqs = pd.merge(
            df_allele_freq,
            merged_ref_var[['chrom', 'position']],
            left_on=['CHROM', 'POS'],
            right_on=['chrom', 'position'],
            how='inner'
        )

        filtered_freqs.reset_index(drop=True, inplace=True)
        allele_freq_matrix = filtered_freqs.drop(columns=['CHROM', 'POS', 'chrom', 'position'])
        allele_freq_matrix = allele_freq_matrix.to_numpy()


        # === CVXPY optimization ===
        A = allele_freq_matrix
        b = observed_var_vector
        c = observed_ref_vector

        x = cvxpy.Variable(strain_count, nonneg=True)
        ps = cvxpy.matmul(A, x)
        ps += np.nextafter(0, 1)  # avoid log(0)

        # Apply entropy weights to log-likelihood
        if args.weight_by_entropy:
            print("Applying weighted log-likelihood model...")
            w = compute_entropy_weights(A)
            weighted_loglik = cvxpy.multiply(w, cvxpy.multiply(b, cvxpy.log(ps)) + cvxpy.multiply(c, cvxpy.log(1 - ps)))
            nllh = -cvxpy.sum(weighted_loglik)

        else:
            print("Applying unweighted log-likelihood model...")
            nllh = -1 * cvxpy.sum(cvxpy.multiply(b, cvxpy.log(ps)) + cvxpy.multiply(c, cvxpy.log(1 - ps)))

        objective = cvxpy.Minimize(nllh)
        constraints = [0 <= x, x <= 1, sum(x) == 1]
        prob = cvxpy.Problem(objective, constraints)

        try:
            print("Trying ECOS solver...")
            result = prob.solve(verbose=False, solver='ECOS', max_iters=50000, warm_start=True)
        except cvxpy.error.SolverError:
            print("ECOS solver failed. Trying SCS solver...")
            try:
                result = prob.solve(verbose=False, solver='SCS', max_iters=50000, warm_start=True)
            except cvxpy.error.SolverError:
                print(f"Both ECOS and SCS solvers failed for {read_counts_file}")
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
    abundance_path = os.path.join(os.path.dirname(read_counts_dir), "abundance_estimates_combined.csv")
    abundance_df.to_csv(abundance_path, index=False)

    print("\nAll samples processed. Combined output saved to 'abundance_estimates_combined.csv'.")