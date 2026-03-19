import sys
import os
import glob
import cvxpy
import numpy as np
import pandas as pd
import argparse
import concurrent.futures


def compute_entropy_weights(A):
    f = np.mean(A, axis=1)
    epsilon = 1e-10
    f_clipped = np.clip(f, epsilon, 1 - epsilon)
    entropy = -f_clipped * np.log2(f_clipped) - (1 - f_clipped) * np.log2(1 - f_clipped)
    return entropy


def run_mle(A, b, c, strain_count, read_counts_file, weight_by_entropy):
    x = cvxpy.Variable(strain_count, nonneg=True)
    ps = cvxpy.matmul(A, x)
    ps += np.nextafter(0, 1)  # avoid log(0)

    if weight_by_entropy:
        print("Applying weighted log-likelihood model...")
        w = compute_entropy_weights(A)
        weighted_loglik = cvxpy.multiply(
            w,
            cvxpy.multiply(b, cvxpy.log(ps)) + cvxpy.multiply(c, cvxpy.log(1 - ps))
        )
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
            # continue

    if result is not None:
        x.value = x.value.astype(float)
        return x.value
    else:
        return None



def process_sample(read_counts_path, df_AF, df_allele_freq, strain_names, strain_count,
                   weight_by_entropy, bootstrap, bootstrap_iterations):
    read_counts_file = read_counts_path  # so run_mle can print it unchanged

    print(f"Processing {read_counts_file}...")
    df_read_counts = pd.read_csv(read_counts_file, sep='\t')

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

    merged_ref_var = pd.merge(ref_reads.iloc[:, :5], var_reads.iloc[:, :5], on=['chrom', 'position'], how='inner')

    if merged_ref_var.empty:
        print(f"Warning: No overlap found in {read_counts_file}")
        return None  # <--- was "continue"

    merged_ref_var['depth'] = merged_ref_var['count_x'] + merged_ref_var['count_y']
    merged_ref_var['vaf_x'] = merged_ref_var['count_x'] / merged_ref_var['depth']
    merged_ref_var['vaf_y'] = merged_ref_var['count_y'] / merged_ref_var['depth']

    merged_ref_var['depth'] = merged_ref_var['count_x'] + merged_ref_var['count_y']

    threshold = 0
    merged_ref_var = merged_ref_var[merged_ref_var['depth'] >= threshold]

    if merged_ref_var.empty:
        print(f"All positions filtered out for {read_counts_file} due to low depth")
        return None  # <--- was "continue"

    observed_ref_vector = merged_ref_var['vaf_x']
    observed_var_vector = merged_ref_var['vaf_y']

    df_allele_freq = df_allele_freq.copy()
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

    A = allele_freq_matrix
    b = observed_var_vector
    c = observed_ref_vector


    output = run_mle(A, b, c, strain_count, read_counts_file, weight_by_entropy)

    sample_id = os.path.basename(read_counts_file).replace("_read_counts.tsv", "")

    abund_col = None
    ci_col = None

    if output is not None:
        rounded_arr = np.round(np.array(output) * 100, 4)
        abund_col = rounded_arr
    else:
        print(f"Optimization failed for {read_counts_file}")


    if bootstrap:
        print(f"Performing bootstrap for {sample_id}...")
        n_sites = A.shape[0]
        bootstrap_estimates = []
        for i in range(bootstrap_iterations):
            sampled_indices = np.random.choice(n_sites, n_sites, replace=True)
            A_sampled = A[sampled_indices, :]
            b_sampled = b[sampled_indices]
            c_sampled = c[sampled_indices]

            output_sampled = run_mle(A_sampled, b_sampled, c_sampled, strain_count, read_counts_file, weight_by_entropy)
            if output_sampled is None:
                continue
            bootstrap_estimates.append(output_sampled)

        if len(bootstrap_estimates) == 0:
            print(f"All bootstrap iterations failed for {sample_id}")
        else:
            print(f"Successfully computed bootstrap confidence intervals for {len(bootstrap_estimates)}/{bootstrap_iterations} iterations for {sample_id}")
            bootstrap_estimates = np.vstack(bootstrap_estimates)
            ci_lower = np.percentile(bootstrap_estimates, 2.5, axis=0) * 100
            ci_upper = np.percentile(bootstrap_estimates, 97.5, axis=0) * 100
            ci_col = [f"{lo.round(4)};{hi.round(4)}" for lo, hi in zip(ci_lower, ci_upper)]

    return (sample_id, abund_col, ci_col)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute abundances")
    parser.add_argument('--read_counts_dir', required=True, help="read counts directory")
    parser.add_argument('--filtered_variants', required=True, help="filtered_variant_matrix.csv")
    parser.add_argument('--weight_by_entropy', action='store_true', help="use weighted likelihood model")
    parser.add_argument('--bootstrap', action='store_true', help="perform bootstrap resampling for confidence intervals")
    parser.add_argument('--bootstrap_iterations', type=int, default=100, help="number of bootstrap iterations (default: 100)")
    parser.add_argument('--threads', type=int, default=8, help="number of threads for parallel processing (default:8)")
    args = parser.parse_args()

    read_counts_dir = args.read_counts_dir
    variant_matrix_file = args.filtered_variants

    df_AF = pd.read_csv(variant_matrix_file)
    if 'CHROM' in df_AF.columns and 'POS' in df_AF.columns:
        df_AF = df_AF.sort_values(['CHROM', 'POS'], ascending=[True, True]).reset_index(drop=True)

    df_allele_freq = df_AF.drop(columns=['ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
    strain_names = df_allele_freq.columns[2:]
    strain_count = len(strain_names)

    abundance_df = pd.DataFrame({'strain name': strain_names})
    bootstrap_df = pd.DataFrame({'strain name': strain_names}) if args.bootstrap else None

    read_count_files = sorted(glob.glob(os.path.join(read_counts_dir, "*_read_counts.tsv")))

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = [
            executor.submit(
                process_sample,
                f, df_AF, df_allele_freq, strain_names, strain_count,
                args.weight_by_entropy, args.bootstrap, args.bootstrap_iterations
            )
            for f in read_count_files
        ]

        for fut in concurrent.futures.as_completed(futures):
            res = fut.result()
            if res is None:
                continue
            sample_id, abund_col, ci_col = res

            if abund_col is not None:
                abundance_df[sample_id] = abund_col

            if args.bootstrap and ci_col is not None:
                bootstrap_df[sample_id] = ci_col

    # Save outputs (same as your original)
    abundance_df['strain name'] = abundance_df['strain name'].str.replace(r'\.fna$|\.fasta$', '', regex=True)
    abundance_path = os.path.join(os.path.dirname(read_counts_dir), "abundance_estimates_combined.csv")
    abundance_df.to_csv(abundance_path, index=False)
    print("\nAll samples processed. Combined output saved to 'abundance_estimates_combined.csv'.")

    if bootstrap_df is not None:
        bootstrap_df['strain name'] = bootstrap_df['strain name'].str.replace(r'\.fna$|\.fasta$', '', regex=True)
        bootstrap_path = os.path.join(os.path.dirname(read_counts_dir), "abundance_estimates_bootstrap_CIs.csv")
        bootstrap_df.to_csv(bootstrap_path, index=False)
        print("\nBootstrap confidence intervals saved to 'abundance_estimates_bootstrap_CIs.csv'.")
