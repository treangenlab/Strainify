import pandas as pd

def vcf_to_dataframe(vcf_file):
    with open(vcf_file, 'r') as vcf:
        data = []
        # Read the file line by line
        header = None
        for line in vcf:
            # Skip metadata lines that start with ##
            if line.startswith("##"):
                continue

            # The header line starts with a single #
            if line.startswith("#"):
                # Set header as the columns from the # line
                header = line.strip().lstrip('#').split("\t")
                continue

            # Parse data lines (variant rows)
            row = line.strip().split("\t")
            data.append(row)

        # Create DataFrame from the data and header
        df = pd.DataFrame(data, columns=header)
        df = df[df['ALT'] != '<INV>']  # remove indels that are inversions



        # Only keep rows with exactly one alternative allele
        filtered_df = df[df['ALT'].apply(lambda x: len(set(x.split(','))) == 1)]
        filtered_df['POS'] = filtered_df.apply(lambda x: int(x['POS']) + 1 if 'SVTYPE=DEL' in str(x['INFO']) else x['POS'], axis=1) # add 1 to positions with deletions
        filtered_df['ALT'] = filtered_df.apply(lambda x: f"-{x['REF'][1:]}" if 'SVTYPE=DEL' in str(x['INFO']) else x['ALT'], axis=1) #show deletions in this format '-xxx'
        filtered_df['ALT'] = filtered_df.apply(lambda x: f"+{x['ALT'][1:]}" if 'SVTYPE=INS' in str(x['INFO']) else x['ALT'], axis=1) #show insertions in this format '+xxx'
        filtered_df['REF'] = filtered_df.apply(lambda x: x['REF'][1] if 'SVTYPE=DEL' in str(x['INFO']) else x['REF'], axis=1)  # change ref base to the 2nd base in deleted sequence
        #filtered_df = filtered_df[filtered_df['ALT']!='<INV>'] # remove indels that are inversions
        filtered_df['CHROM'] = filtered_df['CHROM'].apply(lambda x: x.split("#")[-1])

        filtered_df.reset_index(drop=True, inplace=True)

    return filtered_df

def calculate_allele_frequency(allele_list, allele_num_tot):
    allele = []
    for i in allele_list:
        allele_freq = int(i)/allele_num_tot
        allele.append(allele_freq)
        return allele


vcf_file = "./W0008_as_ref.vcf"  # Replace with your VCF file path
filtered_df = vcf_to_dataframe(vcf_file)
filtered_df['POS'] = filtered_df['POS'].astype(int)
filtered_df = filtered_df[~filtered_df['POS'].duplicated(keep=False)]
pd.set_option('display.max_columns', None)
#filtered_df.to_csv('filtered_allele_frequency.csv', index=False)

# include the ref column when the ref is one of the genomes of interest (second range from 9 instead of 10)
df_AF = filtered_df.iloc[:, list(range(0, 5)) + list(range(9, filtered_df.shape[1]))]
print(df_AF)

df_AF = df_AF[~df_AF['POS'].duplicated(keep=False)]

allele_lists = df_AF.iloc[:, list(range(5, df_AF.shape[1]))].values.tolist()

# Calculate allele frequencies
allele_freq_list = []
strain_count = len(allele_lists[0])
for allele_list in allele_lists:
    allele_freqs = []
    for i in allele_list:
        AF = int(i)
        allele_freqs.append(AF)
    allele_freq_list.append(allele_freqs)


AF_table_header=df_AF['POS'].values.tolist()
df_allele_freq = pd.DataFrame(allele_freq_list).T

df_allele_freq.columns = AF_table_header
df_allele_freq = df_allele_freq.T
positions = df_AF['POS'].tolist()
print(len(positions))
#positions.reset_index(drop=True, inplace=True)
print(df_allele_freq.shape)
df_allele_freq.insert(0, 'position', positions)
strains = list(df_AF.columns[5:])
strains.insert(0, 'position')
print(strains)
df_allele_freq.columns = strains
#df_allele_freq.to_csv("allele_freq_all.csv", index=False)

def find_similar_columns(df, threshold=5):
    similar_columns=[]
    for i in range(df.shape[1]):
        for j in range(i + 1, df.shape[1]):
            diff_count = (df.iloc[:, i] != df.iloc[:, j]).sum()
            if diff_count < threshold:
                print(f"Columns {df.columns[i]} and {df.columns[j]} have {diff_count} differences (less than {threshold})")
                similar_columns.append((df.columns[i],df.columns[j],int(diff_count)))

    return similar_columns

# Find similar columns
similar_columns=find_similar_columns(df_allele_freq)
print(f"these are similar columns: {similar_columns}")

# Create sites for bam-readcounts to count allele frequencies in aligned reads (BAM file)
df_sites = df_AF.iloc[:, list(range(0,2))]

#df_sites.to_csv('sites_merged.txt', sep='\t', index=False, header=None)

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
        df_AF.to_csv('df_AF.csv')

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
        #merged_ref_var.to_csv('merged_ref_var.csv')

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
        #ref_reads.to_csv('ref_reads.csv')

        # Filter frequencies
        df_allele_freq['position'] = df_allele_freq['position'].astype(int)
        merged_ref_var['position'] = merged_ref_var['position'].astype(int)
        filtered_freqs = df_allele_freq[df_allele_freq['position'].isin(merged_ref_var['position'])]

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
        #lambda_reg = 0.1  # Regularization parameter
        #nllh = (-1) * cvxpy.sum(
            #(cvxpy.multiply(b, cvxpy.log(ps))) + cvxpy.multiply(c, cvxpy.log(1 - ps))) + lambda_reg * cvxpy.norm(x,
                                                                                                                 #"fro")
        nllh = (-1) * cvxpy.sum((cvxpy.multiply(b, cvxpy.log(ps))) + cvxpy.multiply(c, cvxpy.log(1 - ps)))


        objective = cvxpy.Minimize(nllh)
        constraints = [0 <= x, x <= 1, sum(x) <= 1]
        prob = cvxpy.Problem(objective, constraints)
        try:
            result = prob.solve(verbose=True, solver='ECOS', max_iters=50000,warm_start=True)
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
        strain_names = pd.DataFrame(strains[1:])  # Adjust this if `strains` is a list
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
combined_results.to_csv('W0008_no_IS.csv', index=False)
print(combined_results)


