import subprocess
import os
from collections import defaultdict
import re
import argparse


def process_maf_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    extracted_data = set()
    ordered_pairs = []
    for line in lines:
        if line.startswith('s'):
            parts = line.strip().split('\t')
            if len(parts) > 1:
                sample_contig = parts[1].split('#')
                if len(sample_contig) > 1:
                    sample_name = sample_contig[0]
                    contig_name = sample_contig[1].split('\t')[0]  # Extract up to the second tab
                    current_pair = (sample_name, contig_name)
                    if current_pair not in extracted_data:
                        extracted_data.add(current_pair)
                        ordered_pairs.append(current_pair)

    return ordered_pairs


def execute_wgatools_bcftools(pairs, maf_filename):
    sample_vcf_files = defaultdict(list)

    # Process each sample-contig pair
    for sample_name, contig_name in pairs:
        vcf_filename = f"{sample_name}_{contig_name}.vcf"
        query_name = f"{sample_name}#{contig_name}"
        command = [
            "wgatools", "call", maf_filename, "-s", "-l0",
            "--sample", sample_name,
            "--query-name", query_name
        ]
        with open(vcf_filename, "w") as output_file:
            subprocess.run(command, stdout=output_file, stderr=subprocess.PIPE, text=True)

        if os.path.exists(vcf_filename) and os.path.getsize(vcf_filename) > 0:  # Ensure file exists and is not empty
            sample_vcf_files[sample_name].append(vcf_filename)
            print(f"VCF file stored: {vcf_filename}")
        else:
            print(f"Skipping empty or missing VCF: {vcf_filename}")

    compressed_vcf_files = []
    index_files = []

    # Compress & Merge VCF files per sample
    for sample_name, vcf_files in sample_vcf_files.items():
        compressed_vcf_files_for_sample = []

        for vcf_file in vcf_files:
            subprocess.run(["bgzip", vcf_file], stderr=subprocess.PIPE, text=True)
            compressed_vcf = f"{vcf_file}.gz"
            if os.path.exists(compressed_vcf):  # Ensure compression succeeded
                compressed_vcf_files_for_sample.append(compressed_vcf)
                print(f"Compressed: {compressed_vcf}")
            else:
                print(f"Failed to compress: {vcf_file}")

        # If multiple contigs for a sample, concatenate the VCFs
        if len(compressed_vcf_files_for_sample) > 1:
            merged_vcf = f"{sample_name}.vcf.gz"
            subprocess.run(["bcftools", "concat", "-o", merged_vcf, "-O", "z"] + compressed_vcf_files_for_sample,
                           stderr=subprocess.PIPE, text=True)
            print(f"Merged VCF for {sample_name}: {merged_vcf}")

            # Index the concatenated VCF
            index_file = f"{merged_vcf}.csi"
            subprocess.run(["bcftools", "index", merged_vcf], stderr=subprocess.PIPE, text=True)
            print(f"Indexed merged VCF: {index_file}")

            compressed_vcf_files.append(merged_vcf)
            index_files.append(index_file)
        elif len(compressed_vcf_files_for_sample) == 1:
            single_vcf = compressed_vcf_files_for_sample[0]
            index_file = f"{single_vcf}.csi"
            subprocess.run(["bcftools", "index", single_vcf], stderr=subprocess.PIPE, text=True)
            print(f"Indexed VCF: {index_file}")

            compressed_vcf_files.append(single_vcf)
            index_files.append(index_file)
        else:
            print(f"No valid VCFs for sample: {sample_name}, skipping.")

    # Ensure we have VCF files to merge
    if not compressed_vcf_files:
        print("No valid VCF files to merge. Exiting.")
        return None

    # Merge all VCF files and remove FORMAT/QI
    merged_vcf_output = "merged.vcf"
    merge_command = ["bcftools", "merge"] + compressed_vcf_files + ["-0"]
    annotate_command = ["bcftools", "annotate", "--remove", "FORMAT/QI"]

    with open(merged_vcf_output, "w") as final_vcf:
        merge_process = subprocess.Popen(merge_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        annotate_process = subprocess.Popen(annotate_command, stdin=merge_process.stdout, stdout=final_vcf,
                                            stderr=subprocess.PIPE, text=True)
        merge_process.stdout.close()
        annotate_process.communicate()

    print(f"Final merged VCF file: {merged_vcf_output}")

    # Cleanup intermediate files (compressed VCFs & index files)
    for vcf_file in compressed_vcf_files:
        os.remove(vcf_file)
        print(f"Removed: {vcf_file}")

    for index_file in index_files:
        os.remove(index_file)
        print(f"Removed index file: {index_file}")

    return merged_vcf_output

def modify_vcf(vcf_filename):
    """Modify the final merged VCF by replacing genotype values."""
    with open(vcf_filename, "r") as file:
        lines = file.readlines()

    with open(vcf_filename, "w") as file:
        for line in lines:
            if line.startswith("#"):  # Keep headers unchanged
                file.write(line)
            else:
                # Replace 0/0 with 0
                line = re.sub(r"\b0/0\b", "0", line)
                # Replace x|x with 1 (e.g., 1|1 -> 1, 2|2 -> 1)
                line = re.sub(r"\b(\d)\|\1\b", "1", line)
                file.write(line)

    print(f"Modified VCF file: {vcf_filename} (Replaced 0/0 with 0, x|x with 1)")

def main():
    parser = argparse.ArgumentParser(description="Process a MAF file to generate a merged VCF.")
    parser.add_argument("--maf_file", help="Path to the input MAF file")

    args = parser.parse_args()
    maf_filename = args.maf_file

    result = process_maf_file(maf_filename)

    if result:
        final_vcf = execute_wgatools_bcftools(result,maf_filename)
        if final_vcf:
            modify_vcf(final_vcf)
            print(f"Final output VCF: {final_vcf}")


if __name__ == "__main__":
    main()
