import re
import argparse
import sys
import shutil
from collections import defaultdict
import os
import subprocess
import tempfile

def split_maf_by_contig_combination(maf_filename, output_dir):
    """
    Split MAF file into sub-MAFs based on the unique set of contigs in each block.
    """
    os.makedirs(output_dir, exist_ok=True)

    maf_by_contig_set = defaultdict(list)  # key: frozenset of contigs, value: list of (block_lines, sample_name, reference_contig)

    current_block = []
    contigs_in_block = set()
    sample_name = None
    reference_contig = None

    def process_block(block_lines, contig_set, sample, ref_contig):
        if block_lines and contig_set and sample and ref_contig:
            key = frozenset(contig_set)
            maf_by_contig_set[key].append((block_lines, sample, ref_contig))

    with open(maf_filename, 'r') as maf_file:
        for line in maf_file:
            if line.startswith('a'):
                process_block(current_block, contigs_in_block, sample_name, reference_contig)
                current_block = [line]
                contigs_in_block = set()
                sample_name = None
                reference_contig = None
            elif line.startswith('s'):
                current_block.append(line)
                parts = line.strip().split()
                if len(parts) > 1 and '#' in parts[1]:
                    sample, contig = parts[1].split('#', 1)
                    contigs_in_block.add(contig)
                    if reference_contig is None:
                        reference_contig = contig
                        sample_name = sample
                else:
                    print(f"Warning: Could not parse sample#contig from line:\n{line}", file=sys.stderr)
            else:
                current_block.append(line)

        process_block(current_block, contigs_in_block, sample_name, reference_contig)

    # Write sub-MAFs
    output_mafs = []
    for i, (contig_set, blocks) in enumerate(maf_by_contig_set.items(), 1):
        for j, (block_lines, sample, ref_contig) in enumerate(blocks):
            output_filename = f"contigset{i}--{sample}#{ref_contig}.maf"
            output_path = os.path.join(output_dir, output_filename)
            with open(output_path, 'a') as out_file:  # use append to group blocks with same contig set
                out_file.writelines(block_lines)
            if output_path not in output_mafs:
                output_mafs.append(output_path)

    return output_mafs

def get_sample_contig_pairs(maf_filename):
    """
    Extract unique (sample, contig) pairs from the MAF file.
    """
    pairs = set()
    with open(maf_filename, 'r') as f:
        for line in f:
            if line.startswith('s'):
                parts = line.strip().split()
                if len(parts) > 1 and '#' in parts[1]:
                    sample, contig = parts[1].split('#')
                    pairs.add((sample, contig))
    return sorted(pairs)

def run_wgatools_call(maf_file, sample, contig, output_vcf):
    """
    Run wgatools call on a given sample and contig within a MAF file.
    """
    query_name = f"{sample}#{contig}"
    cmd = [
        "wgatools", "call", maf_file, "-s", "-l0",
        "--sample", sample,
        "--query-name", query_name
    ]
    #print("Running:", ' '.join(cmd))
    with open(output_vcf, "w") as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)


def is_vcf_empty(vcf_file):
    """Check if VCF file contains only headers (lines starting with #)."""
    result = subprocess.run(
        ["bcftools", "view", "-H", vcf_file],
        capture_output=True,
        text=True
    )
    return result.stdout.strip() == ""


def extract_contig_name_from_filename(vcf_path):
    filename = os.path.basename(vcf_path)
    match = re.search(r'contigset\d+--(.+)\.vcf$', filename)
    if not match:
        raise ValueError(f"Could not extract contig name from {filename}")
    contig = match.group(1)

    return contig


def reheader_vcf_inplace_with_contig(vcf_path, contig_length = 1000000):

    contig_name = extract_contig_name_from_filename(vcf_path)

    # Create a temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as tmpdir:
        header_path = os.path.join(tmpdir, "full_header.txt")
        reheadered_vcf_path = os.path.join(tmpdir, "reheadered.vcf")

        # Extract full original header + add contig line
        contig_line = f"##contig=<ID={contig_name},length={contig_length}>\n"
        with open(vcf_path, 'r') as fin, open(header_path, 'w') as fout:
            for line in fin:
                if line.startswith("#CHROM"):
                    # Insert contig line just before #CHROM
                    fout.write(contig_line)
                    fout.write(line)
                    break
                else:
                    fout.write(line)

        # Run bcftools reheader with new header
        cmd = [
            "bcftools", "reheader",
            "-h", header_path,
            "-o", reheadered_vcf_path,
            vcf_path
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"bcftools reheader failed:\n{result.stderr}")

        # Replace original VCF with reheadered VCF (handle cross-device)
        try:
            os.replace(reheadered_vcf_path, vcf_path)
        except OSError:
            # Different devices, fallback to copy+remove
            shutil.copy2(reheadered_vcf_path, vcf_path)
            os.remove(reheadered_vcf_path)

        #print(f"Reheadered VCF inplace: {vcf_path}")

def compress_index_merge_vcfs(vcf_dir):
    """
    Compress VCFs with bgzip (if needed), index each with bcftools, merge by sample with bcftools concat,
    and index the merged output VCFs.
    If all VCFs for a sample are empty, just copy one as {sample}.vcf.gz.
    """
    vcfs_by_sample = defaultdict(list)

    for filename in os.listdir(vcf_dir):
        if not (filename.endswith(".vcf") or filename.endswith(".vcf.gz")):
            continue

        filepath = os.path.join(vcf_dir, filename)
        sample_name = filename.split("--")[0]
        compressed_vcf = filepath
        reheader_vcf_inplace_with_contig(compressed_vcf)

        #Compress if it's a plain .vcf file
        if filename.endswith(".vcf"):
            #print(f"Compressing {filename}...")
            result = subprocess.run(["bgzip", "-f", filepath],
                                    stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                print(f"Error compressing {filename}:\n{result.stderr}")
                continue
            compressed_vcf = filepath + ".gz"

        # Sort the VCF if needed
        sorted_vcf = compressed_vcf.replace(".vcf.gz", ".sorted.vcf.gz")
        #print(f"Sorting {compressed_vcf} to {sorted_vcf}...")
        result = subprocess.run(
            ["bcftools", "sort", "-O", "z", "-o", sorted_vcf, compressed_vcf],
            stderr=subprocess.PIPE, text=True
        )
        if result.returncode != 0:
            print(f"Error sorting {compressed_vcf}:\n{result.stderr}")
            continue

        # Replace original compressed_vcf with sorted version
        compressed_vcf = sorted_vcf

        # Index the compressed VCF
        #print(f"Indexing {compressed_vcf}...")
        result = subprocess.run(["bcftools", "index", "-f", compressed_vcf],
                                stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            print(f"Error indexing {compressed_vcf}:\n{result.stderr}")
            continue

        vcfs_by_sample[sample_name].append(compressed_vcf)

    merged_vcfs = []

    for sample, vcf_list in vcfs_by_sample.items():
        vcf_list.sort()
        merged_vcf = os.path.join(vcf_dir, f"{sample}.vcf.gz")

        #Check if all VCFs are empty (no records, just headers)
        all_empty = all(is_vcf_empty(vcf) for vcf in vcf_list)

        if all_empty:
            print(f"All VCFs for sample {sample} are empty. Copying header-only VCF.")
            shutil.copyfile(vcf_list[0], merged_vcf)
        else:
            #print(f"Merging VCFs for sample {sample} into {merged_vcf}...")
            result = subprocess.run(["bcftools", "concat", "-a", "-o", merged_vcf, "-O", "z"] + vcf_list,
                                    stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                print(f"Error merging VCFs for {sample}:\n{result.stderr}")
                continue

        #print(f"Indexing {merged_vcf}...")
        result = subprocess.run(["bcftools", "index", "-f", merged_vcf],
                                stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            print(f"Error indexing merged VCF {merged_vcf}:\n{result.stderr}")
            continue

        merged_vcfs.append(merged_vcf)

    return merged_vcfs

def merge_final_vcfs(all_sample_vcfs, final_merged_vcf):
    """
    Merge all sample VCFs into one multi-sample VCF with bcftools merge.
    """
    cmd_merge = ["bcftools", "merge"] + all_sample_vcfs + ["-0", "-o", final_merged_vcf, "-O", "v"]
    print(f"Merging all samples into {final_merged_vcf}")
    result = subprocess.run(cmd_merge, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"bcftools merge failed:\n{result.stderr}", file=sys.stderr)

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
    parser = argparse.ArgumentParser(description="Split MAF by reference contig and generate merged VCF.")
    parser.add_argument("--maf_file", required=True, help="Input MAF file")
    parser.add_argument("--temp_dir", default="temp_maf", help="Temporary directory for intermediate files")
    args = parser.parse_args()

    maf_file = args.maf_file
    temp_dir = args.temp_dir

    os.makedirs(temp_dir, exist_ok=True)

    # Step 1: Split MAF by reference contig
    sub_mafs = split_maf_by_contig_combination(maf_file, temp_dir)

    # Step 2: For each sub-MAF, generate per-sample-per-contig VCFs
    sample_to_vcfs = defaultdict(list)
    samples = []
    contig_names = []

    for sub_maf in sub_mafs:
        pairs = get_sample_contig_pairs(sub_maf)
        for sample, contig in pairs:
            sample_name = f"{sample}.vcf.gz"
            if sample_name not in samples:
                samples.append(sample_name)
            maf_basename = os.path.basename(sub_maf).split('.maf')[0]
            out_vcf = os.path.join(temp_dir, f"{sample}--{contig}--{maf_basename}.vcf")

            #out_vcf = os.path.join(temp_dir, f"{sample}.{contig}.vcf")
            run_wgatools_call(sub_maf, sample, contig, out_vcf)
            if os.path.exists(out_vcf) and os.path.getsize(out_vcf) > 0:
                sample_to_vcfs[sample].append(out_vcf)
            else:
                print(f"Skipping empty or missing VCF: {out_vcf}", file=sys.stderr)

    compress_index_merge_vcfs(temp_dir)

# Merge all VCF files and remove FORMAT/QI
    compressed_vcf_files = [os.path.join(temp_dir, sample) for sample in samples]
    directory = os.path.dirname(maf_file)
    merged_vcf_output = os.path.join(directory, "merged.vcf")
    #merged_vcf_output = "merged.vcf"
    #print(compressed_vcf_files)
    intermediate_vcf = "merged_temp.vcf"

    # Run bcftools merge
    merge_command = ["bcftools", "merge"] + compressed_vcf_files + ["-0", "-o", intermediate_vcf, "-O", "v"]
    result = subprocess.run(merge_command, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Merge failed:\n{result.stderr}")
        return

    # Run bcftools annotate
    annotate_command = ["bcftools", "annotate", "--remove", "FORMAT/QI", intermediate_vcf, "-o", merged_vcf_output]
    result = subprocess.run(annotate_command, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Annotate failed:\n{result.stderr}")
        return

    modify_vcf(merged_vcf_output)
    shutil.rmtree(temp_dir)
    os.remove(intermediate_vcf)

if __name__ == "__main__":
    main()
