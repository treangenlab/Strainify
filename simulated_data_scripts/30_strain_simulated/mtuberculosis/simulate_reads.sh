#!/bin/bash

# -------- Configuration --------
input_csv="mtuberculosis_30_ratios.csv"       # CSV with genome,sample1,sample2,...
read_length=250                     # Length per read (paired-end)
fragment_mean=600                   # Mean fragment length
fragment_std=150                    # Fragment length std dev
profile="MSv3"                      # Illumina ART profile (e.g., HiSeq 2500)
genome_folder="mtuberculosis_downloads"
base_output="mtuberculosis"

# List of total metagenome coverages to simulate
coverages=(10 20 50 100 200)

# -------- Parse Header to Get Sample Names --------
header=$(head -n 1 "$input_csv")
IFS=',' read -ra columns <<< "$header"
samples=("${columns[@]:1}")  # skip 'genome'

# -------- Loop Over Each Coverage --------
for total_coverage in "${coverages[@]}"; do
    echo "===================="
    echo "Simulating at ${total_coverage}x total metagenome coverage"
    echo "===================="
    
    output_folder="${base_output}/${total_coverage}x"
    mkdir -p "$output_folder"

    # -------- For Each Sample --------
    for ((i=0; i<${#samples[@]}; i++)); do
        sample="${samples[$i]}"
        echo "Simulating sample: $sample at ${total_coverage}x"

        r1_list="tmp_${sample}_r1_list.txt"
        r2_list="tmp_${sample}_r2_list.txt"
        > "$r1_list"
        > "$r2_list"

        # -------- Process Each Genome --------
        tail -n +2 "$input_csv" | while IFS=',' read -r line; do
            IFS=',' read -ra fields <<< "$line"
            genome="${fields[0]}"
            proportion="${fields[$((i+1))]}"

            [[ -z "$genome" || -z "$proportion" ]] && continue

            base=$(basename "$genome")
            base=${base%.fna}
            indiv_cov=$(echo "$proportion * $total_coverage" | bc -l)

            echo "  - $genome → ${indiv_cov}× coverage"

            /home/Users/rl152/art_bin_MountRainier/art_illumina \
                -ss $profile \
                -i $genome_folder/$genome \
                -p \
                -l $read_length \
                -f $indiv_cov \
                -m $fragment_mean \
                -s $fragment_std \
                -o $output_folder/"${sample}_${base}_sim"

            echo $output_folder/"${sample}_${base}_sim1.fq" >> "$r1_list"
            echo $output_folder/"${sample}_${base}_sim2.fq" >> "$r2_list"
        done

        # -------- Combine Reads --------
        echo "Merging reads for sample: $sample"
        cat $(cat "$r1_list") > $output_folder/"${sample}_r1.fq"
        cat $(cat "$r2_list") > $output_folder/"${sample}_r2.fq"
        rm -f "$r1_list" "$r2_list"

        echo "Output: ${sample}_r1.fq, ${sample}_r2.fq"
    done

    # -------- Clean Intermediate Simulated Reads --------
    echo "Cleaning up intermediate ART files"
    rm -f $output_folder/*_sim*.fq $output_folder/*.aln
    # -------- Clean temp file lists --------
    rm -f tmp_*_r1_list.txt tmp_*_r2_list.txt
done

