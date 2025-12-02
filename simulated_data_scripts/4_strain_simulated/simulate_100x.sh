#!/bin/bash

# -------- Configuration --------
input_csv="/home/Users/rl152/PycharmProjects/tests/4-strain_simulated_ratios.csv"       # CSV with genome,sample1,sample2,...
total_coverage=100                   # Target total metagenome coverage per sample
read_length=250                     # Length per read (paired-end)
fragment_mean=600                   # Mean fragment length
fragment_std=150                    # Fragment length std dev
profile="MSv3"                      # Illumina ART profile (e.g., HiSeq 2500)
genome_folder="/home/Users/rl152/strain_phasing/tests/4-strain-mock-exp/genomes"
output_folder="/home/Users/rl152/PycharmProjects/tests/100x"

# -------- Parse Header to Get Sample Names --------
header=$(head -n 1 "$input_csv")
IFS=',' read -ra columns <<< "$header"
samples=("${columns[@]:1}")  # skip 'genome'

# -------- For Each Sample --------
for ((i=0; i<${#samples[@]}; i++)); do
    sample="${samples[$i]}"
    echo "🔬 Simulating sample: $sample"

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
        indiv_cov=$(echo "$proportion * $total_coverage" | bc)

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

        echo $output_folder/"${sample}_${base}_sim1.fq" >> $output_folder/"$r1_list"
        echo $output_folder/"${sample}_${base}_sim2.fq" >> $output_folder/"$r2_list"
    done

    # -------- Combine Reads --------
    echo "Merging reads for sample: $sample"
    cat $(cat $output_folder/"$r1_list") > $output_folder/"${sample}_r1.fq"
    cat $(cat $output_folder/"$r2_list") > $output_folder/"${sample}_r2.fq"
    rm -f "$r1_list" "$r2_list"
    rm $output_folder/"*_sim*.fq" $output_folder/"*.aln" "*.txt"

    echo "Output: ${sample}_r1.fq, ${sample}_r2.fq"
done

