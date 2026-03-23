#!/usr/bin/env bash
set -euo pipefail

THREADS=8
GENOME_DIR=""
FASTQ_DIR=""
OUTDIR=""
CONFIG_PATH=""
SNAKEMAKE_CONFIG_ARGS=()

usage() {
    cat <<EOF
Usage:
  ./strainify filter-run --genomes <genome_dir> --fastqs <fastq_dir> --out <output_dir> --configfile <strainify_config.yaml> [--threads <n>] [--config key=value ...]

Required arguments:
  -g, --genomes   Path to folder of input genomes
  -f, --fastqs    Path to folder of FASTQ files
  -o, --out       Output directory
  -c, --configfile    Path to Strainify config.yaml

Optional arguments:
  -t, --threads   Number of threads (default: 8)
      --config    Additional Snakemake config overrides, e.g. --config weight_by_entropy=true
  -h, --help      Show this help message

Example:
  bash filter_genomes.sh \\
    --genomes /path/to/genomes \\
    --fastqs /path/to/fastqs \\
    --out /path/to/output \\
    --configfile /path/to/config.yaml \\
    --threads 12
EOF
    exit 1
}

# ----------------------------
# Parse flags
# ----------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -g|--genomes)
            [[ $# -ge 2 ]] || usage
            GENOME_DIR="$2"
            shift 2
            ;;
        -f|--fastqs)
            [[ $# -ge 2 ]] || usage
            FASTQ_DIR="$2"
            shift 2
            ;;
        -o|--out)
            [[ $# -ge 2 ]] || usage
            OUTDIR="$2"
            shift 2
            ;;
        -c|--configfile)
            [[ $# -ge 2 ]] || usage
            CONFIG_PATH="$2"
            shift 2
            ;;
        -t|--threads)
            [[ $# -ge 2 ]] || usage
            THREADS="$2"
            shift 2
            ;;
        --config)
            shift
            while [[ $# -gt 0 && "$1" != -* ]]; do
                SNAKEMAKE_CONFIG_ARGS+=("$1")
                shift
            done
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            usage
            ;;
    esac
done

# ----------------------------
# Validate inputs
# ----------------------------
if [[ -z "$GENOME_DIR" || -z "$FASTQ_DIR" || -z "$OUTDIR" || -z "$CONFIG_PATH" ]]; then
    echo "ERROR: --genomes, --fastqs, --out, and --configfile are required." >&2
    usage
fi

if [[ ! -d "$GENOME_DIR" ]]; then
    echo "ERROR: Genome directory does not exist: $GENOME_DIR" >&2
    exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "ERROR: FASTQ directory does not exist: $FASTQ_DIR" >&2
    exit 1
fi

if [[ ! -f "$CONFIG_PATH" ]]; then
    echo "ERROR: Config file does not exist: $CONFIG_PATH" >&2
    exit 1
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
    echo "ERROR: --threads must be a positive integer." >&2
    exit 1
fi

for cmd in bwa samtools snakemake awk sort comm cp ls; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: Required command not found: $cmd" >&2
        exit 1
    }
done

mkdir -p "$OUTDIR"
RENAMED_DIR="$OUTDIR/renamed_genomes"
MAP_DIR="$OUTDIR/mapping"
COV_DIR="$OUTDIR/coverage"
FILTERED_DIR="$OUTDIR/filtered_genomes"

mkdir -p "$RENAMED_DIR" "$MAP_DIR" "$COV_DIR" "$FILTERED_DIR"

CONCAT_FASTA="$OUTDIR/concat_genome.fna"
CONTIG_MAP="$OUTDIR/contig_to_genome.tsv"
GENOME_LIST="$OUTDIR/all_input_genomes.txt"
ANY_COVERED="$OUTDIR/.any_covered.tmp"
ZERO_OUT="$OUTDIR/genomes_with_zero_coverage_across_all_samples.txt"

: > "$CONTIG_MAP"
: > "$GENOME_LIST"
: > "$ANY_COVERED"

echo "Preparing genomes..."

found_genome=0
for genome in "$GENOME_DIR"/*.fa "$GENOME_DIR"/*.fna "$GENOME_DIR"/*.fasta; do
    [[ -e "$genome" ]] || continue
    found_genome=1

    genome_base=$(basename "$genome")
    genome_name="${genome_base%.*}"

    echo "$genome_name" >> "$GENOME_LIST"

    renamed="$RENAMED_DIR/${genome_name}.renamed.fna"

    awk -v g="$genome_name" -v map="$CONTIG_MAP" '
        /^>/{
            name=substr($0,2)
            split(name,a,/ /)
            contig=a[1]
            new=g "|" contig
            print new "\t" g >> map
            print ">" new
            next
        }
        {print}
    ' "$genome" > "$renamed"
done

if [[ "$found_genome" -eq 0 ]]; then
    echo "ERROR: No genome FASTA files found in $GENOME_DIR" >&2
    exit 1
fi

sort -u "$GENOME_LIST" -o "$GENOME_LIST"

echo "Concatenating genomes..."
cat "$RENAMED_DIR"/*.fna > "$CONCAT_FASTA"

echo "Indexing reference..."
bwa index "$CONCAT_FASTA"
samtools faidx "$CONCAT_FASTA"

process_sample() {
    local sample="$1"
    local mode="$2"
    local r1="$3"
    local r2="${4:-}"

    local sam="$MAP_DIR/${sample}.sam"
    local bam="$MAP_DIR/${sample}.bam"
    local cov="$COV_DIR/${sample}.coverage.tsv"

    echo "Mapping sample $sample"

    if [[ "$mode" == "paired" ]]; then
        bwa mem -t "$THREADS" "$CONCAT_FASTA" "$r1" "$r2" > "$sam"
    else
        bwa mem -t "$THREADS" "$CONCAT_FASTA" "$r1" > "$sam"
    fi

    samtools view "$sam" -b -F256 -F2048 -F4 -q60 \
        | samtools sort -o "$bam" --write-index -@ "$THREADS"

    rm -f "$sam"

    samtools coverage "$bam" > "$cov"

    awk '
        NR>1 && $5>0 {
            split($1,a,"|")
            print a[1]
        }
    ' "$cov" >> "$ANY_COVERED"
}

echo "Processing FASTQ files..."

found_fastq=0

for r1 in "$FASTQ_DIR"/*_r1.fq "$FASTQ_DIR"/*_R1.fq "$FASTQ_DIR"/*_r1.fastq "$FASTQ_DIR"/*_R1.fastq; do
    [[ -e "$r1" ]] || continue
    found_fastq=1

    r2="${r1/_r1/_r2}"
    r2="${r2/_R1/_R2}"

    if [[ ! -f "$r2" ]]; then
        echo "ERROR: Missing mate pair for $r1" >&2
        exit 1
    fi

    base=$(basename "$r1")

    sample="$base"
    sample="${sample%_r1.fq}"
    sample="${sample%_r2.fq}"
    sample="${sample%_R1.fq}"
    sample="${sample%_R2.fq}"
    sample="${sample%_r1.fastq}"
    sample="${sample%_r2.fastq}"
    sample="${sample%_R1.fastq}"
    sample="${sample%_R2.fastq}"

    process_sample "$sample" paired "$r1" "$r2"
done

for fq in "$FASTQ_DIR"/*.fq "$FASTQ_DIR"/*.fastq; do
    [[ -e "$fq" ]] || continue

    base=$(basename "$fq")

    if [[ "$base" == *_r1* || "$base" == *_R1* || "$base" == *_r2* || "$base" == *_R2* ]]; then
        continue
    fi

    found_fastq=1
    sample="${base%.*}"

    process_sample "$sample" single "$fq"
done

if [[ "$found_fastq" -eq 0 ]]; then
    echo "ERROR: No FASTQ files found in $FASTQ_DIR" >&2
    exit 1
fi

echo "Summarizing coverage..."

sort -u "$ANY_COVERED" -o "$ANY_COVERED"
comm -23 "$GENOME_LIST" "$ANY_COVERED" > "$ZERO_OUT" || true

echo "Copying genomes with coverage..."

while read -r g; do
    [[ -n "$g" ]] || continue
    src=$(ls "$GENOME_DIR"/"$g".* 2>/dev/null | head -n1 || true)
    if [[ -n "$src" ]]; then
        cp "$src" "$FILTERED_DIR/"
    fi
done < "$ANY_COVERED"

n_filtered=$(ls "$FILTERED_DIR"/*.fna "$FILTERED_DIR"/*.fa "$FILTERED_DIR"/*.fasta 2>/dev/null | wc -l || true)
if [[ "$n_filtered" -eq 0 ]]; then
    echo "ERROR: genome filtering removed all genomes. None of the input genomes are likely to be in the input metagenomic data." >&2
    exit 1
fi

echo
echo "Finished genome filtering."
echo
echo "Genomes with coverage copied to:"
echo "  $FILTERED_DIR"
echo
echo "Genomes with zero coverage listed in:"
echo "  $ZERO_OUT"
echo
echo "Running Strainify..."
#snakemake --cores "$THREADS" --configfile "$CONFIG_PATH" --config genome_folder="$FILTERED_DIR"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SNAKEFILE="$SCRIPT_DIR/Snakefile"

snakemake --snakefile "$SNAKEFILE" --cores "$THREADS" --configfile "$CONFIG_PATH" --config genome_folder="$FILTERED_DIR" "${SNAKEMAKE_CONFIG_ARGS[@]}"