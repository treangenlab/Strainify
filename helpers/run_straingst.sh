#!/usr/bin/env bash
#
# run_straingst.sh
#
# Run StrainGST on FASTQ sample(s) to discover which reference strains are
# present, using a pan-genome DB from build_straingst_db.sh. For each sample it
# writes the StrainGST results, the list of discovered strains, and (optionally)
# collects those strains' genome FASTAs into a folder ready to hand to Strainify.
#
# Per sample:
#   straingst kmerize -> straingst run -> parse strains -> (optional) collect genomes
#
# Output (under --outdir, default: straingst_results/):
#   <sample>/<sample>.hdf5            kmerized sample
#   <sample>/<sample>.straingst.tsv   StrainGST results table
#   <sample>/<sample>.strains.txt     discovered strain names (one per line)
#   <sample>/genomes/                 (if --ref-dir) FASTAs of discovered strains
#   discovered_strains.tsv            summary across all samples (sample,strain,score)

set -euo pipefail

usage() {
    cat <<'EOF'
Run StrainGST on FASTQ sample(s) and report the strains discovered in each.

USAGE:
  run_straingst.sh --db pan-genome-db.hdf5 (--fastq-folder DIR | --fastq1 R1 [...]) [options]

REQUIRED:
  --db FILE              pan-genome-db.hdf5 from build_straingst_db.sh

INPUT (choose one):
  --fastq-folder DIR     a folder of fastq sets (paired or single)
  --fastq1 FILE          a single set's reads (add --fastq2 for paired)
  --fastq2 FILE          mate file for a single paired set
  --sample-id NAME       sample name for the single-set form (default: from filename)
  --read-type T          paired | single  (folder mode; default: paired)

OUTPUT:
  --outdir DIR           output directory (default: straingst_results)
  --ref-dir DIR          the DB's genome folder (build_straingst_db.sh's db/);
                         if given, discovered strains' FASTAs are copied to
                         <outdir>/<sample>/genomes/ for Strainify
  --no-collect           skip collecting genome FASTAs even if --ref-dir is set

STRAINGST TUNING (full flexibility):
  --kmer-size K          kmerize -k (MUST match the DB's k; default: 23)
  --kmerize-args "..."   extra flags passed verbatim to `straingst kmerize`
  --run-args "..."       extra flags passed verbatim to `straingst run`, e.g.:
                           -i N           max iterations / strains to find
                           --minscore X   minimum score to report a strain
                           --fraction X   minimum fraction of k-mers
                         (check `straingst run --help` for your version)
  --threads N            threads (default: 8)

ENVIRONMENT (optional; otherwise straingst must be on PATH):
  --strainge-env NAME    conda env with straingst (uses `conda run -n`)

  -h, --help             this help

EXAMPLES:
  # whole folder of single-end samples, collect genomes for Strainify
  run_straingst.sh --db straingst_db/pan-genome-db.hdf5 \
      --ref-dir straingst_db/db \
      --fastq-folder fastqs/single --read-type single \
      --strainge-env strainge --run-args "-i 5 --minscore 0.1"

  # then feed one sample's genomes to Strainify:
  strainify run --genome_folder straingst_results/<sample>/genomes \
      --fastq_folder fastqs/single --read_type single -profile conda
EOF
}

# ---- defaults ----
DB=""
FASTQ_FOLDER=""
FASTQ1=""
FASTQ2=""
SAMPLE_ID=""
READ_TYPE="paired"
OUTDIR="straingst_results"
REF_DIR=""
COLLECT=1
THREADS=8
KMER_SIZE=23
KMERIZE_ARGS=""
RUN_ARGS=""
STRAINGE_ENV=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --db)            DB="$2"; shift 2;;
        --fastq-folder)  FASTQ_FOLDER="$2"; shift 2;;
        --fastq1)        FASTQ1="$2"; shift 2;;
        --fastq2)        FASTQ2="$2"; shift 2;;
        --sample-id)     SAMPLE_ID="$2"; shift 2;;
        --read-type)     READ_TYPE="$2"; shift 2;;
        --outdir)        OUTDIR="$2"; shift 2;;
        --ref-dir)       REF_DIR="$2"; shift 2;;
        --no-collect)    COLLECT=0; shift;;
        --threads)       THREADS="$2"; shift 2;;
        --kmer-size)     KMER_SIZE="$2"; shift 2;;
        --kmerize-args)  KMERIZE_ARGS="$2"; shift 2;;
        --run-args)      RUN_ARGS="$2"; shift 2;;
        --strainge-env)  STRAINGE_ENV="$2"; shift 2;;
        -h|--help)       usage; exit 0;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2;;
    esac
done

_sg() { if [[ -n "$STRAINGE_ENV" ]]; then conda run -n "$STRAINGE_ENV" "$@"; else "$@"; fi; }

[[ -n "$DB" ]]   || { echo "ERROR: --db is required" >&2; usage >&2; exit 2; }
[[ -f "$DB" ]]   || { echo "ERROR: DB not found: $DB" >&2; exit 2; }
mkdir -p "$OUTDIR"
printf 'sample\tstrain\tscore\n' > "$OUTDIR/discovered_strains.tsv"

# Where to look for a discovered strain's genome FASTA. The build script puts
# genomes in <build-outdir>/db, so accept either that dir or its db/ subfolder --
# this way --ref-dir works whether you pass the build outdir or .../db.
REF_SEARCH_DIRS=()
if [[ -n "$REF_DIR" ]]; then
    [[ -d "$REF_DIR/db" ]] && REF_SEARCH_DIRS+=("$REF_DIR/db")
    REF_SEARCH_DIRS+=("$REF_DIR")
fi

# Process one sample: name followed by its read file(s).
run_one() {
    local sample="$1"; shift
    local sdir="$OUTDIR/$sample"
    mkdir -p "$sdir"

    echo ">>> [$sample] kmerize (k=$KMER_SIZE)"
    _sg straingst kmerize -k "$KMER_SIZE" $KMERIZE_ARGS -o "$sdir/$sample.hdf5" "$@"

    echo ">>> [$sample] straingst run"
    _sg straingst run $RUN_ARGS "$DB" "$sdir/$sample.hdf5" -o "$sdir/$sample.straingst.tsv"

    # The results have two blocks; the per-strain block starts at a header line
    # whose first two fields are 'i' and 'strain'. Strain name is column 2,
    # score is the last column.
    awk -F'\t' '
        $1=="i" && $2=="strain" { seen=1; next }
        seen && NF>=2 { print $2 }
    ' "$sdir/$sample.straingst.tsv" > "$sdir/$sample.strains.txt"

    awk -F'\t' -v s="$sample" '
        $1=="i" && $2=="strain" { seen=1; next }
        seen && NF>=2 { print s"\t"$2"\t"$NF }
    ' "$sdir/$sample.straingst.tsv" >> "$OUTDIR/discovered_strains.tsv"

    local n; n=$(wc -l < "$sdir/$sample.strains.txt"); n="$(echo "$n" | tr -d ' ')"
    echo ">>> [$sample] discovered $n strain(s)"

    # Optionally collect the discovered strains' genome FASTAs for Strainify.
    if [[ "$COLLECT" -eq 1 && -n "$REF_DIR" ]]; then
        mkdir -p "$sdir/genomes"
        local missing=0 strain g ext d
        while IFS= read -r strain; do
            [[ -n "$strain" ]] || continue
            g=""
            for d in "${REF_SEARCH_DIRS[@]}"; do
                for ext in fa.gz fasta.gz fna.gz fa fasta fna; do
                    if [[ -e "$d/$strain.$ext" ]]; then g="$d/$strain.$ext"; break 2; fi
                done
            done
            if [[ -n "$g" ]]; then
                cp "$g" "$sdir/genomes/"
            else
                echo "WARN: [$sample] no genome FASTA for '$strain' under $REF_DIR (or $REF_DIR/db)" >&2
                missing=$((missing+1))
            fi
        done < "$sdir/$sample.strains.txt"
        echo ">>> [$sample] genomes for Strainify -> $sdir/genomes/  (missing: $missing)"
    fi
}

# ---- iterate samples (single set or folder), mirroring strainify's globbing ----
found_any=0
shopt -s nullglob

if [[ -n "$FASTQ1" ]]; then
    if [[ -n "$SAMPLE_ID" ]]; then
        sid="$SAMPLE_ID"
    else
        sid="$(basename "$FASTQ1")"; sid="${sid%%.*}"
    fi
    if [[ -n "$FASTQ2" ]]; then run_one "$sid" "$FASTQ1" "$FASTQ2"; else run_one "$sid" "$FASTQ1"; fi
    found_any=1

elif [[ -n "$FASTQ_FOLDER" ]]; then
    if [[ "$READ_TYPE" == "paired" ]]; then
        for r1 in "$FASTQ_FOLDER"/*_[rR]1.f*q "$FASTQ_FOLDER"/*_[rR]1.f*q.gz; do
            [[ -e "$r1" ]] || continue
            r2="${r1/_R1./_R2.}"; r2="${r2/_r1./_r2.}"
            if [[ ! -e "$r2" ]]; then
                echo "WARN: no mate found for $r1 (expected $r2); skipping" >&2; continue
            fi
            sample="$(basename "$r1")"; sample="${sample%%_[rR]1.*}"
            run_one "$sample" "$r1" "$r2"
            found_any=1
        done
        [[ "$found_any" -eq 1 ]] || { echo "No paired fastq sets found in $FASTQ_FOLDER" >&2; exit 1; }
    else
        for r1 in "$FASTQ_FOLDER"/*.f*q "$FASTQ_FOLDER"/*.f*q.gz; do
            [[ -e "$r1" ]] || continue
            case "$(basename "$r1")" in *_[rR][12].*) continue;; esac
            sample="$(basename "$r1")"; sample="${sample%.f*q}"; sample="${sample%.f*q.gz}"
            run_one "$sample" "$r1"
            found_any=1
        done
        [[ "$found_any" -eq 1 ]] || { echo "No single-end fastq files found in $FASTQ_FOLDER" >&2; exit 1; }
    fi
else
    echo "ERROR: provide --fastq-folder DIR or --fastq1 FILE [--fastq2 FILE]" >&2
    usage >&2; exit 2
fi

echo ""
echo ">>> Done. Per-sample results under: $OUTDIR/<sample>/"
echo ">>> Summary of all discovered strains: $OUTDIR/discovered_strains.tsv"
if [[ -n "$REF_DIR" && "$COLLECT" -eq 1 ]]; then
    echo ">>> Hand a sample's genomes to Strainify, e.g.:"
    echo "      strainify run --genome_folder $OUTDIR/<sample>/genomes --fastq_folder <fastqs> --read_type $READ_TYPE -profile conda"
fi