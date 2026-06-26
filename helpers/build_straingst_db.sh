#!/usr/bin/env bash
#
# build_straingst_db.sh
#
# Build a StrainGST (StrainGE) pan-genome k-mer database, either from genomes
# downloaded from NCBI (by genus / "Genus species") OR from a custom folder of
# genome FASTAs you already have. This is a one-time, network- and compute-heavy
# prep step; the resulting DB is reused for every run_straingst.sh run.
#
# Pipeline (mirrors the StrainGE docs / the reference build commands):
#   genomes -> prepare_strainge_db -> kmerize each -> kmersim all-vs-all
#           -> cluster (de-replicate) -> createdb -> pan-genome-db.hdf5
#
# Output (under --outdir, default: straingst_db/):
#   db/                     prepared reference genomes (<strain>.fa.gz) + .hdf5
#   references_meta.tsv     metadata from prepare_strainge_db
#   similarities.tsv        all-vs-all k-mer similarities
#   clusters.tsv            cluster assignments
#   references_to_keep.txt  cluster representatives
#   pan-genome-db.hdf5      THE DATABASE -> pass to run_straingst.sh --db
#
# The db/ folder is ALSO what run_straingst.sh uses (--ref-dir) to pull the
# actual genome FASTAs for the strains it discovers, to hand to Strainify.

set -euo pipefail

# Directory this script lives in, so we can find a prepare_strainge_db.py copied
# alongside it even when the script is invoked from elsewhere.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

usage() {
    cat <<'EOF'
Build a StrainGST pan-genome database.

USAGE:
  build_straingst_db.sh (--taxa NAME | --genome-folder DIR) [options]

GENOME SOURCE (choose one):
  --taxa NAME            NCBI genera or organisms to download, comma-separated
                         (e.g. "Escherichia,Shigella" or "Escherichia coli").
  --genome-folder DIR    Use your own folder of genome FASTAs instead of NCBI.

OUTPUT:
  --outdir DIR           Output directory (default: straingst_db)

NCBI DOWNLOAD (only with --taxa):
  --assembly-level LVL   ncbi-genome-download -l (default: complete)
  --ncbi-format FMT      ncbi-genome-download -F (default: fasta)
  --ncbi-section SEC     refseq | genbank (default: refseq)
  --ncbi-extra "..."     extra flags passed verbatim to ncbi-genome-download

STRAINGST TUNING:
  --threads N            threads for download/kmersim (default: 8)
  --kmer-size K          k-mer size for kmerize (default: 23)
  --cluster-jaccard C    straingst cluster -C  (default: 0.99)
  --cluster-subset c     straingst cluster -c  (default: 0.947)

ENVIRONMENTS (optional; otherwise tools must be on PATH):
  --download-env NAME    conda env with ncbi-genome-download (uses `conda run -n`)
  --strainge-env NAME    conda env with straingst / prepare_strainge_db
  --prepare-db-cmd CMD   prepare-db command. Default: auto-detect a
                         prepare_strainge_db.py next to this script (run via
                         python3), else `prepare_strainge_db.py` on PATH.
                         Override e.g. "python3 /path/StrainGE/bin/prepare_strainge_db.py"
  --prepare-db           force running prepare_strainge_db (default for --taxa)
  --no-prepare-db        skip prepare_strainge_db and use the genomes as-is
                         (default for --genome-folder; clustering still de-reps)

  -h, --help             this help

EXAMPLES:
  # E. coli + Shigella from NCBI, using two conda envs
  build_straingst_db.sh --taxa "Escherichia,Shigella" \
      --download-env ncbi --strainge-env strainge --threads 16

  # From a custom genome folder
  build_straingst_db.sh --genome-folder my_refs/ --strainge-env strainge
EOF
}

# ---- defaults ----
OUTDIR="straingst_db"
GENOME_FOLDER=""
TAXA=""
ASSEMBLY_LEVEL="complete"
NCBI_FORMAT="fasta"
NCBI_SECTION="refseq"
NCBI_EXTRA=""
THREADS=8
KMER_SIZE=23
CLUSTER_JACCARD="0.99"
CLUSTER_SUBSET="0.947"
DOWNLOAD_ENV=""
STRAINGE_ENV=""
PREPARE_DB_CMD=""   # empty -> auto-detect after arg parsing (see below)
RUN_PREPARE=""      # empty -> decide by mode (NCBI: yes, custom folder: no)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --taxa)            TAXA="$2"; shift 2;;
        --genome-folder)   GENOME_FOLDER="$2"; shift 2;;
        --outdir)          OUTDIR="$2"; shift 2;;
        --assembly-level)  ASSEMBLY_LEVEL="$2"; shift 2;;
        --ncbi-format)     NCBI_FORMAT="$2"; shift 2;;
        --ncbi-section)    NCBI_SECTION="$2"; shift 2;;
        --ncbi-extra)      NCBI_EXTRA="$2"; shift 2;;
        --threads)         THREADS="$2"; shift 2;;
        --kmer-size)       KMER_SIZE="$2"; shift 2;;
        --cluster-jaccard) CLUSTER_JACCARD="$2"; shift 2;;
        --cluster-subset)  CLUSTER_SUBSET="$2"; shift 2;;
        --download-env)    DOWNLOAD_ENV="$2"; shift 2;;
        --strainge-env)    STRAINGE_ENV="$2"; shift 2;;
        --prepare-db-cmd)  PREPARE_DB_CMD="$2"; shift 2;;
        --prepare-db)      RUN_PREPARE=1; shift;;
        --no-prepare-db)   RUN_PREPARE=0; shift;;
        -h|--help)         usage; exit 0;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2;;
    esac
done

# Decide whether to run prepare_strainge_db. It cleans up NCBI downloads (using
# their assembly metadata / human_readable layout), so default it ON for --taxa
# downloads and OFF for a custom --genome-folder -- those genomes are used as-is,
# and the later clustering step still de-replicates near-identical assemblies.
# Override either way with --prepare-db / --no-prepare-db.
if [[ -z "$RUN_PREPARE" ]]; then
    if [[ -n "$GENOME_FOLDER" ]]; then RUN_PREPARE=0; else RUN_PREPARE=1; fi
fi

# Resolve prepare_strainge_db.py only if we're going to run it: explicit
# --prepare-db-cmd wins; otherwise a copy sitting next to this script (run via
# python3 so it needn't be on PATH/executable); else `prepare_strainge_db.py`.
if [[ "$RUN_PREPARE" -eq 1 ]]; then
    if [[ -z "$PREPARE_DB_CMD" ]]; then
        if [[ -f "$SCRIPT_DIR/prepare_strainge_db.py" ]]; then
            PREPARE_DB_CMD="python3 $SCRIPT_DIR/prepare_strainge_db.py"
        else
            PREPARE_DB_CMD="prepare_strainge_db.py"
        fi
    fi
    echo ">>> Using prepare-db command: $PREPARE_DB_CMD"
fi

# run a command in a conda env if one was given, else on PATH
_dl() { if [[ -n "$DOWNLOAD_ENV" ]]; then conda run -n "$DOWNLOAD_ENV" "$@"; else "$@"; fi; }
_sg() { if [[ -n "$STRAINGE_ENV" ]]; then conda run -n "$STRAINGE_ENV" "$@"; else "$@"; fi; }

mkdir -p "$OUTDIR"

# --- Step 1: obtain genomes ---
if [[ -n "$GENOME_FOLDER" ]]; then
    [[ -d "$GENOME_FOLDER" ]] || { echo "ERROR: --genome-folder not found: $GENOME_FOLDER" >&2; exit 2; }
    GENOMES_DIR="$GENOME_FOLDER"
    echo ">>> Using custom genome folder: $GENOMES_DIR"
else
    [[ -n "$TAXA" ]] || { echo "ERROR: provide --taxa NAME or --genome-folder DIR" >&2; usage >&2; exit 2; }
    echo ">>> Downloading '$TAXA' (level: $ASSEMBLY_LEVEL, section: $NCBI_SECTION) from NCBI ..."
    _dl ncbi-genome-download bacteria \
        -l "$ASSEMBLY_LEVEL" \
        -g "$TAXA" \
        -H -F "$NCBI_FORMAT" \
        -s "$NCBI_SECTION" \
        -p "$THREADS" \
        -o "$OUTDIR/ref_genomes" \
        $NCBI_EXTRA
    GENOMES_DIR="$OUTDIR/ref_genomes/human_readable"
fi

# --- Step 2: assemble the reference set into $OUTDIR/db ---
mkdir -p "$OUTDIR/db"
if [[ "$RUN_PREPARE" -eq 1 ]]; then
    echo ">>> Preparing StrainGE reference set -> $OUTDIR/db ..."
    _sg $PREPARE_DB_CMD "$GENOMES_DIR" -o "$OUTDIR/db" > "$OUTDIR/references_meta.tsv"
else
    echo ">>> Copying genomes into $OUTDIR/db (skipping prepare_strainge_db) ..."
    _copied=0
    for f in "$GENOMES_DIR"/*.fa    "$GENOMES_DIR"/*.fasta    "$GENOMES_DIR"/*.fna \
             "$GENOMES_DIR"/*.fa.gz "$GENOMES_DIR"/*.fasta.gz "$GENOMES_DIR"/*.fna.gz; do
        [[ -e "$f" ]] || continue
        cp "$f" "$OUTDIR/db/"; _copied=$((_copied+1))
    done
    [[ "$_copied" -gt 0 ]] || { echo "ERROR: no genome FASTAs (.fa/.fasta/.fna, optionally .gz) found in $GENOMES_DIR" >&2; exit 1; }
    echo ">>> Copied $_copied genome(s). (Each filename becomes a strain name.)"
fi

# --- Step 3: kmerize each reference (accepts gzipped OR unzipped FASTAs) ---
echo ">>> Kmerizing references (k=$KMER_SIZE) ..."
_have_genome=0
for f in "$OUTDIR"/db/*.fa.gz "$OUTDIR"/db/*.fasta.gz "$OUTDIR"/db/*.fna.gz \
         "$OUTDIR"/db/*.fa    "$OUTDIR"/db/*.fasta    "$OUTDIR"/db/*.fna; do
    [[ -e "$f" ]] || continue
    _have_genome=1
    # strip the FASTA extension (.gz, then .fa/.fasta/.fna) to name the .hdf5
    base="$f"; base="${base%.gz}"; base="${base%.fa}"; base="${base%.fasta}"; base="${base%.fna}"
    _sg straingst kmerize -k "$KMER_SIZE" -o "${base}.hdf5" "$f"
done
[[ "$_have_genome" -eq 1 ]] || {
    echo "ERROR: no genome FASTAs (.fa/.fasta/.fna, optionally .gz) found in $OUTDIR/db" >&2
    echo "       db/ currently contains:" >&2
    ls -A "$OUTDIR/db" >&2 || true
    echo "       If you used --genome-folder, make sure it holds FASTAs; if prepare_strainge_db" >&2
    echo "       produced nothing usable from a custom folder, re-run with --no-prepare-db." >&2
    exit 1
}

# --- Step 4: all-vs-all k-mer similarities ---
echo ">>> Computing all-vs-all similarities ..."
_sg straingst kmersim --all-vs-all -t "$THREADS" -S jaccard -S subset \
    "$OUTDIR"/db/*.hdf5 > "$OUTDIR/similarities.tsv"

# --- Step 5: cluster / de-replicate ---
echo ">>> Clustering references (-C $CLUSTER_JACCARD -c $CLUSTER_SUBSET) ..."
_sg straingst cluster -i "$OUTDIR/similarities.tsv" -d \
    -C "$CLUSTER_JACCARD" -c "$CLUSTER_SUBSET" \
    --clusters-out "$OUTDIR/clusters.tsv" \
    "$OUTDIR"/db/*.hdf5 > "$OUTDIR/references_to_keep.txt"

# --- Step 6: build the pan-genome database ---
echo ">>> Building pan-genome database ..."
_sg straingst createdb -f "$OUTDIR/references_to_keep.txt" -o "$OUTDIR/pan-genome-db.hdf5"

echo ""
echo ">>> Done."
echo "    Database  : $OUTDIR/pan-genome-db.hdf5   (-> run_straingst.sh --db)"
echo "    Genomes   : $OUTDIR/db                    (-> run_straingst.sh --ref-dir)"