"""
genome_list_input.py

Adds a "bring your own genomes" input mode to MAGnet.

Instead of starting from a taxonomic-classifier report (taxids -> NCBI Datasets
download -> ANI clustering), this module lets the user hand MAGnet an explicit
list of strain genomes (local FASTA files). It produces exactly the same
`downloaded_assemblies` DataFrame schema and the same on-disk layout
(`<output>/reference_genomes/{Assembly Accession ID}.fasta`) that the rest of the
MAGnet pipeline already consumes, so nothing downstream has to change:
fastANI clustering, merge_reference_fasta, get_seq2assembly_dict, alignment,
coverage, ANI and present/absence calling all keep working as-is.

Genome-list file format (TSV by default; CSV if the path ends in .csv).
A header row is required. Columns (case-insensitive):

    genome_path   (REQUIRED) path to a FASTA file, optionally .gz
    genome_id     (optional) unique id used as the "accession"/filename.
                  Defaults to the FASTA basename. Sanitized + de-duplicated.
    species       (optional) species / organism name for reporting
    strain        (optional) strain label for reporting
    taxid         (optional) NCBI taxid for reporting (not required)
    assembly_level(optional) e.g. "Complete Genome", "Contig". Default
                  "Complete Genome". Only affects representative choice when
                  clustering is enabled.
    abundance     (optional) carried through for reporting / pre-filtering

Only `genome_path` is mandatory. Everything else is best-effort metadata.
"""

import os
import gzip
import shutil

import pandas as pd
from Bio import SeqIO

# The column set that the rest of magnet.py expects on `downloaded_assemblies`
# before clustering. Matches what prepare_reference_genomes() returns.
METADATA_COLUMNS = [
    "Taxonomy ID",
    "Assembly Accession ID",
    "Source Database",
    "Is Representative",
    "Assembly Level",
    "Organism of Assembly",
    "Strain",
    "Total Length",
    "Downloaded",
    "Species",
]

# Map of accepted (lowercased) header names -> canonical key
_COLUMN_ALIASES = {
    "genome_path": "genome_path",
    "path": "genome_path",
    "fasta": "genome_path",
    "file": "genome_path",
    "genome_id": "genome_id",
    "id": "genome_id",
    "accession": "genome_id",
    "label": "genome_id",
    "name": "genome_id",
    "species": "species",
    "organism": "species",
    "strain": "strain",
    "taxid": "taxid",
    "tax_id": "taxid",
    "taxonomy_id": "taxid",
    "assembly_level": "assembly_level",
    "level": "assembly_level",
    "abundance": "abundance",
}


def _sanitize_id(raw_id):
    """Make an id safe for use as a filename and a fastANI matrix label."""
    safe = "".join(c if (c.isalnum() or c in "._-") else "_" for c in str(raw_id).strip())
    return safe or "genome"


def _open_maybe_gzip(path):
    """Open a (optionally gzipped) text FASTA for reading."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _normalize_columns(df):
    """Rename incoming columns to canonical keys; ignore unknown columns."""
    rename = {}
    for col in df.columns:
        key = _COLUMN_ALIASES.get(str(col).strip().lower())
        if key is not None:
            rename[col] = key
    return df.rename(columns=rename)


def parse_genome_list(genome_list_path):
    """Read the genome-list file into a normalized DataFrame with a genome_path column."""
    sep = "," if str(genome_list_path).lower().endswith(".csv") else "\t"
    df = pd.read_csv(genome_list_path, sep=sep, dtype=str, comment="#")
    df = _normalize_columns(df)

    if "genome_path" not in df.columns:
        raise ValueError(
            "Genome list must contain a 'genome_path' column "
            "(aliases: path, fasta, file). Got columns: %s" % list(df.columns)
        )

    df["genome_path"] = df["genome_path"].astype(str).str.strip()
    df = df[df["genome_path"].astype(bool)]  # drop blank rows
    if df.empty:
        raise ValueError("Genome list contained no usable rows.")
    return df.reset_index(drop=True)


def _write_reference_fasta(src_path, dst_path, accession, prefix_contigs=True):
    """
    Copy a genome FASTA into the reference_genomes directory, optionally
    prefixing every contig header with the accession id so that contig names
    are globally unique across the concatenated reference. Returns total length.
    """
    total_length = 0
    with _open_maybe_gzip(src_path) as in_f, open(dst_path, "w") as out_f:
        records = SeqIO.parse(in_f, "fasta")
        wrote_any = False
        for record in records:
            wrote_any = True
            total_length += len(record.seq)
            if prefix_contigs:
                # Keep original id readable but guarantee uniqueness.
                record.id = f"{accession}|{record.id}"
                record.name = record.id
                record.description = ""
            SeqIO.write(record, out_f, "fasta")
    if not wrote_any:
        raise ValueError(f"No FASTA records found in {src_path}")
    return total_length


def prepare_reference_genomes_from_list(genome_list_path, working_directory,
                                        prefix_contigs=True):
    """
    Drop-in replacement for the taxid -> NCBI download front-half of MAGnet.

    Returns a `downloaded_assemblies`-shaped DataFrame and writes each genome to
    <working_directory>/reference_genomes/{Assembly Accession ID}.fasta.
    """
    df = parse_genome_list(genome_list_path)

    reference_genome_path = os.path.join(working_directory, "reference_genomes")
    os.makedirs(reference_genome_path, exist_ok=True)

    seen_ids = {}
    rows = []

    for _, row in df.iterrows():
        src = row["genome_path"]
        if not os.path.isfile(src):
            raise FileNotFoundError(f"Genome file not found: {src}")

        # Derive a unique accession-like id.
        raw_id = row.get("genome_id")
        if raw_id is None or (isinstance(raw_id, float)) or str(raw_id).strip() == "" or str(raw_id) == "nan":
            base = os.path.basename(src)
            for ext in (".gz", ".fasta", ".fa", ".fna", ".fas"):
                if base.lower().endswith(ext):
                    base = base[: -len(ext)]
            raw_id = base
        accession = _sanitize_id(raw_id)

        # De-duplicate ids so two files never collide on the same filename.
        if accession in seen_ids:
            seen_ids[accession] += 1
            accession = f"{accession}_{seen_ids[accession]}"
        else:
            seen_ids[accession] = 0

        dst = os.path.join(reference_genome_path, f"{accession}.fasta")
        total_length = _write_reference_fasta(src, dst, accession,
                                              prefix_contigs=prefix_contigs)

        species = row.get("species")
        strain = row.get("strain")
        taxid = row.get("taxid")
        level = row.get("assembly_level")
        species = None if (species is None or str(species) == "nan") else species
        strain = None if (strain is None or str(strain) == "nan") else strain
        taxid = None if (taxid is None or str(taxid) == "nan") else taxid
        level = None if (level is None or str(level) == "nan") else level

        rows.append({
            "Taxonomy ID": taxid if taxid is not None else "",
            "Assembly Accession ID": accession,
            "Source Database": "USER_PROVIDED",
            "Is Representative": False,
            # Default to Complete Genome so that, if clustering is enabled, the
            # representative-selection heuristic behaves sensibly. Override per
            # row with an assembly_level column if you have better information.
            "Assembly Level": level if level is not None else "Complete Genome",
            "Organism of Assembly": species or strain or accession,
            "Strain": strain or "",
            "Total Length": float(total_length),
            "Downloaded": True,
            "Species": species or accession,
        })

    reference_metadata = pd.DataFrame(rows, columns=METADATA_COLUMNS)
    # Mirror magnet.py's convention: `downloaded_assemblies` is the rows that
    # were successfully obtained. Here that is all of them.
    return reference_metadata[reference_metadata["Downloaded"]].reset_index(drop=True)