#!/usr/bin/env python3
import argparse
import sys
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Convert a sylph-tax merged table into a species-level MAGnet "
            "classification TSV using an accession->NCBI taxid map from TaxonKit."
        )
    )
    p.add_argument(
        "-s", "--sylph-merged", required=True,
        help="sylph-tax merge output (e.g. *_sylph_mpa_merged_abundance_table.tsv)"
    )
    p.add_argument(
        "--acc-taxmap", required=True,
        help="TSV with accession, species_name, ncbi_taxid (from TaxonKit name2taxid)"
    )
    p.add_argument(
        "--sample", required=True,
        help=(
            "Sample identifier to pull from the merged table. "
            "Can be the full column name or a substring like '10_1.chop.fastq.gz'."
        )
    )
    p.add_argument(
        "-o", "--out", required=True,
        help="Output MAGnet classification TSV (e.g. 10_1_magnet_species.tsv)"
    )
    p.add_argument(
        "--min-abundance", type=float, default=0.0,
        help="Minimum abundance (in Sylph units) to keep at genome level."
    )
    p.add_argument(
        "--scale-to-fraction", action="store_true",
        help=(
            "If set, divide abundances by 100 so values are in [0,1] "
            "for MAGnet's --min-abundance convention."
        )
    )
    return p.parse_args()


def main():
    args = parse_args()

    # Load Sylph merged abundance table
    df = pd.read_csv(args.sylph_merged, sep="\t")

    if "clade_name" not in df.columns:
        sys.stderr.write("[ERROR] Expected a 'clade_name' column in sylph merged table.\n")
        sys.exit(1)

    # Identify the sample column
    if args.sample in df.columns:
        sample_col = args.sample
    else:
        matches = [c for c in df.columns if args.sample in c]
        if not matches:
            sys.stderr.write(
                f"[ERROR] No columns in {args.sylph_merged} contain '{args.sample}'.\n"
            )
            sys.exit(1)
        if len(matches) > 1:
            sys.stderr.write(
                "[WARN] Multiple columns matched sample string; using first:\n"
                f"       {matches}\n"
                f"       -> {matches[0]}\n"
            )
        sample_col = matches[0]

    sys.stderr.write(f"[INFO] Using sample column: {sample_col}\n")

    # Keep only strain-level rows with t__GCF/GCA...
    df["accession"] = df["clade_name"].str.extract(r"t__([^|]+)$")[0]
    df = df.dropna(subset=["accession"])

    before = len(df)
    if args.min_abundance > 0:
        df = df[df[sample_col] > args.min_abundance]
    after = len(df)
    if args.min_abundance > 0:
        sys.stderr.write(
            f"[INFO] Filtered genomes by min_abundance {args.min_abundance}: "
            f"{before} -> {after} rows\n"
        )

    if df.empty:
        sys.stderr.write("[ERROR] No genome-level rows left after filtering.\n")
        sys.exit(1)

    # Load accession -> species_name -> ncbi_taxid map
    acc_map = pd.read_csv(
        args.acc_taxmap,
        sep="\t",
        header=None,
        names=["accession", "species_name", "ncbi_taxid"],
        dtype={"accession": str, "species_name": str, "ncbi_taxid": "Int64"},
    )

    acc_map = acc_map.dropna(subset=["ncbi_taxid"])
    if acc_map.empty:
        sys.stderr.write("[ERROR] acc-taxmap has no valid ncbi_taxid entries.\n")
        sys.exit(1)

    # Make sure each accession appears only once
    dup_mask = acc_map.duplicated(subset=["accession"])
    if dup_mask.any():
        n_dup = dup_mask.sum()
        sys.stderr.write(
            f"[WARN] {n_dup} duplicate accession rows in acc-taxmap; "
            "keeping the first occurrence per accession.\n"
        )
        acc_map = acc_map[~dup_mask]


    # Merge Sylph table with taxid map
    merged = df.merge(
        acc_map[["accession", "ncbi_taxid"]],
        on="accession",
        how="left",
        validate="m:1",
    )

    missing = merged["ncbi_taxid"].isna().sum()
    if missing:
        sys.stderr.write(
            f"[WARN] {missing} genome rows had no ncbi_taxid and will be dropped.\n"
        )

    merged = merged.dropna(subset=["ncbi_taxid"])

    if merged.empty:
        sys.stderr.write(
            "[ERROR] No rows remain after joining with accession->taxid map.\n"
        )
        sys.exit(1)

    # --- Collapse to species level (one row per taxid) ---
    species_abund = (
        merged
        .groupby("ncbi_taxid", as_index=False)[sample_col]
        .sum()
    )

    if args.scale_to_fraction:
        species_abund[sample_col] = species_abund[sample_col] / 100.0
        sys.stderr.write("[INFO] Scaled abundances by 1/100 to get fractions in [0,1].\n")

    out_df = species_abund.rename(
        columns={"ncbi_taxid": "taxid", sample_col: "abundance"}
    )

    out_df.to_csv(args.out, sep="\t", index=False)
    sys.stderr.write(
        f"[INFO] Wrote {len(out_df)} species-level rows to {args.out}\n"
        "[INFO] Columns: taxid (NCBI species taxid), abundance\n"
    )


if __name__ == "__main__":
    main()

