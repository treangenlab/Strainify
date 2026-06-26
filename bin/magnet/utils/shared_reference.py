# """
# shared_reference.py

# Build ONE shared reference layout from a genome list and reuse it across every
# sample in a batch, so the reference lay-out, FASTA merge, and (for illumina) the
# `bwa index` happen exactly once instead of once per sample.

# Layout written under <shared_dir>:
#     reference_genomes/{accession}.fasta   per-genome FASTAs (contigs prefixed)
#     merged_reference.fasta                concatenation of the above
#     merged_reference.fasta.{amb,ann,bwt,pac,sa}   bwa index (illumina only)
#     reference_metadata.csv                downloaded_assemblies + cluster columns
#     .shared_reference_complete            completion marker

# Per sample, `load_shared_reference` symlinks <working_dir>/reference_genomes ->
# <shared_dir>/reference_genomes (so get_seq2assembly_dict works unchanged) and
# returns the cached metadata DataFrame plus the prebuilt merged FASTA path.

# Shared mode assumes the candidate genomes are used as-is (i.e. --no-cluster);
# magnet enforces that. The build is marker-guarded and lock-protected so it is
# safe to call from every per-sample invocation (and tolerant of parallel runs).
# """

# import os
# import time
# import shutil
# import subprocess

# import pandas as pd

# MARKER = ".shared_reference_complete"
# LOCKDIR = ".lock"
# _LOCK_WAIT_SECONDS = 2
# _LOCK_TIMEOUT_SECONDS = 60 * 60  # 1h ceiling waiting for another builder


# def _truthy(series):
#     return series.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


# def _acquire_build(shared_dir):
#     """
#     Returns True if THIS caller should perform the build (holds the lock),
#     False if the reference is already built (marker present).
#     Blocks while another process is building.
#     """
#     marker = os.path.join(shared_dir, MARKER)
#     lock = os.path.join(shared_dir, LOCKDIR)
#     waited = 0
#     while True:
#         if os.path.exists(marker):
#             return False
#         try:
#             os.mkdir(lock)  # atomic
#         except FileExistsError:
#             time.sleep(_LOCK_WAIT_SECONDS)
#             waited += _LOCK_WAIT_SECONDS
#             if waited > _LOCK_TIMEOUT_SECONDS:
#                 # Stale lock guard: give up waiting and try to proceed.
#                 try:
#                     os.rmdir(lock)
#                 except OSError:
#                     pass
#             continue
#         # We hold the lock. Re-check the marker in case it appeared meanwhile.
#         if os.path.exists(marker):
#             _release_build(shared_dir)
#             return False
#         return True


# def _release_build(shared_dir):
#     lock = os.path.join(shared_dir, LOCKDIR)
#     if os.path.isdir(lock):
#         try:
#             os.rmdir(lock)
#         except OSError:
#             pass


# def build_shared_reference(genome_list, shared_dir, prefix_contigs=True,
#                            index_for_bwa=False, threads=1):
#     """
#     Build the shared reference layout once. No-op if already built.
#     Returns the metadata DataFrame (also cached to reference_metadata.csv).
#     """
#     shared_dir = str(shared_dir)
#     os.makedirs(shared_dir, exist_ok=True)
#     marker = os.path.join(shared_dir, MARKER)

#     if os.path.exists(marker):
#         return pd.read_csv(os.path.join(shared_dir, "reference_metadata.csv"))

#     if not _acquire_build(shared_dir):
#         # Built by someone else while we waited.
#         return pd.read_csv(os.path.join(shared_dir, "reference_metadata.csv"))

#     try:
#         # Lay out reference_genomes/{acc}.fasta and get the metadata DataFrame.
#         from utils.genome_list_input import prepare_reference_genomes_from_list
#         df = prepare_reference_genomes_from_list(
#             genome_list, shared_dir, prefix_contigs=prefix_contigs).copy()

#         # Shared mode == use-as-is: every genome is its own representative.
#         df["Cluster Representative"] = True
#         df["Cluster Members"] = df["Assembly Accession ID"].astype(str)

#         # Concatenate per-genome FASTAs into one merged reference (deterministic
#         # order = DataFrame order). Contig headers are already globally unique.
#         ref_genomes = os.path.join(shared_dir, "reference_genomes")
#         merged = os.path.join(shared_dir, "merged_reference.fasta")
#         with open(merged, "w") as out_f:
#             for acc in df["Assembly Accession ID"]:
#                 part = os.path.join(ref_genomes, f"{acc}.fasta")
#                 with open(part) as in_f:
#                     shutil.copyfileobj(in_f, out_f)

#         # Build the bwa index ONCE (only needed for the illumina/bwa path).
#         if index_for_bwa:
#             subprocess.run(
#                 ["bwa", "index", merged],
#                 stdout=open(os.path.join(shared_dir, "bwa_index.log"), "a"),
#                 stderr=open(os.path.join(shared_dir, "bwa_index.err"), "a"),
#                 check=True,
#             )

#         df.to_csv(os.path.join(shared_dir, "reference_metadata.csv"), index=False)

#         # Mark complete last, so a partial build is never treated as done.
#         with open(marker, "w") as mf:
#             mf.write("ok\n")
#         return df
#     finally:
#         _release_build(shared_dir)


# def load_shared_reference(shared_dir, working_directory):
#     """
#     Reuse a prebuilt shared reference for one sample.

#     - Symlinks <working_directory>/reference_genomes -> <shared_dir>/reference_genomes
#       (so get_seq2assembly_dict reads the shared per-genome FASTAs unchanged).
#     - Returns (downloaded_assemblies_df, merged_reference_fasta_path).
#     """
#     shared_dir = str(shared_dir)
#     working_directory = str(working_directory)

#     meta_csv = os.path.join(shared_dir, "reference_metadata.csv")
#     df = pd.read_csv(meta_csv)
#     # CSV round-trips bools to strings; coerce the columns magnet relies on.
#     if "Cluster Representative" in df.columns:
#         df["Cluster Representative"] = _truthy(df["Cluster Representative"])
#     if "Downloaded" in df.columns:
#         df["Downloaded"] = _truthy(df["Downloaded"])

#     os.makedirs(working_directory, exist_ok=True)
#     link = os.path.join(working_directory, "reference_genomes")
#     target = os.path.abspath(os.path.join(shared_dir, "reference_genomes"))
#     # os.path.lexists is true for broken symlinks too, so we never clobber.
#     if not os.path.lexists(link):
#         os.symlink(target, link)

#     merged = os.path.join(shared_dir, "merged_reference.fasta")
#     return df, merged




"""
shared_reference.py

Build ONE shared reference layout from a genome list and reuse it across every
sample in a batch, so the reference lay-out, FASTA merge, and (for illumina) the
`bwa index` happen exactly once instead of once per sample.

Layout written under <shared_dir>:
    reference_genomes/{accession}.fasta   per-genome FASTAs (contigs prefixed)
    merged_reference.fasta                concatenation of the above
    merged_reference.fasta.{amb,ann,bwt,pac,sa}   bwa index (illumina only)
    reference_metadata.csv                downloaded_assemblies + cluster columns
    .shared_reference_complete            completion marker

Per sample, `load_shared_reference` symlinks <working_dir>/reference_genomes ->
<shared_dir>/reference_genomes (so get_seq2assembly_dict works unchanged) and
returns the cached metadata DataFrame plus the prebuilt merged FASTA path.

Shared mode assumes the candidate genomes are used as-is (i.e. --no-cluster);
magnet enforces that. The build is marker-guarded and lock-protected so it is
safe to call from every per-sample invocation (and tolerant of parallel runs).
"""

import os
import time
import shutil
import subprocess

import pandas as pd

MARKER = ".shared_reference_complete"
LOCKDIR = ".lock"
_LOCK_WAIT_SECONDS = 2
_LOCK_TIMEOUT_SECONDS = 60 * 60  # 1h ceiling waiting for another builder


def _truthy(series):
    return series.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def _acquire_build(shared_dir):
    """
    Returns True if THIS caller should perform the build (holds the lock),
    False if the reference is already built (marker present).
    Blocks while another process is building.
    """
    marker = os.path.join(shared_dir, MARKER)
    lock = os.path.join(shared_dir, LOCKDIR)
    waited = 0
    while True:
        if os.path.exists(marker):
            return False
        try:
            os.mkdir(lock)  # atomic
        except FileExistsError:
            time.sleep(_LOCK_WAIT_SECONDS)
            waited += _LOCK_WAIT_SECONDS
            if waited > _LOCK_TIMEOUT_SECONDS:
                # Stale lock guard: give up waiting and try to proceed.
                try:
                    os.rmdir(lock)
                except OSError:
                    pass
            continue
        # We hold the lock. Re-check the marker in case it appeared meanwhile.
        if os.path.exists(marker):
            _release_build(shared_dir)
            return False
        return True


def _release_build(shared_dir):
    lock = os.path.join(shared_dir, LOCKDIR)
    if os.path.isdir(lock):
        try:
            os.rmdir(lock)
        except OSError:
            pass


def build_shared_reference(genome_list, shared_dir, prefix_contigs=True,
                           index_for_bwa=False, threads=1):
    """
    Build the shared reference layout once. No-op if already built.
    Returns the metadata DataFrame (also cached to reference_metadata.csv).
    """
    shared_dir = str(shared_dir)
    os.makedirs(shared_dir, exist_ok=True)
    marker = os.path.join(shared_dir, MARKER)

    if os.path.exists(marker):
        return pd.read_csv(os.path.join(shared_dir, "reference_metadata.csv"))

    if not _acquire_build(shared_dir):
        # Built by someone else while we waited.
        return pd.read_csv(os.path.join(shared_dir, "reference_metadata.csv"))

    try:
        # Lay out reference_genomes/{acc}.fasta and get the metadata DataFrame.
        from utils.genome_list_input import prepare_reference_genomes_from_list
        df = prepare_reference_genomes_from_list(
            genome_list, shared_dir, prefix_contigs=prefix_contigs).copy()

        # Shared mode == use-as-is: every genome is its own representative.
        df["Cluster Representative"] = True
        df["Cluster Members"] = df["Assembly Accession ID"].astype(str)

        # Concatenate per-genome FASTAs into one merged reference (deterministic
        # order = DataFrame order). Contig headers are already globally unique.
        ref_genomes = os.path.join(shared_dir, "reference_genomes")
        merged = os.path.join(shared_dir, "merged_reference.fasta")
        with open(merged, "w") as out_f:
            for acc in df["Assembly Accession ID"]:
                part = os.path.join(ref_genomes, f"{acc}.fasta")
                with open(part) as in_f:
                    shutil.copyfileobj(in_f, out_f)

        # Build the bwa index ONCE (only needed for the illumina/bwa path).
        # Drop a completion sentinel after success so run_bwa_index_once trusts
        # this index (and rebuilds it if it is ever found partial/stale).
        if index_for_bwa:
            sentinel = merged + ".bwa_index_complete"
            if os.path.exists(sentinel):
                os.remove(sentinel)
            subprocess.run(
                ["bwa", "index", merged],
                stdout=open(os.path.join(shared_dir, "bwa_index.log"), "a"),
                stderr=open(os.path.join(shared_dir, "bwa_index.err"), "a"),
                check=True,
            )
            with open(sentinel, "w") as sf:
                sf.write("ok\n")

        df.to_csv(os.path.join(shared_dir, "reference_metadata.csv"), index=False)

        # Mark complete last, so a partial build is never treated as done.
        with open(marker, "w") as mf:
            mf.write("ok\n")
        return df
    finally:
        _release_build(shared_dir)


def load_shared_reference(shared_dir, working_directory):
    """
    Reuse a prebuilt shared reference for one sample.

    - Symlinks <working_directory>/reference_genomes -> <shared_dir>/reference_genomes
      (so get_seq2assembly_dict reads the shared per-genome FASTAs unchanged).
    - Returns (downloaded_assemblies_df, merged_reference_fasta_path).
    """
    shared_dir = str(shared_dir)
    working_directory = str(working_directory)

    meta_csv = os.path.join(shared_dir, "reference_metadata.csv")
    df = pd.read_csv(meta_csv)
    # CSV round-trips bools to strings; coerce the columns magnet relies on.
    if "Cluster Representative" in df.columns:
        df["Cluster Representative"] = _truthy(df["Cluster Representative"])
    if "Downloaded" in df.columns:
        df["Downloaded"] = _truthy(df["Downloaded"])

    os.makedirs(working_directory, exist_ok=True)
    link = os.path.join(working_directory, "reference_genomes")
    target = os.path.abspath(os.path.join(shared_dir, "reference_genomes"))
    # os.path.lexists is true for broken symlinks too, so we never clobber.
    if not os.path.lexists(link):
        os.symlink(target, link)

    merged = os.path.join(shared_dir, "merged_reference.fasta")
    return df, merged