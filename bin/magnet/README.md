# MAGNET (modified, vendored in Strainify)

This directory contains a **modified version of MAGNET**
(<https://github.com/treangenlab/magnet>), the whole-genome read-mapping
presence/absence caller from the Treangen Lab, integrated into Strainify for
strain-level reference filtering.

The standalone MAGNET command-line interface and its original README **do not
apply to this vendored copy**. Here, MAGNET is driven by Strainify's `filter-run`
subworkflow rather than invoked directly — see Strainify's main documentation for
usage.

## Attribution

Original tool by the Treangen Lab. Please cite:

> Sapoval, Nicolae, et al. "Lightweight taxonomic profiling of long-read metagenomic datasets with Lemur and Magnet." bioRxiv (2024).

## Main modifications
 
This copy differs from upstream MAGNET in the following ways:
 
- **"Bring your own genomes" input mode** (`utils/genome_list_input.py`): accepts
  an explicit list of local strain genome FASTAs instead of going through
  taxids → NCBI Datasets download, producing the same internal metadata schema
  the rest of the pipeline expects.
- **Shared-reference batch mode** (`utils/shared_reference.py`): lays out the
  per-genome FASTAs, the merged reference, and the `bwa` index **once** and
  reuses them across every sample in a batch (lock-protected and marker-guarded),
  instead of rebuilding the reference per sample.
- **Strain-level filtering integration**: present/absence calls are used as an
  upstream filter feeding Strainify, with the supplied genomes used as-is
  (no ANI clustering) in this mode.
- **Absolute evidence gate for false-positive removal**
  (`apply_absolute_gate`): the present/absent call is a breadth ÷ expected-breadth
  ratio that saturates to ~1.0 at tiny read counts, so a genome covered by a
  handful of reads (primary breadth ~1e-4) can score 1.0 and be called present.
  This gate runs *after* the present/absent call and requires a minimum absolute
  amount of primary evidence — primary breadth and/or primary mapped-read count —
  before a present call is trusted. It only ever downgrades present → absent, so
  lowering `--min-covscore` to recover borderline true positives no longer admits
  these near-empty false positives. **Off by default in standalone MAGNET, but
  Strainify enables it:** the prefilter
  (`subworkflows/local/prefilter_genomes.nf`, `MAGNET_FILTER`) passes
  `--min-breadth 0.005 --min-reads 0` by default, so a present call is downgraded
  to absent unless its primary breadth ≥ 0.005. Tunable from the pipeline via
  `--magnet_min_breadth` and `--magnet_min_reads`.
- **Primary-evidence rescue for low-abundance strains**
  (`apply_primary_evidence_rescue`): the mirror image of the gate, it only ever
  upgrades absent → present. The coverage-score call thresholds on the
  breadth/expected ratio, which normalizes away a genuinely-present but
  low-abundance strain — on a near-identical panel its reads cluster so observed
  breadth ≪ expected, the score drops below `--min-covscore`, and it's called
  absent even though it has thousands of *uniquely-mapped (primary)* reads (reads
  mapping to it and no other genome in the panel). That unique-read signal is
  competition-proof evidence of presence and cleanly separates true low-abundance
  strains (thousands of primary reads) from near-identical false positives (a
  handful). A genome called absent is upgraded to present when it has at least
  `min_reads` primary reads **or** at least `min_breadth` primary breadth, **and**
  a Consensus ANI of at least `min_ani`. It operates on the existing competitive
  (merged) alignment rather than per-genome realignment, so conserved-core
  cross-mapping doesn't inflate calls. **Off by default in standalone MAGNET, but
  Strainify enables it:** the same prefilter passes
  `--rescue-min-reads 100 --rescue-min-ani 0.97` by default, so a genome called
  absent is rescued to present when it has ≥ 100 primary reads **and** Consensus
  ANI ≥ 0.97. (`--rescue-min-breadth` is left unset, so only the read-count branch
  is active by default.) Tunable from the pipeline via `--magnet_rescue_min_reads`,
  `--magnet_rescue_min_ani`, and `--magnet_rescue_min_breadth`; the underlying
  present/absent call is set by `--magnet_min_covscore` (default 0.5).
- _(Add any other notable changes you made here.)_
## Strainify parameters
 
These pipeline parameters (set on `strainify filter-run`) drive MAGNET; each maps to the
underlying MAGNET CLI flag shown. Defaults are the values Strainify's prefilter applies.
 
| Parameter | Default | MAGNET flag | Effect |
|---|---|---|---|
| `--magnet_min_covscore` | `0.5` | `--min-covscore` | Coverage-score threshold for the base present/absent call |
| `--magnet_min_breadth` | `0.005` | `--min-breadth` | Gate: min primary breadth for a present call to stand |
| `--magnet_min_reads` | `0` | `--min-reads` | Gate: min primary reads for a present call to stand |
| `--magnet_rescue_min_reads` | `100` | `--rescue-min-reads` | Rescue: primary reads needed to recover an absent call |
| `--magnet_rescue_min_ani` | `0.97` | `--rescue-min-ani` | Rescue: minimum consensus ANI required to rescue |
| `--magnet_rescue_min_breadth` | _off_ | `--rescue-min-breadth` | Rescue: alternative primary-breadth trigger (unset by default) |
| `--magnet_min_mapq` | _off_ | `--min-mapq` | Minimum MAPQ for primary alignments (unset by default) |
| `--magnet_mode` | `illumina` | `-m` | Sequencing mode (`illumina` or `ont`) |
| `--magnet_threads` | `task.cpus` | `-t` | Threads for MAGNET |
| `--magnet_extra` | _none_ | _(appended verbatim)_ | Any extra MAGNET flags to pass through |
 
The absolute gate is active whenever `--magnet_min_breadth` or `--magnet_min_reads` is > 0; the
rescue is active whenever `--magnet_rescue_min_reads` or `--magnet_rescue_min_breadth` is > 0.
If too few genomes pass (Strainify needs at least 2 references), loosen these — e.g. lower
`--magnet_min_covscore` or `--magnet_rescue_min_reads`, or relax `--magnet_rescue_min_ani`.
