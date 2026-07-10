# CircleSeeker benchmark

This repository contains the simulation description, tool-run scripts,
evaluation code and metric tables used for the platform-matched CircleSeeker
benchmark. It also contains the round-3 analyses that separate the common
UeccDNA detection task from the UMC-aware stress tests and test sensitivity to
the study-derived PacBio HiFi sampling profile.

## Interpretation of the benchmark

The three scenarios answer different questions and should not be pooled into a
single generic ranking.

- **Human-Simple (`human_U_10000`)** contains 10,000 UeccDNAs and no MeccDNA
  or CeccDNA. This is the primary common head-to-head detection task because
  precision, recall and F1 have the same UeccDNA-only interpretation for all
  six tools.
- **Human-Complex (`human_23000`)** contains 20,000 UeccDNAs, 2,000 MeccDNAs
  and 1,000 CeccDNAs.
- **Arabidopsis (`ara_UMC_5200`)** contains 4,000 UeccDNAs, 1,000 MeccDNAs
  and 200 CeccDNAs.
- Human-Complex and Arabidopsis are **UMC-aware stress tests**. Their overall
  and per-type results assess recovery under multi-mapping and chimeric truth;
  they are not presented as the primary generic cross-tool comparison.
- MeccDNA and CeccDNA values for Human-Simple are not applicable, not zero,
  because those truth classes are absent from that scenario.

`CReSIL-HiFi` is an in-house HiFi adaptation of CReSIL created by the authors
for this benchmark. It is not an independently published comparator.

## Tools and native input platforms

| Tool | Archived benchmark version | Input |
|---|---|---|
| [CircleSeeker](https://github.com/leoxqy/CircleSeeker) | v1.1.2.dev0 | PacBio HiFi |
| [CReSIL-HiFi](https://github.com/YxZhang-XHCY/cresil-hifi) | v1.2.0+hifi; author-created adaptation | PacBio HiFi |
| [CReSIL](https://github.com/visanuwan/cresil) | v1.1.0 | Oxford Nanopore |
| [eccDNA_RCA_nanopore](https://github.com/icebert/eccDNA_RCA_nanopore) | commit `3f4b1dd` | Oxford Nanopore |
| [Circle-Map](https://github.com/iprada/Circle-Map) with [Circle-Map-cpp](https://github.com/BGI-Qingdao/Circle-Map-cpp) | v1.1.4 with v1.0.0 | Illumina paired-end |
| [ecc_finder](https://github.com/njaupan/ecc_finder) | v1.1.0 | Illumina paired-end |

Each tool received reads from its intended sequencing platform. Platform
matching is more realistic than feeding one read type to every caller, but it
also means differences reflect the complete platform-plus-caller workflow and
should not be interpreted as an isolated software-only comparison.

## Simulation design

The benchmark used the `ecc simulate` module of
[eccToolkit](https://github.com/YxZhang-XHCY/eccToolkit), following the general
simulation framework of Gao et al.
([Nature Communications, 2024](https://doi.org/10.1038/s41467-024-53496-8)).
Human truth was generated against T2T-CHM13v2.0 and plant truth against
ColCEN.

1. UeccDNA truth entries were sampled as single contiguous genomic intervals.
2. MeccDNA candidates were sampled from repeat-rich regions and retained when
   their simulated sequence mapped to multiple loci under minimap2-based
   criteria.
3. CeccDNA candidates were assembled from two to five non-collinear genomic
   fragments, including intra- and inter-chromosomal structures.
4. RCA-derived long reads contained variable concatemer copy numbers, random
   starting offsets and random linearisation breakpoints.
5. Circular reads were mixed 1:1 with simulated linear genomic-background
   reads.
6. Each scenario was generated at 10x, 30x and 50x template coverage with
   three independent replicates. Matched HiFi, Nanopore and Illumina read sets
   were generated for each condition.

The main HiFi simulation used the eccToolkit `HeLa_HiFi_2k` PBSIM2 sampling
profile, which was derived from 2,000 empirical HiFi reads. This affects the
simulated HiFi read-length, quality and sampling distributions; it does not
define the ground-truth circle coordinates or UMC composition.

### Scope of the simulator

The simulation represents circle structure, RCA concatemers, platform-specific
read errors and a linear-background challenge. It does not model every
biochemical feature of cell lysis, chromatin accessibility, circle extraction,
size-dependent recovery or phi29 amplification. Results therefore support
caller performance under the stated in-silico conditions, not absolute
experimental recovery efficiency.

For CeccDNA, the non-collinear truth fragments were joined into the circular
template before concatemer generation. These are truth-positive multi-fragment
circle reads. The archived benchmark did not enable a separate stochastic
inter-template phi29-chimera injection step, so it should not be interpreted as
modelling the full false-chimera spectrum of an experimental RCA reaction.

The commands used to generate the three truth populations were:

```bash
ecc simulate -r chm13v2.0.fa -o human_U_10000 \
  -u 10000 --replicates 3 --cov 10 --cov 30 --cov 50 -t 24

ecc simulate -r chm13v2.0.fa -o human_23000 \
  -u 20000 -m 2000 -c 1000 \
  --replicates 3 --cov 10 --cov 30 --cov 50 -t 24

ecc simulate -r ColCEN.fasta -o ara_UMC_5200 \
  -u 4000 -m 1000 -c 200 \
  --replicates 3 --cov 10 --cov 30 --cov 50 -t 24
```

## Evaluation

The evaluation uses detection groups rather than counting every output
fragment as an independent call.

The core evaluator uses the Python standard library. The synthetic-profile
utility requires Python 3.10 or later, NumPy and pandas; PBSIM2 and CircleSeeker
are additionally required for its `simulate` subcommand.

1. CeccDNA truth is evaluated first. A CeccDNA is recovered only when all of
   its fragments are present in one detected group with similarity at least
   0.90, where similarity is overlap divided by the shorter interval length.
2. Remaining detections are matched to UeccDNA and MeccDNA truth by at least
   0.90 reciprocal overlap.
3. Overall precision, recall and F1 are calculated from grouped detections.
4. Per-type recall is reported where that truth class exists. Per-type
   precision is reported only when the output contains sufficient structural
   or explicit type information; an absent value is not treated as zero.
5. Read-level eccDNA_RCA output is collapsed at 99% reciprocal overlap before
   evaluation; the residual redundancy ratio is retained in the result table.

Run the evaluator on a populated benchmark collection with:

```bash
python scripts/benchmark_analysis.py /path/to/benchmark_collect
```

The expected directory names are `human_U_10000`, `human_23000` and
`ara_UMC_5200`, each containing `rep1` to `rep3` and
`sequencing_10X`, `sequencing_30X` and `sequencing_50X` subdirectories.

## Common-task result

Mean performance across the three Human-Simple 30x replicates is shown below.
The full replicate-level data and all depths are in
`results/uecc_only/`.

| Tool | Precision | Recall | F1 |
|---|---:|---:|---:|
| CircleSeeker | 99.60% | 98.54% | 99.07% |
| CReSIL-HiFi (author-created adaptation) | 99.76% | 89.53% | 94.37% |
| CReSIL | 98.72% | 80.61% | 88.75% |
| Circle-Map | 98.28% | 78.77% | 87.45% |
| eccDNA_RCA_nanopore | 72.86% | 85.28% | 78.58% |
| ecc_finder | 94.01% | 35.60% | 51.64% |

Mixed-scenario overall and per-type values remain available in
`results/benchmark_results.csv` for the explicitly labelled UMC-aware stress
tests.

## Synthetic HiFi-profile sensitivity

To test whether the Human-Simple result depended on the study-derived
`HeLa_HiFi_2k` profile, CircleSeeker was rerun in a paired sensitivity analysis
with a fully synthetic, prespecified alternative profile. The alternative used
2,000 random sequences, a truncated log-normal length distribution (3-50 kb;
median parameter 12 kb; log-scale sigma 0.65; seed 20260710) and HiFi-range
base qualities. Truth subsets, coverage, background ratio, random seed and
CircleSeeker settings were held constant within each of three replicate pairs.

| Profile | Precision | Recall | F1 |
|---|---:|---:|---:|
| Built-in `HeLa_HiFi_2k` | 99.46% | 98.87% | 99.16% |
| Fully synthetic | 99.32% | 99.08% | 99.20% |

The paired synthetic-minus-built-in F1 difference was +0.035 percentage
points. This supports the narrow conclusion that the observed CircleSeeker
performance is not contingent on the study-derived read-length profile. It
does not make the simulator a complete model of experimental recovery bias.
The deterministic generator, selected truth coordinates and replicate-level
outputs are in `scripts/round3_hifi_profile_sensitivity.py` and
`results/profile_sensitivity/`.

## Repository contents

```text
CircleSeeker-benchmark/
|-- README.md
|-- scripts/
|   |-- benchmark_analysis.py
|   |-- round3_hifi_profile_sensitivity.py
|   |-- CircleMap_Enhanced.py
|   |-- collect_benchmark_results.sh
|   `-- sbatch_*.sh
`-- results/
    |-- benchmark_results.csv
    |-- benchmark_rca_stats.csv
    |-- length_stratified_results.csv
    |-- runtime_memory_stats.csv
    |-- uecc_only/
    `-- profile_sensitivity/
```

- `results/benchmark_results.csv` contains 540 rows: three scenarios, three
  depths, three replicates, six tools and the applicable overall/per-type
  evaluations.
- `results/length_stratified_results.csv` contains replicate-level recall in
  four truth-size bins.
- `results/runtime_memory_stats.csv` contains runtime and peak-memory records.
- `results/uecc_only/` contains the round-3 common-task summaries and
  mixed-scenario UeccDNA recall exports.
- `results/profile_sensitivity/` contains the paired synthetic-profile inputs,
  provenance records, read summaries and metrics.

## Data availability and reproducibility boundary

This repository includes the scripts and metric-level source data needed to
audit the reported benchmark summaries. The large simulated read files and
complete intermediate output directories from all six tools are not stored in
Git. Rerunning the benchmark from reads therefore requires regenerating the
datasets with the commands and versions above, then running the tool-specific
scripts. Cluster paths and scheduler settings in `sbatch_*.sh` must be adapted
to the local environment.
