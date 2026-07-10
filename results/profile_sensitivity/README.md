# Synthetic HiFi profile sensitivity analysis

## Question

Does CircleSeeker's Human-Simple performance depend on the study-derived
`HeLa_HiFi_2k` PBSIM sampling profile?

## Paired design

- Scenario: Human-Simple at nominal 30x template coverage.
- Replicates: three independently generated benchmark replicates.
- Per replicate: 2,000 UeccDNA templates and 2,000 linear-background
  templates sampled from the existing truth library.
- Paired controls: truth coordinates, template coverage, background ratio,
  simulation seed, CircleSeeker settings and total output bases were held
  constant within each replicate.
- Built-in profile: `HeLa_HiFi_2k`, generated from 2,000 empirical HiFi
  reads (mean length 13,425.6 bp; maximum 51,185 bp; mean accuracy 0.997883).
- Synthetic profile: 2,000 random DNA reads with lengths drawn from a
  prespecified truncated log-normal distribution (3-50 kb; median parameter
  12 kb; log-scale sigma 0.65; seed 20260710). Base qualities were generated
  in the same HiFi range as the built-in condition.
- Caller: CircleSeeker v1.1.2, public tag commit
  `b22217ace545dbb85c6b17de7b8e29c6d0998f08`.
- Evaluation: 90% reciprocal overlap against the selected UeccDNA truth.

This is a sensitivity analysis using a fully synthetic, prespecified
alternative distribution. It is independent of study reads, but it is not
presented as an independent empirical sequencing dataset.

## Output matching

The simulator generated nearly identical mean total output for the two paired
conditions: 881,063,638 bases with the built-in profile and 881,063,897 bases
with the synthetic profile. Mean output base qualities were Q38.86 and Q39.01,
respectively. The synthetic condition had a broader output length distribution
(mean length SD 8,687 bp versus 6,732 bp; N50 15,551 bp versus 13,260 bp).

## Results

| Profile | Precision | Recall | F1 |
|---|---:|---:|---:|
| Built-in `HeLa_HiFi_2k` | 99.46% | 98.87% | 99.16% |
| Fully synthetic | 99.32% | 99.08% | 99.20% |

The mean paired synthetic-minus-built-in F1 difference was +0.035 percentage
points. Length-stratified recall was also stable; differences were +0.319,
-0.138, +1.431 and +1.487 percentage points for <1 kb, 1-5 kb, 5-50 kb and
>=50 kb, respectively. The >=50-kb estimate is based on few truth molecules
and was variable across replicates.

The result supports the narrow conclusion that CircleSeeker's high
Human-Simple performance is not contingent on the study-derived read-length
profile. It does not convert the idealised simulator into a model of every
phi29 amplification or experimental recovery bias.

## Files

- `circleseeker_profile_sensitivity_replicates.csv`: replicate-level calls and
  metrics.
- `circleseeker_profile_sensitivity_summary.csv`: means and standard
  deviations.
- `circleseeker_profile_sensitivity_paired_differences.csv`: paired effect
  sizes.
- `hifi_output_profile_replicate_stats.csv` and
  `hifi_output_profile_summary.csv`: generated-read matching checks.
- `builtin_profile_fastq_stats.json` and
  `synthetic_lognormal_q39_profile_stats.json`: input profile summaries.
- `rep*_subset_truth.bed` and `rep*_subset_metadata.json`: selected truth
  coordinates and subset provenance.

Paths beginning with `benchmark_source/`, `raw_reads_not_in_repository/` or
`raw_outputs_not_in_repository/` are logical provenance paths. The large
source libraries, simulated FASTQ files and complete caller-output directories
are not stored in Git.

The deterministic generator, simulator wrapper, evaluator and aggregators are
implemented in `../../scripts/round3_hifi_profile_sensitivity.py`.
