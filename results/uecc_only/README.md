# Round-3 UeccDNA-only benchmark outputs

## Purpose

Reviewer 5 objected that the previous benchmark presentation mixed a common
UeccDNA detection task with CircleSeeker-specific MeccDNA/CeccDNA classification.
These tables separate the shared task using the existing benchmark source data.

## Source files

- `../benchmark_results.csv`
- `../length_stratified_results.csv`

## Analysis definition

- `human_U_10000` contains 10,000 simulated human UeccDNAs and no MeccDNA/CeccDNA.
- For this dataset, the `Overall` rows are therefore UeccDNA-only detection metrics
  with precision, recall and F1 available for every tool.
- For the mixed UMC datasets (`human_23000` and `ara_UMC_5200`), per-type UeccDNA
  recall is exported separately. Per-type precision/F1 should not be used as a
  cross-tool comparison for tools that do not type their output calls.
- `CReSIL_HiFi` is labelled in the output as an author-created HiFi adaptation.

## Key 30X common-task result

In the 30X human UeccDNA-only benchmark, the top-ranked method by mean F1 was
CircleSeeker (F1 99.1%, precision 99.6%,
recall 98.5%). The second row was CReSIL-HiFi (author-created HiFi adaptation)
(F1 94.4%).

CircleSeeker's 30X UeccDNA-only common-task values are precision
99.6%, recall 98.5%, and F1
99.1%. CReSIL-HiFi's corresponding values are precision
99.8%, recall 89.5%,
and F1 94.4%. Among external comparator tools excluding
CircleSeeker and the author-created CReSIL-HiFi adaptation, the highest 30X
UeccDNA-only F1 is CReSIL with F1
88.8%.

## Output files

- `human_simple_uecc_only_replicate_metrics.csv`: replicate-level common-task metrics.
- `human_simple_uecc_only_summary_by_depth_tool.csv`: mean/s.d. across three replicates.
- `human_simple_30x_uecc_only_summary.csv`: 30X summary ranked by F1.
- `mixed_umc_datasets_uecc_recall_replicate_metrics.csv`: UeccDNA recall rows in mixed UMC datasets.
- `mixed_umc_datasets_uecc_recall_summary_by_depth_tool.csv`: mixed-dataset UeccDNA recall means/s.d.
- `human_simple_uecc_only_length_stratified_replicates.csv`: Uecc-only length-bin recall replicates.
- `human_simple_uecc_only_length_stratified_recall_summary.csv`: Uecc-only length-bin recall means/s.d.

## Source-of-truth note

The tables in this directory were regenerated directly from the two CSV files
listed above. They are the intended source for the round-3 common-task summary.
