# CircleSeeker-benchmark

Comprehensive benchmarking of eccDNA detection tools on simulated long-read and short-read sequencing data, covering three types of eccDNA: unique-mapping (Uecc), multi-mapping (Mecc), and chimeric (Cecc).

## Overview

This repository contains the benchmark framework used to evaluate four eccDNA detection tools:

| Tool | Version | Input Data | Read Type |
|------|---------|-----------|-----------|
| [CircleSeeker](https://github.com/YaoXinZH/CircleSeeker) | v1.0 | HiFi long reads | PacBio HiFi |
| [CircleMap Enhanced](https://github.com/iprada/Circle-Map) | v1.1.4 | NGS paired-end reads | Illumina |
| [CReSIL](https://github.com/visanuwan/cresil) | v1.1 | ONT long reads | Oxford Nanopore |
| [eccDNA_RCA_nanopore](https://github.com/yourrepo/eccDNA_RCA_nanopore) | latest | ONT long reads | Oxford Nanopore |

## Simulated Datasets

Benchmark data was generated on two T2T reference genomes with three datasets of varying complexity:

### Reference Genomes

| Genome | Species | Description |
|--------|---------|-------------|
| Human T2T (CHM13v2.0) | *Homo sapiens* | T2T human reference genome |
| ColCEN | *Arabidopsis thaliana* | T2T Arabidopsis genome with centromere regions |

### Dataset Composition

| Dataset | Genome | Uecc | Mecc | Cecc | Total |
|---------|--------|------|------|------|-------|
| human_10000_simple | Human T2T | 10,000 | - | - | 10,000 |
| human_23000 | Human T2T | 20,000 | 2,000 | 1,000 | 23,000 |
| ColCEN_5200 | Arabidopsis T2T | 4,000 | 1,000 | 200 | 5,200 |

### eccDNA Types

- **Uecc (Unique eccDNA)**: eccDNA derived from unique genomic regions with a single mapping site
- **Mecc (Multi-mapping eccDNA)**: eccDNA originating from repetitive regions with multiple mapping sites in the genome
- **Cecc (Chimeric eccDNA)**: eccDNA composed of multiple non-contiguous genomic fragments joined together

### Sequencing Simulation

Each dataset was simulated at three sequencing depths with three biological replicates:

- **Depths**: 10X, 30X, 50X
- **Replicates**: 3 per condition
- **Total samples**: 3 datasets x 3 depths x 3 replicates = **27 samples**

Each sample includes matched simulated reads for all applicable sequencing platforms (HiFi, ONT, Illumina paired-end).

## Evaluation Methodology

### Unified Matching Strategy

All tools are evaluated using a single unified matching pass to avoid double-counting or misattribution of detections:

1. **Phase 1 — Cecc (strict matching)**: Each detected eccDNA group is first tested against chimeric truth entries. A Cecc truth is matched only if **ALL** of its fragments are found in the detected group (similarity >= 90%, defined as overlap / min(length1, length2)).

2. **Phase 2 — Uecc/Mecc (reciprocal overlap)**: Remaining unmatched detected groups are matched against Uecc and Mecc truth entries using reciprocal overlap >= 90% (overlap must cover >= 90% of both the truth fragment and the detected region).

### Metrics

- **Overall**: Precision, Recall, and F1 are computed across all eccDNA types combined. This is the primary evaluation metric.
- **Per-type (Uecc/Mecc/Cecc)**: Only **Recall** is reported, as Precision cannot be meaningfully attributed to individual types when tools do not classify their detections by type (unmatched detections cannot be assigned to a specific type).

### Detection Grouping

Detected regions are grouped by eccDNA ID (name). Each group counts as **one detection**, regardless of how many individual regions it contains. This correctly handles:
- **Mecc**: Multiple mapping sites for the same eccDNA count as one detection
- **Cecc**: Multiple fragments of the same chimeric eccDNA count as one detection

### eccDNA_RCA Deduplication

eccDNA_RCA_nanopore outputs read-level results with substantial redundancy. Before evaluation, its output is deduplicated using **99% reciprocal overlap** to collapse near-identical detections. Even after this processing, significant redundancy remains (1.7x at 10X, 3.2x at 30X, 4.3x at 50X).

## Results

### Overall Performance (3-replicate average)

#### human_10000_simple (10,000 Uecc)

| Depth | Tool | Precision | Recall | F1 |
|-------|------|-----------|--------|------|
| 10X | CircleSeeker | **99.7%** | **93.5%** | **96.5%** |
| 10X | CircleMap | 99.7% | 25.3% | 36.6% |
| 10X | CReSIL | 94.3% | 51.9% | 66.9% |
| 10X | eccDNA_RCA | 87.4% | 73.1% | 79.6% |
| 30X | CircleSeeker | **99.6%** | **98.6%** | **99.1%** |
| 30X | CircleMap | 98.5% | 78.7% | 87.5% |
| 30X | CReSIL | 98.7% | 80.5% | 88.7% |
| 30X | eccDNA_RCA | 72.9% | 85.4% | 78.7% |
| 50X | CircleSeeker | **99.4%** | **99.1%** | **99.2%** |
| 50X | CircleMap | 97.4% | 84.4% | 90.5% |
| 50X | CReSIL | 99.0% | 85.8% | 91.9% |
| 50X | eccDNA_RCA | 63.7% | 88.4% | 74.0% |

#### human_23000 (20,000 Uecc + 2,000 Mecc + 1,000 Cecc)

| Depth | Tool | Precision | Recall | F1 |
|-------|------|-----------|--------|------|
| 10X | CircleSeeker | **98.2%** | **91.0%** | **94.4%** |
| 10X | CircleMap | 94.4% | 46.1% | 62.0% |
| 10X | CReSIL | 89.8% | 43.2% | 58.4% |
| 10X | eccDNA_RCA | 80.9% | 63.3% | 71.1% |
| 30X | CircleSeeker | **98.1%** | **96.4%** | **97.2%** |
| 30X | CircleMap | 93.6% | 68.5% | 79.1% |
| 30X | CReSIL | 95.3% | 67.6% | 79.0% |
| 30X | eccDNA_RCA | 65.2% | 74.2% | 69.4% |
| 50X | CircleSeeker | **97.9%** | **97.1%** | **97.5%** |
| 50X | CircleMap | 92.8% | 74.0% | 82.3% |
| 50X | CReSIL | 96.4% | 72.4% | 82.7% |
| 50X | eccDNA_RCA | 56.5% | 76.5% | 65.0% |

#### ColCEN_5200 (4,000 Uecc + 1,000 Mecc + 200 Cecc)

| Depth | Tool | Precision | Recall | F1 |
|-------|------|-----------|--------|------|
| 10X | CircleSeeker | **98.9%** | **90.8%** | **94.7%** |
| 10X | CircleMap | 88.2% | 19.5% | 29.6% |
| 10X | CReSIL | 91.6% | 29.1% | 44.2% |
| 10X | eccDNA_RCA | 82.7% | 57.3% | 67.7% |
| 30X | CircleSeeker | **98.9%** | **96.1%** | **97.5%** |
| 30X | CircleMap | 94.8% | 60.4% | 73.8% |
| 30X | CReSIL | 96.1% | 43.5% | 59.8% |
| 30X | eccDNA_RCA | 66.6% | 66.9% | 66.7% |
| 50X | CircleSeeker | **98.8%** | **97.0%** | **97.9%** |
| 50X | CircleMap | 94.2% | 64.5% | 76.5% |
| 50X | CReSIL | 96.8% | 46.9% | 63.2% |
| 50X | eccDNA_RCA | 57.4% | 69.3% | 62.8% |

### Per-type Recall (3-replicate average)

#### Uecc Recall

| Dataset | Depth | CircleSeeker | CircleMap | CReSIL | eccDNA_RCA |
|---------|-------|-------------|-----------|--------|------------|
| human_10000_simple | 10X | **93.5%** | 25.3% | 51.9% | 73.1% |
| human_10000_simple | 30X | **98.6%** | 78.7% | 80.5% | 85.4% |
| human_10000_simple | 50X | **99.1%** | 84.4% | 85.8% | 88.4% |
| human_23000 | 10X | **93.1%** | 53.0% | 49.0% | 72.8% |
| human_23000 | 30X | **98.3%** | 78.5% | 76.1% | 85.3% |
| human_23000 | 50X | **98.8%** | 84.2% | 81.0% | 88.0% |
| ColCEN_5200 | 10X | **91.9%** | 25.3% | 37.2% | 74.5% |
| ColCEN_5200 | 30X | **96.9%** | 77.4% | 55.1% | 86.9% |
| ColCEN_5200 | 50X | **97.5%** | 82.0% | 59.0% | 90.1% |

#### Mecc Recall

| Dataset | Depth | CircleSeeker | CircleMap | CReSIL | eccDNA_RCA |
|---------|-------|-------------|-----------|--------|------------|
| human_23000 | 10X | **77.2%** | 0.4% | 1.6% | 0.0% |
| human_23000 | 30X | **82.2%** | 3.0% | 4.5% | 0.0% |
| human_23000 | 50X | **83.7%** | 9.2% | 6.4% | 0.0% |
| ColCEN_5200 | 10X | **90.1%** | 0.2% | 1.8% | 0.0% |
| ColCEN_5200 | 30X | **95.2%** | 4.5% | 4.5% | 0.0% |
| ColCEN_5200 | 50X | **96.5%** | 7.2% | 6.0% | 0.0% |

#### Cecc Recall

| Dataset | Depth | CircleSeeker | CircleMap | CReSIL | eccDNA_RCA |
|---------|-------|-------------|-----------|--------|------------|
| human_23000 | 10X | **76.9%** | 0.0% | 10.8% | 0.0% |
| human_23000 | 30X | **87.7%** | 0.0% | 22.6% | 0.0% |
| human_23000 | 50X | **90.8%** | 0.0% | 33.2% | 0.0% |
| ColCEN_5200 | 10X | **72.7%** | 0.0% | 4.0% | 0.0% |
| ColCEN_5200 | 30X | **85.0%** | 0.0% | 6.8% | 0.0% |
| ColCEN_5200 | 50X | **89.7%** | 0.0% | 9.5% | 0.0% |

### eccDNA_RCA Redundancy (after 99% deduplication)

| Dataset | 10X | 30X | 50X |
|---------|-----|-----|-----|
| human_10000_simple | 1.68x | 3.27x | 4.50x |
| human_23000 | 1.66x | 3.14x | 4.31x |
| ColCEN_5200 | 1.69x | 3.16x | 4.31x |

## Repository Structure

```
CircleSeeker-benchmark/
├── README.md
├── .gitignore
├── scripts/
│   ├── benchmark_analysis.py             # Evaluation script (unified matching)
│   ├── CircleMap_Enhanced.py             # CircleMap Enhanced implementation
│   ├── collect_benchmark_results.sh      # Collect tool outputs into unified directory
│   ├── sbatch_CircleMap_Enhanced.sh     # SLURM job: CircleMap Enhanced
│   ├── sbatch_CircleSeeker.sh           # SLURM job: CircleSeeker
│   ├── sbatch_CReSIL.sh                # SLURM job: CReSIL
│   └── sbatch_eccDNA_RCA_nanopore.sh   # SLURM job: eccDNA_RCA
└── results/
    ├── benchmark_results.csv             # Full benchmark results (all conditions)
    └── benchmark_rca_stats.csv           # eccDNA_RCA redundancy statistics
```

### Scripts

- **`benchmark_analysis.py`**: Core evaluation script. Parses truth BED files and tool outputs, performs unified matching (Cecc strict → Uecc/Mecc reciprocal overlap), and outputs per-sample and per-type metrics.
- **`collect_benchmark_results.sh`**: Aggregates outputs from all tools into the `benchmark_collect/` directory structure used by the analysis script.
- **`sbatch_*.sh`**: SLURM array job scripts for HPC cluster execution. Processes all 27 samples in parallel. Modify `BASE_DIR` to match your data directory before use.
- **`CircleMap_Enhanced.py`**: Enhanced CircleMap pipeline integrating fastp, BWA, samtools, and Circle-Map with automated filtering.

### Results CSV Format

**benchmark_results.csv** columns:
- `Genome`, `Rep`, `Depth`, `Tool`, `EvalType`: Sample and evaluation identifiers
- `Truth`: Number of ground-truth eccDNA entries
- `Detected`: Number of detected eccDNA groups (Overall only)
- `TP`, `FP`, `FN`: True positives, false positives, false negatives (Overall only)
- `Precision`, `Recall`, `F1`: Performance metrics (Precision/F1 for Overall only; Recall for all types)
- `Redundancy`: Detection redundancy ratio (eccDNA_RCA only)

## Reproducing the Benchmark

### Prerequisites

- Python 3.7+
- conda/mamba with bioconda channel
- Tool-specific environments (see individual run scripts for conda environment setup)

### Step 1: Download Benchmark Data

Download the simulated benchmark dataset from figshare:

```bash
# Download and extract benchmark data
# DOI: [TO BE ADDED]
tar -xzf benchmark_collect.tar.gz
```

The extracted `benchmark_collect/` directory contains:
```
benchmark_collect/
├── ColCEN_5200/
│   ├── rep1/
│   │   ├── sequencing_10X/
│   │   │   ├── truth_all.bed          # Ground truth
│   │   │   ├── CircleSeeker_summary.csv
│   │   │   ├── CircleMap_filtered.bed
│   │   │   ├── CReSIL_eccDNA_final.txt
│   │   │   └── eccDNA_RCA_info.tsv
│   │   ├── sequencing_30X/
│   │   └── sequencing_50X/
│   ├── rep2/
│   └── rep3/
├── human_10000_simple/
└── human_23000/
```

### Step 2: Run Evaluation

```bash
python scripts/benchmark_analysis.py benchmark_collect
```

This will output:
- Console summary of all results
- `results/benchmark_results.csv` — full benchmark results
- `results/benchmark_rca_stats.csv` — eccDNA_RCA redundancy statistics

### Step 3 (Optional): Rerun Tools from Raw Data

To rerun individual tools from simulated sequencing reads (not included in `benchmark_collect`), use the corresponding `run_*.sh` or `sbatch_*.sh` scripts. Refer to each script's header for required conda environments and dependencies.

## Data Availability

- **Benchmark results**: Included in this repository (`results/`)
- **Tool outputs and ground truth**: Available on figshare (DOI: [TO BE ADDED])

## Citation

If you use this benchmark in your research, please cite:

```
[TO BE ADDED]
```

## License

[TO BE ADDED]
