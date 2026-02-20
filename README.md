# CircleSeeker-benchmark

Comprehensive benchmarking of eccDNA detection tools on simulated long-read and short-read sequencing data, covering three types of eccDNA: unique-mapping (Uecc), multi-mapping (Mecc), and chimeric (Cecc).

## Overview

This repository contains the benchmark framework used to evaluate six eccDNA detection tools across three sequencing platforms:

| Tool | Version | Input Data | Read Type |
|------|---------|-----------|-----------|
| [CircleSeeker](https://github.com/YaoXinZH/CircleSeeker) | v1.0 | HiFi long reads | PacBio HiFi |
| [CReSIL-HiFi](https://github.com/YxZhang-XHCY/cresil-hifi) | v1.2.0+hifi | HiFi long reads | PacBio HiFi |
| [CReSIL](https://github.com/visanuwan/cresil) | v1.1.0 | ONT long reads | Oxford Nanopore |
| [eccDNA_RCA_nanopore](https://github.com/icebert/eccDNA_RCA_nanopore) | commit [3f4b1dd](https://github.com/icebert/eccDNA_RCA_nanopore/commit/3f4b1dd) | ONT long reads | Oxford Nanopore |
| [Circle-Map](https://github.com/iprada/Circle-Map) + [Circle-Map-cpp](https://github.com/BGI-Qingdao/Circle-Map-cpp) | v1.1.4 + v1.0.0 | NGS paired-end reads | Illumina |
| [ecc_finder](https://github.com/njaupan/ecc_finder) | v1.1.0 | NGS paired-end reads | Illumina |

### Excluded Tool

[CIDER-Seq2](https://github.com/devang-mehta/ciderseq2) (v2.0) was initially considered but excluded from the benchmark. Its MUSCLE-based de-concatenation (O(n·L²) per read for read length L) and exhaustive BLAST with O(n²) hit reordering are intractable for large eukaryotic genomes. Even after replacing the DeConcat step with TideHunter, the downstream BLAST-based eccDNA detection remained infeasible for human genome-scale data (>300K reads, 3.1 Gb reference).

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

Since existing eccDNA detection tools output only genomic coordinate intervals without eccDNA type classification (unique / multi-mapping / chimeric), we designed a unified two-phase evaluation framework for fair cross-tool comparison:

1. **Phase 1 — Cecc (strict matching)**: Each detected eccDNA group is first tested against chimeric truth entries. A Cecc truth is matched only if **ALL** of its fragments are found in the detected group (similarity >= 90%, defined as overlap / min(length1, length2)).

2. **Phase 2 — Uecc/Mecc (reciprocal overlap)**: Remaining unmatched detected groups are matched against Uecc and Mecc truth entries using reciprocal overlap >= 90% (overlap must cover >= 90% of both the truth fragment and the detected region).

### Metrics

- **Overall**: Precision, Recall, and F1 are computed across all eccDNA types combined. This is the primary cross-tool comparison metric.
- **Per-type Recall (Uecc/Mecc/Cecc)**: Leveraging the type labels in simulated ground truth, per-type Recall reveals how well each tool detects eccDNA of different structural complexities. Per-type Precision is **not** reported for most tools, because their outputs do not carry type labels — false positive detections cannot be attributed to a specific eccDNA type.
- **Per-type Precision (selected tools)**: Per-type Precision is reported for tools whose outputs carry structural or explicit type information:
  - **CircleSeeker**: Full per-type Precision for Uecc, Mecc, and Cecc — the only tool that explicitly classifies each detection by eccDNA type.
  - **CReSIL / CReSIL-HiFi**: Cecc Precision only — multi-fragment detection groups (multiple regions under one eccDNA ID) are structurally identifiable as chimeric eccDNA candidates.
  - **Other tools**: No per-type Precision (outputs are untyped single-region detections).

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
| 10X | CReSIL-HiFi | 98.4% | 66.1% | 79.1% |
| 10X | CReSIL | 94.3% | 51.9% | 66.9% |
| 10X | eccDNA_RCA | 87.4% | 73.1% | 79.6% |
| 10X | CircleMap | 99.7% | 25.3% | 36.6% |
| 10X | ecc_finder | 95.7% | 13.9% | 23.3% |
| 30X | CircleSeeker | **99.6%** | **98.6%** | **99.1%** |
| 30X | CReSIL-HiFi | 99.7% | 89.5% | 94.3% |
| 30X | CReSIL | 98.7% | 80.5% | 88.7% |
| 30X | eccDNA_RCA | 72.9% | 85.4% | 78.7% |
| 30X | CircleMap | 98.5% | 78.7% | 87.5% |
| 30X | ecc_finder | 94.0% | 35.6% | 51.7% |
| 50X | CircleSeeker | **99.4%** | **99.1%** | **99.2%** |
| 50X | CReSIL-HiFi | 99.8% | 91.8% | 95.6% |
| 50X | CReSIL | 99.0% | 85.8% | 91.9% |
| 50X | eccDNA_RCA | 63.7% | 88.4% | 74.0% |
| 50X | CircleMap | 97.4% | 84.4% | 90.5% |
| 50X | ecc_finder | 94.5% | 42.2% | 58.3% |

#### human_23000 (20,000 Uecc + 2,000 Mecc + 1,000 Cecc)

| Depth | Tool | Precision | Recall | F1 |
|-------|------|-----------|--------|------|
| 10X | CircleSeeker | **98.2%** | **91.0%** | **94.4%** |
| 10X | CReSIL-HiFi | 94.0% | 55.2% | 69.5% |
| 10X | CReSIL | 89.8% | 43.2% | 58.4% |
| 10X | eccDNA_RCA | 80.9% | 63.3% | 71.1% |
| 10X | CircleMap | 94.4% | 46.1% | 62.0% |
| 10X | ecc_finder | 71.8% | 8.1% | 14.6% |
| 30X | CircleSeeker | **98.1%** | **96.4%** | **97.2%** |
| 30X | CReSIL-HiFi | 95.6% | 74.7% | 83.9% |
| 30X | CReSIL | 95.3% | 67.6% | 79.1% |
| 30X | eccDNA_RCA | 65.2% | 74.2% | 69.4% |
| 30X | CircleMap | 93.6% | 68.5% | 79.1% |
| 30X | ecc_finder | 68.6% | 10.0% | 17.4% |
| 50X | CircleSeeker | **97.9%** | **97.1%** | **97.5%** |
| 50X | CReSIL-HiFi | 95.7% | 76.9% | 85.3% |
| 50X | CReSIL | 96.4% | 72.4% | 82.7% |
| 50X | eccDNA_RCA | 56.5% | 76.6% | 65.0% |
| 50X | CircleMap | 92.8% | 74.0% | 82.3% |
| 50X | ecc_finder | 67.3% | 9.8% | 17.2% |

#### ColCEN_5200 (4,000 Uecc + 1,000 Mecc + 200 Cecc)

| Depth | Tool | Precision | Recall | F1 |
|-------|------|-----------|--------|------|
| 10X | CircleSeeker | **98.9%** | **90.8%** | **94.7%** |
| 10X | CReSIL-HiFi | 95.0% | 35.6% | 51.8% |
| 10X | CReSIL | 91.6% | 29.1% | 44.2% |
| 10X | eccDNA_RCA | 82.7% | 57.3% | 67.7% |
| 10X | CircleMap | 88.2% | 19.5% | 29.6% |
| 10X | ecc_finder | 18.9% | 0.6% | 1.1% |
| 30X | CircleSeeker | **98.9%** | **96.1%** | **97.5%** |
| 30X | CReSIL-HiFi | 96.3% | 46.8% | 63.0% |
| 30X | CReSIL | 96.1% | 43.5% | 59.8% |
| 30X | eccDNA_RCA | 66.6% | 66.9% | 66.7% |
| 30X | CircleMap | 94.8% | 60.4% | 73.8% |
| 30X | ecc_finder | 17.0% | 1.5% | 2.7% |
| 50X | CircleSeeker | **98.8%** | **97.0%** | **97.9%** |
| 50X | CReSIL-HiFi | 96.1% | 48.7% | 64.6% |
| 50X | CReSIL | 96.8% | 46.9% | 63.2% |
| 50X | eccDNA_RCA | 57.4% | 69.3% | 62.8% |
| 50X | CircleMap | 94.2% | 64.5% | 76.5% |
| 50X | ecc_finder | 19.6% | 2.0% | 3.7% |

### Per-type Recall (3-replicate average)

#### Uecc Recall

| Dataset | Depth | CircleSeeker | CReSIL-HiFi | CReSIL | eccDNA_RCA | CircleMap | ecc_finder |
|---------|-------|:-----------:|:-----------:|:------:|:----------:|:---------:|:----------:|
| human_10000_simple | 10X | **93.5%** | 66.1% | 51.9% | 73.1% | 25.3% | 13.9% |
| human_10000_simple | 30X | **98.6%** | 89.5% | 80.5% | 85.4% | 78.7% | 35.6% |
| human_10000_simple | 50X | **99.1%** | 91.8% | 85.8% | 88.4% | 84.4% | 42.2% |
| human_23000 | 10X | **93.1%** | 62.6% | 49.0% | 72.8% | 53.0% | 9.3% |
| human_23000 | 30X | **98.3%** | 84.5% | 76.1% | 85.3% | 78.5% | 11.4% |
| human_23000 | 50X | **98.8%** | 86.9% | 81.0% | 88.0% | 84.2% | 11.3% |
| ColCEN_5200 | 10X | **91.9%** | 44.7% | 37.2% | 74.5% | 25.3% | 0.7% |
| ColCEN_5200 | 30X | **96.9%** | 58.9% | 55.1% | 86.9% | 77.4% | 1.9% |
| ColCEN_5200 | 50X | **97.5%** | 60.8% | 59.0% | 90.1% | 82.0% | 2.6% |

#### Mecc Recall

| Dataset | Depth | CircleSeeker | CReSIL-HiFi | CReSIL | eccDNA_RCA | CircleMap | ecc_finder |
|---------|-------|:-----------:|:-----------:|:------:|:----------:|:---------:|:----------:|
| human_23000 | 10X | **77.2%** | 5.8% | 1.6% | 0.0% | 0.4% | 0.3% |
| human_23000 | 30X | **82.2%** | 9.0% | 4.5% | 0.0% | 3.0% | 0.4% |
| human_23000 | 50X | **83.7%** | 10.3% | 6.4% | 0.0% | 9.2% | 0.4% |
| ColCEN_5200 | 10X | **90.1%** | 5.6% | 1.8% | 0.0% | 0.2% | 0.0% |
| ColCEN_5200 | 30X | **95.2%** | 7.4% | 4.5% | 0.0% | 4.5% | 0.3% |
| ColCEN_5200 | 50X | **96.5%** | 9.7% | 6.0% | 0.0% | 7.2% | 0.3% |

#### Cecc Recall

| Dataset | Depth | CircleSeeker | CReSIL-HiFi | CReSIL | eccDNA_RCA | CircleMap | ecc_finder |
|---------|-------|:-----------:|:-----------:|:------:|:----------:|:---------:|:----------:|
| human_23000 | 10X | **76.9%** | 5.4% | 10.8% | 0.0% | 0.0% | 0.0% |
| human_23000 | 30X | **87.7%** | 9.4% | 22.6% | 0.0% | 0.0% | 0.0% |
| human_23000 | 50X | **90.8%** | 10.9% | 33.2% | 0.0% | 0.0% | 0.0% |
| ColCEN_5200 | 10X | **72.7%** | 1.8% | 4.0% | 0.0% | 0.0% | 0.0% |
| ColCEN_5200 | 30X | **85.0%** | 2.2% | 6.8% | 0.0% | 0.0% | 0.0% |
| ColCEN_5200 | 50X | **89.7%** | 1.8% | 9.5% | 0.0% | 0.0% | 0.0% |

### CircleSeeker Per-type Precision (3-replicate average)

CircleSeeker is the only tool that classifies each detection by eccDNA type, enabling full per-type Precision evaluation.

| Dataset | Depth | Uecc Precision | Mecc Precision | Cecc Precision |
|---------|-------|:--------------:|:--------------:|:--------------:|
| human_10000_simple | 10X | **99.8%** | - | - |
| human_10000_simple | 30X | **99.8%** | - | - |
| human_10000_simple | 50X | **99.6%** | - | - |
| human_23000 | 10X | **99.2%** | **86.0%** | **95.4%** |
| human_23000 | 30X | **99.2%** | **85.0%** | **93.8%** |
| human_23000 | 50X | **99.2%** | **84.5%** | **92.8%** |
| ColCEN_5200 | 10X | **99.2%** | **96.8%** | **98.2%** |
| ColCEN_5200 | 30X | **99.3%** | **96.7%** | **96.8%** |
| ColCEN_5200 | 50X | **99.3%** | **96.5%** | **96.0%** |

### Cecc Precision Comparison (3-replicate average)

Tools that produce multi-fragment detection groups can be evaluated for Cecc Precision. CircleSeeker explicitly classifies detections; CReSIL/CReSIL-HiFi use structural inference (multi-region groups = Cecc candidates).

| Dataset | Depth | CircleSeeker | CReSIL | CReSIL-HiFi |
|---------|-------|:-----------:|:------:|:-----------:|
| human_23000 | 10X | **95.4%** | 42.7% | 28.5% |
| human_23000 | 30X | **93.8%** | 67.8% | 40.5% |
| human_23000 | 50X | **92.8%** | 79.9% | 44.0% |
| ColCEN_5200 | 10X | **98.2%** | 33.8% | 21.6% |
| ColCEN_5200 | 30X | **96.8%** | 52.4% | 31.1% |
| ColCEN_5200 | 50X | **96.0%** | 59.6% | 26.8% |

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
