#!/usr/bin/env python3
"""Build and evaluate the round-3 synthetic HiFi profile sensitivity test."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


DNA = np.array(list("ACGT"))


def truncated_lognormal(
    rng: np.random.Generator,
    n: int,
    median: float,
    sigma: float,
    minimum: int,
    maximum: int,
) -> np.ndarray:
    values: list[int] = []
    while len(values) < n:
        draw = rng.lognormal(mean=math.log(median), sigma=sigma, size=max(n, 1024))
        keep = draw[(draw >= minimum) & (draw <= maximum)].astype(int)
        values.extend(keep.tolist())
    return np.asarray(values[:n], dtype=np.int64)


def n50(lengths: np.ndarray) -> int:
    ordered = np.sort(lengths)[::-1]
    return int(ordered[np.searchsorted(np.cumsum(ordered), ordered.sum() / 2)])


def generate_profile(args: argparse.Namespace) -> None:
    rng = np.random.default_rng(args.seed)
    lengths = truncated_lognormal(
        rng,
        args.reads,
        args.median_length,
        args.sigma,
        args.min_length,
        args.max_length,
    )
    qscores = np.clip(
        np.rint(rng.normal(args.qscore_mean, args.qscore_sd, size=args.reads)),
        args.qscore_min,
        args.qscore_max,
    ).astype(int)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        for idx, (length, qscore) in enumerate(zip(lengths, qscores, strict=True), start=1):
            sequence = "".join(rng.choice(DNA, size=int(length)))
            quality = chr(int(qscore) + 33) * int(length)
            handle.write(f"@synthetic_profile_{idx:06d}\n{sequence}\n+\n{quality}\n")

    stats = {
        "profile": args.profile_name,
        "parameter_source": "prespecified synthetic truncated log-normal distribution",
        "seed": args.seed,
        "n_reads": int(len(lengths)),
        "length_min": int(lengths.min()),
        "length_max": int(lengths.max()),
        "length_mean": float(lengths.mean()),
        "length_sd": float(lengths.std(ddof=1)),
        "length_median": float(np.median(lengths)),
        "length_n50": n50(lengths),
        "qscore_mean": float(qscores.mean()),
        "qscore_sd": float(qscores.std(ddof=1)),
        "distribution_median_parameter": args.median_length,
        "distribution_sigma_parameter": args.sigma,
        "distribution_min_parameter": args.min_length,
        "distribution_max_parameter": args.max_length,
    }
    args.stats.write_text(json.dumps(stats, indent=2) + "\n")


def parse_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def prepare_subset(args: argparse.Namespace) -> None:
    args.output_lib.parent.mkdir(parents=True, exist_ok=True)
    lib = pd.read_csv(args.source_lib, sep="\t")
    background = parse_bool(lib["is_background"])
    circular = lib.loc[~background]
    linear = lib.loc[background]
    if len(circular) < args.n_circular or len(linear) < args.n_background:
        raise ValueError("Requested subset is larger than the available molecule pool")

    selected_circular = circular.sample(n=args.n_circular, random_state=args.seed)
    selected_linear = linear.sample(n=args.n_background, random_state=args.seed + 1)
    subset = pd.concat([selected_circular, selected_linear], ignore_index=True)
    subset.to_csv(args.output_lib, sep="\t", index=False)

    truth_names = set(selected_circular["ecc_ids"].astype(str))
    with args.source_truth.open() as source, args.output_truth.open("w") as output:
        for line in source:
            if line.startswith("#"):
                output.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4 and fields[3] in truth_names:
                output.write(line)

    metadata = {
        "source_lib": str(args.source_lib),
        "source_truth": str(args.source_truth),
        "seed": args.seed,
        "n_circular": args.n_circular,
        "n_background": args.n_background,
        "selected_total": int(len(subset)),
        "tempcov_mean_circular": float(selected_circular["tempcov"].mean()),
        "tempcov_mean_background": float(selected_linear["tempcov"].mean()),
    }
    args.metadata.write_text(json.dumps(metadata, indent=2) + "\n")


def simulate(args: argparse.Namespace) -> None:
    from ecctoolkit.simulate.readsim import fqsim

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fqsim(
        sample=args.sample,
        csv=str(args.input_lib),
        path=str(args.output_dir),
        seed=args.seed,
        thread=args.threads,
        skip_sr=True,
        skip_hifi=False,
        skip_ont=True,
        hifi_sample_fastq=str(args.sample_fastq) if args.sample_fastq else None,
        hifi_profile_id=args.profile_id,
        hifi_profile_root=str(args.profile_root) if args.profile_root else None,
        generate_truth=True,
    )


def detect_profile_format(path: Path) -> str:
    with path.open() as handle:
        first_four = [handle.readline().rstrip("\n") for _ in range(4)]
    if (
        len(first_four) == 4
        and first_four[0].startswith("@")
        and first_four[2].startswith("+")
        and len(first_four[1]) == len(first_four[3])
    ):
        return "fastq"
    return "pbsim-profile"


def read_fastq_stats(path: Path, input_format: str = "auto") -> dict[str, float | int | str]:
    if input_format == "auto":
        input_format = detect_profile_format(path)

    lengths: list[int] = []
    qsum = 0
    qcount = 0
    if input_format == "fastq":
        with path.open() as handle:
            for line_number, line in enumerate(handle):
                offset = line_number % 4
                if offset == 1:
                    lengths.append(len(line.rstrip("\n")))
                elif offset == 3:
                    qualities = line.rstrip("\n")
                    qsum += sum(ord(char) - 33 for char in qualities)
                    qcount += len(qualities)
    elif input_format == "pbsim-profile":
        # PBSIM sampling profiles contain one quality string per source read,
        # not four-line FASTQ records, despite the conventional .fastq suffix.
        with path.open() as handle:
            for line in handle:
                qualities = line.rstrip("\n")
                if not qualities:
                    continue
                lengths.append(len(qualities))
                qsum += sum(ord(char) - 33 for char in qualities)
                qcount += len(qualities)
    else:
        raise ValueError(f"Unsupported input format: {input_format}")

    arr = np.asarray(lengths, dtype=np.int64)
    return {
        "fastq": str(path),
        "input_format": input_format,
        "n_reads": int(len(arr)),
        "total_bases": int(arr.sum()),
        "length_min": int(arr.min()),
        "length_max": int(arr.max()),
        "length_mean": float(arr.mean()),
        "length_sd": float(arr.std(ddof=1)),
        "length_median": float(np.median(arr)),
        "length_n50": n50(arr),
        "mean_base_qscore": qsum / qcount,
    }


def summarize_fastq(args: argparse.Namespace) -> None:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(read_fastq_stats(args.fastq, args.input_format), indent=2) + "\n"
    )


def read_truth(path: Path) -> list[tuple[str, int, int, str]]:
    truth = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, name = line.rstrip("\n").split("\t")[:4]
            truth.append((chrom, int(start), int(end), name))
    return truth


def parse_circleseeker(path: Path) -> dict[str, list[tuple[str, int, int]]]:
    detected: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    with path.open() as handle:
        for row in csv.DictReader(handle):
            ecc_id = row["eccDNA_id"]
            if row.get("type") == "Cecc":
                locations = row.get("location", "").split(";")
            elif row.get("type") == "Mecc":
                locations = row.get("location", "").split("|")
            else:
                detected[ecc_id].append((row["chr"], int(row["start"]), int(row["end"])))
                continue
            for location in locations:
                match = re.match(r"([^:]+):(\d+)-(\d+)", location)
                if match:
                    detected[ecc_id].append(
                        (match.group(1), int(match.group(2)), int(match.group(3)))
                    )
    return detected


def parse_cresil(path: Path) -> dict[str, list[tuple[str, int, int]]]:
    detected: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    with path.open() as handle:
        next(handle, None)
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            for location in fields[1].split(","):
                match = re.match(r"([^:]+):(\d+)-(\d+)", location)
                if match:
                    detected[fields[0]].append(
                        (match.group(1), int(match.group(2)), int(match.group(3)))
                    )
    return detected


def reciprocal_overlap(a_start: int, a_end: int, b_start: int, b_end: int, threshold: float) -> bool:
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start))
    if overlap == 0:
        return False
    return overlap / (a_end - a_start) >= threshold and overlap / (b_end - b_start) >= threshold


def evaluate_calls(
    truth: list[tuple[str, int, int, str]],
    detected: dict[str, list[tuple[str, int, int]]],
    threshold: float,
) -> tuple[dict[str, float | int], set[str]]:
    truth_by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    truth_length: dict[str, int] = {}
    for chrom, start, end, name in truth:
        truth_by_chrom[chrom].append((start, end, name))
        truth_length[name] = end - start
    for chrom in truth_by_chrom:
        truth_by_chrom[chrom].sort()

    matched_truth: set[str] = set()
    matched_detection: set[str] = set()
    for detected_id, regions in detected.items():
        for chrom, start, end in regions:
            for truth_start, truth_end, truth_name in truth_by_chrom.get(chrom, []):
                if truth_name in matched_truth:
                    continue
                if truth_start >= end:
                    break
                if truth_end <= start:
                    continue
                if reciprocal_overlap(start, end, truth_start, truth_end, threshold):
                    matched_truth.add(truth_name)
                    matched_detection.add(detected_id)
                    break
            if detected_id in matched_detection:
                break

    tp = len(matched_truth)
    detected_count = len(detected)
    truth_count = len(truth)
    precision = tp / detected_count if detected_count else 0.0
    recall = tp / truth_count if truth_count else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    result: dict[str, float | int] = {
        "truth": truth_count,
        "detected": detected_count,
        "tp": tp,
        "fp": detected_count - tp,
        "fn": truth_count - tp,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }
    for label, lower, upper in [
        ("lt_1kb", 0, 1_000),
        ("1_5kb", 1_000, 5_000),
        ("5_50kb", 5_000, 50_000),
        ("ge_50kb", 50_000, math.inf),
    ]:
        names = {name for name, length in truth_length.items() if lower <= length < upper}
        result[f"truth_{label}"] = len(names)
        result[f"recall_{label}"] = len(names & matched_truth) / len(names) if names else math.nan
    return result, matched_truth


def evaluate(args: argparse.Namespace) -> None:
    truth = read_truth(args.truth)
    if args.tool == "CircleSeeker":
        detected = parse_circleseeker(args.calls)
    else:
        detected = parse_cresil(args.calls)
    metrics, _ = evaluate_calls(truth, detected, args.threshold)
    metrics.update(
        {
            "replicate": args.replicate,
            "profile": args.profile,
            "tool": args.tool,
            "truth_file": str(args.truth),
            "calls_file": str(args.calls),
            "overlap_threshold": args.threshold,
        }
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([metrics]).to_csv(args.output, index=False)


def aggregate(args: argparse.Namespace) -> None:
    frames = [pd.read_csv(path) for path in args.metrics]
    replicate_metrics = pd.concat(frames, ignore_index=True)
    args.output_replicates.parent.mkdir(parents=True, exist_ok=True)
    replicate_metrics.to_csv(args.output_replicates, index=False)

    metric_columns = [
        "precision",
        "recall",
        "f1",
        "recall_lt_1kb",
        "recall_1_5kb",
        "recall_5_50kb",
        "recall_ge_50kb",
    ]
    summary_rows = []
    for (tool, profile), group in replicate_metrics.groupby(["tool", "profile"]):
        row: dict[str, float | int | str] = {
            "tool": tool,
            "profile": profile,
            "n_replicates": len(group),
        }
        for metric in metric_columns:
            row[f"{metric}_mean"] = group[metric].mean()
            row[f"{metric}_sd"] = group[metric].std(ddof=1)
        summary_rows.append(row)
    pd.DataFrame(summary_rows).to_csv(args.output_summary, index=False)

    profiles = sorted(replicate_metrics["profile"].unique())
    paired_rows = []
    if len(profiles) == 2:
        for tool, group in replicate_metrics.groupby("tool"):
            for metric in metric_columns:
                pivot = group.pivot(index="replicate", columns="profile", values=metric).dropna()
                differences = pivot[profiles[1]] - pivot[profiles[0]]
                paired_rows.append(
                    {
                        "tool": tool,
                        "metric": metric,
                        "profile_a": profiles[0],
                        "profile_b": profiles[1],
                        "difference_definition": "profile_b_minus_profile_a",
                        "n_paired_replicates": len(differences),
                        "mean_difference": differences.mean(),
                        "sd_difference": differences.std(ddof=1),
                        "min_difference": differences.min(),
                        "max_difference": differences.max(),
                    }
                )
    pd.DataFrame(paired_rows).to_csv(args.output_paired, index=False)


def aggregate_fastq(args: argparse.Namespace) -> None:
    rows = []
    for path in args.stats:
        row = json.loads(path.read_text())
        row["profile"] = path.parent.name
        row["replicate"] = path.parent.parent.name
        rows.append(row)
    replicate_stats = pd.DataFrame(rows)
    args.output_replicates.parent.mkdir(parents=True, exist_ok=True)
    replicate_stats.to_csv(args.output_replicates, index=False)

    numeric = [
        "n_reads",
        "total_bases",
        "length_mean",
        "length_sd",
        "length_median",
        "length_n50",
        "mean_base_qscore",
    ]
    summary_rows = []
    for profile, group in replicate_stats.groupby("profile"):
        row: dict[str, float | int | str] = {
            "profile": profile,
            "n_replicates": len(group),
        }
        for metric in numeric:
            row[f"{metric}_mean"] = group[metric].mean()
            row[f"{metric}_sd"] = group[metric].std(ddof=1)
        summary_rows.append(row)
    pd.DataFrame(summary_rows).to_csv(args.output_summary, index=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    profile = subparsers.add_parser("generate-profile")
    profile.add_argument("--output", required=True, type=Path)
    profile.add_argument("--stats", required=True, type=Path)
    profile.add_argument("--profile-name", default="synthetic_lognormal_q39")
    profile.add_argument("--reads", default=2_000, type=int)
    profile.add_argument("--median-length", default=12_000.0, type=float)
    profile.add_argument("--sigma", default=0.65, type=float)
    profile.add_argument("--min-length", default=3_000, type=int)
    profile.add_argument("--max-length", default=50_000, type=int)
    profile.add_argument("--qscore-mean", default=39.0, type=float)
    profile.add_argument("--qscore-sd", default=2.0, type=float)
    profile.add_argument("--qscore-min", default=30, type=int)
    profile.add_argument("--qscore-max", default=45, type=int)
    profile.add_argument("--seed", default=20260710, type=int)
    profile.set_defaults(func=generate_profile)

    subset = subparsers.add_parser("prepare-subset")
    subset.add_argument("--source-lib", required=True, type=Path)
    subset.add_argument("--source-truth", required=True, type=Path)
    subset.add_argument("--output-lib", required=True, type=Path)
    subset.add_argument("--output-truth", required=True, type=Path)
    subset.add_argument("--metadata", required=True, type=Path)
    subset.add_argument("--n-circular", default=2_000, type=int)
    subset.add_argument("--n-background", default=2_000, type=int)
    subset.add_argument("--seed", required=True, type=int)
    subset.set_defaults(func=prepare_subset)

    simulation = subparsers.add_parser("simulate")
    simulation.add_argument("--input-lib", required=True, type=Path)
    simulation.add_argument("--output-dir", required=True, type=Path)
    simulation.add_argument("--sample", required=True)
    simulation.add_argument("--sample-fastq", type=Path)
    simulation.add_argument("--profile-id")
    simulation.add_argument("--profile-root", type=Path)
    simulation.add_argument("--threads", default=4, type=int)
    simulation.add_argument("--seed", required=True, type=int)
    simulation.set_defaults(func=simulate)

    fastq = subparsers.add_parser("summarize-fastq")
    fastq.add_argument("--fastq", required=True, type=Path)
    fastq.add_argument("--output", required=True, type=Path)
    fastq.add_argument(
        "--input-format",
        default="auto",
        choices=["auto", "fastq", "pbsim-profile"],
    )
    fastq.set_defaults(func=summarize_fastq)

    evaluation = subparsers.add_parser("evaluate")
    evaluation.add_argument("--truth", required=True, type=Path)
    evaluation.add_argument("--calls", required=True, type=Path)
    evaluation.add_argument("--output", required=True, type=Path)
    evaluation.add_argument("--tool", required=True, choices=["CircleSeeker", "CReSIL-HiFi"])
    evaluation.add_argument("--replicate", required=True)
    evaluation.add_argument("--profile", required=True)
    evaluation.add_argument("--threshold", default=0.9, type=float)
    evaluation.set_defaults(func=evaluate)

    aggregation = subparsers.add_parser("aggregate")
    aggregation.add_argument("--metrics", required=True, nargs="+", type=Path)
    aggregation.add_argument("--output-replicates", required=True, type=Path)
    aggregation.add_argument("--output-summary", required=True, type=Path)
    aggregation.add_argument("--output-paired", required=True, type=Path)
    aggregation.set_defaults(func=aggregate)

    fastq_aggregation = subparsers.add_parser("aggregate-fastq")
    fastq_aggregation.add_argument("--stats", required=True, nargs="+", type=Path)
    fastq_aggregation.add_argument("--output-replicates", required=True, type=Path)
    fastq_aggregation.add_argument("--output-summary", required=True, type=Path)
    fastq_aggregation.set_defaults(func=aggregate_fastq)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
