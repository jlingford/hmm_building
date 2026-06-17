#!/usr/bin/env python3
"""Kay-style coarse MSA occupancy row filter.

This script implements a transparent Python interpretation of the coarse
alignment filter described in Kay et al. supplementary information. It removes
poorly aligned rows/sequences, not columns.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np


GAP_BYTES = {ord("-"), ord(".")}


def read_fasta(path):
    records = []
    name = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def write_fasta(records, path, width=60):
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i : i + width] + "\n")


def strip_alignment_gaps(records):
    stripped = []
    for name, seq in records:
        stripped.append((name, "".join(ch for ch in seq if ord(ch) not in GAP_BYTES)))
    return stripped


def load_alignment_matrix(records):
    if not records:
        raise ValueError("No FASTA records found")
    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise ValueError(f"Input is not aligned: found lengths {sorted(lengths)[:10]}")
    aln_len = len(records[0][1])
    joined = "".join(seq.upper() for _, seq in records).encode("ascii")
    return np.frombuffer(joined, dtype=np.uint8).reshape(len(records), aln_len)


def quantile(values, q):
    return float(np.quantile(values, q)) if len(values) else None


def run_filter(
    alignment,
    outdir,
    prefix=None,
    row_sum_fraction=0.80,
    column_score_fraction=0.80,
    low_occupancy_fraction=0.10,
    high_occupancy_fraction=0.90,
    max_low_occupancy_columns=0.10,
    max_high_occupancy_breaks=0.10,
    verbose=True,
):
    start = time.perf_counter()
    alignment = Path(alignment)
    outdir = Path(outdir)
    prefix = prefix or alignment.name.removesuffix(".faa")
    outdir.mkdir(parents=True, exist_ok=True)

    records = read_fasta(alignment)
    names = [name for name, _ in records]
    matrix = load_alignment_matrix(records)
    nseq, aln_len = matrix.shape

    is_gap = np.isin(matrix, list(GAP_BYTES))
    non_gap = ~is_gap
    residue_counts = non_gap.sum(axis=0).astype(np.float64)
    occupancy_fraction = residue_counts / nseq

    # Kay-style occupancy score: an occupied cell receives the occupancy count
    # of its column; gap cells receive zero.
    occupancy_scores = non_gap * residue_counts
    row_sum = occupancy_scores.sum(axis=1)
    non_gap_count = non_gap.sum(axis=1).astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        row_mean_column_score = np.where(non_gap_count > 0, row_sum / non_gap_count, 0.0)

    low_occ_columns = occupancy_fraction < low_occupancy_fraction
    high_occ_columns = occupancy_fraction > high_occupancy_fraction
    low_occ_total = int(low_occ_columns.sum())
    high_occ_total = int(high_occ_columns.sum())

    # Interpreted per row:
    # - low occupancy columns: fraction of the row's occupied residues that sit
    #   in columns occupied by <10% of all sequences.
    # - high occupancy breaks: fraction of highly occupied columns where this
    #   row has a gap, i.e. the row breaks otherwise well-occupied columns.
    low_occ_hits = (non_gap[:, low_occ_columns]).sum(axis=1) if low_occ_total else np.zeros(nseq)
    with np.errstate(divide="ignore", invalid="ignore"):
        low_occ_column_fraction = np.where(non_gap_count > 0, low_occ_hits / non_gap_count, 0.0)

    high_occ_breaks = (is_gap[:, high_occ_columns]).sum(axis=1) if high_occ_total else np.zeros(nseq)
    high_occ_break_fraction = (
        high_occ_breaks / high_occ_total if high_occ_total else np.zeros(nseq, dtype=np.float64)
    )

    row_sum_mean = float(row_sum.mean())
    row_mean_column_score_mean = float(row_mean_column_score.mean())
    row_sum_cutoff = row_sum_fraction * row_sum_mean
    row_mean_column_score_cutoff = column_score_fraction * row_mean_column_score_mean

    fail_row_sum = row_sum < row_sum_cutoff
    fail_column_score = row_mean_column_score < row_mean_column_score_cutoff
    fail_low_occ = low_occ_column_fraction > max_low_occupancy_columns
    fail_high_break = high_occ_break_fraction > max_high_occupancy_breaks
    remove = fail_row_sum | fail_column_score | fail_low_occ | fail_high_break

    fields = [
        "safe_id",
        "non_gap_count",
        "gap_fraction",
        "row_sum",
        "row_sum_cutoff",
        "row_mean_column_score",
        "row_mean_column_score_cutoff",
        "low_occupancy_column_fraction",
        "high_occupancy_break_fraction",
        "fail_row_sum",
        "fail_column_score",
        "fail_low_occupancy_columns",
        "fail_high_occupancy_breaks",
        "remove_candidate",
        "flags",
    ]

    qc_path = outdir / f"{prefix}.coarse_occupancy_qc.tsv"
    remove_path = outdir / f"{prefix}.coarse_occupancy_remove.tsv"
    kept_aligned_path = outdir / f"{prefix}.coarse_occupancy.filtered.aln.faa"
    kept_unaligned_path = outdir / f"{prefix}.coarse_occupancy.filtered.unaligned.faa"
    summary_path = outdir / f"{prefix}.coarse_occupancy_summary.json"

    kept_records = []
    removed_records = []
    with open(qc_path, "w") as qc, open(remove_path, "w") as rem:
        qc.write("\t".join(fields) + "\n")
        rem.write("\t".join(fields) + "\n")
        for i, (name, seq) in enumerate(records):
            flags = []
            if fail_row_sum[i]:
                flags.append("low_row_sum")
            if fail_column_score[i]:
                flags.append("low_mean_column_score")
            if fail_low_occ[i]:
                flags.append("too_many_low_occupancy_residue_columns")
            if fail_high_break[i]:
                flags.append("too_many_high_occupancy_gap_breaks")

            row = {
                "safe_id": names[i],
                "non_gap_count": int(non_gap_count[i]),
                "gap_fraction": float(is_gap[i].mean()),
                "row_sum": float(row_sum[i]),
                "row_sum_cutoff": row_sum_cutoff,
                "row_mean_column_score": float(row_mean_column_score[i]),
                "row_mean_column_score_cutoff": row_mean_column_score_cutoff,
                "low_occupancy_column_fraction": float(low_occ_column_fraction[i]),
                "high_occupancy_break_fraction": float(high_occ_break_fraction[i]),
                "fail_row_sum": bool(fail_row_sum[i]),
                "fail_column_score": bool(fail_column_score[i]),
                "fail_low_occupancy_columns": bool(fail_low_occ[i]),
                "fail_high_occupancy_breaks": bool(fail_high_break[i]),
                "remove_candidate": bool(remove[i]),
                "flags": ",".join(flags),
            }
            line = "\t".join(str(row[field]) for field in fields) + "\n"
            qc.write(line)
            if remove[i]:
                rem.write(line)
                removed_records.append((name, seq))
            else:
                kept_records.append((name, seq))

    write_fasta(kept_records, kept_aligned_path)
    write_fasta(strip_alignment_gaps(kept_records), kept_unaligned_path)

    elapsed = time.perf_counter() - start
    summary = {
        "alignment": str(alignment),
        "num_sequences": nseq,
        "alignment_length": aln_len,
        "kept_sequences": len(kept_records),
        "removed_sequences": len(removed_records),
        "removed_fraction": len(removed_records) / nseq if nseq else 0,
        "elapsed_seconds": elapsed,
        "thresholds": {
            "row_sum_remove_if_lt": row_sum_cutoff,
            "row_sum_fraction_of_mean": row_sum_fraction,
            "row_mean_column_score_remove_if_lt": row_mean_column_score_cutoff,
            "column_score_fraction_of_mean": column_score_fraction,
            "low_occupancy_column_definition_lt_fraction": low_occupancy_fraction,
            "remove_if_low_occupancy_column_fraction_gt": max_low_occupancy_columns,
            "high_occupancy_column_definition_gt_fraction": high_occupancy_fraction,
            "remove_if_high_occupancy_break_fraction_gt": max_high_occupancy_breaks,
        },
        "column_counts": {
            "low_occupancy_columns": low_occ_total,
            "high_occupancy_columns": high_occ_total,
        },
        "metric_summaries": {
            "row_sum": {
                "min": float(row_sum.min()),
                "q1": quantile(row_sum, 0.25),
                "median": quantile(row_sum, 0.5),
                "q3": quantile(row_sum, 0.75),
                "max": float(row_sum.max()),
                "mean": row_sum_mean,
            },
            "row_mean_column_score": {
                "min": float(row_mean_column_score.min()),
                "q1": quantile(row_mean_column_score, 0.25),
                "median": quantile(row_mean_column_score, 0.5),
                "q3": quantile(row_mean_column_score, 0.75),
                "max": float(row_mean_column_score.max()),
                "mean": row_mean_column_score_mean,
            },
        },
        "flag_counts": {
            "low_row_sum": int(fail_row_sum.sum()),
            "low_mean_column_score": int(fail_column_score.sum()),
            "too_many_low_occupancy_residue_columns": int(fail_low_occ.sum()),
            "too_many_high_occupancy_gap_breaks": int(fail_high_break.sum()),
        },
        "interpretation_notes": [
            "This filter removes rows/sequences from an existing MSA.",
            "It is a Python interpretation of Kay et al.'s prose description, not the original awk script.",
            "The filtered unaligned FASTA is intended for realignment after sequence removal.",
        ],
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    if verbose:
        print(f"wrote {qc_path}")
        print(f"wrote {remove_path}")
        print(f"wrote {kept_aligned_path}")
        print(f"wrote {kept_unaligned_path}")
        print(f"wrote {summary_path}")
        print(json.dumps(summary, indent=2))

    return {
        "qc_path": qc_path,
        "remove_path": remove_path,
        "kept_aligned_path": kept_aligned_path,
        "kept_unaligned_path": kept_unaligned_path,
        "summary_path": summary_path,
        "summary": summary,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=(
            "Remove poorly aligned MSA rows using a Kay-style coarse occupancy filter. "
            "Rows are removed if they fail any of the four occupancy criteria."
        )
    )
    ap.add_argument("--alignment", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--prefix")
    ap.add_argument("--row-sum-fraction", type=float, default=0.80)
    ap.add_argument("--column-score-fraction", type=float, default=0.80)
    ap.add_argument("--low-occupancy-fraction", type=float, default=0.10)
    ap.add_argument("--high-occupancy-fraction", type=float, default=0.90)
    ap.add_argument("--max-low-occupancy-columns", type=float, default=0.10)
    ap.add_argument("--max-high-occupancy-breaks", type=float, default=0.10)
    args = ap.parse_args(argv)
    run_filter(
        alignment=args.alignment,
        outdir=args.outdir,
        prefix=args.prefix,
        row_sum_fraction=args.row_sum_fraction,
        column_score_fraction=args.column_score_fraction,
        low_occupancy_fraction=args.low_occupancy_fraction,
        high_occupancy_fraction=args.high_occupancy_fraction,
        max_low_occupancy_columns=args.max_low_occupancy_columns,
        max_high_occupancy_breaks=args.max_high_occupancy_breaks,
    )


if __name__ == "__main__":
    main()
