#!/usr/bin/env python3
"""Column QC and gap-fraction trimming for aligned FASTA files."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path


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
                name = line[1:]
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


def quantile(values, q):
    vals = sorted(values)
    if not vals:
        return None
    pos = (len(vals) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(vals) - 1)
    if lo == hi:
        return vals[lo]
    return vals[lo] * (hi - pos) + vals[hi] * (pos - lo)


def column_qc(records, gap_chars):
    nseq = len(records)
    aln_len = len(records[0][1]) if records else 0
    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise ValueError(f"Input is not aligned: found sequence lengths {sorted(lengths)[:10]}")

    rows = []
    for pos in range(aln_len):
        counts = Counter(seq[pos].upper() for _, seq in records)
        gap_count = sum(counts.get(g, 0) for g in gap_chars)
        non_gap_count = nseq - gap_count
        residues = {k: v for k, v in counts.items() if k not in gap_chars}
        if residues:
            dominant_residue, dominant_count = max(residues.items(), key=lambda kv: kv[1])
            dominant_fraction = dominant_count / non_gap_count
        else:
            dominant_residue = ""
            dominant_fraction = 0.0
        rows.append(
            {
                "column_1based": pos + 1,
                "gap_count": gap_count,
                "gap_fraction": gap_count / nseq,
                "non_gap_count": non_gap_count,
                "non_gap_fraction": non_gap_count / nseq,
                "dominant_residue": dominant_residue,
                "dominant_residue_fraction_of_non_gap": dominant_fraction,
            }
        )
    return rows


def trim_records(records, keep_positions):
    keep = set(keep_positions)
    trimmed = []
    for name, seq in records:
        trimmed.append((name, "".join(base for i, base in enumerate(seq) if i in keep)))
    return trimmed


def write_qc(rows, path):
    fields = [
        "column_1based",
        "gap_count",
        "gap_fraction",
        "non_gap_count",
        "non_gap_fraction",
        "dominant_residue",
        "dominant_residue_fraction_of_non_gap",
    ]
    with open(path, "w") as out:
        out.write("\t".join(fields) + "\n")
        for row in rows:
            out.write("\t".join(str(row[f]) for f in fields) + "\n")


def summarize(rows, nseq, aln_len, gap_threshold, keep_positions):
    gaps = [row["gap_fraction"] for row in rows]
    return {
        "num_sequences": nseq,
        "input_alignment_length": aln_len,
        "gap_threshold_remove_if_gt": gap_threshold,
        "kept_columns": len(keep_positions),
        "removed_columns": aln_len - len(keep_positions),
        "removed_fraction": (aln_len - len(keep_positions)) / aln_len if aln_len else 0,
        "gap_fraction_min": min(gaps) if gaps else None,
        "gap_fraction_q1": quantile(gaps, 0.25),
        "gap_fraction_median": quantile(gaps, 0.5),
        "gap_fraction_q3": quantile(gaps, 0.75),
        "gap_fraction_max": max(gaps) if gaps else None,
        "columns_gt_50pct_gap": sum(1 for g in gaps if g > 0.5),
        "columns_gt_70pct_gap": sum(1 for g in gaps if g > 0.7),
        "columns_gt_80pct_gap": sum(1 for g in gaps if g > 0.8),
        "columns_gt_90pct_gap": sum(1 for g in gaps if g > 0.9),
        "columns_gt_95pct_gap": sum(1 for g in gaps if g > 0.95),
    }


def run_msa_column_qc_trim(alignment, outdir, prefix=None, gap_threshold=0.8, verbose=True):
    alignment = Path(alignment)
    outdir = Path(outdir)
    prefix = prefix or alignment.name.removesuffix(".faa")
    outdir.mkdir(parents=True, exist_ok=True)

    records = read_fasta(alignment)
    if not records:
        raise ValueError(f"No FASTA records found in {alignment}")
    nseq = len(records)
    aln_len = len(records[0][1])
    gap_chars = {"-", "."}

    rows = column_qc(records, gap_chars)
    keep_positions = [i for i, row in enumerate(rows) if row["gap_fraction"] <= gap_threshold]
    trimmed = trim_records(records, keep_positions)

    qc_path = outdir / f"{prefix}.column_qc.tsv"
    trimmed_path = outdir / f"{prefix}.gap{int(gap_threshold * 100)}.trimmed.aln.faa"
    summary_path = outdir / f"{prefix}.gap{int(gap_threshold * 100)}.trim_summary.json"

    write_qc(rows, qc_path)
    write_fasta(trimmed, trimmed_path)
    summary = summarize(rows, nseq, aln_len, gap_threshold, keep_positions)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    result = {
        "qc_path": qc_path,
        "trimmed_path": trimmed_path,
        "summary_path": summary_path,
        "summary": summary,
    }
    if verbose:
        print(f"wrote {qc_path}")
        print(f"wrote {trimmed_path}")
        print(f"wrote {summary_path}")
        print(json.dumps(summary, indent=2))
    return result


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--alignment", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--prefix")
    ap.add_argument("--gap-threshold", type=float, default=0.8)
    args = ap.parse_args(argv)
    run_msa_column_qc_trim(
        alignment=args.alignment,
        outdir=args.outdir,
        prefix=args.prefix,
        gap_threshold=args.gap_threshold,
    )


if __name__ == "__main__":
    main()
