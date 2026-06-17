#!/usr/bin/env python3
"""Plot a binned global gap heatmap for a protein MSA.

The script writes a binary PPM image directly, so it does not require
matplotlib/Pillow. Use ImageMagick's `convert` to make PNGs if desired.
"""

from __future__ import annotations

import argparse
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


def load_gap_matrix(records):
    if not records:
        raise ValueError("No FASTA records found")
    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise ValueError(f"Input is not aligned: found lengths {sorted(lengths)[:10]}")
    joined = "".join(seq.upper() for _, seq in records).encode("ascii")
    matrix = np.frombuffer(joined, dtype=np.uint8).reshape(len(records), len(records[0][1]))
    return np.isin(matrix, list(GAP_BYTES)).astype(np.float32)


def bin_matrix(matrix, row_bins, col_bins):
    nrow, ncol = matrix.shape
    row_edges = np.linspace(0, nrow, min(row_bins, nrow) + 1, dtype=int)
    col_edges = np.linspace(0, ncol, min(col_bins, ncol) + 1, dtype=int)
    out = np.zeros((len(row_edges) - 1, len(col_edges) - 1), dtype=np.float32)
    for i in range(len(row_edges) - 1):
        block_rows = matrix[row_edges[i] : row_edges[i + 1], :]
        for j in range(len(col_edges) - 1):
            out[i, j] = block_rows[:, col_edges[j] : col_edges[j + 1]].mean()
    return out


def heatmap_rgb(values):
    """Simple blue-to-yellow heatmap for values in [0, 1]."""
    x = np.clip(values, 0, 1)
    r = (255 * x).astype(np.uint8)
    g = (255 * np.sqrt(x)).astype(np.uint8)
    b = (255 * (1 - x)).astype(np.uint8)
    return np.dstack([r, g, b])


def write_ppm(rgb, path):
    height, width, _ = rgb.shape
    with open(path, "wb") as out:
        out.write(f"P6\n{width} {height}\n255\n".encode("ascii"))
        out.write(rgb.tobytes())


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--alignment", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path, help="Output .ppm path")
    ap.add_argument("--row-bins", type=int, default=800)
    ap.add_argument("--col-bins", type=int, default=800)
    args = ap.parse_args(argv)

    records = read_fasta(args.alignment)
    gap = load_gap_matrix(records)
    binned = bin_matrix(gap, args.row_bins, args.col_bins)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_ppm(heatmap_rgb(binned), args.out)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
