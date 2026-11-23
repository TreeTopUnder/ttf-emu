#!/usr/bin/env python3
"""Normalize sample metadata TSV to ampliseq / ttf-rstudio-ampliseq conventions.

The first column must be ID or sample. It is rewritten as ID.
Research use only. Not for clinical or diagnostic use.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Normalize metadata TSV first column to ID.")
    parser.add_argument("--input", required=True, help="Input metadata TSV.")
    parser.add_argument("--output", required=True, help="Normalized sample_metadata.tsv.")
    return parser.parse_args()


def normalize_metadata(in_path: Path, out_path: Path) -> None:
    with in_path.open(newline="") as handle:
        # Skip comment lines while keeping the header.
        lines = [line for line in handle if line.strip() and not line.startswith("#")]
    if not lines:
        raise ValueError("Metadata TSV is empty")

    reader = csv.reader(lines, delimiter="\t")
    rows = list(reader)
    header = rows[0]
    if not header or not header[0].strip():
        raise ValueError("Metadata TSV is missing a header")

    first = header[0].strip()
    if first.lower() not in {"id", "sample"}:
        raise ValueError("Metadata first column must be ID or sample (ampliseq convention)")

    header[0] = "ID"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows[1:])


def main() -> None:
    args = parse_args()
    normalize_metadata(Path(args.input), Path(args.output))


if __name__ == "__main__":
    main()
