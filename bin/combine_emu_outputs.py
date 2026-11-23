#!/usr/bin/env python3
"""Combine Emu *_rel-abundance.tsv reports into species/genus tables.

Mirrors `emu combine-outputs --split-tables` filenames used by EMU_COMBINE_OUTPUTS:
  emu-combined-abundance-species.tsv
  emu-combined-abundance-species-counts.tsv
  emu-combined-abundance-genus.tsv

Research use only. Not for clinical or diagnostic use.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Combine Emu rel-abundance TSVs.")
    parser.add_argument(
        "--reports",
        nargs="+",
        required=True,
        help="Per-sample *_rel-abundance.tsv files.",
    )
    parser.add_argument("--outdir", required=True, help="Directory for combined tables.")
    return parser.parse_args()


def sample_name(path: Path) -> str:
    stem = path.name
    for suffix in ("_rel-abundance.tsv", "-rel-abundance.tsv"):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return path.stem


def _col(row: dict[str, str], name: str) -> str:
    lowered = {key.strip().lower(): key for key in row}
    key = lowered.get(name.lower())
    if key is None:
        return ""
    return (row.get(key) or "").strip()


def _metric(row: dict[str, str], name: str) -> float:
    value = _col(row, name)
    if not value:
        return 0.0
    try:
        return float(value)
    except ValueError:
        return 0.0


def load_reports(paths: list[Path]) -> list[tuple[str, list[dict[str, str]]]]:
    loaded: list[tuple[str, list[dict[str, str]]]] = []
    for path in paths:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            loaded.append((sample_name(path), list(reader)))
    return loaded


def combine_rank(
    reports: list[tuple[str, list[dict[str, str]]]],
    rank: str,
    metric: str,
) -> tuple[list[str], list[list[str]]]:
    samples = [name for name, _ in reports]
    values: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for sample, rows in reports:
        for row in rows:
            taxon = _col(row, rank)
            if not taxon:
                taxon = _col(row, "tax_id") or "unmapped"
            values[taxon][sample] += _metric(row, metric)
    header = [rank] + samples
    table_rows = []
    for taxon in sorted(values):
        table_rows.append([taxon] + [f"{values[taxon].get(sample, 0.0):g}" for sample in samples])
    return header, table_rows


def write_table(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def combine_outputs(report_paths: list[Path], outdir: Path) -> dict[str, Path]:
    reports = load_reports(report_paths)
    outdir.mkdir(parents=True, exist_ok=True)
    species_ab_header, species_ab_rows = combine_rank(reports, "species", "abundance")
    species_ct_header, species_ct_rows = combine_rank(reports, "species", "estimated counts")
    genus_ab_header, genus_ab_rows = combine_rank(reports, "genus", "abundance")
    outputs = {
        "species": outdir / "emu-combined-abundance-species.tsv",
        "species_counts": outdir / "emu-combined-abundance-species-counts.tsv",
        "genus": outdir / "emu-combined-abundance-genus.tsv",
    }
    write_table(outputs["species"], species_ab_header, species_ab_rows)
    write_table(outputs["species_counts"], species_ct_header, species_ct_rows)
    write_table(outputs["genus"], genus_ab_header, genus_ab_rows)
    return outputs


def main() -> None:
    args = parse_args()
    combine_outputs([Path(path) for path in args.reports], Path(args.outdir))


if __name__ == "__main__":
    main()
