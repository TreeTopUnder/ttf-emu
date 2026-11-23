#!/usr/bin/env python3
"""Convert an Emu relative-abundance TSV into Krona ktImportText input.

Research use only. Not for clinical or diagnostic use.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

RANKS = (
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert Emu rel-abundance TSV to Krona text.")
    parser.add_argument("--report", required=True, help="Emu *_rel-abundance.tsv.")
    parser.add_argument("--output", required=True, help="Krona tab-delimited text output.")
    return parser.parse_args()


def _is_unassigned(tax_id: str, species: str) -> bool:
    tokens = {tax_id.strip().lower(), species.strip().lower()}
    return bool(tokens & {"", "unassigned", "unmapped", "unclassified", "unknown", "nan"})


def detect_metric(header: list[str]) -> str:
    lowered = {col.strip().lower(): col for col in header}
    if "estimated counts" in lowered:
        return lowered["estimated counts"]
    if "abundance" in lowered:
        return lowered["abundance"]
    raise ValueError("Emu report is missing 'estimated counts' and 'abundance' columns")


def read_report(path: Path) -> list[tuple[float, list[str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        header = list(reader.fieldnames)
        metric_col = detect_metric(header)
        rank_cols = []
        for rank in RANKS:
            for col in header:
                if col.strip().lower() == rank:
                    rank_cols.append(col)
                    break
        rows: list[tuple[float, list[str]]] = []
        for raw in reader:
            tax_id = (raw.get("tax_id") or "").strip()
            species = ""
            for col in header:
                if col.strip().lower() == "species":
                    species = (raw.get(col) or "").strip()
                    break
            try:
                abundance = float(raw.get(metric_col) or 0)
            except ValueError:
                continue
            if abundance <= 0:
                continue
            if _is_unassigned(tax_id, species):
                rows.append((abundance, ["UNCLASSIFIED"]))
                continue
            lineage = []
            for col in rank_cols:
                value = (raw.get(col) or "").strip()
                if value and value.lower() not in {"nan", "none", "null"}:
                    lineage.append(value)
            rows.append((abundance, lineage or ["UNCLASSIFIED"]))
        return rows


def write_krona(path: Path, rows: list[tuple[float, list[str]]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for abundance, lineage in rows:
            writer.writerow([f"{abundance:g}"] + lineage)


def main() -> None:
    args = parse_args()
    write_krona(Path(args.output), read_report(Path(args.report)))


if __name__ == "__main__":
    main()
