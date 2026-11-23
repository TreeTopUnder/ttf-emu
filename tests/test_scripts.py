#!/usr/bin/env python3
"""Unit tests for ttf-emu helper scripts and JSON schemas."""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
FIXTURES = ROOT / "tests" / "fixtures"


class KronaConverterTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmpdir = Path(tempfile.mkdtemp(prefix="ttf-emu-krona-"))
        self.addCleanup(shutil.rmtree, self.tmpdir, True)

    def test_uses_estimated_counts_and_full_lineage(self) -> None:
        out = self.tmpdir / "sample1.krona.txt"
        subprocess.check_call(
            [
                sys.executable,
                str(BIN / "emu_to_krona.py"),
                "--report",
                str(FIXTURES / "emu" / "sample1_rel-abundance.tsv"),
                "--output",
                str(out),
            ]
        )
        rows = [line.split("\t") for line in out.read_text().splitlines() if line.strip()]
        taxa = [row[-1] for row in rows]
        self.assertIn("Escherichia coli", taxa)
        self.assertIn("UNCLASSIFIED", taxa)
        ecoli = [row for row in rows if row[-1] == "Escherichia coli"][0]
        self.assertEqual(ecoli[0], "100.4")
        self.assertEqual(ecoli[1], "Bacteria")
        self.assertEqual(ecoli[2], "Proteobacteria")
        self.assertEqual(ecoli[7], "Escherichia coli")
        unmapped = [row for row in rows if row[-1] == "UNCLASSIFIED"][0]
        self.assertEqual(unmapped[0], "20")


class CombineOutputsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmpdir = Path(tempfile.mkdtemp(prefix="ttf-emu-combine-"))
        self.addCleanup(shutil.rmtree, self.tmpdir, True)

    def test_combine_species_counts_and_genus(self) -> None:
        subprocess.check_call(
            [
                sys.executable,
                str(BIN / "combine_emu_outputs.py"),
                "--reports",
                str(FIXTURES / "emu" / "sample1_rel-abundance.tsv"),
                str(FIXTURES / "emu" / "sample2_rel-abundance.tsv"),
                "--outdir",
                str(self.tmpdir),
            ]
        )
        species = self.tmpdir / "emu-combined-abundance-species.tsv"
        counts = self.tmpdir / "emu-combined-abundance-species-counts.tsv"
        genus = self.tmpdir / "emu-combined-abundance-genus.tsv"
        self.assertTrue(species.exists())
        self.assertTrue(counts.exists())
        self.assertTrue(genus.exists())

        species_rows = [line.split("\t") for line in species.read_text().splitlines()]
        self.assertEqual(species_rows[0], ["species", "sample1", "sample2"])
        by_species = {row[0]: row for row in species_rows[1:]}
        self.assertEqual(by_species["Escherichia coli"][1], "0.5")
        self.assertEqual(by_species["Escherichia coli"][2], "0.1")
        self.assertEqual(by_species["Bacteroides fragilis"][1], "0.4")
        self.assertEqual(by_species["Bacteroides fragilis"][2], "0.7")

        count_rows = [line.split("\t") for line in counts.read_text().splitlines()]
        by_counts = {row[0]: row for row in count_rows[1:]}
        self.assertEqual(by_counts["Escherichia coli"][1], "100.4")
        self.assertEqual(by_counts["Escherichia coli"][2], "20.2")
        self.assertEqual(by_counts["Bacteroides fragilis"][1], "80.1")
        self.assertEqual(by_counts["Bacteroides fragilis"][2], "140.6")

        genus_rows = [line.split("\t") for line in genus.read_text().splitlines()]
        by_genus = {row[0]: row for row in genus_rows[1:]}
        self.assertAlmostEqual(float(by_genus["Escherichia"][1]), 0.5)
        self.assertAlmostEqual(float(by_genus["Bacteroides"][1]), 0.4)
        self.assertAlmostEqual(float(by_genus["Bacteroides"][2]), 0.7)


class MetadataNormalizationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmpdir = Path(tempfile.mkdtemp(prefix="ttf-emu-md-"))
        self.addCleanup(shutil.rmtree, self.tmpdir, True)

    def test_renames_sample_column_to_id(self) -> None:
        out = self.tmpdir / "sample_metadata.tsv"
        subprocess.check_call(
            [
                sys.executable,
                str(BIN / "normalize_metadata.py"),
                "--input",
                str(FIXTURES / "metadata.tsv"),
                "--output",
                str(out),
            ]
        )
        rows = [line.split("\t") for line in out.read_text().splitlines()]
        self.assertEqual(rows[0][0], "ID")
        self.assertEqual(rows[1][0], "sample1")
        self.assertEqual(rows[2][1], "treatment")

    def test_keeps_existing_id_column(self) -> None:
        out = self.tmpdir / "sample_metadata.tsv"
        subprocess.check_call(
            [
                sys.executable,
                str(BIN / "normalize_metadata.py"),
                "--input",
                str(FIXTURES / "metadata_id.tsv"),
                "--output",
                str(out),
            ]
        )
        rows = [line.split("\t") for line in out.read_text().splitlines()]
        self.assertEqual(rows[0][0], "ID")
        self.assertEqual(rows[1], ["sample1", "control", "1"])

    def test_rejects_unknown_first_column(self) -> None:
        bad = self.tmpdir / "bad.tsv"
        bad.write_text("name\tgroup\nfoo\tbar\n")
        with self.assertRaises(subprocess.CalledProcessError):
            subprocess.check_call(
                [
                    sys.executable,
                    str(BIN / "normalize_metadata.py"),
                    "--input",
                    str(bad),
                    "--output",
                    str(self.tmpdir / "out.tsv"),
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )


class PhyloseqScriptTests(unittest.TestCase):
    def test_build_phyloseq_if_rscript_available(self) -> None:
        rscript = shutil.which("Rscript")
        if not rscript:
            self.skipTest("Rscript is not available in this environment")
        tmpdir = Path(tempfile.mkdtemp(prefix="ttf-emu-phyloseq-"))
        self.addCleanup(shutil.rmtree, tmpdir, True)
        rds = tmpdir / "emu_phyloseq.rds"
        try:
            subprocess.check_call(
                [
                    rscript,
                    str(BIN / "build_phyloseq.R"),
                    "--abundances",
                    str(FIXTURES / "emu" / "sample1_rel-abundance.tsv"),
                    str(FIXTURES / "emu" / "sample2_rel-abundance.tsv"),
                    "--metadata",
                    str(FIXTURES / "metadata.tsv"),
                    "--output",
                    str(rds),
                ]
            )
        except subprocess.CalledProcessError as exc:
            self.skipTest(f"phyloseq is not installed: {exc}")
        self.assertTrue(rds.exists())
        self.assertGreater(rds.stat().st_size, 0)

        verify = tmpdir / "verify.R"
        verify.write_text(
            f"""
            ps <- readRDS("{rds}")
            stopifnot(phyloseq::taxa_are_rows(ps))
            otu <- as(phyloseq::otu_table(ps), "matrix")
            stopifnot(all(otu == round(otu)))
            stopifnot("562" %in% phyloseq::taxa_names(ps))
            stopifnot(all(c("sample1", "sample2") %in% phyloseq::sample_names(ps)))
            stopifnot(identical(phyloseq::rank_names(ps), c("Kingdom","Phylum","Class","Order","Family","Genus","Species")))
            tax <- as(phyloseq::tax_table(ps), "matrix")
            stopifnot(tax["562", "Kingdom"] == "Bacteria")
            stopifnot(tax["562", "Species"] == "Escherichia coli")
            stopifnot(otu["562", "sample1"] == 100L)
            stopifnot(otu["562", "sample2"] == 20L)
            stopifnot(otu["817", "sample1"] == 80L)
            stopifnot(otu["817", "sample2"] == 141L)
            stopifnot(!is.null(phyloseq::sample_data(ps, errorIfNULL = FALSE)))
            stopifnot(is.null(phyloseq::phy_tree(ps, errorIfNULL = FALSE)))
            cat("ok\\n")
            """
        )
        subprocess.check_call([rscript, str(verify)])


class SchemaTests(unittest.TestCase):
    def test_schemas_are_valid_json_schema(self) -> None:
        try:
            from jsonschema.validators import validator_for
        except ImportError as exc:  # pragma: no cover
            self.skipTest(f"jsonschema is not installed: {exc}")

        for rel in ("nextflow_schema.json", "assets/schema_input.json"):
            path = ROOT / rel
            with path.open() as handle:
                schema = json.load(handle)
            validator_for(schema).check_schema(schema)
            self.assertIn("$schema", schema)
            self.assertEqual(schema["type"] in {"object", "array"}, True)

    def test_nextflow_schema_exposes_required_params(self) -> None:
        schema = json.loads((ROOT / "nextflow_schema.json").read_text())
        defs = schema.get("$defs") or schema.get("definitions")
        self.assertIsNotNone(defs)
        props = {}
        for group in defs.values():
            props.update(group.get("properties", {}))
        for key in [
            "input",
            "outdir",
            "metadata",
            "emu_db",
            "skip_chopper",
            "min_length",
            "max_length",
            "min_quality",
            "emu_min_abundance",
            "skip_krona",
            "skip_phyloseq",
        ]:
            self.assertIn(key, props)
        self.assertEqual(props["min_length"]["default"], 1000)
        self.assertEqual(props["max_length"]["default"], 1800)
        self.assertEqual(props["min_quality"]["default"], 10)
        self.assertEqual(props["emu_min_abundance"]["default"], 0.0001)
        self.assertFalse(props["skip_chopper"]["default"])
        self.assertEqual(props["emu_db"]["default"], "s3://codonyat/refdb/public/emu/default")
        self.assertIn("refdb.treetopunder.com/emu/default", props["emu_db"]["description"])
        self.assertNotIn("pending", json.dumps(props["emu_db"]).lower())
        self.assertNotIn("will be wired", json.dumps(props["emu_db"]).lower())

    def test_nextflow_config_s3_default_and_offline_test_profile(self) -> None:
        config = (ROOT / "nextflow.config").read_text()
        self.assertIn("s3://codonyat/refdb/public/emu/default", config)
        self.assertIn("refdb.treetopunder.com/emu/default", config)
        test_config = (ROOT / "conf" / "test.config").read_text()
        self.assertIn("tests/fixtures/emu_db", test_config)
        self.assertNotIn("s3://codonyat", test_config)


if __name__ == "__main__":
    os.chdir(ROOT)
    unittest.main(verbosity=2)
