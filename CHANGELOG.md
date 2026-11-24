# Changelog

All notable changes to this project are documented in this file.

## 0.1.1 - 2026-11-24

### Fixed

- **MultiQC** — `assets/multiqc_config.yml` listed `run_modules: [nanoplot]`. MultiQC 1.29 (container `biocontainers/multiqc:1.29--pyhdfd78af_0`) has no module named `nanoplot`; NanoPlot `NanoStats.txt` is parsed by **`nanostat`** ([MultiQC v1.29 `modules/nanostat`](https://github.com/MultiQC/MultiQC/blob/v1.29/multiqc/modules/nanostat/nanostat.py), [`search_patterns.yaml`](https://github.com/MultiQC/MultiQC/blob/v1.29/multiqc/search_patterns.yaml), [docs](https://github.com/MultiQC/MultiQC/blob/v1.29/docs/markdown/modules/nanostat.md)). The unknown module caused `No analysis results found` and a missing `*multiqc_report.html`, which failed the pipeline. There is no Chopper MultiQC module (`pychopper` / `porechop` are different tools).
- **MultiQC inputs** — NanoPlot `NanoStats.txt` is still passed into `MULTIQC`; Chopper `*.log` is staged into the MultiQC work directory. Search pattern `*NanoStats.txt` is set explicitly.
- **MultiQC empty runs** — if MultiQC finds nothing it now writes a placeholder HTML report and `multiqc_data/` instead of failing the whole pipeline.
- **Emu threads / resources** — `EMU_ABUNDANCE` already passed `--threads ${task.cpus}`, but `withName` overrides in `conf/base.config`, `conf/modules.config`, `conf/test.config`, and `conf/awsbatch.config` pinned the process to 4 vCPU / 8 GB (so production ran `--threads 4`). Production pins are removed; the process uses `label: process_high` (16 CPUs, 32 GB, 12 h). `-profile test` still caps local stub-runs at 2 CPUs / 4 GB.
- **Chopper defaults** — filter flags were already passed via `ext.args2` (`--minlength` / `--maxlength` / `--quality`). The previous `--minlength 1000` default is too permissive for typical ONT 16S libraries (production kept 100% of ~59k reads). Defaults are now **Q10, 1200–1800 bp** for full-length 16S. Chopper stderr (including `Kept N reads out of M reads`) is published under `qc/chopper/<sample>/`.

### Added

- `--emu_max_secondary` (default `50`) maps to Emu `--N` (max minimap2 secondary alignments per read). Lowering it (e.g. 20) trades some EM re-estimation evidence for speed and RAM; 50 is Emu's own default.
- `--emu_keep_read_assignments` (default `false`). Phyloseq and Krona only need `*_rel-abundance.tsv` from `--keep-counts`; the per-read assignment TSV is large and slow (especially on diverse mocks).

## 0.1.0 - 2026-11-23

Initial TreeTopFlow pipeline: NanoPlot, Chopper, Emu abundance, combined tables, phyloseq, Krona, MultiQC.
