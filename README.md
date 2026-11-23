# ttf-emu

Full-Length 16S Profiling (Emu) — a lean Nextflow DSL2 pipeline for TreeTopFlow.

**Research use only. Not for clinical or diagnostic use.**

Version `0.1.0` profiles full-length ONT / PacBio HiFi 16S reads with Emu only: NanoPlot QC, optional Chopper length/quality filtering, Emu abundance, combined tables, optional phyloseq RDS (ampliseq v2.0.0 notebook-compatible) and Krona plots, and MultiQC.

## Inputs

`--input` samplesheet CSV (required):

```csv
sample,fastq_1,fastq_2
SAMPLE_ONT,/path/to/SAMPLE_ONT.fastq.gz,
SAMPLE_HIFI,/path/to/SAMPLE_HIFI.fastq.gz,
```

- Single-end long reads only (ONT 16S or PacBio HiFi 16S)
- `fastq_2` is optional and ignored; it is accepted so the samplesheet schema matches other TreeTopFlow pipelines
- Paths may be absolute or relative to the launch directory / pipeline root

Optional `--metadata` TSV for phyloseq `sample_data` and a published `sample_metadata.tsv`. The first column must be `ID` or `sample` and must match samplesheet sample IDs:

```tsv
sample	condition	batch
SAMPLE_ONT	control	1
SAMPLE_HIFI	treatment	1
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | Samplesheet CSV |
| `--outdir` | required | Output directory (local or `s3://`) |
| `--metadata` | `null` | Optional sample metadata TSV |
| `--emu_db` | `s3://codonyat/refdb/public/emu/default` | Unpacked Emu default 16S directory (or `.tar`/`.tar.gz`). HTTPS: `https://refdb.treetopunder.com/emu/default/` |
| `--skip_chopper` | `false` | Skip Chopper length/quality filtering |
| `--min_length` | `1000` | Chopper `--minlength` |
| `--max_length` | `1800` | Chopper `--maxlength` |
| `--min_quality` | `10` | Chopper `--quality` |
| `--emu_min_abundance` | `0.0001` | Emu `--min-abundance` |
| `--skip_krona` | `false` | Skip Krona plots |
| `--skip_phyloseq` | `false` | Skip phyloseq RDS |

`--emu_db` may be an unpacked directory (local or `s3://`, trailing slash optional) or a `.tar` / `.tar.gz` archive (unpacked with `nf-core/untar`). Unpacked directories skip `UNTAR` and are staged as directory channels so Emu sees `species_taxid.fasta` and `taxonomy.tsv`.

## Outputs

```
outdir/
  qc/nanoplot/<sample>/*.txt
  qc/nanoplot/<sample>/*.html
  emu/*_rel-abundance.tsv
  emu/emu-combined-abundance-species.tsv
  emu/emu-combined-abundance-species-counts.tsv
  emu/emu-combined-abundance-genus.tsv
  phyloseq/emu_phyloseq.rds
  sample_metadata.tsv
  krona/krona_combined.html
  krona/<sample>.krona.html
  multiqc/multiqc_report.html
  multiqc/multiqc_data/multiqc_general_stats.txt
  pipeline_info/software_versions.yml
```

`phyloseq/emu_phyloseq.rds` is compatible with the `nf-core/ampliseq` v2.0.0 hosted RStudio notebook (`ttf-rstudio-ampliseq`):

- `otu_table`: integer estimated read counts (`taxa_are_rows = TRUE`, rounded from Emu `estimated counts`, rows = `tax_id`, columns = sample ID)
- `tax_table`: `Kingdom, Phylum, Class, Order, Family, Genus, Species` (Emu `superkingdom` mapped to `Kingdom`)
- `sample_data`: attached from `--metadata` when provided
- `phy_tree`: omitted (`ttfampliseq` v2.0.0 falls back from weighted UniFrac to Bray-Curtis)

## Reference database

The default is an unpacked prefix on S3 (`s3://codonyat/refdb/public/emu/default`), MD5-verified and mirrored at `https://refdb.treetopunder.com`. Override with `--emu_db` for a local copy or a `.tar`/`.tar.gz` archive. `-profile test` points `--emu_db` at `tests/fixtures/emu_db` so stub-runs stay offline.

| Resource | Default S3 prefix | HTTPS mirror | Contents |
|---|---|---|---|
| Emu default 16S | `s3://codonyat/refdb/public/emu/default` | https://refdb.treetopunder.com/emu/default/ | Unpacked dir: `species_taxid.fasta`, `taxonomy.tsv` |

## Run

```bash
nextflow run . \
  --input samplesheet.csv \
  --outdir results \
  --metadata metadata.tsv \
  -profile docker
```

The S3 default above is used unless you pass `--emu_db`.

AWS Batch (TreeTopFlow controllers pass the worker config with `-c`; do not require `--aws_queue`):

```bash
nextflow run . \
  -profile awsbatch \
  -c /path/to/worker/awsbatch.config \
  -work-dir s3://bucket/work \
  --input samplesheet.csv \
  --outdir s3://bucket/results
```

Optional CLI overrides (`--aws_queue`, `--aws_region`, `--aws_cli_path`) are applied only when set, so they do not clobber the worker `-c` file.

Stub-run (DAG compile, no containers required):

```bash
nextflow run . -profile test -stub-run --outdir results_chopper
nextflow run . -profile test -stub-run --skip_chopper true --outdir results_nochopper
```

Helper-script tests:

```bash
python tests/test_scripts.py
```

## Modules

Standard `nf-core/modules` (see `modules.json`):

- `nf-core/nanoplot`
- `nf-core/chopper`
- `nf-core/untar`
- `nf-core/emu/abundance` (container `emu=3.5.1`)
- `nf-core/krona/ktimporttext`
- `nf-core/multiqc`

Local modules: `EMU_COMBINE_OUTPUTS` (same Emu BioContainer), `EMU_TO_KRONA`, `BUILD_PHYLOSEQ`, `NORMALIZE_METADATA`.

## License

MIT. nf-core modules retain their upstream MIT license. Emu, Chopper, NanoPlot, Krona, MultiQC, and phyloseq retain their respective licenses.
