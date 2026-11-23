#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TTF_EMU } from './workflows/emu'

def helpMessage() {
    log.info(
        """
        ttf-emu ${workflow.manifest.version}
        Full-Length 16S Profiling (Emu) for TreeTopFlow

        Research use only. Not for clinical or diagnostic use.

        Required:
          --input             Samplesheet CSV with columns sample,fastq_1 (fastq_2 optional/ignored)
          --outdir            Output directory
          --emu_db            Emu database dir or .tar/.tar.gz (default: s3://codonyat/refdb/public/emu/default)

        Optional:
          --metadata          Sample metadata TSV (first column ID or sample)
          --skip_chopper      Skip Chopper length/quality filtering (default: false)
          --min_length        Chopper minimum read length (default: 1000)
          --max_length        Chopper maximum read length (default: 1800)
          --min_quality       Chopper minimum average quality (default: 10)
          --emu_min_abundance Emu --min-abundance (default: 0.0001)
          --skip_krona        Skip Krona plots (default: false)
          --skip_phyloseq     Skip phyloseq RDS (default: false)
        """.stripIndent()
    )
}

def isArchive(path) {
    return RefPaths.isArchive(path)
}

def asRefFile(path) {
    def value = RefPaths.normalize(path)
    if (!value) {
        return null
    }
    // Unpacked dirs (local or s3://, trailing slash optional) stage as directories.
    return isArchive(value) ? file(value) : file(value, type: 'dir')
}

def resolveInputPath(raw) {
    if (raw == null) {
        return null
    }
    def value = raw.toString().trim()
    if (!value) {
        return null
    }
    def candidate = file(value)
    if (candidate.exists()) {
        return candidate
    }
    def fromProject = file("${workflow.projectDir}/${value}")
    if (fromProject.exists()) {
        return fromProject
    }
    return candidate
}

def createFastqChannel(LinkedHashMap row) {
    def meta = [:]
    meta.id = row.sample?.toString()?.trim()
    if (!meta.id) {
        error "Samplesheet rows must include a non-empty 'sample' value"
    }

    def fastq_1 = resolveInputPath(row.fastq_1)
    if (!fastq_1) {
        error "Sample ${meta.id} is missing fastq_1"
    }

    // Optional fastq_2 is accepted for schema uniformity and ignored:
    // Emu profiles single-end full-length ONT / PacBio HiFi 16S reads.
    meta.single_end = true
    return [meta, fastq_1]
}

workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }

    if (!params.input) {
        error "Missing required parameter --input (samplesheet CSV)"
    }
    if (!params.outdir) {
        error "Missing required parameter --outdir"
    }
    if (!params.emu_db) {
        error "Missing required parameter --emu_db"
    }

    ch_input = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> createFastqChannel(row) }

    ch_metadata = params.metadata ? Channel.value(file(params.metadata, checkIfExists: true)) : Channel.value([])

    def emu_db_file = asRefFile(params.emu_db)
    ch_emu_db_archive = isArchive(params.emu_db) ? Channel.of([[id: 'emu_db'], emu_db_file]) : Channel.empty()
    ch_emu_db_dir     = isArchive(params.emu_db) ? Channel.empty() : Channel.value(emu_db_file)

    TTF_EMU(
        ch_input,
        ch_metadata,
        ch_emu_db_archive,
        ch_emu_db_dir
    )
}

workflow.onComplete {
    log.info(workflow.success ? "ttf-emu completed successfully (research use only)" : "ttf-emu failed")
}
