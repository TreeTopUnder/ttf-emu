include { NANOPLOT } from '../modules/nf-core/nanoplot/main'
include { CHOPPER } from '../modules/nf-core/chopper/main'
include { UNTAR } from '../modules/nf-core/untar/main'
include { EMU_ABUNDANCE } from '../modules/nf-core/emu/abundance/main'
include { KRONA_KTIMPORTTEXT as KRONA_PER_SAMPLE } from '../modules/nf-core/krona/ktimporttext/main'
include { KRONA_KTIMPORTTEXT as KRONA_COMBINED } from '../modules/nf-core/krona/ktimporttext/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { EMU_COMBINE_OUTPUTS } from '../modules/local/emu_combine_outputs/main'
include { EMU_TO_KRONA } from '../modules/local/emu_to_krona/main'
include { BUILD_PHYLOSEQ } from '../modules/local/build_phyloseq/main'
include { NORMALIZE_METADATA } from '../modules/local/normalize_metadata/main'

def isArchivePath(path) {
    // Unpacked dirs (local or s3://, trailing slash optional) skip UNTAR.
    return RefPaths.isArchive(path)
}

workflow TTF_EMU {
    take:
    ch_input
    ch_metadata
    ch_emu_db_archive
    ch_emu_db_dir

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    NANOPLOT(ch_input)
    ch_versions      = ch_versions.mix(NANOPLOT.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.txt.collect { it[1] })

    if (!params.skip_chopper) {
        CHOPPER(ch_input, [])
        ch_versions = ch_versions.mix(CHOPPER.out.versions)
        ch_reads    = CHOPPER.out.fastq
    } else {
        ch_reads = ch_input
    }

    if (isArchivePath(params.emu_db)) {
        UNTAR(ch_emu_db_archive)
        ch_emu_db   = UNTAR.out.untar.map { _meta, dir -> dir }
        ch_versions = ch_versions.mix(UNTAR.out.versions)
    } else {
        ch_emu_db = ch_emu_db_dir
    }

    ch_emu_db_value = ch_emu_db.collect().map { dirs -> dirs[0] }

    EMU_ABUNDANCE(
        ch_reads,
        ch_emu_db_value
    )
    ch_versions = ch_versions.mix(EMU_ABUNDANCE.out.versions)

    ch_abundance_files = EMU_ABUNDANCE.out.report
        .map { _meta, report -> report }
        .collect()

    EMU_COMBINE_OUTPUTS(ch_abundance_files)
    ch_versions = ch_versions.mix(EMU_COMBINE_OUTPUTS.out.versions)

    if (params.metadata) {
        NORMALIZE_METADATA(ch_metadata)
        ch_versions = ch_versions.mix(NORMALIZE_METADATA.out.versions)
    }

    if (!params.skip_phyloseq) {
        BUILD_PHYLOSEQ(
            ch_abundance_files,
            ch_metadata
        )
        ch_versions = ch_versions.mix(BUILD_PHYLOSEQ.out.versions)
    }

    if (!params.skip_krona) {
        EMU_TO_KRONA(EMU_ABUNDANCE.out.report)
        ch_versions = ch_versions.mix(EMU_TO_KRONA.out.versions)

        KRONA_PER_SAMPLE(EMU_TO_KRONA.out.txt)
        ch_versions = ch_versions.mix(KRONA_PER_SAMPLE.out.versions)

        ch_krona_combined = EMU_TO_KRONA.out.txt
            .map { _meta, txt -> txt }
            .collect()
            .map { files -> [[id: 'krona_combined'], files] }
        KRONA_COMBINED(ch_krona_combined)
        ch_versions = ch_versions.mix(KRONA_COMBINED.out.versions)
    }

    ch_collated_versions = ch_versions
        .map { path -> path.text }
        .unique()
        .mix(
            Channel.of(
                """
                Workflow:
                    ttf-emu: ${workflow.manifest.version}
                    nextflow: ${workflow.nextflow.version}
                """.stripIndent()
            )
        )
        .collectFile(
            name: 'software_versions.yml',
            storeDir: "${params.outdir}/pipeline_info",
            sort: true,
            newLine: true
        )
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    MULTIQC(
        ch_multiqc_files.collect(),
        file("${projectDir}/assets/multiqc_config.yml"),
        [],
        [],
        [],
        []
    )

    emit:
    versions = ch_versions
}
