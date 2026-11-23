process EMU_COMBINE_OUTPUTS {
    tag "combine"
    label 'process_low'

    conda "${projectDir}/modules/nf-core/emu/abundance/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/emu:3.5.1--hdfd78af_0'
        : 'biocontainers/emu:3.5.1--hdfd78af_0'}"

    input:
    path reports

    output:
    path "emu-combined-abundance-species.tsv"       , emit: species
    path "emu-combined-abundance-species-counts.tsv", emit: species_counts
    path "emu-combined-abundance-genus.tsv"         , emit: genus
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p emu_reports
    for report in ${reports}; do
        cp "\${report}" emu_reports/
    done

    emu combine-outputs emu_reports species --split-tables
    emu combine-outputs emu_reports species --split-tables --counts
    emu combine-outputs emu_reports genus --split-tables

    cp emu_reports/emu-combined-abundance-species.tsv .
    cp emu_reports/emu-combined-abundance-species-counts.tsv .
    cp emu_reports/emu-combined-abundance-genus.tsv .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emu: \$(echo \$(emu --version 2>&1) | sed 's/^.*emu //; s/Using.*\$//' )
    END_VERSIONS
    """

    stub:
    """
    touch emu-combined-abundance-species.tsv
    touch emu-combined-abundance-species-counts.tsv
    touch emu-combined-abundance-genus.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emu: \$(echo \$(emu --version 2>&1) | sed 's/^.*emu //; s/Using.*\$//' )
    END_VERSIONS
    """
}
