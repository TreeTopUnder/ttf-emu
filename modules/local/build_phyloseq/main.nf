process BUILD_PHYLOSEQ {
    tag "phyloseq"
    label 'process_low'

    conda "bioconda::bioconductor-phyloseq=1.46.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-phyloseq:1.46.0--r43hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-phyloseq:1.46.0--r43hdfd78af_0' }"

    input:
    path abundances
    path metadata

    output:
    path "emu_phyloseq.rds", emit: rds
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def metadata_arg = metadata ? "--metadata ${metadata}" : ""
    """
    build_phyloseq.R \\
        --abundances ${abundances} \\
        ${metadata_arg} \\
        --output emu_phyloseq.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -n 1 | sed 's/R version //; s/ .*//')
        phyloseq: \$(Rscript -e 'cat(as.character(packageVersion("phyloseq")))')
    END_VERSIONS
    """

    stub:
    """
    touch emu_phyloseq.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: 4.3.0
        phyloseq: 1.46.0
    END_VERSIONS
    """
}
