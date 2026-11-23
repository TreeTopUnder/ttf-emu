process NORMALIZE_METADATA {
    tag "metadata"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    path metadata

    output:
    path "sample_metadata.tsv", emit: tsv
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    normalize_metadata.py \\
        --input ${metadata} \\
        --output sample_metadata.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    if [ -f "${metadata}" ]; then
        normalize_metadata.py --input ${metadata} --output sample_metadata.tsv
    else
        printf "ID\\tcondition\\n" > sample_metadata.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
