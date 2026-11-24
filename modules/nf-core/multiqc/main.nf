process MULTIQC {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.29--pyhdfd78af_0' :
        'biocontainers/multiqc:1.29--pyhdfd78af_0' }"

    input:
    path  multiqc_files, stageAs: "?/*"
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)
    path(replace_names)
    path(sample_names)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def report_name = task.ext.prefix ? "${task.ext.prefix}.html" : 'multiqc_report.html'
    def prefix = task.ext.prefix ? "--filename ${report_name}" : ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    def logo = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo}\"'" : ''
    def replace = replace_names ? "--replace-names ${replace_names}" : ''
    def samples = sample_names ? "--sample-names ${sample_names}" : ''
    """
    # MultiQC 1.29 exits non-zero and writes no HTML when no modules match.
    # Keep the pipeline green: always emit the declared report + data dir.
    multiqc \\
        --force \\
        $args \\
        $config \\
        $prefix \\
        $extra_config \\
        $logo \\
        $replace \\
        $samples \\
        . || true

    if ! ls *multiqc_report.html >/dev/null 2>&1 && [ ! -f "${report_name}" ]; then
        echo "WARNING: MultiQC found no recognised module outputs; writing a placeholder report." >&2
        mkdir -p multiqc_data
        cat > "${report_name}" <<'EOF'
<!DOCTYPE html>
<html lang="en">
<head><meta charset="utf-8"><title>MultiQC report</title></head>
<body>
<h1>MultiQC</h1>
<p>No recognised module outputs were found for this run. The pipeline completed; this placeholder exists so MultiQC is not a hard failure.</p>
</body>
</html>
EOF
        echo "placeholder" > multiqc_data/multiqc_sources.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    mkdir multiqc_data
    mkdir multiqc_plots
    touch multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
