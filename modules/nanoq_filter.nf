process NANOQ_FILTER {
    tag   "${meta.id}"
    label 'process_low'
    conda "bioconda::nanoq=0.10.0"
    publishDir "${params.outdir}/qc/nanoq", mode: params.publish_mode, pattern: "*.json"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.filt.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.nanoq.json"),   emit: stats
    path "versions.yml",                              emit: versions

    script:
    """
    # Pre-filter statistics, so the effect of filtering is auditable rather than assumed.
    nanoq -i ${reads} -s -j > ${meta.id}.raw.json

    nanoq \\
        -i ${reads} \\
        -l ${params.min_read_len} \\
        -q ${params.min_read_q} \\
        -o ${meta.id}.filt.fastq.gz \\
        -s -j > ${meta.id}.nanoq.json

    printf '"${task.process}"\\n    nanoq: %s\\n' "\$(nanoq --version 2>&1 | sed 's/nanoq //')" > versions.yml
    """

    stub:
    """
    echo -e "@r1\\nACGT\\n+\\nIIII" | gzip > ${meta.id}.filt.fastq.gz
    echo '{}' > ${meta.id}.nanoq.json
    touch versions.yml
    """
}
