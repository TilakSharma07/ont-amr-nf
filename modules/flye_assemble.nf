process FLYE_ASSEMBLE {
    tag   "${meta.id}"
    label 'process_high'
    conda "bioconda::flye=2.9.5"
    publishDir "${params.outdir}/assembly", mode: params.publish_mode,
               saveAs: { fn -> fn.endsWith('.fasta') || fn.endsWith('.txt') ? fn : null }

    input:
    tuple val(meta), path(reads)

    // Deliberately NOT optional. An assembly that fails must still emit a file, because
    // the negative control is EXPECTED to fail here: if this output were optional, the
    // decoy would silently drop out of the pipeline and the "negative control is empty"
    // check would never run — the pipeline would look clean by omitting its own control.
    // A zero-sequence FASTA carries that failure forward to the validation gate instead.
    output:
    tuple val(meta), path("${meta.id}.assembly.fasta"),    emit: assembly
    tuple val(meta), path("${meta.id}.flye_status.txt"),   emit: status
    path "${meta.id}.assembly_info.txt",                   emit: info, optional: true
    path "versions.yml",                                   emit: versions

    script:
    """
    # --nano-hq is the correct Flye mode for R10.4.1 chemistry with SUP basecalling
    # (<=5% error). Using --nano-raw here would be a silent accuracy regression.
    #
    # A failed assembly must not kill the run: the negative control is EXPECTED to
    # assemble poorly or not at all, and that outcome is itself a result. We capture
    # the status and let the validation gate decide what it means.
    set +e
    flye \\
        ${params.flye_mode} ${reads} \\
        --genome-size ${params.genome_size} \\
        --iterations ${params.flye_iterations} \\
        --threads ${task.cpus} \\
        --out-dir flye_out
    FLYE_RC=\$?
    set -e

    if [ -s flye_out/assembly.fasta ]; then
        cp flye_out/assembly.fasta      ${meta.id}.assembly.fasta
        cp flye_out/assembly_info.txt   ${meta.id}.assembly_info.txt
        printf 'sample_id\\tflye_status\\tflye_exit\\n${meta.id}\\tassembled\\t%s\\n' "\$FLYE_RC" \\
            > ${meta.id}.flye_status.txt
    else
        # Empty-but-present FASTA: keeps this sample on the single shared code path so
        # its failure is reported by the validation gate rather than hidden by a
        # missing channel item.
        : > ${meta.id}.assembly.fasta
        printf 'sample_id\\tflye_status\\tflye_exit\\n${meta.id}\\tno_assembly\\t%s\\n' "\$FLYE_RC" \\
            > ${meta.id}.flye_status.txt
        echo "flye produced no assembly (exit \$FLYE_RC) for ${meta.id}" >&2
    fi

    printf '"${task.process}"\\n    flye: %s\\n' "\$(flye --version 2>&1)" > versions.yml
    """

    stub:
    """
    if [ "${meta.role}" = "negative_control" ]; then
        : > ${meta.id}.assembly.fasta
        printf 'sample_id\\tflye_status\\tflye_exit\\n${meta.id}\\tno_assembly\\t1\\n' > ${meta.id}.flye_status.txt
    else
        printf '>ctg1\\nACGTACGTACGT\\n' > ${meta.id}.assembly.fasta
        printf 'sample_id\\tflye_status\\tflye_exit\\n${meta.id}\\tassembled\\t0\\n' > ${meta.id}.flye_status.txt
    fi
    touch ${meta.id}.assembly_info.txt versions.yml
    """
}
