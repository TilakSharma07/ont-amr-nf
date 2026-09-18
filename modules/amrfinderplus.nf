process AMRFINDERPLUS {
    tag   "${meta.id}"
    label 'process_medium'
    conda "bioconda::ncbi-amrfinderplus=4.0.19"
    publishDir "${params.outdir}/amr", mode: params.publish_mode

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}.amrfinder.tsv"), emit: calls
    path "${meta.id}.amrfinder.log",                   emit: log
    path "versions.yml",                               emit: versions

    script:
    // Organism-aware calling matters: with -O, AMRFinderPlus applies species-specific
    // point-mutation rules and suppresses core genes that are intrinsic to that
    // species rather than acquired resistance. Without it, intrinsic chromosomal
    // genes get reported as findings and the call set is inflated.
    def org = [
        'Klebsiella pneumoniae': 'Klebsiella_pneumoniae',
        'Escherichia coli'     : 'Escherichia'
    ].get(meta.species, null)
    def org_arg = (org && meta.role != 'negative_control') ? "--organism ${org}" : ""
    def db_arg  = params.amrfinder_db ? "--database ${params.amrfinder_db}" : ""
    """
    # An empty assembly reaches here when Flye produced nothing (expected for the
    # negative control). Emit a header-only result rather than failing: "no assembly,
    # therefore no calls" is a real outcome that the validation gate must see, and it
    # is different from "assembly present, no calls found".
    if [ ! -s ${assembly} ]; then
        printf 'Name\\tElement symbol\\tElement name\\tClass\\tSubclass\\t%% Coverage of reference\\t%% Identity to reference\\n' \\
            > ${meta.id}.amrfinder.tsv
        echo "no assembly for ${meta.id}: emitted header-only call set" > ${meta.id}.amrfinder.log
        printf '"%s":\\n    amrfinderplus: %s\\n' "${task.process}" "\$(amrfinder --version)" > versions.yml
        exit 0
    fi

    # The negative control is deliberately run WITHOUT --organism: no organism claim
    # can be made about shuffled sequence, and we want the least-filtered, most
    # permissive call mode applied to it. If anything survives there, we want to see it.
    amrfinder \\
        --nucleotide ${assembly} \\
        ${org_arg} \\
        ${db_arg} \\
        --ident_min ${params.amr_min_ident} \\
        --coverage_min ${params.amr_min_cov} \\
        --threads ${task.cpus} \\
        --name ${meta.id} \\
        --plus \\
        > ${meta.id}.amrfinder.tsv 2> ${meta.id}.amrfinder.log

    # Fail loudly on an empty file (no header at all) — that is a tool error.
    # Zero DATA rows is a legitimate scientific result and must NOT fail.
    if [ ! -s ${meta.id}.amrfinder.tsv ]; then
        echo "ERROR: amrfinder produced no output for ${meta.id}" >&2
        exit 1
    fi

    printf '"${task.process}"\\n    amrfinderplus: %s\\n    amrfinderplus_db: %s\\n' "\$(amrfinder --version)" "\$(amrfinder --database_version 2>&1 | grep -oP 'Database version: \\K.*' || echo 'bundled')" > versions.yml
    """

    stub:
    """
    printf 'Name\\tElement symbol\\tElement name\\tClass\\tSubclass\\n' > ${meta.id}.amrfinder.tsv
    touch ${meta.id}.amrfinder.log versions.yml
    """
}
