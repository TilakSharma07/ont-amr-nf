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
    // Quote the database path. It is user-supplied and routinely contains spaces —
    // an external volume ("/media/user/MY DRIVE/db"), a macOS "Macintosh HD" path, a
    // Windows "Program Files" path. Unquoted, the shell splits it on whitespace and
    // amrfinder rejects the fragment as a positional parameter. This is invisible on
    // a developer machine whose paths happen to have no spaces.
    def db_arg  = params.amrfinder_db ? "--database '${params.amrfinder_db}'" : ""
    // --ident_min is only passed when the user asks for a non-default floor. This is not
    // cosmetic. In amrfinder.cpp the flag reaches the reporting engine only when the
    // value is not -1:
    //     (ident == -1 ? noString : "  -ident_min " + toString (ident))
    // and the option's own help reads "-1 means use a curated threshold if it exists and
    // 0.9 otherwise". So ANY explicit value -- 0.9 included -- replaces the per-gene
    // curated cutoffs with one flat number for every reference in the database. Those
    // cutoffs are what separate closely-related alleles in families like blaOXA, tet and
    // qnr, where a few percent of identity is the difference between two enzymes with
    // different substrate spectra. Passing 0.9 looks like "keep the default" and is not.
    def ident_arg = (params.amr_min_ident == null || params.amr_min_ident < 0)
                    ? "" : "--ident_min ${params.amr_min_ident}"
    """
    # Resolve the database version FIRST, so it is recorded on every path through this
    # process — including the empty-assembly path below. It must be queried WITH
    # --database, or amrfinder looks in its default location, finds nothing, and the
    # version silently becomes a placeholder. An AMR call is only interpretable against
    # a known database version, so an unidentifiable database is a hard failure rather
    # than a note in a log.
    db_version=\$(amrfinder --database_version ${db_arg} 2>&1 | grep -oP 'Database version: \\K.*' || true)
    if [ -z "\$db_version" ]; then
        echo "ERROR: could not determine the AMRFinderPlus database version for ${meta.id}" >&2
        echo "  (queried with: amrfinder --database_version ${db_arg})" >&2
        amrfinder --database_version ${db_arg} >&2 2>&1 || true
        exit 1
    fi

    write_versions() {
        printf '"${task.process}":\\n    amrfinderplus: %s\\n    amrfinderplus_db: %s\\n' \\
            "\$(amrfinder --version)" "\$db_version" > versions.yml
    }

    # An empty assembly reaches here when Flye produced nothing (expected for the
    # negative control). Emit a header-only result rather than failing: "no assembly,
    # therefore no calls" is a real outcome that the validation gate must see, and it
    # is different from "assembly present, no calls found".
    if [ ! -s "${assembly}" ]; then
        printf 'Name\\tElement symbol\\tElement name\\tClass\\tSubclass\\t%% Coverage of reference\\t%% Identity to reference\\n' \\
            > ${meta.id}.amrfinder.tsv
        echo "no assembly for ${meta.id}: emitted header-only call set" > ${meta.id}.amrfinder.log
        write_versions
        exit 0
    fi

    # The negative control is deliberately run WITHOUT --organism: no organism claim
    # can be made about shuffled sequence, and we want the least-filtered, most
    # permissive call mode applied to it. If anything survives there, we want to see it.
    amrfinder \\
        --nucleotide "${assembly}" \\
        ${org_arg} \\
        ${db_arg} \\
        ${ident_arg} \\
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

    write_versions
    """

    stub:
    """
    printf 'Name\\tElement symbol\\tElement name\\tClass\\tSubclass\\n' > ${meta.id}.amrfinder.tsv
    touch ${meta.id}.amrfinder.log versions.yml
    """
}
