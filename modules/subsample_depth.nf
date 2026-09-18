process SUBSAMPLE_DEPTH {
    tag   "${meta.id}_${target_depth}x_rep${replicate}"
    label 'process_low'
    conda "bioconda::seqkit=2.8.2"

    input:
    tuple val(meta), path(reads), val(target_depth), val(replicate)

    // The output name embeds depth and replicate so it cannot collide with the staged
    // input file, which is named after the parent sample.
    output:
    tuple val(meta), val(target_depth), val(replicate),
          path("${meta.id}_d${target_depth}_r${replicate}.fastq.gz"), emit: reads

    script:
    def seed = (params.titration_seed as int) + (replicate as int) * 1000 + (target_depth as int)
    def out  = "${meta.id}_d${target_depth}_r${replicate}.fastq.gz"
    """
    # Subsample to a TARGET DEPTH, not a target read count: read-length distributions
    # differ between isolates, so a fixed read count would mean different coverage per
    # sample and the titration curve would not be comparable.
    #
    # The seed is derived from (replicate, depth) so every point is independently
    # resampled but the whole experiment re-runs identically.
    #
    # Sampling is proportional (Bernoulli per read) rather than "shuffle and take a
    # prefix": shuffling holds the whole file in memory, and closing the pipe early
    # kills the upstream process with SIGPIPE. Proportional sampling is single-pass
    # and preserves the read-length distribution in expectation. The depth this
    # actually achieves is not assumed — ASSEMBLY_STATS measures it by remapping, and
    # that realised value is what the titration table reports.
    GENOME_BP=\$(python3 -c "
    s = '${params.genome_size}'.lower().strip()
    mult = {'k': 1e3, 'm': 1e6, 'g': 1e9}
    print(int(float(s[:-1]) * mult[s[-1]]) if s[-1] in mult else int(float(s)))
    ")

    TARGET_BASES=\$(( ${target_depth} * GENOME_BP ))
    # -F'\\t' is load-bearing: seqkit's table is tab-separated, and a staged path
    # containing a space would shift the columns under awk's default splitting.
    TOTAL_BASES=\$(seqkit stats -T ${reads} | awk -F'\\t' 'NR==2{print \$5}')

    PROP=\$(python3 -c "
    t, a = ${target_depth} * \$GENOME_BP, \$TOTAL_BASES
    print(f'{min(1.0, t / a):.6f}' if a > 0 else '1.0')
    ")

    if [ "\$PROP" = "1.000000" ]; then
        # Requested depth is at or above what the run actually contains. Copying rather
        # than silently returning a shallower set keeps this point on the curve, and the
        # measured depth will show it sitting below its nominal target.
        echo "WARNING: ${meta.id} has only \$TOTAL_BASES bp; ${target_depth}x needs \$TARGET_BASES bp. Using all reads." >&2
        cp ${reads} ${out}
    else
        seqkit sample -s ${seed} -p \$PROP ${reads} 2> sample.log | gzip > ${out}
    fi

    ACTUAL=\$(seqkit stats -T ${out} | awk -F'\\t' 'NR==2{print \$5}')
    echo "${meta.id} d${target_depth} r${replicate}: target \$TARGET_BASES bp, prop \$PROP, got \$ACTUAL bp" >&2
    """

    stub:
    """
    printf '@s1\\nACGT\\n+\\nIIII\\n' | gzip > ${meta.id}_d${target_depth}_r${replicate}.fastq.gz
    """
}
