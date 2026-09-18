process ASSEMBLY_STATS {
    tag   "${meta.id}"
    label 'process_low'
    conda "bioconda::seqkit=2.8.2 bioconda::minimap2=2.28 bioconda::samtools=1.21"
    publishDir "${params.outdir}/assembly_stats", mode: params.publish_mode

    input:
    tuple val(meta), path(assembly), path(reads)

    output:
    tuple val(meta), path("${meta.id}.stats.tsv"), emit: stats
    path "versions.yml",                          emit: versions

    script:
    """
    # No assembly (expected for the negative control) still produces a stats row, so
    # the sample stays visible to the validation gate with zeros rather than vanishing.
    if [ ! -s ${assembly} ]; then
        printf 'sample_id\\trole\\tspecies\\tn_contigs\\ttotal_len\\tn50\\tlargest\\tgc_percent\\tmean_depth\\tbreadth_1x\\n' \\
            > ${meta.id}.stats.tsv
        printf '${meta.id}\\t${meta.role}\\t${meta.species}\\t0\\t0\\t0\\t0\\tNA\\t0.00\\t0.0000\\n' \\
            >> ${meta.id}.stats.tsv
        touch versions.yml
        exit 0
    fi

    # Contiguity from the assembly itself ...
    seqkit stats -a -T ${assembly} > seqkit.tsv

    # ... and realised depth by mapping the reads BACK onto their own assembly.
    # Depth must be measured, not inferred from input yield / expected genome size:
    # filtered-out reads and unassembled content both break that assumption.
    minimap2 -ax map-ont -t ${task.cpus} ${assembly} ${reads} 2> minimap2.log \\
      | samtools sort -@ ${task.cpus} -o aln.bam -
    samtools index aln.bam
    samtools depth -a aln.bam > depth.txt

    python3 - <<'PY'
    import csv

    with open("seqkit.tsv") as fh:
        s = list(csv.DictReader(fh, delimiter="\\t"))[0]

    tot = cov = 0
    n = 0
    with open("depth.txt") as fh:
        for line in fh:
            d = int(line.rsplit("\\t", 1)[1])
            tot += d
            n += 1
            if d > 0:
                cov += 1

    mean_depth = tot / n if n else 0.0
    breadth    = cov / n if n else 0.0

    out = {
        "sample_id":   "${meta.id}",
        "role":        "${meta.role}",
        "species":     "${meta.species}",
        "n_contigs":   s["num_seqs"],
        "total_len":   s["sum_len"],
        "n50":         s["N50"],
        "largest":     s["max_len"],
        "gc_percent":  s.get("GC(%)", "NA"),
        "mean_depth":  f"{mean_depth:.2f}",
        "breadth_1x":  f"{breadth:.4f}",
    }
    with open("${meta.id}.stats.tsv", "w") as fh:
        fh.write("\\t".join(out) + "\\n")
        fh.write("\\t".join(str(v) for v in out.values()) + "\\n")
    PY

    printf '"${task.process}"\\n    seqkit: %s\\n    minimap2: %s\\n    samtools: %s\\n' "\$(seqkit version | sed 's/seqkit v//')" "\$(minimap2 --version)" "\$(samtools --version | head -1 | sed 's/samtools //')" > versions.yml
    """

    stub:
    """
    printf 'sample_id\\trole\\n${meta.id}\\t${meta.role}\\n' > ${meta.id}.stats.tsv
    touch versions.yml
    """
}
