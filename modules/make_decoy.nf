process MAKE_DECOY {
    tag   "decoy_from:${meta.id}"
    label 'process_low'
    conda "conda-forge::python=3.11"
    publishDir "${params.outdir}/controls", mode: params.publish_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("NEG_DECOY.fastq.gz"), emit: reads
    path "decoy_composition.tsv",                emit: report

    script:
    """
    #!/usr/bin/env python3
    import gzip, random, collections

    random.seed(1234)

    # A negative control is only meaningful if it is hard to distinguish from real input
    # on any axis EXCEPT gene content. So: keep the real read-length distribution, keep
    # per-read base composition, and destroy only the sequence ORDER. A caller that
    # reports AMR genes here is matching on composition, not on genes.
    n_in = n_out = 0
    comp_in  = collections.Counter()
    comp_out = collections.Counter()

    with gzip.open("${reads}", "rt") as fh, gzip.open("NEG_DECOY.fastq.gz", "wt") as out:
        while True:
            h = fh.readline()
            if not h: break
            s = fh.readline().strip(); p = fh.readline(); q = fh.readline().strip()
            n_in += 1
            comp_in.update(s)
            if n_in > 40000:      # decoy does not need to be deeper than the real thing
                break
            chars = list(s)
            random.shuffle(chars)
            sh = "".join(chars)
            comp_out.update(sh)
            n_out += 1
            out.write(f"@decoy_{n_out} shuffled_from=${meta.id}\\n{sh}\\n+\\n{q}\\n")

    with open("decoy_composition.tsv", "w") as r:
        r.write("metric\\tsource\\tdecoy\\n")
        r.write(f"reads\\t{n_in}\\t{n_out}\\n")
        tot_i = sum(comp_in[b] for b in "ACGT") or 1
        tot_o = sum(comp_out[b] for b in "ACGT") or 1
        for b in "ACGT":
            r.write(f"frac_{b}\\t{comp_in[b]/tot_i:.4f}\\t{comp_out[b]/tot_o:.4f}\\n")
    """

    stub:
    """
    printf '@d1\\nACGT\\n+\\nIIII\\n' | gzip > NEG_DECOY.fastq.gz
    printf 'metric\\tsource\\tdecoy\\nreads\\t1\\t1\\n' > decoy_composition.tsv
    """
}
