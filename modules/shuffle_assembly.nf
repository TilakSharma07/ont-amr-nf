process SHUFFLE_ASSEMBLY {
    tag   "caller_control_from:${meta.id}"
    label 'process_low'
    conda "conda-forge::python=3.11"
    publishDir "${params.outdir}/controls", mode: params.publish_mode

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("CALLER_CONTROL.fasta"), emit: assembly
    path "caller_control_composition.tsv",         emit: report

    script:
    """
    #!/usr/bin/env python3
    import random, collections

    # WHY THIS EXISTS, separately from the read-level decoy:
    #
    # The read-level decoy (MAKE_DECOY) does not assemble — which is the expected and
    # correct outcome for shuffled reads. But that means the resistance-gene caller
    # never sees any sequence from it: it receives an empty assembly and returns an
    # empty call set trivially. So that control tests the ASSEMBLER, not the CALLER.
    #
    # This control closes that gap. It takes a real, good assembly and shuffles the
    # bases WITHIN each contig. The result preserves contig count, contig lengths, GC
    # content and k-mer-free base composition, and destroys gene content. It is then
    # passed to the same AMRFINDERPLUS process, with the same parameters, as every real
    # isolate. Any call made here is a call driven by composition rather than by gene
    # identity — a false positive, by construction.
    random.seed(20240117)

    def read_fasta(path):
        name, buf = None, []
        with open(path) as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if name:
                        yield name, "".join(buf)
                    name, buf = line[1:], []
                else:
                    buf.append(line)
        if name:
            yield name, "".join(buf)

    comp_in, comp_out = collections.Counter(), collections.Counter()
    n_contigs = 0
    total_len = 0

    with open("CALLER_CONTROL.fasta", "w") as out:
        for name, seq in read_fasta("${assembly}"):
            n_contigs += 1
            total_len += len(seq)
            comp_in.update(seq.upper())
            chars = list(seq)
            random.shuffle(chars)
            shuf = "".join(chars)
            comp_out.update(shuf.upper())
            out.write(f">caller_control_{n_contigs} shuffled_from={name.split()[0]} len={len(shuf)}\\n")
            for i in range(0, len(shuf), 60):
                out.write(shuf[i:i + 60] + "\\n")

    with open("caller_control_composition.tsv", "w") as r:
        r.write("metric\\tsource\\tshuffled\\n")
        r.write(f"contigs\\t{n_contigs}\\t{n_contigs}\\n")
        r.write(f"total_bp\\t{total_len}\\t{total_len}\\n")
        ti = sum(comp_in[b] for b in "ACGT") or 1
        to = sum(comp_out[b] for b in "ACGT") or 1
        for b in "ACGT":
            r.write(f"frac_{b}\\t{comp_in[b] / ti:.6f}\\t{comp_out[b] / to:.6f}\\n")
        gc_i = (comp_in["G"] + comp_in["C"]) / ti
        gc_o = (comp_out["G"] + comp_out["C"]) / to
        r.write(f"gc_fraction\\t{gc_i:.6f}\\t{gc_o:.6f}\\n")

    print(f"[caller_control] {n_contigs} contigs, {total_len:,} bp, "
          f"GC {gc_i:.4f} -> {gc_o:.4f} (must be identical)")
    """

    stub:
    """
    printf '>caller_control_1 shuffled_from=ctg1\\nACGTACGTACGT\\n' > CALLER_CONTROL.fasta
    printf 'metric\\tsource\\tshuffled\\ncontigs\\t1\\t1\\n' > caller_control_composition.tsv
    """
}
