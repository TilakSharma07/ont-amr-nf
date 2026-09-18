process VALIDATION_GATE {
    tag   "${meta.id}"
    label 'process_low'
    conda "conda-forge::python=3.11"
    publishDir "${params.outdir}/validation", mode: params.publish_mode

    input:
    tuple val(meta), path(stats), path(amr_calls), path(flye_status)

    output:
    tuple val(meta), path("${meta.id}.verdict.tsv"), emit: verdict
    path "${meta.id}.verdict.json",                  emit: json

    script:
    """
    #!/usr/bin/env python3
    import csv, json

    # The point of this process: an AMR result is only reportable if the assembly it
    # came from is good enough to support it. A clean pipeline run is not evidence.
    # Absence of a resistance gene in a 6x assembly is not absence of the gene.

    with open("${stats}") as fh:
        s = list(csv.DictReader(fh, delimiter="\\t"))[0]

    with open("${amr_calls}") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\\t")]

    with open("${flye_status}") as fh:
        flye = list(csv.DictReader(fh, delimiter="\\t"))[0]
    assembled = flye["flye_status"] == "assembled"

    role       = "${meta.role}"
    depth      = float(s["mean_depth"])
    n50        = int(s["n50"])
    total_len  = int(s["total_len"])
    n_contigs  = int(s["n_contigs"])
    breadth    = float(s["breadth_1x"])
    n_amr      = len(rows)

    checks = []
    def check(name, ok, detail):
        checks.append({"check": name, "status": "PASS" if ok else "FAIL", "detail": detail})

    check("assembly_produced", assembled,
          f"flye status: {flye['flye_status']} (exit {flye['flye_exit']})")
    check("depth",       depth     >= ${params.min_depth_x},
          f"mean realised depth {depth:.1f}x (floor {${params.min_depth_x}}x)")
    check("n50",         n50       >= ${params.min_n50},
          f"N50 {n50:,} bp (floor {${params.min_n50}:,} bp)")
    check("genome_size", ${params.min_total_len} <= total_len <= ${params.max_total_len},
          f"assembly {total_len:,} bp (expected {${params.min_total_len}:,}-{${params.max_total_len}:,})")
    check("fragmentation", n_contigs <= ${params.max_contigs},
          f"{n_contigs} contigs (ceiling {${params.max_contigs}})")
    check("breadth",     breadth   >= 0.95,
          f"{breadth:.1%} of assembly covered at >=1x")

    # Control-specific expectations. These are the checks that would catch a pipeline
    # that silently produces plausible-looking garbage.
    if role == "negative_control":
        check("negative_control_is_empty", n_amr == 0,
              f"{n_amr} AMR calls on shuffled sequence (must be 0)")
    else:
        check("positive_expectation_met", n_amr > 0,
              f"{n_amr} AMR determinants found; these are clinical carbapenemase/ESBL "
              f"isolates, so zero calls would indicate pipeline failure")

    hard_fail = [c for c in checks if c["status"] == "FAIL"]
    # A negative control that assembles badly has not failed — that is the expected
    # behaviour of shuffled input. Only its emptiness check is binding.
    if role == "negative_control":
        binding = [c for c in hard_fail if c["check"] == "negative_control_is_empty"]
    else:
        binding = hard_fail

    verdict = "PASS" if not binding else "FAIL"
    interpretable = verdict == "PASS"

    out = {
        "sample_id": "${meta.id}", "role": role, "verdict": verdict,
        "flye_status": flye["flye_status"],
        "amr_calls": n_amr, "mean_depth": depth, "n50": n50,
        "total_len": total_len, "n_contigs": n_contigs, "breadth_1x": breadth,
        "amr_result_interpretable": interpretable,
        "checks": checks,
    }
    json.dump(out, open("${meta.id}.verdict.json", "w"), indent=2)

    with open("${meta.id}.verdict.tsv", "w") as fh:
        fh.write("sample_id\\trole\\tcheck\\tstatus\\tdetail\\n")
        for c in checks:
            fh.write(f"${meta.id}\\t{role}\\t{c['check']}\\t{c['status']}\\t{c['detail']}\\n")
        fh.write(f"${meta.id}\\t{role}\\tOVERALL\\t{verdict}\\t"
                 f"amr_interpretable={interpretable}\\n")

    print(f"[validation] ${meta.id} ({role}): {verdict} — {n_amr} AMR calls, {depth:.1f}x, N50 {n50:,}")
    for c in checks:
        if c["status"] == "FAIL":
            print(f"  FAIL {c['check']}: {c['detail']}")
    """

    stub:
    """
    printf 'sample_id\\trole\\tcheck\\tstatus\\tdetail\\n' > ${meta.id}.verdict.tsv
    echo '{}' > ${meta.id}.verdict.json
    """
}
