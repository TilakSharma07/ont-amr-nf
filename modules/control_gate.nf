process CONTROL_GATE {
    tag   "${meta.id}"
    label 'process_low'
    conda "conda-forge::python=3.11"
    publishDir "${params.outdir}/validation", mode: params.publish_mode

    input:
    tuple val(meta), path(amr_calls)

    output:
    tuple val(meta), path("${meta.id}.verdict.tsv"), emit: verdict
    path "${meta.id}.verdict.json",                  emit: json

    script:
    // Why this process exists separately from VALIDATION_GATE:
    //
    // The caller-level control has no reads, therefore no depth, N50 or breadth, so it
    // cannot pass through the assembly gate — that gate's whole premise is that an AMR
    // result is only reportable if the assembly behind it is good enough. The control's
    // assembly is a real one with its bases shuffled, so assembly quality is not the
    // question being asked of it.
    //
    // But leaving it ungated meant the single most important control in this pipeline
    // was the one result nobody checked: its calls went straight to the aggregator, and
    // "the control came back empty" was a claim in the README rather than a check in the
    // code. If the shuffle silently stopped shuffling, nothing would have failed.
    //
    // So it gets its own gate with exactly one binding check, and its verdict lands in
    // the same validation_summary.tsv as every other sample.
    """
    #!/usr/bin/env python3
    import csv, json

    def pick(row, *names):
        # AMRFinderPlus has renamed these columns across versions ("Element type" ->
        # "Type"). Accept either rather than pinning one spelling.
        for n in names:
            if n in row and row[n] != "":
                return row[n]
        return ""

    rows = []
    with open("${amr_calls}", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\\t"):
            rows.append(row)

    n_elements  = len(rows)
    etypes      = [pick(r, "element_type", "Element type", "Type").upper() for r in rows]
    n_amr       = sum(1 for e in etypes if e == "AMR")
    n_stress    = sum(1 for e in etypes if e == "STRESS")
    n_virulence = sum(1 for e in etypes if e == "VIRULENCE")

    checks = []
    def check(name, ok, detail):
        checks.append({"check": name, "status": "PASS" if ok else "FAIL", "detail": detail})

    # A control must return nothing at all, not merely nothing under one label: an
    # element of any class surviving on shuffled sequence is a false positive, and
    # whether the caller labelled it AMR, STRESS or VIRULENCE does not change that.
    check("caller_control_is_empty", n_elements == 0,
          f"{n_elements} elements of any type on shuffled assembly (must be 0; "
          f"of which {n_amr} AMR, {n_stress} stress, {n_virulence} virulence)")

    verdict = "PASS" if not [c for c in checks if c["status"] == "FAIL"] else "FAIL"

    # Assembly-quality fields are null rather than zero. Zero would read as a measured
    # value; this control has no reads, so those quantities do not exist for it.
    out = {
        "sample_id": "${meta.id}", "role": "${meta.role}", "verdict": verdict,
        "flye_status": "not_applicable",
        "amr_calls": n_amr, "total_elements": n_elements,
        "stress_calls": n_stress, "virulence_calls": n_virulence,
        "mean_depth": None, "n50": None,
        "total_len": None, "n_contigs": None, "breadth_1x": None,
        "amr_result_interpretable": verdict == "PASS",
        "shuffled_from": "${meta.shuffled_from ?: 'unknown'}",
        "checks": checks,
    }
    json.dump(out, open("${meta.id}.verdict.json", "w"), indent=2)

    with open("${meta.id}.verdict.tsv", "w") as fh:
        fh.write("sample_id\\trole\\tcheck\\tstatus\\tdetail\\n")
        for c in checks:
            fh.write(f"${meta.id}\\t${meta.role}\\t{c['check']}\\t{c['status']}\\t{c['detail']}\\n")
        fh.write(f"${meta.id}\\t${meta.role}\\tOVERALL\\t{verdict}\\tcaller-level control\\n")

    print(f"[validation] ${meta.id} (${meta.role}): {verdict} — {n_elements} elements "
          f"on shuffled sequence derived from ${meta.shuffled_from ?: 'unknown'}")
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
