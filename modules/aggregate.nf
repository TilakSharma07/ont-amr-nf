process AGGREGATE_RESULTS {
    label 'process_low'
    conda "conda-forge::python=3.11"
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path amr_tsvs,    stageAs: "amr/*"
    path stats_tsvs,  stageAs: "stats/*"
    path verdicts,    stageAs: "verdict/*"

    output:
    path "amr_calls.tsv",          emit: amr
    path "assembly_metrics.tsv",   emit: stats
    path "validation_summary.tsv", emit: validation
    path "run_summary.md",         emit: summary
    // Written only when titration points are present in the run.
    path "depth_titration.tsv",    emit: titration, optional: true

    script:
    """
    #!/usr/bin/env python3
    import csv, glob, json, os

    def num(x):
        # The caller-level control has no reads, so depth/N50/length are null rather
        # than zero -- zero would read as a measured value. Render them as NA.
        return "NA" if x is None else x

    # ---- AMR calls: one tidy long table, every sample, one row per determinant ----
    amr_rows = []
    for fn in sorted(glob.glob("amr/*.tsv")):
        sample = os.path.basename(fn).replace(".amrfinder.tsv", "")
        with open(fn) as fh:
            for r in csv.DictReader(fh, delimiter="\\t"):
                # AMRFinderPlus column names have shifted across major versions;
                # resolve them defensively rather than assuming one schema.
                def pick(*names, default="NA"):
                    for n in names:
                        if n in r and r[n] not in (None, ""):
                            return r[n]
                    return default
                amr_rows.append({
                    "sample_id":    sample,
                    "gene_symbol":  pick("Element symbol", "Gene symbol"),
                    "gene_name":    pick("Element name", "Sequence name"),
                    "element_type": pick("Element type", "Type"),
                    "subtype":      pick("Element subtype", "Subtype"),
                    "drug_class":   pick("Class"),
                    "subclass":     pick("Subclass"),
                    "method":       pick("Method"),
                    "pct_identity": pick("% Identity to reference", "% Identity to reference sequence"),
                    "pct_coverage": pick("% Coverage of reference", "% Coverage of reference sequence"),
                    "contig":       pick("Contig id", "Contig"),
                })

    amr_cols = ["sample_id","gene_symbol","gene_name","element_type","subtype",
                "drug_class","subclass","method","pct_identity","pct_coverage","contig"]
    with open("amr_calls.tsv","w",newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=amr_cols, delimiter="\\t")
        w.writeheader(); w.writerows(amr_rows)

    # ---- assembly metrics ----
    stat_rows = []
    for fn in sorted(glob.glob("stats/*.tsv")):
        with open(fn) as fh:
            stat_rows += list(csv.DictReader(fh, delimiter="\\t"))
    if stat_rows:
        with open("assembly_metrics.tsv","w",newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(stat_rows[0]), delimiter="\\t")
            w.writeheader(); w.writerows(stat_rows)

    # ---- validation verdicts ----
    vrows = []
    for fn in sorted(glob.glob("verdict/*.json")):
        with open(fn) as fh:
            v = json.load(fh)
        if v:
            vrows.append(v)
    with open("validation_summary.tsv","w",newline="") as fh:
        fh.write("sample_id\\trole\\tverdict\\tflye_status\\tamr_calls\\ttotal_elements\\t"
                 "mean_depth\\tn50\\t"
                 "total_len\\tn_contigs\\tinterpretable\\tfailed_checks\\n")
        for v in sorted(vrows, key=lambda x: x.get("sample_id","")):
            failed = ";".join(c["check"] for c in v.get("checks",[]) if c["status"]=="FAIL") or "-"
            fh.write(f"{v['sample_id']}\\t{v['role']}\\t{v['verdict']}\\t"
                     f"{v.get('flye_status','NA')}\\t{v['amr_calls']}\\t"
                     f"{v.get('total_elements','NA')}\\t"
                     f"{num(v['mean_depth'])}\\t{num(v['n50'])}\\t"
                     f"{num(v['total_len'])}\\t{num(v['n_contigs'])}\\t"
                     f"{v['amr_result_interpretable']}\\t{failed}\\n")

    # ---- depth titration: determinants recovered as a function of depth ----
    # Only written when titration points are present. Reports the RECOVERY FRACTION
    # against the full-depth call set of the same parent sample, which is what makes
    # the curve interpretable: "how much of the truth do I still see at this depth".
    # Built here rather than lower down: the titration block below reads it, and it
    # was previously defined after its first use (a NameError that only fired with
    # --run_titration true, which no completed run had used).
    # AMR only. AMRFinderPlus returns AMR, STRESS and VIRULENCE elements in one
    # table, and counting all three under "determinants recovered" both inflates the
    # denominator and flatters the curve: the stress-tolerance genes are numerous and
    # largely depth-insensitive here, so mixing them in dilutes the loss of the calls
    # that matter. It also silently disagreed with every other number in the repo --
    # the gate's positive expectation, the figures and the results table all count AMR
    # only. Measured on this run, the mix reported 0.8837 recovery at 10x where the
    # AMR-only figure is 0.7500, and the inflated version is the one that would have
    # been quoted.
    genes_per_sample = {}
    elements_per_sample = {}
    for r in amr_rows:
        elements_per_sample.setdefault(r["sample_id"], set()).add(r["gene_symbol"])
        if str(r.get("element_type", "")).upper() == "AMR":
            genes_per_sample.setdefault(r["sample_id"], set()).add(r["gene_symbol"])

    tit = [v for v in vrows if v["role"] == "titration"]
    full_sets = {}
    if tit:
        for v in vrows:
            if v["role"] != "titration":
                full_sets[v["sample_id"]] = genes_per_sample.get(v["sample_id"], set())

        with open("depth_titration.tsv","w",newline="") as fh:
            fh.write("sample_id\\tparent_id\\ttarget_depth\\treplicate\\trealised_depth\\t"
                     "n50\\tn_genes\\tn_genes_full_depth\\trecovery_fraction\\tgenes_missed\\tverdict\\n")
            for v in sorted(tit, key=lambda x: x["sample_id"]):
                sid = v["sample_id"]
                # sample_id encodes parent/depth/replicate: <parent>_d<depth>_r<rep>
                base = sid.rsplit("_d", 1)
                parent = base[0]
                tail = base[1] if len(base) > 1 else ""
                depth_s, _, rep_s = tail.partition("_r")
                found = genes_per_sample.get(sid, set())
                truth = full_sets.get(parent, set())
                missed = sorted(truth - found)
                frac = (len(found & truth) / len(truth)) if truth else float("nan")
                fh.write(f"{sid}\\t{parent}\\t{depth_s}\\t{rep_s}\\t{v['mean_depth']:.2f}\\t"
                         f"{v['n50']}\\t{len(found)}\\t{len(truth)}\\t{frac:.4f}\\t"
                         f"{';'.join(missed) or '-'}\\t{v['verdict']}\\n")

    # ---- human-readable run summary ----
    n_pass = sum(1 for v in vrows if v["verdict"] == "PASS")
    neg = [v for v in vrows if v["role"] == "negative_control"]

    L = ["# Run summary\\n",
         f"- samples evaluated: **{len(vrows)}**",
         f"- passed validation: **{n_pass}/{len(vrows)}**",
         f"- total AMR determinant calls: **{len(amr_rows)}**\\n",
         "## Per-sample\\n",
         "| Sample | Role | Verdict | Depth | N50 | AMR genes |",
         "|---|---|---|---|---|---|"]
    for v in sorted(vrows, key=lambda x: (x["role"], x["sample_id"])):
        g = len(genes_per_sample.get(v["sample_id"], []))
        d  = f"{v['mean_depth']:.1f}x" if v.get("mean_depth") is not None else "n/a"
        n5 = f"{v['n50']:,}"            if v.get("n50")        is not None else "n/a"
        L.append(f"| `{v['sample_id']}` | {v['role']} | **{v['verdict']}** | "
                 f"{d} | {n5} | {g} |")
    if tit:
        L.append("\\n## Depth titration\\n")
        L.append("Recovery fraction is measured against the full-depth call set of the "
                 "same parent isolate.\\n")
        L.append("| Point | Target | Realised | Genes | Recovery |")
        L.append("|---|---|---|---|---|")
        for v in sorted(tit, key=lambda x: x["sample_id"]):
            sid = v["sample_id"]
            base = sid.rsplit("_d", 1)
            parent = base[0]
            depth_s = base[1].partition("_r")[0] if len(base) > 1 else "?"
            found = genes_per_sample.get(sid, set())
            truth = full_sets.get(parent, set())
            frac = (len(found & truth) / len(truth)) if truth else 0.0
            L.append(f"| `{sid}` | {depth_s}x | {v['mean_depth']:.1f}x | "
                     f"{len(found)}/{len(truth)} | {frac:.0%} |")

    L.append("\\n## Control outcomes\\n")
    if neg:
        n = neg[0]
        ok = "as required" if n.get("total_elements", n["amr_calls"]) == 0 else "UNEXPECTED — investigate"
        L.append(f"**Read-level decoy** (tests the assembler) — `{n['sample_id']}` returned "
                 f"**{n.get('total_elements', n['amr_calls'])}** elements of any type ({ok}); flye status "
                 f"`{n.get('flye_status','NA')}`, which is the expected outcome for "
                 f"shuffled reads.")
    else:
        L.append("No read-level decoy was run (`--make_decoy false`).")

    # The caller-level control is gated by CONTROL_GATE (it bypasses the assembly gate
    # by design -- no reads, so no depth or N50 -- but its emptiness is a binding check
    # there, and its verdict appears in vrows above). The narrative below reads its call
    # set directly so the prose can name the genes if any survived.
    #
    # Presence is read from the verdict rows, NOT from the call set. Inferring it from
    # the calls made a control that returned zero calls indistinguishable from one that
    # never ran -- and zero calls is the PASS case, so the pipeline reported its own
    # best result as "No caller-level control was run". That happened: results_main has
    # CALLER_CONTROL.amrfinder.tsv with zero calls and a run summary denying the control
    # existed. A reviewer reading that summary would conclude the control was skipped.
    cc_rows = [v for v in vrows if v["role"] == "caller_control"]
    if not cc_rows:
        L.append("\\nNo caller-level control was run (`--caller_control false`).")
    else:
        cc_genes = elements_per_sample.get("CALLER_CONTROL", set())
        # Elements of ANY class, not AMR only: a control must return nothing at all,
        # not merely nothing under one label.
        n_cc = sum(1 for r in amr_rows if r["sample_id"] == "CALLER_CONTROL")
        if n_cc == 0:
            L.append(f"\\n**Caller-level control** (tests the gene caller) — a real assembly "
                     f"with its bases shuffled within each contig, identical in contig "
                     f"count, length and GC, returned **0** calls. The caller is matching "
                     f"on gene identity, not composition.")
        else:
            L.append(f"\\n**Caller-level control** — returned **{n_cc}** calls "
                     f"({', '.join(sorted(cc_genes))}) on shuffled sequence. This is a "
                     f"FALSE POSITIVE by construction: no genes are present. "
                     f"**Investigate before trusting any call in this run.**")
    open("run_summary.md","w").write("\\n".join(L) + "\\n")

    print(f"[aggregate] {len(vrows)} samples, {len(amr_rows)} AMR calls, {n_pass} passed validation")
    """

    stub:
    """
    touch amr_calls.tsv assembly_metrics.tsv validation_summary.tsv run_summary.md
    """
}
