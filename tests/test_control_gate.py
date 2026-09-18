#!/usr/bin/env python3
"""Tests for modules/control_gate.nf and the aggregator's handling of its verdict.

Run:  python3 tests/test_control_gate.py <results_dir>

Why this file exists
--------------------
The caller-level control is the sharpest control in this pipeline: a real assembly
with its bases shuffled inside each contig, identical in contig count, length and
GC, re-called with identical parameters. If it returns determinants, the caller is
matching on composition rather than gene identity.

It was also, until CONTROL_GATE existed, the one result nothing checked. It has no
reads, so it cannot pass through the assembly gate, and its calls went straight to
the aggregator. "The control came back empty" was a sentence in the README, not a
check in the code — if the shuffle had silently stopped shuffling, every isolate's
call set would have been reproduced on the control and the run would still have
reported success.

These tests run the gate's real body (extracted from the .nf, as with the
validation gate) against the published control call set and against a contaminated
one, and they pin the two things the aggregator must get right about a sample with
no reads: null depth/N50 must not crash it, and must not print as "None".
"""
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GATE = os.path.join(ROOT, "modules", "control_gate.nf")
AGG = os.path.join(ROOT, "modules", "aggregate.nf")

AMRFINDER_COLS = ["Name", "Element symbol", "Element name", "Type", "Subtype",
                  "Class", "Subclass", "Method", "% Identity to reference",
                  "% Coverage of reference", "Contig id"]


def gate_source(amr_calls, sample_id="CALLER_CONTROL", role="caller_control",
                shuffled_from="KP_ES_7636"):
    nf = open(GATE).read()
    m = re.search(r'script:.*?"""(.*?)"""', nf, re.S)
    assert m, "script block not found in control_gate.nf"
    src = textwrap.dedent(m.group(1)).replace("\\\\", "\\")
    for key, val in [("${amr_calls}", amr_calls), ("${meta.id}", sample_id),
                     ("${meta.role}", role),
                     ("${meta.shuffled_from ?: 'unknown'}", shuffled_from)]:
        src = src.replace(key, val)
    left = re.findall(r"\$\{[^}]+\}", src)
    assert not left, f"unresolved interpolations: {sorted(set(left))}"
    return src


def run_gate(tmp, amr_calls, sample_id="CALLER_CONTROL"):
    d = tempfile.mkdtemp(dir=tmp)
    script = os.path.join(d, "gate.py")
    open(script, "w").write(gate_source(amr_calls, sample_id=sample_id))
    r = subprocess.run([sys.executable, script], cwd=d, capture_output=True, text=True)
    assert r.returncode == 0, f"gate crashed:\n{r.stderr[-1200:]}"
    return json.load(open(os.path.join(d, f"{sample_id}.verdict.json")))


def write_calls(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=AMRFINDER_COLS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "NA") for c in AMRFINDER_COLS})


def aggregator_source():
    """Extract the aggregator's Python body with its interpolations resolved."""
    nf = open(AGG).read()
    m = re.search(r'script:.*?"""(.*?)"""', nf, re.S)
    assert m, "script block not found in aggregate.nf"
    src = textwrap.dedent(m.group(1)).replace("\\\\", "\\")
    src = src.replace("${params.outdir}", "results")
    left = re.findall(r"\$\{[^}]+\}", src)
    assert not left, f"unresolved interpolations: {sorted(set(left))}"
    return src


def aggregator_handles_null_depth(tmp, results):
    """Run the real aggregator over one normal verdict plus a null-depth control.

    This is the failure the caller control introduces: formatting None with :.1f
    raises TypeError, so an unguarded aggregator crashes the run at the last step,
    after every assembly and every AMR call has already been computed.
    """
    d = tempfile.mkdtemp(dir=tmp)
    for sub in ("verdict", "amr", "stats"):
        os.makedirs(os.path.join(d, sub), exist_ok=True)

    normal = {"sample_id": "KP_ES_7636", "role": "test", "verdict": "PASS",
              "flye_status": "assembled", "amr_calls": 20, "total_elements": 43,
              "stress_calls": 23, "virulence_calls": 0, "mean_depth": 29.9,
              "n50": 5425447, "total_len": 5887566, "n_contigs": 11,
              "breadth_1x": 1.0, "amr_result_interpretable": True, "checks": []}
    control = {"sample_id": "CALLER_CONTROL", "role": "caller_control",
               "verdict": "PASS", "flye_status": "not_applicable",
               "amr_calls": 0, "total_elements": 0, "stress_calls": 0,
               "virulence_calls": 0, "mean_depth": None, "n50": None,
               "total_len": None, "n_contigs": None, "breadth_1x": None,
               "amr_result_interpretable": True, "shuffled_from": "KP_ES_7636",
               "checks": [{"check": "caller_control_is_empty", "status": "PASS",
                           "detail": "0 elements"}]}
    for v in (normal, control):
        json.dump(v, open(os.path.join(d, "verdict", f"{v['sample_id']}.json"), "w"))

    write_calls(os.path.join(d, "amr", "KP_ES_7636.amrfinder.tsv"),
                [{"Name": "KP_ES_7636", "Element symbol": "blaOXA-48", "Type": "AMR"}])
    write_calls(os.path.join(d, "amr", "CALLER_CONTROL.amrfinder.tsv"), [])
    with open(os.path.join(d, "stats", "KP_ES_7636.stats.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "mean_depth", "n50", "total_len", "n_contigs", "breadth_1x"])
        w.writerow(["KP_ES_7636", 29.9, 5425447, 5887566, 11, 1.0])

    script = os.path.join(d, "agg.py")
    open(script, "w").write(aggregator_source())
    r = subprocess.run([sys.executable, script], cwd=d, capture_output=True, text=True)
    if r.returncode != 0:
        tail = (r.stderr or "").strip().splitlines()[-1:] or ["(no stderr)"]
        return False, f"aggregator exited {r.returncode}: {tail[0][:120]}"

    md = open(os.path.join(d, "run_summary.md")).read()
    tsv = open(os.path.join(d, "validation_summary.tsv")).read()
    cc_line = [l for l in md.splitlines() if "CALLER_CONTROL" in l and l.startswith("|")]
    if not cc_line:
        return False, "the control is absent from the per-sample table"
    if "None" in cc_line[0]:
        return False, f"null rendered as 'None': {cc_line[0].strip()}"
    if "NA" not in [f.strip() for f in
                    next(l for l in tsv.splitlines() if l.startswith("CALLER_CONTROL")).split("\t")]:
        return False, "null numerics not rendered as NA in validation_summary.tsv"
    return True, f"control row: {cc_line[0].strip()}"


def mutation_breaks_null_handling(tmp, results):
    """Restore the pre-fix formatting and confirm the aggregator crashes.

    Without this, the check above would pass on an aggregator that never touched
    depth at all — a green test that proves nothing about the guard.
    """
    target = """    d  = f"{v['mean_depth']:.1f}x" if v.get("mean_depth") is not None else "n/a"
    n5 = f"{v['n50']:,}"            if v.get("n50")        is not None else "n/a"
    L.append(f"| `{v['sample_id']}` | {v['role']} | **{v['verdict']}** | "
             f"{d} | {n5} | {g} |")"""
    replacement = """    L.append(f"| `{v['sample_id']}` | {v['role']} | **{v['verdict']}** | "
             f"{v['mean_depth']:.1f}x | {v['n50']:,} | {g} |")"""

    real = aggregator_source
    src = real()
    if target not in src:
        return False, "mutation target not found — the guard was rewritten; update this test"

    globals()["aggregator_source"] = lambda: src.replace(target, replacement)
    try:
        ok, detail = aggregator_handles_null_depth(tmp, results)
    finally:
        globals()["aggregator_source"] = real
    return (not ok), ("the mutant still passed — the check is vacuous"
                      if ok else f"mutant correctly rejected ({detail[:80]})")


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"\n        {detail}" if detail else ""))
    return bool(cond)


def main(results):
    tmp = tempfile.mkdtemp(prefix="ccgate-")
    ok = []

    # --- 1. The real published control -------------------------------------
    real = os.path.join(results, "amr", "CALLER_CONTROL.amrfinder.tsv")
    if os.path.exists(real):
        v = run_gate(tmp, real)
        ok.append(check("published caller control PASSES with zero elements",
                        v["verdict"] == "PASS" and v["total_elements"] == 0,
                        f"verdict={v['verdict']}, elements={v['total_elements']}"))
        ok.append(check("assembly-quality fields are null, not zero",
                        v["mean_depth"] is None and v["n50"] is None,
                        "a control with no reads has no depth; zero would read as measured"))
        # NOT "shuffled_from is set" — this test substitutes that value itself, so
        # such a check would only confirm its own substitution. What needs checking
        # is the wiring: that the workflow puts the source sample id into the
        # control's meta, and that the gate reads it back out.
        main_nf = open(os.path.join(ROOT, "main.nf")).read()
        i_cc = main_nf.find("id       : 'CALLER_CONTROL'")
        cc_meta = main_nf[i_cc:main_nf.find("]", i_cc)] if i_cc != -1 else ""
        gate_nf = open(GATE).read()
        ok.append(check("the workflow records which sample the control was shuffled from",
                        "shuffled_from" in cc_meta and "shuffled_from" in v
                        and "${meta.shuffled_from" in gate_nf,
                        f"meta sets it={('shuffled_from' in cc_meta)}, "
                        f"gate reads it={('${meta.shuffled_from' in gate_nf)}, "
                        f"verdict field present={('shuffled_from' in v)}"))
    else:
        print(f"  SKIP  no published control call set at {real}")

    # --- 2. A contaminated control must FAIL -------------------------------
    # This is the case the gate exists for. A single surviving element of any class
    # means the caller matched something on shuffled sequence.
    for etype in ("AMR", "STRESS", "VIRULENCE"):
        p = os.path.join(tmp, f"contaminated_{etype}.tsv")
        write_calls(p, [{"Name": "CALLER_CONTROL", "Element symbol": "blaOXA-48",
                         "Type": etype, "Subtype": "AMR"}])
        v = run_gate(tmp, p)
        ok.append(check(f"a surviving {etype} element FAILS the control",
                        v["verdict"] == "FAIL" and v["amr_result_interpretable"] is False,
                        f"verdict={v['verdict']}"))

    # --- 3. The aggregator must survive null depth -------------------------
    # Formatting None with :.1f raises TypeError; printing it bare writes "None"
    # into a numeric column. Both are defects, so pin the NA rendering.
    agg = open(AGG).read()
    ok.append(check("aggregator renders null numerics as NA",
                    "def num(x)" in agg and 'return "NA" if x is None else x' in agg,
                    "num() helper present"))
    # The per-sample markdown table iterates over ALL verdicts including the control,
    # so its depth/N50 cells must tolerate null. Test the BEHAVIOUR, not the source
    # text: a guarded expression still contains a :.1f format spec, so pattern-matching
    # the source cannot tell a guarded line from an unguarded one.
    ok.append(check("the per-sample table renders a null-depth sample without crashing",
                    *aggregator_handles_null_depth(tmp, results)))

    # And the mutation: with the guard removed the aggregator must actually break,
    # otherwise the check above proves nothing. The aggregator runs last, so this
    # failure would land after every assembly and every AMR call had been computed.
    ok.append(check("...and fails if that guard is removed (mutation test)",
                    *mutation_breaks_null_handling(tmp, results)))

    # --- 4. The NameError that titration would have hit --------------------
    # genes_per_sample was read in the titration branch 27 lines before it was
    # defined. Never fired, because no completed run had titration enabled.
    i_def = agg.find("genes_per_sample = {}")
    i_use = agg.find('tit = [v for v in vrows if v["role"] == "titration"]')
    ok.append(check("genes_per_sample is defined before the titration block reads it",
                    i_def != -1 and i_use != -1 and i_def < i_use,
                    f"defined at offset {i_def}, titration block at {i_use}"))

    print()
    if all(ok):
        print(f"all {len(ok)} checks passed")
        return 0
    print(f"{sum(1 for x in ok if not x)} of {len(ok)} checks FAILED")
    return 1


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: test_control_gate.py <results_dir>")
    sys.exit(main(os.path.abspath(sys.argv[1])))
