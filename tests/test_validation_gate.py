#!/usr/bin/env python3
"""Tests for the counting logic in modules/validation_gate.nf.

Run:  python3 tests/test_validation_gate.py <results_dir>

Why this file exists
--------------------
AMRFinderPlus returns three classes of element in one table: AMR (acquired and
mutational resistance determinants), STRESS (biocide, metal and heat tolerance)
and VIRULENCE. An earlier version of this gate counted every row as an "AMR
call", which did two things:

  * inflated the reported determinant count by roughly 2x on these isolates
    (KP_ES_7636: 43 rows, of which 20 are resistance determinants); and
  * quietly weakened `positive_expectation_met`, because an isolate with zero
    resistance genes and twenty metal-tolerance genes satisfied a check that
    claims to prove the resistance caller works.

The second is the dangerous one: the gate would have passed a run in which
resistance calling had failed completely. These tests pin the distinction, using
the real published call table plus synthetic tables for the cases a real run does
not happen to contain.

The gate's body is Nextflow-interpolated Python. Rather than duplicating it here
(which would test the copy, not the pipeline), the test extracts the script block
from the .nf file and substitutes the interpolations, so it runs the same source
the pipeline runs.
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
GATE = os.path.join(ROOT, "modules", "validation_gate.nf")

# Mirrors the defaults in nextflow.config. Kept explicit so a threshold change
# that breaks these tests is a deliberate decision, not a silent drift.
PARAMS = {
    "params.min_depth_x": "20",
    "params.min_n50": "50000",
    "params.min_total_len": "4000000",
    "params.max_total_len": "7000000",
    "params.max_contigs": "200",
    "params.outdir": "results",
}

AMR_COLS = ["sample_id", "gene_symbol", "gene_name", "element_type", "subtype",
            "drug_class", "subclass", "method", "pct_identity", "pct_coverage", "contig"]


def gate_source(sample_id, role, stats, amr_calls, flye_status):
    """Extract the gate's Python body and resolve Nextflow's interpolations."""
    nf = open(GATE).read()
    m = re.search(r'script:\s*"""(.*?)"""', nf, re.S)
    assert m, "could not find the script block in validation_gate.nf"
    src = m.group(1)

    # Nextflow renders a Groovy GString: an escaped backslash becomes one
    # backslash, so the \\t written in the .nf reaches Python as \t.
    src = src.replace("\\\\", "\\")

    # The block is indented inside the process definition; Nextflow does not care
    # but Python does. Strip the common leading whitespace.
    src = textwrap.dedent(src)

    subs = dict(PARAMS)
    subs.update({"meta.id": sample_id, "meta.role": role, "stats": stats,
                 "amr_calls": amr_calls, "flye_status": flye_status})
    for key, val in subs.items():
        src = src.replace("${" + key + "}", str(val))

    left = re.findall(r"\$\{[^}]+\}", src)
    assert not left, f"unresolved interpolations in extracted gate: {sorted(set(left))}"
    return src


def write_amr(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=AMR_COLS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "NA") for c in AMR_COLS})


def write_stats(path, sample_id, depth=30.0, n50=5_000_000, total_len=5_400_000,
                n_contigs=11, breadth=1.0):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "mean_depth", "n50", "total_len", "n_contigs", "breadth_1x"])
        w.writerow([sample_id, depth, n50, total_len, n_contigs, breadth])


def write_flye(path, status="assembled", exit_code=0):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["flye_status", "flye_exit"])
        w.writerow([status, exit_code])


def run_gate(tmp, sample_id, role, amr_rows, **stats_kw):
    """Run the extracted gate in a scratch dir; return its verdict dict."""
    d = tempfile.mkdtemp(dir=tmp)
    amr = os.path.join(d, "calls.tsv")
    st = os.path.join(d, "stats.tsv")
    fl = os.path.join(d, "flye.tsv")
    if isinstance(amr_rows, str):          # a path to a real table
        amr = amr_rows
    else:
        write_amr(amr, amr_rows)
    write_stats(st, sample_id, **stats_kw)
    write_flye(fl)

    src = gate_source(sample_id, role, st, amr, fl)
    script = os.path.join(d, "gate.py")
    open(script, "w").write(src)
    r = subprocess.run([sys.executable, script], cwd=d, capture_output=True, text=True)
    assert r.returncode == 0, f"gate crashed:\n{r.stderr[-1500:]}"
    return json.load(open(os.path.join(d, f"{sample_id}.verdict.json")))


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"\n        {detail}" if detail else ""))
    return bool(cond)


def main(results):
    tmp = tempfile.mkdtemp(prefix="gatetest-")
    ok = []

    # --- 1. Real data: the counts must separate AMR from the rest -------------
    real_calls = os.path.join(results, "amr_calls.tsv")
    rows = [r for r in csv.DictReader(open(real_calls), delimiter="\t")
            if r["sample_id"] == "KP_ES_7636"]
    n_amr_expected = sum(1 for r in rows if r["element_type"] == "AMR")
    one = os.path.join(tmp, "one.tsv")
    write_amr(one, rows)

    v = run_gate(tmp, "KP_ES_7636", "test", rows)
    ok.append(check("resistance determinants counted separately from all elements",
                    v["amr_calls"] == n_amr_expected and v["total_elements"] == len(rows),
                    f"amr_calls={v['amr_calls']} (expected {n_amr_expected}), "
                    f"total_elements={v['total_elements']} (expected {len(rows)})"))
    ok.append(check("stress and virulence reported, not folded into the AMR count",
                    v["amr_calls"] + v["stress_calls"] + v["virulence_calls"]
                    == v["total_elements"],
                    f"{v['amr_calls']} AMR + {v['stress_calls']} stress + "
                    f"{v['virulence_calls']} virulence = {v['total_elements']}"))

    # --- 2. The bug this gate was rewritten to catch -------------------------
    # An isolate whose resistance calling produced nothing, but which carries
    # plenty of metal-tolerance genes. The old count made this a PASS.
    stress_only = [{"sample_id": "FAKE", "gene_symbol": f"mer{i}",
                    "element_type": "STRESS", "subtype": "METAL"} for i in range(20)]
    v = run_gate(tmp, "FAKE", "test", stress_only)
    pos = [c for c in v["checks"] if c["check"] == "positive_expectation_met"][0]
    ok.append(check("isolate with 0 AMR but 20 stress genes FAILS positive expectation",
                    pos["status"] == "FAIL" and v["verdict"] == "FAIL",
                    f"{pos['status']}: {pos['detail'][:90]}"))

    # --- 3. Negative control must be empty of everything, not just of AMR ----
    v = run_gate(tmp, "NEG", "negative_control", stress_only,
                 depth=2.0, n50=1000, total_len=50_000, n_contigs=3, breadth=0.4)
    neg = [c for c in v["checks"] if c["check"] == "negative_control_is_empty"][0]
    ok.append(check("negative control with stress-only hits FAILS emptiness",
                    neg["status"] == "FAIL" and v["verdict"] == "FAIL",
                    f"{neg['status']}: {neg['detail'][:90]}"))

    # A genuinely empty control still passes, and its poor assembly is not fatal.
    v = run_gate(tmp, "NEG", "negative_control", [],
                 depth=2.0, n50=1000, total_len=50_000, n_contigs=3, breadth=0.4)
    ok.append(check("empty negative control PASSES despite failing assembly checks",
                    v["verdict"] == "PASS" and v["amr_calls"] == 0,
                    f"verdict={v['verdict']}, failed checks are non-binding for this role"))

    # --- 4. Depth floor still gates interpretability -------------------------
    v = run_gate(tmp, "LOW", "test", rows, depth=6.0)
    ok.append(check("below-floor depth makes the AMR result uninterpretable",
                    v["amr_result_interpretable"] is False,
                    f"6.0x < {PARAMS['params.min_depth_x']}x floor -> "
                    f"interpretable={v['amr_result_interpretable']}"))

    print()
    if all(ok):
        print(f"all {len(ok)} checks passed")
        return 0
    print(f"{sum(1 for x in ok if not x)} of {len(ok)} checks FAILED")
    return 1


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: test_validation_gate.py <results_dir>")
    sys.exit(main(os.path.abspath(sys.argv[1])))
