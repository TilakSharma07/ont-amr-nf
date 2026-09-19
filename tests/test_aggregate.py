#!/usr/bin/env python3
"""Tests for modules/aggregate.nf -- the process that turns per-sample outputs into
the tables and the run summary a reader actually sees.

Run:  python3 tests/test_aggregate.py

Why this file exists
--------------------
The aggregator had no tests, and two real defects were sitting in it. Both were
found by reading its output against its inputs rather than by anything failing:

1. `recovery_fraction` in depth_titration.tsv counted AMR, STRESS and VIRULENCE
   elements together, while every other number in the repo -- the gate's positive
   expectation, the figures, the results table -- counts resistance determinants
   only. On the real titration that reported 0.8837 recovery at 10x where the
   AMR-only figure is 0.7500. The inflated number is the one that would have been
   quoted, and it flatters the method: the stress-tolerance genes are numerous and
   largely depth-insensitive, so mixing them in dilutes the loss of the calls that
   matter.

2. The run summary inferred whether the caller-level control had run from whether
   it had produced any calls. Zero calls is the PASS case, so a control that
   behaved perfectly was reported as "No caller-level control was run". That is not
   hypothetical: results_main contains CALLER_CONTROL.amrfinder.tsv with zero calls
   and a run summary denying the control existed. A reviewer reading it would
   conclude the sharpest control in the pipeline had been skipped.

The aggregator's script block takes no Nextflow interpolations and reads its inputs
from amr/, stats/ and verdict/ directories, so it runs standalone in a tmpdir.
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
AGG = os.path.join(ROOT, "modules", "aggregate.nf")

AMR_COLS = ["Name", "Element symbol", "Element name", "Type", "Subtype", "Class",
            "Subclass", "Method", "% Identity to reference",
            "% Coverage of reference", "Contig id"]


def agg_source():
    nf = open(AGG).read()
    m = re.search(r'script:.*?"""(.*?)"""', nf, re.S)
    assert m, "script block not found in aggregate.nf"
    src = textwrap.dedent(m.group(1)).replace("\\\\", "\\")
    left = re.findall(r"\$\{[^}]+\}", src)
    assert not left, f"unresolved interpolations: {sorted(set(left))}"
    return src


def write_amr(path, sample_id, calls):
    """calls: [(gene_symbol, element_type)]"""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(AMR_COLS)
        for gene, etype in calls:
            w.writerow([sample_id, gene, f"{gene} product", etype, "AMR", "BETA-LACTAM",
                        "-", "EXACTX", "100.00", "100.00", "contig_1"])


def write_stats(path, sample_id, role, n50, depth, total_len=5_400_000):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "role", "species", "n_contigs", "total_len", "n50",
                    "largest", "gc_percent", "mean_depth", "breadth_1x"])
        w.writerow([sample_id, role, "Klebsiella pneumoniae", 1, total_len, n50,
                    n50, 57.0, depth, 0.99])


def write_verdict(path, sample_id, role, verdict="PASS", mean_depth=30.0,
                  n50=5_400_000, amr_calls=0, total_elements=0):
    doc = {
        "sample_id": sample_id, "role": role, "verdict": verdict,
        "flye_status": "ok", "amr_calls": amr_calls,
        "total_elements": total_elements, "mean_depth": mean_depth, "n50": n50,
        "total_len": 5_400_000, "n_contigs": 1,
        "amr_result_interpretable": "yes", "checks": [],
    }
    json.dump(doc, open(path, "w"))


def run_aggregator(workdir):
    src = os.path.join(workdir, "agg.py")
    open(src, "w").write(agg_source())
    r = subprocess.run([sys.executable, src], cwd=workdir,
                       capture_output=True, text=True)
    return r


def scenario_dir(samples, titration=()):
    """samples: [(sample_id, role, [(gene, etype)], n50, depth)]
    titration: same shape; these get role 'titration'."""
    d = tempfile.mkdtemp(prefix="agg-")
    for sub in ("amr", "stats", "verdict"):
        os.makedirs(os.path.join(d, sub))
    for sid, role, calls, n50, depth in list(samples) + list(titration):
        write_amr(os.path.join(d, "amr", f"{sid}.amrfinder.tsv"), sid, calls)
        n_amr = sum(1 for _g, t in calls if t.upper() == "AMR")
        if role != "caller_control":
            write_stats(os.path.join(d, "stats", f"{sid}.stats.tsv"), sid, role, n50, depth)
            write_verdict(os.path.join(d, "verdict", f"{sid}.verdict.json"), sid, role,
                          mean_depth=depth, n50=n50, amr_calls=n_amr,
                          total_elements=len(calls))
        else:
            # No reads: depth and N50 do not exist for this sample.
            write_verdict(os.path.join(d, "verdict", f"{sid}.verdict.json"), sid, role,
                          mean_depth=None, n50=None, amr_calls=n_amr,
                          total_elements=len(calls))
    return d


def check_titration_counts_amr_only():
    """recovery_fraction must be measured over resistance determinants only."""
    full = [(f"amr{i}", "AMR") for i in range(4)] + [(f"pco{i}", "STRESS") for i in range(16)]
    # At the low-depth point: 2 of 4 AMR genes survive, but every STRESS gene does.
    # Mixed counting gives 18/20 = 0.90; AMR-only gives 2/4 = 0.50.
    low = [("amr0", "AMR"), ("amr1", "AMR")] + [(f"pco{i}", "STRESS") for i in range(16)]
    d = scenario_dir(
        samples=[("KP_X", "test", full, 5_400_000, 30.0)],
        titration=[("KP_X_d5_r1", "titration", low, 53_000, 4.7)])
    r = run_aggregator(d)
    if r.returncode != 0:
        print("  FAIL  titration recovery counts resistance determinants only")
        print(f"        aggregator exited {r.returncode}: {r.stderr.strip()[-200:]}")
        return False
    rows = list(csv.DictReader(open(os.path.join(d, "depth_titration.tsv")), delimiter="\t"))
    assert len(rows) == 1, rows
    frac = float(rows[0]["recovery_fraction"])
    denom = int(rows[0]["n_genes_full_depth"])
    missed = rows[0]["genes_missed"]
    problems = []
    if denom != 4:
        problems.append(f"denominator is {denom}, expected 4 (AMR only, not 20 with STRESS)")
    if abs(frac - 0.5) > 1e-6:
        problems.append(f"recovery_fraction is {frac:.4f}, expected 0.5000 "
                        f"(mixed counting would give 0.9000)")
    if "pco0" in missed:
        problems.append(f"genes_missed names a STRESS gene: {missed}")
    if problems:
        print("  FAIL  titration recovery counts resistance determinants only")
        for p in problems:
            print(f"        {p}")
        return False
    print("  PASS  titration recovery counts resistance determinants only")
    print(f"        4 AMR + 16 STRESS at full depth -> denominator 4, recovery {frac:.4f}")
    return True


def check_clean_control_is_reported():
    """A control that returned zero calls must not be reported as absent."""
    d = scenario_dir(samples=[
        ("KP_X", "test", [("amr0", "AMR")], 5_400_000, 30.0),
        ("CALLER_CONTROL", "caller_control", [], None, None),
    ])
    r = run_aggregator(d)
    if r.returncode != 0:
        print("  FAIL  a clean caller-level control is reported, not erased")
        print(f"        aggregator exited {r.returncode}: {r.stderr.strip()[-200:]}")
        return False
    summary = open(os.path.join(d, "run_summary.md")).read()
    problems = []
    if "No caller-level control was run" in summary:
        problems.append("summary says the control was not run, but it ran and "
                        "returned zero calls (this is the PASS case)")
    if "returned **0** calls" not in summary:
        problems.append("summary does not report the zero-call result")
    if problems:
        print("  FAIL  a clean caller-level control is reported, not erased")
        for p in problems:
            print(f"        {p}")
        return False
    print("  PASS  a clean caller-level control is reported, not erased")
    print("        zero calls reads as a result, not as an absent control")
    return True


def check_absent_control_is_reported_absent():
    """The other direction: no control in the run must still say so."""
    d = scenario_dir(samples=[("KP_X", "test", [("amr0", "AMR")], 5_400_000, 30.0)])
    r = run_aggregator(d)
    if r.returncode != 0:
        print("  FAIL  an absent caller-level control is reported absent")
        print(f"        aggregator exited {r.returncode}: {r.stderr.strip()[-200:]}")
        return False
    summary = open(os.path.join(d, "run_summary.md")).read()
    if "No caller-level control was run" not in summary:
        print("  FAIL  an absent caller-level control is reported absent")
        print(f"        summary does not say the control was skipped")
        return False
    print("  PASS  an absent caller-level control is reported absent")
    return True


def check_contaminated_control_is_loud():
    """A control with calls must be named as a false positive, not summarised away."""
    d = scenario_dir(samples=[
        ("KP_X", "test", [("amr0", "AMR")], 5_400_000, 30.0),
        ("CALLER_CONTROL", "caller_control", [("blaOXA-48", "AMR")], None, None),
    ])
    r = run_aggregator(d)
    if r.returncode != 0:
        print("  FAIL  a contaminated caller-level control is loud")
        print(f"        aggregator exited {r.returncode}: {r.stderr.strip()[-200:]}")
        return False
    summary = open(os.path.join(d, "run_summary.md")).read()
    problems = []
    for want in ("FALSE POSITIVE", "blaOXA-48", "Investigate"):
        if want not in summary:
            problems.append(f"summary does not contain {want!r}")
    if problems:
        print("  FAIL  a contaminated caller-level control is loud")
        for p in problems:
            print(f"        {p}")
        return False
    print("  PASS  a contaminated caller-level control is loud")
    print("        names the surviving gene and says to investigate")
    return True


def check_null_depth_does_not_print_as_none():
    """The caller control has no depth or N50; those must render as n/a, not None."""
    d = scenario_dir(samples=[
        ("KP_X", "test", [("amr0", "AMR")], 5_400_000, 30.0),
        ("CALLER_CONTROL", "caller_control", [], None, None),
    ])
    r = run_aggregator(d)
    if r.returncode != 0:
        print("  FAIL  null depth and N50 render as n/a, not None")
        print(f"        aggregator exited {r.returncode}: {r.stderr.strip()[-200:]}")
        return False
    summary = open(os.path.join(d, "run_summary.md")).read()
    val = open(os.path.join(d, "validation_summary.tsv")).read()
    problems = []
    if "None" in summary:
        bad = [ln for ln in summary.splitlines() if "None" in ln][:2]
        problems.append(f"run_summary.md prints None: {bad}")
    if "None" in val:
        bad = [ln for ln in val.splitlines() if "None" in ln][:2]
        problems.append(f"validation_summary.tsv prints None: {bad}")
    if problems:
        print("  FAIL  null depth and N50 render as n/a, not None")
        for p in problems:
            print(f"        {p}")
        return False
    print("  PASS  null depth and N50 render as n/a, not None")
    return True


def check_published_counts_agree(results):
    """Every count the summaries report must be recomputable from the call table.

    The three published files can disagree without anything crashing, and that is
    exactly what shipped: validation_summary.tsv reported amr_calls=55 for an isolate
    whose call table holds 29 resistance determinants (55 is every element, including
    20 metal-tolerance genes), while run_summary.md reported 26 (distinct AMR symbols)
    and the new total_elements column was NA. Three files, three numbers, no error
    anywhere -- the verdicts were published by the pre-fix gate and the aggregator
    faithfully copied them forward.

    The other checks in this file build their own fixtures, so none of them can see a
    stale input. This one reads what is actually shipped.
    """
    import csv as _csv
    import os as _os

    label = "published counts agree with the call table"
    calls = _os.path.join(results, "amr_calls.tsv")
    summary = _os.path.join(results, "validation_summary.tsv")
    if not (_os.path.exists(calls) and _os.path.exists(summary)):
        print(f"  SKIP  {label} (no published amr_calls.tsv/validation_summary.tsv)")
        return None

    with open(calls) as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    with open(summary) as fh:
        vrows = list(_csv.DictReader(fh, delimiter="\t"))

    def etype(r):
        return (r.get("element_type") or r.get("Element type") or "").upper()

    per_amr, per_all = {}, {}
    for r in rows:
        s = r["sample_id"]
        per_all[s] = per_all.get(s, 0) + 1
        if etype(r) == "AMR":
            per_amr[s] = per_amr.get(s, 0) + 1

    bad = []
    for v in vrows:
        s = v["sample_id"]
        exp_amr = per_amr.get(s, 0)
        exp_all = per_all.get(s, 0)
        got_amr = v.get("amr_calls", "NA")
        got_all = v.get("total_elements", "NA")
        if str(got_amr) != str(exp_amr):
            bad.append(f"{s}: amr_calls={got_amr}, call table has {exp_amr} AMR rows")
        # total_elements may legitimately be absent for samples predating the column,
        # but if it is populated it has to be right.
        if str(got_all) not in ("NA", "") and str(got_all) != str(exp_all):
            bad.append(f"{s}: total_elements={got_all}, call table has {exp_all} rows")

    if bad:
        print(f"  FAIL  {label}")
        for b in bad[:8]:
            print(f"        {b}")
        return False
    print(f"  PASS  {label}")
    print(f"        {len(vrows)} samples reconcile against {len(rows)} call rows")
    return True


def main(results=None):
    print("aggregate.nf\n")
    ok = [
        check_titration_counts_amr_only(),
        check_clean_control_is_reported(),
        check_absent_control_is_reported_absent(),
        check_contaminated_control_is_loud(),
        check_null_depth_does_not_print_as_none(),
    ]
    if results:
        r = check_published_counts_agree(results)
        if r is not None:
            ok.append(r)
    print()
    if all(ok):
        print(f"all {len(ok)} checks passed")
        return 0
    print(f"{sum(1 for x in ok if not x)} of {len(ok)} checks FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else None))