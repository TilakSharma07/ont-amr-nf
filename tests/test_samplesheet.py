#!/usr/bin/env python3
"""Tests for samplesheet validation in main.nf.

Run:  python3 tests/test_samplesheet.py

Why this file exists
--------------------
An unvalidated samplesheet does not fail — it produces a wrong result that looks
like a right one. Two cases were found by hand-feeding a malformed sheet:

  * A row with an empty sample_id produced publish files named ".verdict.json" and
    ".amrfinder.tsv" — dotfiles, invisible to `ls`, and the run still printed
    "completed : OK".
  * A role typo ('negatve_control') does not match the control branch of the
    validation gate, so the decoy would be held to a POSITIVE expectation: the one
    sample that must come back empty would be scored as though it had to come back
    full, and its emptiness would be reported as a failure while a contaminated
    decoy passed.

Both now abort at parse time. These tests drive the real pipeline with malformed
sheets, because the validation lives in Groovy inside the workflow and the thing
worth testing is whether the run actually refuses.

Each case is a stub run, so no assembly or AMR calling happens.
"""
import copy
import csv
import io
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SHEET = os.path.join(ROOT, "assets", "samplesheet.csv")


def load_sheet():
    raw = open(SHEET, newline="").read()
    rows = list(csv.DictReader(io.StringIO(raw)))
    return list(rows[0].keys()), rows


def write_sheet(path, fields, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def run_pipeline(sheet, tmp, tag):
    """Stub-run the workflow against `sheet`; return (completed_ok, message)."""
    outdir = os.path.join(tmp, f"out_{tag}")
    work = os.path.join(tmp, f"work_{tag}")
    cmd = ["nextflow", "run", "main.nf", "-profile", "local_env", "-stub-run",
           "--samplesheet", sheet, "--make_decoy", "false",
           "--caller_control", "false", "--outdir", outdir, "-work-dir", work]
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=600)
    txt = (r.stdout or "") + (r.stderr or "")
    complaint = next((l.strip() for l in txt.splitlines() if "Samplesheet" in l), "")
    return ("completed : OK" in txt), complaint


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"\n        {detail}" if detail else ""))
    return bool(cond)


def main():
    if not shutil.which("nextflow"):
        print("  SKIP  nextflow not on PATH — samplesheet validation not exercised")
        return 0

    fields, rows = load_sheet()
    tmp = tempfile.mkdtemp(prefix="sheet-")
    ok = []

    # The unmodified sheet must still run, or the validation is too strict.
    good = os.path.join(tmp, "good.csv")
    write_sheet(good, fields, rows)
    done, _ = run_pipeline(good, tmp, "good")
    ok.append(check("the shipped samplesheet still runs", done,
                    "validation must not reject valid input"))

    # A row with no sample_id: the dotfile case.
    r_noid = copy.deepcopy(rows)
    blank = {f: "" for f in fields}
    blank.update({"accession": "SRR_FAKE", "species": "Klebsiella pneumoniae",
                  "role": "test"})
    r_noid.append(blank)
    p = os.path.join(tmp, "no_id.csv")
    write_sheet(p, fields, r_noid)
    done, msg = run_pipeline(p, tmp, "no_id")
    ok.append(check("a row with no sample_id aborts the run",
                    not done and "no sample_id" in msg, msg[:110] or "no complaint emitted"))

    # A mistyped role: the silently-inverted-control case.
    r_typo = copy.deepcopy(rows)
    ctrl = next((i for i, r in enumerate(r_typo) if r["role"] == "negative_control"), None)
    assert ctrl is not None, "shipped sheet has no negative_control row to mistype"
    r_typo[ctrl]["role"] = "negatve_control"
    p = os.path.join(tmp, "typo_role.csv")
    write_sheet(p, fields, r_typo)
    done, msg = run_pipeline(p, tmp, "typo")
    ok.append(check("a mistyped role aborts the run",
                    not done and "not recognised" in msg, msg[:110] or "no complaint emitted"))

    # A duplicate id: two samples would overwrite each other's outputs.
    r_dup = copy.deepcopy(rows) + [copy.deepcopy(rows[0])]
    p = os.path.join(tmp, "dup.csv")
    write_sheet(p, fields, r_dup)
    done, msg = run_pipeline(p, tmp, "dup")
    ok.append(check("a duplicate sample_id aborts the run",
                    not done and "more than once" in msg, msg[:110] or "no complaint emitted"))

    # An id that is not filename-safe.
    r_bad = copy.deepcopy(rows)
    bad = copy.deepcopy(rows[0])
    bad["sample_id"] = "KP/../etc"
    r_bad.append(bad)
    p = os.path.join(tmp, "bad_chars.csv")
    write_sheet(p, fields, r_bad)
    done, msg = run_pipeline(p, tmp, "bad")
    ok.append(check("a sample_id that is not filename-safe aborts the run",
                    not done and "not usable as" in msg, msg[:110] or "no complaint emitted"))

    print()
    if all(ok):
        print(f"all {len(ok)} checks passed")
        return 0
    print(f"{sum(1 for x in ok if not x)} of {len(ok)} checks FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
