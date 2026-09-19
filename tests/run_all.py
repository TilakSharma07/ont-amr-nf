#!/usr/bin/env python3
"""Run every test in this directory and summarise.

Usage:  python3 tests/run_all.py <results_dir>

Some suites need a results directory to check the shipped tables and figures against;
the rest are self-contained and build their own fixtures. Both kinds are listed here
so that "did the tests pass" is one command rather than six, and so a newly added
suite is not silently left out of the README's list.

Exit status is non-zero if any suite fails. A suite that SKIPs (a missing tool, not a
defect) does not fail the run, but the skip is reported -- an unnoticed skip is how a
test stops testing anything.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# (filename, needs_results_dir)
SUITES = [
    ("test_module_paths.py", False),
    ("test_samplesheet.py", False),
    ("test_validation_gate.py", True),
    ("test_control_gate.py", True),
    # Mostly self-contained fixtures, but one check reconciles the published
    # summaries against the published call table, so it needs the results dir.
    ("test_aggregate.py", True),
    ("test_figures.py", True),
]


def main():
    if len(sys.argv) < 2:
        sys.exit(f"usage: {os.path.basename(__file__)} <results_dir>")
    results = os.path.abspath(sys.argv[1])
    if not os.path.isdir(results):
        sys.exit(f"not a directory: {results}")

    # Guard against the failure mode the figure tests hit: an empty or wrong directory
    # let the suite report success while checking nothing.
    required = ["validation_summary.tsv", "amr_calls.tsv", "assembly_metrics.tsv"]
    missing = [f for f in required if not os.path.exists(os.path.join(results, f))]
    if missing:
        sys.exit(f"{results} is missing {', '.join(missing)} — "
                 "pass the directory the pipeline wrote, not an empty one")

    on_disk = {f for f in os.listdir(HERE)
               if f.startswith("test_") and f.endswith(".py")}
    listed = {name for name, _ in SUITES}
    unlisted = sorted(on_disk - listed)
    ghosts = sorted(listed - on_disk)

    failed, skipped = [], []
    for name, needs_results in SUITES:
        path = os.path.join(HERE, name)
        cmd = [sys.executable, path] + ([results] if needs_results else [])
        r = subprocess.run(cmd, capture_output=True, text=True)
        out = (r.stdout or "") + (r.stderr or "")
        tag = "ok  " if r.returncode == 0 else "FAIL"
        # Distinguish a suite that ran nothing from one where a single check skipped.
        # Matching "SKIP" anywhere labelled test_control_gate.py -- 7 checks passed,
        # one skipped for want of a published control set -- as a skipped suite.
        n_skip = sum(1 for l in out.splitlines() if l.strip().startswith("SKIP"))
        n_ran = sum(1 for l in out.splitlines()
                    if l.strip().startswith(("PASS", "FAIL")))
        if n_skip and not n_ran:
            tag = "skip"
            skipped.append(name)
        elif n_skip:
            tag = "ok/-"
            skipped.append(f"{name} ({n_skip} check{'s' if n_skip > 1 else ''})")
        if r.returncode != 0:
            failed.append(name)
        summary = next((l.strip() for l in reversed(out.splitlines())
                        if "checks passed" in l or "FAILED" in l or "SKIP" in l), "")
        print(f"  [{tag}] {name:28s} {summary[:70]}")
        if r.returncode != 0:
            for line in out.splitlines():
                if line.startswith("  FAIL") or "Error" in line:
                    print(f"          {line.strip()[:100]}")

    print()
    if unlisted:
        print(f"  WARNING: {len(unlisted)} suite(s) on disk but not run by this "
              f"script: {', '.join(unlisted)}")
    if ghosts:
        print(f"  WARNING: {len(ghosts)} suite(s) listed but missing: {', '.join(ghosts)}")
    if skipped:
        print(f"  {len(skipped)} suite(s) skipped: {', '.join(skipped)}")
    if failed or unlisted or ghosts:
        print(f"{len(failed)} of {len(SUITES)} suites FAILED")
        return 1
    print(f"all {len(SUITES)} suites passed"
          + (f" ({len(skipped)} skipped)" if skipped else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())