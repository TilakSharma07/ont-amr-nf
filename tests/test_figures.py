#!/usr/bin/env python3
"""Tests for bin/make_figures.py.

Run:  python3 tests/test_figures.py <results_dir>

Why this file exists
--------------------
The figure script carries its own self-checks (no overlapping text, nothing off
canvas, every leader line ends at the point belonging to the sample it names). A
self-check is only worth the lines it occupies if it actually fails when the thing
it describes is broken — and the first version of the leader check did not. It
asked "does this leader end on *a* marker?", which is true for every possible
mis-pairing of labels to points, so it would have passed a figure in which every
sample was labelled with its neighbour's name.

These tests are mutation tests: each one breaks the figure in a specific,
plausible way and asserts that the script refuses to write it. If a check here
stops failing, the corresponding self-check has gone vacuous.

No test framework is required, so this runs in the pipeline's own environment.
"""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "bin", "make_figures.py")


def load(mutate=None, results=None, figdir=None):
    """Exec make_figures.py with an optional source mutation, without running main."""
    src = open(SRC).read()
    if mutate:
        old, new = mutate
        assert old in src, f"mutation target not found in make_figures.py: {old[:60]!r}"
        src = src.replace(old, new)
    src = src.replace('if __name__ == "__main__":', "if False:")
    g = {"__name__": "make_figures_under_test"}
    exec(compile(src, "make_figures.py", "exec"), g)
    if results:
        g["RESULTS"], g["FIGDIR"] = results, figdir
    return g


def expect_failure(name, mutate, results, figdir):
    """The mutated figure must raise, not render."""
    g = load(mutate, results, figdir)
    try:
        g["figure_controls"]()
    except (AssertionError, KeyError) as e:
        print(f"  PASS  {name}\n        caught: {str(e)[:100]}")
        return True
    print(f"  FAIL  {name}: the broken figure rendered without complaint")
    return False


def main(results):
    figdir = tempfile.mkdtemp(prefix="figtest-")
    print(f"results={results}\nscratch={figdir}\n")

    ok = []

    # 1. Baseline: the real figure must render clean, or the mutation tests below
    #    prove nothing — a script that always raised would "pass" every one of them.
    g = load(results=results, figdir=figdir)
    g["figure_controls"]()
    print("  PASS  baseline figure renders with all self-checks clean")
    ok.append(True)

    # 2. Mislabelled leader lines: label text from one sample, leader endpoint from
    #    another. This is the defect the leader check exists for. Note that merely
    #    reordering the loop does NOT reproduce it — label and coordinate travel in
    #    the same tuple — which is why the mutation moves the endpoint directly.
    ok.append(expect_failure(
        "mislabelled leader is caught",
        ('ax2.annotate(sid.replace("_", " "), xy=(d, n50), xytext=(lx, ly),',
         'ax2.annotate(sid.replace("_", " "), '
         'xy=iso_sorted[(k + 1) % len(iso_sorted)][:2], xytext=(lx, ly),'),
        results, figdir))

    # 3. A sample silently losing its label.
    ok.append(expect_failure(
        "dropped leader label is caught",
        ("        for k, (d, n50, sid, role) in enumerate(iso_sorted):",
         "        for k, (d, n50, sid, role) in enumerate(iso_sorted[:-1]):"),
        results, figdir))

    # 4. An unknown samplesheet role must stop the figure rather than render grey.
    ok.append(expect_failure(
        "unknown role is caught",
        ('    "test":             (C_TEST,   "clinical isolate",     False),',
         "    # 'test' role deliberately removed by the test suite"),
        results, figdir))

    print()
    if all(ok):
        print(f"all {len(ok)} checks passed")
        return 0
    print(f"{sum(1 for x in ok if not x)} of {len(ok)} checks FAILED")
    return 1


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: test_figures.py <results_dir>")
    sys.exit(main(os.path.abspath(sys.argv[1])))
