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


def titration_titles_match_data(results, figdir):
    """Render fig3 against synthetic titration tables and check the titles.

    The figure is given three shapes of data in turn — recovery complete at every
    depth, recovery collapsing at low depth, and recovery high but not perfect —
    and each rendered title must describe the data it was given. A title that says
    "falls off below the depth floor" over a flat curve is worse than no title.
    """
    import csv as _csv
    import shutil as _shutil

    scenarios = {
        # (name, rows) -> substrings the title must and must not contain
        "complete recovery": (
            [(5, 1.0, 5_400_000), (10, 1.0, 5_410_000),
             (20, 1.0, 5_420_000), (40, 1.0, 5_425_000)],
            ["Every determinant is recovered"], ["falls"]),
        "collapse at low depth": (
            [(5, 0.35, 900_000), (10, 0.70, 2_100_000),
             (20, 0.95, 4_800_000), (40, 1.0, 5_425_000)],
            ["Recovery falls to 35%"], ["Every determinant"]),
        "high but imperfect": (
            [(5, 0.97, 5_100_000), (10, 0.99, 5_300_000),
             (20, 1.0, 5_400_000), (40, 1.0, 5_425_000)],
            ["97%"], ["Every determinant is recovered"]),
    }

    failures = []
    for label, (rows_, must, must_not) in scenarios.items():
        scratch = tempfile.mkdtemp(prefix="tit-")
        with open(os.path.join(scratch, "depth_titration.tsv"), "w", newline="") as fh:
            w = _csv.writer(fh, delimiter="\t")
            w.writerow(["sample_id", "parent_id", "target_depth", "replicate",
                        "realised_depth", "n50", "n_genes", "n_genes_full_depth",
                        "recovery_fraction", "genes_missed", "verdict"])
            for depth, frac, n50 in rows_:
                w.writerow([f"KP_X_d{depth}_r1", "KP_X", depth, 1, depth * 0.95,
                            n50, int(20 * frac), 20, f"{frac:.4f}", "-", "PASS"])

        g = load(results=scratch, figdir=scratch)

        # figure_titration() closes its figure before returning, so the titles have
        # to be captured at draw time rather than read off a surviving figure. Wrap
        # savefig -- it is called while the figure is still open, and it is the point
        # at which the titles are final.
        captured = []
        real_savefig = g["plt"].Figure.savefig

        def spy(self, *a, **kw):
            # get_title() reads the *centre* slot, but make_figures.py sets
            # rcParams["axes.titlelocation"] = "left", so set_title() writes the left
            # slot and get_title() returns "". Read all three slots; a check that
            # looked only at the default would report every title as missing.
            for ax in self.get_axes():
                for loc in ("left", "center", "right"):
                    txt = ax.get_title(loc=loc)
                    if txt:
                        captured.append(txt)
            return real_savefig(self, *a, **kw)

        g["plt"].Figure.savefig = spy
        try:
            g["figure_titration"]()
        finally:
            g["plt"].Figure.savefig = real_savefig

        assert captured, f"{label}: figure_titration produced no titled axes"
        joined = " | ".join(captured)

        for s in must:
            if s not in joined:
                failures.append(f"{label}: expected {s!r} in titles, got {joined!r}")
        for s in must_not:
            if s in joined:
                failures.append(f"{label}: title wrongly claims {s!r}: {joined!r}")
        _shutil.rmtree(scratch, ignore_errors=True)

    if failures:
        print("  FAIL  titration titles describe the data they were given")
        for f in failures[:3]:
            print(f"        {f[:150]}")
        return False
    print("  PASS  titration titles describe the data they were given"
          f"\n        {len(scenarios)} data shapes, each title checked against its numbers")
    return True


def main(results):
    figdir = tempfile.mkdtemp(prefix="figtest-")
    print(f"results={results}\nscratch={figdir}\n")

    # 0. Refuse to run against a results directory that has no tables in it.
    #
    #    This guard exists because the suite passed without it. make_figures.py
    #    returns early when validation_summary.tsv is absent (load_tsv returns []
    #    for a missing file), so pointing the suite at a wrong path made the
    #    baseline "render clean" and every mutation "render without complaint" --
    #    four green checks over a figure that was never drawn. A test suite that
    #    reports success when handed the wrong directory is worse than no suite,
    #    because it is trusted.
    required = ["validation_summary.tsv", "amr_calls.tsv", "assembly_metrics.tsv"]
    missing = [f for f in required if not os.path.exists(os.path.join(results, f))]
    if missing:
        print(f"  FAIL  results directory is not usable: missing {', '.join(missing)}")
        print(f"        looked in {results}")
        print("        (the figures would be skipped, and every check below would "
              "pass vacuously)")
        return 1
    print(f"  PASS  results directory has all {len(required)} required tables")

    ok = [True]

    # 1. Baseline: the real figure must render clean, or the mutation tests below
    #    prove nothing — a script that always raised would "pass" every one of them.
    g = load(results=results, figdir=figdir)
    g["figure_controls"]()
    fig1 = os.path.join(figdir, "fig1_controls_and_quality.png")
    # Not just "did not raise" -- a skipped figure also does not raise.
    assert os.path.exists(fig1) and os.path.getsize(fig1) > 20_000, (
        f"baseline produced no usable figure: exists={os.path.exists(fig1)}, "
        f"size={os.path.getsize(fig1) if os.path.exists(fig1) else 0}")
    print("  PASS  baseline figure renders with all self-checks clean"
          f"\n        wrote {os.path.getsize(fig1) // 1024} KB")
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

    # 5. The titration figure's titles must agree with the titration data.
    #    They were originally hard-coded ("Recovery ... falls off below the depth
    #    floor"), which would have survived a run showing complete recovery at every
    #    depth. A title is a claim; this checks the claim against the numbers.
    ok.append(titration_titles_match_data(results, figdir))

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
