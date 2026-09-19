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


def control_titles_match_data():
    """Figure 1's titles were hard-coded too: "Both controls return zero calls" was
    asserted regardless of what the controls returned. That is the one claim in the
    whole repo that must never be made on faith -- a contaminated control means every
    call in the run is suspect, and a title that says "zero calls" over a dirty control
    actively hides it. Check all three branches, contamination included.
    """
    import csv as _csv
    import shutil as _shutil

    scenarios = {
        # (samples: [(sample_id, role, n_amr, n50_mb, depth)], must, must_not)
        "clean controls": (
            [("KP_A", "test", 20, 5.4, 30.0), ("NEG_DECOY", "negative_control", 0, 0.0, 0.0),
             ("CALLER_CONTROL", "caller_control", 0, 0.0, 0.0)],
            ["controls return zero calls"], ["contamination"]),
        "contaminated read-level control": (
            [("KP_A", "test", 20, 5.4, 30.0), ("NEG_DECOY", "negative_control", 3, 0.0, 0.0)],
            ["Control contamination", "NEG DECOY"], ["return zero calls"]),
        "no controls in the run": (
            [("KP_A", "test", 20, 5.4, 30.0)],
            ["no controls in this run"], ["controls return zero calls", "contamination"]),
    }

    failures = []
    for label, (samples, must, must_not) in scenarios.items():
        scratch = tempfile.mkdtemp(prefix="ctl-")
        with open(os.path.join(scratch, "validation_summary.tsv"), "w", newline="") as fh:
            w = _csv.writer(fh, delimiter="\t")
            w.writerow(["sample_id", "role", "verdict", "flye_status", "amr_calls",
                        "mean_depth", "n50", "total_len", "n_contigs",
                        "interpretable", "failed_checks"])
            for sid, role, n_amr, n50, depth in samples:
                w.writerow([sid, role, "PASS", "ok", n_amr, depth, int(n50 * 1e6),
                            int(n50 * 1e6), 1, "yes", "-"])
        with open(os.path.join(scratch, "amr_calls.tsv"), "w", newline="") as fh:
            w = _csv.writer(fh, delimiter="\t")
            w.writerow(["sample_id", "gene_symbol", "gene_name", "element_type",
                        "subtype", "drug_class", "subclass", "method",
                        "pct_identity", "pct_coverage", "contig"])
            for sid, role, n_amr, n50, depth in samples:
                for k in range(n_amr):
                    w.writerow([sid, f"gene{k}", f"gene {k}", "AMR", "POINT",
                                "BETA-LACTAM", "-", "EXACTX", "100.0", "100.0", "c1"])
        with open(os.path.join(scratch, "assembly_metrics.tsv"), "w", newline="") as fh:
            w = _csv.writer(fh, delimiter="\t")
            w.writerow(["sample_id", "role", "species", "n_contigs", "total_len",
                        "n50", "largest", "gc_percent", "mean_depth", "breadth_1x"])
            for sid, role, n_amr, n50, depth in samples:
                w.writerow([sid, role, "Klebsiella pneumoniae", 1, int(n50 * 1e6),
                            int(n50 * 1e6), int(n50 * 1e6), 57.0, depth, 0.99])

        # A real results dir also holds per-sample caller output under amr/, and
        # figure_controls() reads amr/CALLER_CONTROL.amrfinder.tsv directly: the
        # caller-level control once had no validation row, so it was appended from
        # the presence of that file. Without the file this fixture never reached that
        # branch, which is how "All 3 controls return zero calls" over two controls --
        # the row appended twice, drawn twice, counted twice -- got past this check
        # and into a committed figure. Write the file so the branch is exercised.
        if any(role == "caller_control" for _s, role, _n, _m, _d in samples):
            os.makedirs(os.path.join(scratch, "amr"), exist_ok=True)
            with open(os.path.join(scratch, "amr", "CALLER_CONTROL.amrfinder.tsv"),
                      "w", newline="") as fh:
                w = _csv.writer(fh, delimiter="\t")
                w.writerow(["Name", "Element symbol", "Type", "% Identity to reference"])

        g = load(results=scratch, figdir=scratch)
        captured = []
        ytick_labels = []
        _drawn = {s for s, _r, _n, _m, _d in samples}
        expected_labels = {s.replace("_", " ") for s in _drawn} | _drawn
        real_savefig = g["plt"].Figure.savefig

        def spy(self, *a, **kw):
            for ax in self.get_axes():
                for loc in ("left", "center", "right"):
                    txt = ax.get_title(loc=loc)
                    if txt:
                        captured.append(txt)
                # Read the drawn rows off the figure being saved. Only panel a puts
                # samples on its y axis; panel b's is assembly N50 in Mb, whose
                # repeated numeric ticks are not duplicate samples. Reading every
                # open figure instead would pick up the figures earlier checks in
                # this process left open, and report every sample as a duplicate.
                labs = [t.get_text() for t in ax.get_yticklabels() if t.get_text()]
                if any(lab in expected_labels for lab in labs):
                    ytick_labels.extend(labs)
            return real_savefig(self, *a, **kw)

        g["plt"].Figure.savefig = spy
        try:
            g["figure_controls"]()
        finally:
            g["plt"].Figure.savefig = real_savefig

        dupes = {lab for lab in ytick_labels if ytick_labels.count(lab) > 1}
        if dupes:
            failures.append(f"{label}: sample drawn more than once: {sorted(dupes)}")

        if not captured:
            failures.append(f"{label}: figure_controls produced no titled axes")
            _shutil.rmtree(scratch, ignore_errors=True)
            continue
        joined = " | ".join(captured)
        for s in must:
            if s not in joined:
                failures.append(f"{label}: expected {s!r} in titles, got {joined!r}")
        for s in must_not:
            if s in joined:
                failures.append(f"{label}: title wrongly claims {s!r}: {joined!r}")
        _shutil.rmtree(scratch, ignore_errors=True)

    if failures:
        print("  FAIL  control titles describe the data they were given")
        for f in failures[:3]:
            print(f"        {f[:150]}")
        return False
    print("  PASS  control titles describe the data they were given"
          f"\n        {len(scenarios)} data shapes incl. a contaminated control")
    return True


def _write_profile_inputs(scratch, samples, genes_per_sample):
    """Minimal amr_calls.tsv + validation_summary.tsv that figure_profile() accepts."""
    import csv as _csv
    with open(os.path.join(scratch, "amr_calls.tsv"), "w", newline="") as fh:
        w = _csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "gene_symbol", "element_type", "element_subtype",
                    "class_", "pct_identity", "pct_coverage", "method"])
        for s in samples:
            for gene in genes_per_sample[s]:
                w.writerow([s, gene, "AMR", "AMR", "BETA-LACTAM",
                            "99.50", "100.00", "BLASTX"])
    with open(os.path.join(scratch, "validation_summary.tsv"), "w", newline="") as fh:
        w = _csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "role", "verdict", "checks_failed",
                    "amr_calls", "total_elements", "mean_depth", "n50"])
        for s in samples:
            w.writerow([s, "test", "PASS", "", len(genes_per_sample[s]),
                        len(genes_per_sample[s]), "30.0", "5000000"])


def profile_title_matches_data():
    """Figure 2's title must describe the matrix it drew, not a fixed narrative."""
    import shutil as _shutil
    scenarios = {
        # one sample: nothing is shared, so the title must not claim sharing
        "single sample": (
            ["KP_X"], {"KP_X": ["blaCTX-M-15", "qnrB1", "tet(A)"]},
            ["3", "KP X"], ["shared", "lineage"]),
        # every gene in every sample: sharing is total
        "fully shared": (
            ["KP_X", "KP_Y"], {"KP_X": ["blaCTX-M-15", "qnrB1"],
                               "KP_Y": ["blaCTX-M-15", "qnrB1"]},
            ["All 2", "shared", "2 isolates"], ["of 2 resistance"]),
        # partial overlap: the count of shared genes must be stated, not implied
        "partly shared": (
            ["KP_X", "KP_Y"], {"KP_X": ["blaCTX-M-15", "qnrB1"],
                               "KP_Y": ["blaCTX-M-15", "tet(A)"]},
            ["1 of 3", "shared"], ["All 3"]),
    }
    failures = []
    for label, (samples, genes, must, must_not) in scenarios.items():
        scratch = tempfile.mkdtemp(prefix="prof-")
        _write_profile_inputs(scratch, samples, genes)
        g = load(results=scratch, figdir=scratch)

        captured = []
        real_savefig = g["plt"].Figure.savefig

        def spy(self, *a, **kw):
            for ax in self.get_axes():
                for loc in ("left", "center", "right"):
                    txt = ax.get_title(loc=loc)
                    if txt:
                        captured.append(txt)
            return real_savefig(self, *a, **kw)

        g["plt"].Figure.savefig = spy
        try:
            g["figure_profile"]()
        finally:
            g["plt"].Figure.savefig = real_savefig

        if not captured:
            failures.append(f"{label}: figure_profile produced no titled axes")
            _shutil.rmtree(scratch, ignore_errors=True)
            continue
        joined = " | ".join(captured)
        for s in must:
            if s not in joined:
                failures.append(f"{label}: expected {s!r} in title, got {joined!r}")
        for s in must_not:
            if s.lower() in joined.lower():
                failures.append(f"{label}: title wrongly claims {s!r}: {joined!r}")
        _shutil.rmtree(scratch, ignore_errors=True)

    if failures:
        print("  FAIL  figure 2 title describes the data it drew")
        for f in failures[:3]:
            print(f"        {f[:150]}")
        return False
    print("  PASS  figure 2 title describes the data it drew")
    return True


def titration_dir_does_not_overwrite_main_figures():
    """fig1 and fig2 must skip on a titration results dir, not overwrite the real ones.

    The README documents rendering twice into the same figure directory -- once from
    the main results, once from the titration results, because fig3 needs data the
    main run does not have. Every figure function used to run unconditionally, so the
    second command re-wrote fig1 and fig2 from the titration directory: one full-depth
    isolate and its subsamples, no controls. The result was a figure captioned
    "1 isolate(s) ... (no controls in this run)" sitting at the filename of the
    six-sample figure, and the run printed nothing but success. This was shipped.

    The subject of both figures is absent from that directory, so both must decline.
    """
    import shutil as _shutil
    import csv as _csv
    import os as _os

    label = "fig1/fig2 skip on a titration dir instead of overwriting"
    scratch = tempfile.mkdtemp(prefix="titdir-")
    try:
        res = _os.path.join(scratch, "results_titration")
        figdir = _os.path.join(scratch, "figs")
        _os.makedirs(res)
        _os.makedirs(figdir)

        # a titration directory: one full-depth isolate + subsample rows, no controls
        subs = ["KP_X_d5_r1", "KP_X_d10_r1", "KP_X_d20_r1", "KP_X_d40_r1"]
        with open(_os.path.join(res, "validation_summary.tsv"), "w", newline="") as fh:
            w = _csv.writer(fh, delimiter="\t")
            w.writerow(["sample_id", "role", "verdict", "checks_failed",
                        "amr_calls", "total_elements", "mean_depth", "n50"])
            w.writerow(["KP_X", "test", "PASS", "", "3", "3", "30.0", "5000000"])
            for s in subs:
                w.writerow([s, "titration", "PASS", "", "2", "2", "10.0", "200000"])
        with open(_os.path.join(res, "amr_calls.tsv"), "w", newline="") as fh:
            w = _csv.writer(fh, delimiter="\t")
            w.writerow(["sample_id", "gene_symbol", "element_type", "element_subtype",
                        "class_", "pct_identity", "pct_coverage", "method"])
            for s in ["KP_X"] + subs:
                for gene in ("blaCTX-M-15", "qnrB1"):
                    w.writerow([s, gene, "AMR", "AMR", "BETA-LACTAM",
                                "99.50", "100.00", "BLASTX"])
        with open(_os.path.join(res, "depth_titration.tsv"), "w", newline="") as fh:
            w = _csv.writer(fh, delimiter="\t")
            w.writerow(["sample_id", "target_depth", "replicate", "realised_depth",
                        "n50", "n_genes_full_depth", "n_genes_recovered",
                        "recovery_fraction", "genes_missed"])
            for d, rd, n50 in ((5, 4.8, 50000), (10, 9.7, 200000),
                               (20, 19.4, 5000000), (40, 38.9, 5100000)):
                w.writerow(["KP_X", d, 1, rd, n50, 2, 2, "1.0000", ""])

        # sentinel files standing in for real, correct figures already rendered
        sentinel = b"SENTINEL-NOT-A-PNG"
        for name in ("fig1_controls_and_quality.png", "fig2_determinant_profile.png"):
            with open(_os.path.join(figdir, name), "wb") as fh:
                fh.write(sentinel)

        g = load(results=res, figdir=figdir)
        g["figure_controls"]()
        g["figure_profile"]()
        g["figure_titration"]()

        survived = []
        for name in ("fig1_controls_and_quality.png", "fig2_determinant_profile.png"):
            with open(_os.path.join(figdir, name), "rb") as fh:
                if fh.read() != sentinel:
                    survived.append(name)
        fig3 = _os.path.join(figdir, "fig3_depth_titration.png")
        wrote_fig3 = _os.path.exists(fig3) and _os.path.getsize(fig3) > 20_000

        if survived:
            print(f"  FAIL  {label}")
            for s in survived:
                print(f"        {s} was overwritten from titration data")
            return False
        if not wrote_fig3:
            print(f"  FAIL  {label}")
            print("        fig3 was not written — the guards are too broad")
            return False
        print(f"  PASS  {label}")
        print("        fig1/fig2 untouched, fig3 rendered from the titration table")
        return True
    finally:
        _shutil.rmtree(scratch, ignore_errors=True)


def geometry_check_is_enforced():
    """verify() must raise on a violation, and must ignore undrawn tick labels.

    Two separate defects lived here. It only printed, so a real fig3 violation was
    written and reported as a success; and it counted tick labels outside the axis
    view interval, which matplotlib never paints, so it invented violations that
    could not be fixed by moving anything.
    """
    import shutil as _shutil
    scratch = tempfile.mkdtemp(prefix="geom-")
    g = load(results=scratch, figdir=scratch)
    plt = g["plt"]
    failures = []

    # (a) deliberately overlapping text must raise
    fig, ax = plt.subplots(figsize=(3, 2), dpi=100)
    ax.text(0.5, 0.5, "AAAAAAAAAAAAAAAA", ha="center", transform=ax.transAxes)
    ax.text(0.5, 0.5, "BBBBBBBBBBBBBBBB", ha="center", transform=ax.transAxes)
    try:
        g["verify"](fig, "overlap-probe")
        failures.append("overlapping text did not raise")
    except g["FigureGeometryError"]:
        pass
    plt.close(fig)

    # (b) text pushed off the canvas must raise
    fig, ax = plt.subplots(figsize=(3, 2), dpi=100)
    ax.text(-3.0, 0.5, "far off to the left", transform=ax.transAxes)
    try:
        g["verify"](fig, "offcanvas-probe")
        failures.append("off-canvas text did not raise")
    except g["FigureGeometryError"]:
        pass
    plt.close(fig)

    # (c) out-of-view tick labels must NOT raise: matplotlib does not draw them.
    #     This is the fig3 geometry -- an x range starting at 1.5 keeps the "0"
    #     tick label out of view, and a y range starting at 0 does the same to "-1".
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7, 3), dpi=300)
    a1.plot([5, 10, 20, 40], [0.70, 0.75, 1.0, 1.0])
    a1.set_xlim(1.5, 43.5)
    a1.set_ylim(-0.04, 1.12)
    a2.plot([5, 10, 20, 40], [0.05, 0.20, 5.42, 5.42])
    a2.set_xlim(1.5, 43.5)
    fig.tight_layout()
    try:
        g["verify"](fig, "undrawn-ticks-probe")
    except g["FigureGeometryError"] as e:
        failures.append(f"undrawn tick labels counted as a violation: {e}")
    plt.close(fig)

    _shutil.rmtree(scratch, ignore_errors=True)
    if failures:
        print("  FAIL  geometry check is enforced and counts only drawn text")
        for f in failures:
            print(f"        {f[:150]}")
        return False
    print("  PASS  geometry check is enforced and counts only drawn text")
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

    # 6. Same rule for figure 1, where the stakes are higher: "Both controls return
    #    zero calls" was also hard-coded, so it would have been printed over a
    #    contaminated control -- hiding the one result that invalidates the whole run.
    ok.append(control_titles_match_data())

    # 7. Same rule for figure 2, which claimed "Shared core resistance genes plus
    #    lineage-specific determinants" unconditionally -- false on a single-sample
    #    panel, where nothing is shared and there are no lineages to differ.
    ok.append(profile_title_matches_data())

    # 8. verify() must raise, not print. It spent several runs detecting a real fig3
    #    violation while the script wrote the figure and exited 0.
    ok.append(geometry_check_is_enforced())

    # 9. Rendering twice into one figure directory -- which the README documents,
    #    because fig3 needs the titration data -- must not let the second command
    #    replace fig1 and fig2 with one-isolate, no-control versions of themselves.
    ok.append(titration_dir_does_not_overwrite_main_figures())

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
