#!/usr/bin/env python3
"""Render the result figures for an ont-amr-nf run.

Usage:  make_figures.py <results_dir> <figure_dir>

Reads only the published result tables, so it can be re-run on any completed run
without re-executing the pipeline.
"""
import csv
import os
import sys
from collections import defaultdict

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

RESULTS = sys.argv[1] if len(sys.argv) > 1 else "results"
FIGDIR = sys.argv[2] if len(sys.argv) > 2 else "figures"
os.makedirs(FIGDIR, exist_ok=True)

# --- house style -------------------------------------------------------------
BASE, MID, SMALL = 8, 7, 6
mpl.rcParams.update({
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.size": BASE, "axes.titlesize": BASE, "axes.labelsize": BASE,
    "legend.fontsize": MID, "xtick.labelsize": SMALL, "ytick.labelsize": SMALL,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlelocation": "left", "axes.titleweight": "regular",
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "legend.frameon": False, "figure.constrained_layout.use": True,
})

# One hue per role, threaded through every panel in every figure.
C_TEST = "#2166ac"      # real clinical isolate
C_DECOY = "#b2182b"     # read-level control  (alarm hue, controls only)
C_CALLER = "#d6604d"    # caller-level control
C_TIT = "#4d9221"       # titration series
C_GREY = "#9e9e9e"

# Roles come from the samplesheet, and several distinct roles are all real clinical
# isolates: `test`, plus `clonal_replicate` (second isolate of the same outbreak clone)
# and `cross_species` (a second species). Resolving role -> style in ONE place is the
# point: a per-panel colour lookup keyed only on 'test' silently rendered two genuine
# isolates in the "unknown" grey, which reads as a third category that does not exist.
ROLE_STYLE = {
    "test":             (C_TEST,   "clinical isolate",     False),
    "clonal_replicate": (C_TEST,   "clinical isolate",     False),
    "cross_species":    (C_TEST,   "clinical isolate",     False),
    "negative_control": (C_DECOY,  "read-level control",   True),
    "caller_control":   (C_CALLER, "caller-level control", True),
    "titration":        (C_TIT,    "depth titration",      False),
}


def role_style(role):
    """Style for a samplesheet role. Unknown roles are loud, never silently grey:
    a new role added to the samplesheet must be given a deliberate style here."""
    if role not in ROLE_STYLE:
        raise KeyError(
            f"role {role!r} has no style; add it to ROLE_STYLE in make_figures.py "
            f"(known: {sorted(ROLE_STYLE)})")
    return ROLE_STYLE[role]


def is_isolate(role):
    return role in ROLE_STYLE and not ROLE_STYLE[role][2] and role != "titration"


def load_tsv(name):
    path = os.path.join(RESULTS, name)
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def place_point_labels(ax, points, fig, max_iter=60):
    """Direct-label scatter points, then push labels apart until no two overlap.

    Alternating fixed offsets are not enough: two points at similar y and moderate
    x-separation still collide, because the label of the left point extends into the
    right point's label. So place, render, measure, and nudge vertically until the
    rendered boxes are disjoint (§9.1) — measured, not assumed.
    """
    anns = []
    for x, y, label, colour in points:
        a = ax.annotate(label, (x, y), xytext=(7, 4), textcoords="offset points",
                        fontsize=SMALL, color=colour, zorder=4,
                        annotation_clip=False)
        anns.append(a)

    for _ in range(max_iter):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        boxes = [a.get_window_extent(r) for a in anns]
        moved = False
        for i in range(len(anns)):
            for j in range(i + 1, len(anns)):
                if not boxes[i].overlaps(boxes[j]):
                    continue
                # Push the lower-y label down and the higher-y label up, so a label
                # always stays on the side of its own marker it started on.
                lo, hi = (i, j) if anns[i].xy[1] <= anns[j].xy[1] else (j, i)
                for idx, step in ((lo, -6), (hi, +6)):
                    dx, dy = anns[idx].xyann
                    anns[idx].xyann = (dx, dy + step)
                moved = True
        if not moved:
            break
    return anns


def assert_leaders_label_their_own_point(ax, fig, expected, tol_px=1.5):
    """Every leader must end at the point belonging to the sample it names.

    A weaker version of this check — "the leader ends on *a* marker" — is vacuous:
    if the label and point lists fall out of step, each leader still lands on some
    marker, just the wrong one, and the figure looks tidy while mislabelling every
    sample. So the check is against `expected`: a label -> data-coordinate mapping
    built from the source rows, independent of draw order.

    `expected` keys are the label strings as drawn.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    marks = []
    for ln in ax.get_lines():
        if ln.get_linestyle() in ("None", "none", "") and ln.get_marker() not in ("", "None"):
            marks.extend(ax.transData.transform(ln.get_xydata()).tolist())

    checked = 0
    for a in ax.texts:
        # ax.texts mixes plain Text (captions) with Annotation; only the latter has a
        # leader, and only annotations given arrowprops actually draw one.
        if not isinstance(a, mpl.text.Annotation) or a.arrowprops is None:
            continue
        label = a.get_text()
        assert label in expected, f"leader label {label!r} is not in the expected mapping"
        want = ax.transData.transform(expected[label])
        got = ax.transData.transform(a.xy)
        off = ((want[0] - got[0]) ** 2 + (want[1] - got[1]) ** 2) ** 0.5
        assert off <= tol_px, (
            f"leader for {label!r} ends {off:.1f}px from that sample's own point "
            f"(label/point lists are out of step)")
        bb = mpl.text.Text.get_window_extent(a, r)
        on_top = [(mx, my) for mx, my in marks if bb.contains(mx, my)]
        assert not on_top, f"label {label!r} is drawn on top of a marker"
        checked += 1
    assert checked == len(expected), (
        f"{checked} leader(s) drawn but {len(expected)} expected — a sample lost its label")
    return checked


def panel_letter(ax, letter):
    ax.text(-0.10, 1.12, letter, transform=ax.transAxes,
            fontsize=BASE + 2, fontweight="bold", va="bottom", ha="left")


def text_extent(t, r):
    """The glyph box only.

    `Annotation.get_window_extent()` returns the union of the text box and the
    leader-line patch, so a leader-labelled point reports a box hundreds of pixels
    wide, and every leader converging on a cluster "overlaps" every other one. That
    is a permanent false positive, and a check that always fires cannot tell a real
    collision from its own noise. Measure what the reader sees: the text.
    """
    try:
        return mpl.text.Text.get_window_extent(t, r)
    except Exception:
        return t.get_window_extent(r)


class FigureGeometryError(RuntimeError):
    """A rendered figure has overlapping or off-canvas text."""


def verify(fig, tag):
    """Geometric check: no text-text or text-spine overlap, nothing off-canvas.

    Raises FigureGeometryError on a violation, so a broken figure fails the run
    instead of being written and reported as fine.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()

    # A tick label whose tick falls outside its axis view interval still reports
    # get_visible() == True and still has a window extent, but matplotlib never paints
    # it — hiding those labels leaves the rendered PNG byte-for-byte identical. Counting
    # them produced phantom violations: fig3 reported an overlap between an x-tick "0"
    # and a y-tick "-1" on a panel whose x range starts at 1.5 and whose y range starts
    # at 0, so neither label was ever on the canvas. Drop them before measuring.
    undrawn = set()
    for ax in fig.axes:
        for axis, lim in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
            lo, hi = sorted(lim)
            for locs, labels in ((axis.get_majorticklocs(), axis.get_majorticklabels()),
                                 (axis.get_minorticklocs(), axis.get_minorticklabels())):
                undrawn.update(lab for loc, lab in zip(locs, labels)
                               if not (lo <= loc <= hi))

    texts = [(t, text_extent(t, r)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible() and t not in undrawn]
    spines = [(s, s.get_window_extent(r)) for ax in fig.axes
              for s in ax.spines.values() if s.get_visible()]
    ticks = {ax: set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"))
             for ax in fig.axes}
    bad = [(a.get_text(), b.get_text())
           for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:] if ba.overlaps(bb)]
    bad += [(t.get_text(), "spine") for t, bt in texts for s, bs in spines
            if bt.overlaps(bs) and t not in ticks.get(s.axes, set())]
    fb = fig.bbox
    off = [t.get_text() for t, bt in texts
           if bt.x0 < fb.x0 - 1 or bt.y0 < fb.y0 - 1
           or bt.x1 > fb.x1 + 1 or bt.y1 > fb.y1 + 1]
    print(f"[verify:{tag}] overlaps={len(bad)} offcanvas={len(off)}"
          + (f" {bad[:4]}" if bad else "") + (f" {off[:4]}" if off else ""))
    if bad or off:
        for t_, bt in texts:
            if any(t_.get_text() in pair for pair in bad) or t_.get_text() in off:
                ax_ = getattr(t_, "axes", None)
                panel = fig.axes.index(ax_) if ax_ in fig.axes else "fig"
                print(f"    {t_.get_text()!r:22s} panel={panel} "
                      f"bbox=({bt.x0:.0f},{bt.y0:.0f})-({bt.x1:.0f},{bt.y1:.0f})")
        # Raise rather than print. This check existed for several runs while only
        # printing, and a real fig3 violation went unnoticed because the script still
        # wrote the figure and exited 0 — a self-check with no failure path is a comment.
        raise FigureGeometryError(
            f"{tag}: {len(bad)} overlapping text pair(s), {len(off)} off-canvas label(s)")


# =============================================================================
# Figure 1 — what the pipeline produced, and whether the controls behaved
# =============================================================================
def figure_controls():
    val = load_tsv("validation_summary.tsv")
    amr = load_tsv("amr_calls.tsv")
    if not val:
        return

    # Skip on a titration results directory. The README documents rendering twice into
    # one figure directory -- once from the main results, once from the titration
    # results, because fig3 needs data the main run does not have -- and this function
    # writes a fixed filename, so the second command silently replaced the real
    # six-sample figure with the titration parent isolate alone, captioned "no controls
    # in this run". figure_titration() already skips when handed no titration table;
    # this is the same courtesy in the other direction.
    #
    # The test is deliberately narrower than "has no controls": --make_decoy false
    # --caller_control false is a legitimate way to run the pipeline, and in that case
    # the absence has to be visible in the figure rather than hidden by a missing one.
    # A titration directory is the specific case of titration-role rows with no
    # controls; a combined run has controls and still renders.
    have_controls = any(v["role"] in ("negative_control", "caller_control") for v in val) \
        or os.path.exists(os.path.join(RESULTS, "amr", "CALLER_CONTROL.amrfinder.tsv"))
    if any(v["role"] == "titration" for v in val) and not have_controls:
        print("titration results dir, not a controls panel; skipping fig1")
        return

    # Count from amr_calls.tsv by element type rather than trusting a single summary
    # column: AMRFinderPlus reports resistance determinants alongside STRESS
    # (biocide/metal/heat) and VIRULENCE elements, and conflating them roughly doubles
    # the apparent number of resistance genes. The axis says "resistance determinants",
    # so only element_type == AMR may be counted under it.
    n_amr = defaultdict(int)
    n_other = defaultdict(int)
    for r in amr:
        if r["element_type"].upper() == "AMR":
            n_amr[r["sample_id"]] += 1
        else:
            n_other[r["sample_id"]] += 1

    rows = [(v["sample_id"], v["role"], n_amr[v["sample_id"]]) for v in val
            if v["role"] != "titration"]

    # The caller-level control has no validation row: it bypasses the assembly gate by
    # design (it IS an assembly, so assembly checks are meaningless for it). It must
    # still appear here — a control that is absent from the figure proves nothing.
    if os.path.exists(os.path.join(RESULTS, "amr", "CALLER_CONTROL.amrfinder.tsv")):
        rows.append(("CALLER_CONTROL", "caller_control", n_amr["CALLER_CONTROL"]))

    order = {"test": 0, "clonal_replicate": 1, "cross_species": 2,
             "negative_control": 3, "caller_control": 4}
    rows.sort(key=lambda r: (order.get(r[1], 9), r[0]))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.1, 3.0))

    # -- panel a: determinants called, per sample -----------------------------
    for i, (sid, role, n) in enumerate(rows):
        col, _lab, _ctl = role_style(role)
        ax1.plot([0, n], [i, i], color=col, lw=0.8, zorder=1, alpha=0.55)
        # A zero-valued point still needs a visible mark at the baseline, otherwise the
        # most important result in the figure — an empty control — renders as nothing.
        ax1.plot([n], [i], "o", color=col, ms=5.5, zorder=3)
        ax1.annotate(f"{n}", (n, i), xytext=(6, 0), textcoords="offset points",
                     va="center", fontsize=MID, color=col, zorder=4)
        # Open marker: elements the caller reported that are NOT resistance
        # determinants. Drawn so "excluded" is visible rather than silently dropped.
        tot = n + n_other[sid]
        if n_other[sid]:
            ax1.plot([tot], [i], "o", mfc="none", mec=col, mew=0.9, ms=5.5, zorder=2)
    ax1.set_yticks(range(len(rows)))
    ax1.set_yticklabels([r[0].replace("_", " ") for r in rows])
    ax1.invert_yaxis()
    ax1.set_xlabel("elements called (filled = resistance determinants)")
    # Same rule as the titration titles: state what these rows actually show. This
    # figure is rendered for any run, including ones launched with --make_decoy false
    # and --caller_control false, where "both controls return zero calls" would be a
    # claim about samples that are not on the axis.
    ctl_rows = [r for r in rows if r[1] in ("negative_control", "caller_control")]
    iso_rows = [r for r in rows if is_isolate(r[1])]
    ctl_clean = [r for r in ctl_rows if r[2] == 0]
    if ctl_rows and len(ctl_clean) == len(ctl_rows) and iso_rows:
        word = "Both controls" if len(ctl_rows) == 2 else (
            "The control" if len(ctl_rows) == 1 else f"All {len(ctl_rows)} controls")
        t1 = (f"{word} return zero calls;\n"
              f"every isolate returns a full determinant set")
    elif ctl_rows and ctl_clean != ctl_rows:
        dirty = [r[0] for r in ctl_rows if r[2] > 0]
        t1 = ("Control contamination: "
              + ", ".join(s.replace("_", " ") for s in dirty)
              + "\nreturned calls and must not be reported")
    else:
        lo = min(r[2] for r in iso_rows) if iso_rows else 0
        hi = max(r[2] for r in iso_rows) if iso_rows else 0
        span = f"{lo}" if lo == hi else f"{lo}-{hi}"
        t1 = (f"{len(iso_rows)} isolate(s), {span} resistance determinants each\n"
              "(no controls in this run)")
    ax1.set_title(t1)
    xmax = max(r[2] + n_other[r[0]] for r in rows)
    ax1.set_xlim(-3, xmax * 1.22)
    ax1.set_xticks([x for x in range(0, int(xmax) + 1, 10)])
    ax1.margins(y=0.13)
    ax1.spines["left"].set_visible(False)
    ax1.tick_params(axis="y", length=0)

    # Legend keyed on the roles actually drawn, deduplicated by label so the two
    # isolate roles share one entry.
    seen, handles = set(), []
    for _sid, role, _n in rows:
        col, lab, _ctl = role_style(role)
        if lab not in seen:
            seen.add(lab)
            handles.append(plt.Line2D([], [], marker="o", ls="", ms=5.5,
                                      color=col, label=lab))
    handles.append(plt.Line2D([], [], marker="o", ls="", ms=5.5, mfc="none",
                              mec=C_GREY, mew=0.9, color=C_GREY,
                              label="+ stress / virulence\n(not resistance)"))
    ax1.legend(handles=handles, loc="lower right", handletextpad=0.4,
               borderaxespad=0.1)
    panel_letter(ax1, "a")

    # -- panel b: assembly quality vs. the interpretability floor -------------
    stats = [s for s in load_tsv("assembly_metrics.tsv")
             if s["role"] != "titration"]
    if stats:
        pts = sorted(((float(s["mean_depth"]), int(s["n50"]) / 1e6,
                       s["sample_id"], s["role"]) for s in stats),
                     key=lambda q: q[0])
        for d, n50, sid, role in pts:
            col, _lab, _ctl = role_style(role)
            ax2.plot([d], [n50], "o", color=col, ms=6, mec="white", mew=0.5, zorder=3)

        # The read-level control sits at the origin and the isolates cluster in the
        # top-right, so ~60% of the panel is empty. Rather than cram four labels into
        # the cluster (where they collide with each other and with the markers), stack
        # them in that empty region and connect each with a thin leader (§6.9). The
        # leader endpoints are asserted to land on their own marker after rendering.
        iso = [q for q in pts if not role_style(q[3])[2]]
        ctl = [q for q in pts if role_style(q[3])[2]]
        xlo, xhi = ax2.get_xlim()
        ylo, yhi = ax2.get_ylim()
        lx = xlo + 0.30 * (xhi - xlo)
        iso_sorted = sorted(iso, key=lambda q: -q[1])
        # Built from the source rows, so it is not derived from the drawing loop it checks.
        expect_leaders = {sid.replace("_", " "): (d, n50) for d, n50, sid, _r in iso}
        for k, (d, n50, sid, role) in enumerate(iso_sorted):
            col = role_style(role)[0]
            ly = yhi - (0.09 + 0.115 * k) * (yhi - ylo)
            ax2.annotate(sid.replace("_", " "), xy=(d, n50), xytext=(lx, ly),
                         textcoords="data", ha="right", va="center",
                         fontsize=SMALL, color=col, zorder=4,
                         arrowprops=dict(arrowstyle="-", lw=0.5, color=col,
                                         shrinkA=1.0, shrinkB=3.5,
                                         connectionstyle="arc3,rad=0.0"))
        # Controls are labelled at a fixed offset from their marker, but a run where
        # both controls behave correctly puts them at the same point -- zero calls,
        # zero depth, zero N50 -- so a single offset stacks the two labels on top of
        # each other. That is the expected case, not an edge case. Offset each label
        # within a coincident group so the healthy run is the readable one.
        ctl_groups = defaultdict(list)
        for d, n50, sid, role in ctl:
            ctl_groups[(round(d, 3), round(n50, 3))].append((sid, role))
        for (d, n50), members in ctl_groups.items():
            for k, (sid, role) in enumerate(sorted(members)):
                col = role_style(role)[0]
                ax2.annotate(sid.replace("_", " "), (d, n50),
                             xytext=(8, 2 + 9 * k), textcoords="offset points",
                             fontsize=SMALL, color=col, zorder=4)
        floor = float(os.environ.get("MIN_DEPTH_X", 20))
        ax2.axvline(floor, color=C_GREY, lw=0.7, ls="--", zorder=1)
        ax2.annotate(f"depth floor ({floor:.0f}x)\nabsence not reported below this",
                     (floor, 0.02), xytext=(5, 0), textcoords="offset points",
                     fontsize=SMALL, color=C_GREY, va="bottom", zorder=2)
        ax2.set_xlabel("realised depth (x), measured by remapping")
        ax2.set_ylabel("assembly N50 (Mb)")
        # Derived, not asserted. "near-complete chromosomes" is a claim about N50 and
        # "does not assemble" is a claim about the decoy, and a run without a decoy
        # (--make_decoy false) has neither on the axis.
        n50s_iso = [q[1] for q in iso]
        lo_n50, hi_n50 = (min(n50s_iso), max(n50s_iso)) if n50s_iso else (0, 0)
        span = f"{lo_n50:.1f}" if abs(hi_n50 - lo_n50) < 0.05 else f"{lo_n50:.1f}-{hi_n50:.1f}"
        if iso and ctl:
            t2 = (f"Isolates assemble to {span} Mb N50;\n"
                  "the read-level control does not assemble")
        elif iso:
            t2 = (f"Isolates assemble to {span} Mb N50\n"
                  "(no read-level control in this run)")
        else:
            t2 = "No isolate assemblies in this run"
        ax2.set_title(t2)
        dmax = max(q[0] for q in pts)
        ax2.set_xlim(-4, dmax * 1.30)
        ax2.set_xticks([x for x in range(0, int(dmax) + 5, 10)])
        ax2.set_ylim(-0.6, max(q[1] for q in pts) * 1.18)
        ax2.set_yticks([y for y in range(0, int(max(q[1] for q in pts)) + 2)])
        ax2.margins(y=0.20)
        ax2.text(1.0, -0.20, "higher = better", transform=ax2.transAxes,
                 ha="right", va="top", fontsize=SMALL, color=C_GREY)
    n_leaders = assert_leaders_label_their_own_point(ax2, fig, expect_leaders)
    panel_letter(ax2, "b")

    out = os.path.join(FIGDIR, "fig1_controls_and_quality.png")
    fig.savefig(out)
    verify(fig, "fig1")
    print(f"[verify:fig1] {n_leaders} leader label(s) land on their own marker")
    plt.close(fig)
    print("wrote", out)


# =============================================================================
# Figure 2 — determinant profile across isolates
# =============================================================================
def figure_profile():
    amr = load_tsv("amr_calls.tsv")
    val = load_tsv("validation_summary.tsv")
    if not amr:
        return
    keep_roles = {v["sample_id"]: v["role"] for v in val}
    # AMR determinants only: the STRESS/metal-tolerance calls are real output but
    # they are not antimicrobial resistance and would pad the panel misleadingly.
    rows = [r for r in amr
            if r["element_type"] == "AMR"
            and is_isolate(keep_roles.get(r["sample_id"], ""))]
    if not rows:
        return

    # Same overwrite hazard as fig1: a titration directory contains one full-depth
    # isolate plus its subsamples, so this would draw a single-column "profile" over
    # the real multi-isolate one. The profile is a comparison; one isolate is not one.
    if len({r["sample_id"] for r in rows}) < 2 and any(
            v["role"] == "titration" for v in val):
        print("titration results dir, not an isolate panel; skipping fig2")
        return

    samples = sorted({r["sample_id"] for r in rows})
    genes = sorted({r["gene_symbol"] for r in rows})
    ident = {(r["sample_id"], r["gene_symbol"]): float(r["pct_identity"]) for r in rows}

    M = np.full((len(genes), len(samples)), np.nan)
    for i, g in enumerate(genes):
        for j, s in enumerate(samples):
            if (s, g) in ident:
                M[i, j] = ident[(s, g)]

    h = max(3.0, 0.155 * len(genes) + 1.0)
    fig, ax = plt.subplots(figsize=(0.95 * len(samples) + 3.1, h))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#f2f2f2")
    im = ax.imshow(M, cmap=cmap, vmin=90, vmax=100, aspect="auto")

    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels([s.replace("_", " ") for s in samples], rotation=30, ha="right")
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=SMALL, style="italic")
    # Derive the title from the matrix. This said "Shared core resistance genes plus
    # lineage-specific determinants" unconditionally, which is false whenever the panel
    # holds a single sample — as it does on a titration run, where there is nothing to
    # share and no lineages to differ. Count what is actually in M instead.
    n_core = int(sum(1 for i in range(len(genes)) if not np.isnan(M[i]).any()))
    if len(samples) == 1:
        title = (f"{len(genes)} resistance determinants in "
                 f"{samples[0].replace('_', ' ')}")
    elif n_core == len(genes):
        title = (f"All {len(genes)} resistance determinants shared "
                 f"across {len(samples)} isolates")
    else:
        title = (f"{n_core} of {len(genes)} resistance determinants shared "
                 f"across {len(samples)} isolates")
    # A long title on a narrow figure (width scales with sample count) ran off the
    # canvas; wrapping keeps it inside the figure instead of relying on it being short.
    ax.set_title(title, wrap=True)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)

    # Small enough to read: print the value in every occupied cell.
    if M.size <= 220:
        for i in range(len(genes)):
            for j in range(len(samples)):
                if not np.isnan(M[i, j]):
                    v = M[i, j]
                    ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                            fontsize=SMALL - 0.5,
                            color="white" if v < 97 else "#1a1a1a")

    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("% identity to reference\n(grey = not detected)", fontsize=MID)
    cb.ax.tick_params(labelsize=SMALL)
    cb.outline.set_visible(False)

    out = os.path.join(FIGDIR, "fig2_determinant_profile.png")
    fig.savefig(out)
    verify(fig, "fig2")
    plt.close(fig)
    print("wrote", out)


# =============================================================================
# Figure 3 — depth titration: where the method stops working
# =============================================================================
def figure_titration():
    tit = load_tsv("depth_titration.tsv")
    if not tit:
        print("no titration table; skipping fig3")
        return

    by_depth = defaultdict(list)
    for r in tit:
        by_depth[int(r["target_depth"])].append(r)
    depths = sorted(by_depth)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.1, 2.9), sharex=True)

    # -- panel a: recovery fraction vs depth ----------------------------------
    med = [np.median([float(r["recovery_fraction"]) for r in by_depth[d]]) for d in depths]
    for d in depths:
        for r in by_depth[d]:
            ax1.plot([d], [float(r["recovery_fraction"])], "o", color=C_TIT,
                     ms=4, alpha=0.45, zorder=2)
    ax1.plot(depths, med, "-", color=C_TIT, lw=1.3, zorder=3)
    ax1.plot(depths, med, "s", color=C_TIT, ms=5, mec="white", mew=0.6, zorder=4)
    ax1.axhline(1.0, color=C_GREY, lw=0.7, ls="--", zorder=1)
    # Offset the reference-line label away from the data, whichever side that is.
    dy = -11 if min(med) > 0.9 else 4
    ax1.annotate("complete recovery", (depths[0], 1.0), xytext=(2, dy),
                 textcoords="offset points", fontsize=SMALL, color=C_GREY)
    ax1.set_ylim(-0.04, 1.12)
    ax1.set_ylabel("determinants recovered\n(fraction of full-depth set)")
    ax1.set_xlabel("target depth (x)")
    # Title states what the data shows, decided from the data. A pre-written
    # conclusion would have survived a result that contradicted it.
    lowest, highest = med[0], med[-1]
    if lowest >= 0.999:
        t1 = ("Every determinant is recovered at all depths tested\n"
              f"(down to {depths[0]}x)")
    elif lowest >= 0.95:
        t1 = (f"Recovery stays above {lowest:.0%} down to {depths[0]}x,\n"
              "so the depth floor is conservative here")
    else:
        t1 = (f"Recovery falls to {lowest:.0%} at {depths[0]}x\n"
              f"from {highest:.0%} at {depths[-1]}x")
    ax1.set_title(t1)
    ax1.margins(x=0.10)
    panel_letter(ax1, "a")

    # -- panel b: contiguity vs depth -----------------------------------------
    n50s = [np.median([int(r["n50"]) / 1e6 for r in by_depth[d]]) for d in depths]
    for d in depths:
        for r in by_depth[d]:
            ax2.plot([d], [int(r["n50"]) / 1e6], "o", color=C_TIT, ms=4, alpha=0.45, zorder=2)
    ax2.plot(depths, n50s, "-", color=C_TIT, lw=1.3, zorder=3)
    ax2.plot(depths, n50s, "s", color=C_TIT, ms=5, mec="white", mew=0.6, zorder=4)
    ax2.set_ylabel("assembly N50 (Mb)")
    ax2.set_xlabel("target depth (x)")
    # Same rule, and the causal claim is only made when both quantities move.
    n50_drop = (n50s[-1] - n50s[0]) / n50s[-1] if n50s[-1] else 0.0
    if n50_drop >= 0.25 and med[0] < 0.999:
        t2 = (f"N50 falls {n50_drop:.0%} from {depths[-1]}x to {depths[0]}x,\n"
              "which is what drives the missed calls")
    elif n50_drop >= 0.25:
        t2 = (f"N50 falls {n50_drop:.0%} from {depths[-1]}x to {depths[0]}x\n"
              "without costing any determinant calls")
    else:
        t2 = (f"Contiguity is stable across the range\n"
              f"({n50s[0]:.1f}-{n50s[-1]:.1f} Mb N50)")
    ax2.set_title(t2)
    ax2.margins(0.10)
    ax2.text(1.0, -0.22, "higher = better", transform=ax2.transAxes,
             ha="right", va="top", fontsize=SMALL, color=C_GREY)
    panel_letter(ax2, "b")

    out = os.path.join(FIGDIR, "fig3_depth_titration.png")
    fig.savefig(out)
    verify(fig, "fig3")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    figure_controls()
    figure_profile()
    figure_titration()
