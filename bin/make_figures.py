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


def load_tsv(name):
    path = os.path.join(RESULTS, name)
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def panel_letter(ax, letter):
    ax.text(-0.02, 1.06, letter, transform=ax.transAxes,
            fontsize=BASE + 2, fontweight="bold", va="bottom", ha="right")


def verify(fig, tag):
    """Geometric check: no text-text or text-spine overlap, nothing off-canvas."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible()]
    spines = [(s, s.get_window_extent(r)) for ax in fig.axes
              for s in ax.spines.values() if s.get_visible()]
    ticks = {ax: set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"))
             for ax in fig.axes}
    bad = [(a.get_text(), b.get_text())
           for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:] if ba.overlaps(bb)]
    bad += [(t.get_text(), "spine") for t, bt in texts for s, bs in spines
            if bt.overlaps(bs) and t not in ticks.get(s.axes, set())]
    off = [t.get_text() for t, bt in texts if not fig.bbox.contains(bt.x0, bt.y0)]
    print(f"[verify:{tag}] overlaps={len(bad)} offcanvas={len(off)}"
          + (f" {bad[:4]}" if bad else "") + (f" {off[:4]}" if off else ""))


# =============================================================================
# Figure 1 — what the pipeline produced, and whether the controls behaved
# =============================================================================
def figure_controls():
    val = load_tsv("validation_summary.tsv")
    amr = load_tsv("amr_calls.tsv")
    if not val:
        return

    genes = defaultdict(set)
    for r in amr:
        genes[r["sample_id"]].add(r["gene_symbol"])

    # Caller control has no verdict row (it bypasses the assembly gate by design),
    # so add it from the call table to keep both controls on one axis.
    rows = [(v["sample_id"], v["role"], int(v["amr_calls"]), float(v["mean_depth"]))
            for v in val]
    if "CALLER_CONTROL" in genes or any(r["sample_id"] == "CALLER_CONTROL" for r in amr):
        n_cc = sum(1 for r in amr if r["sample_id"] == "CALLER_CONTROL")
        rows.append(("CALLER_CONTROL", "caller_control", n_cc, float("nan")))

    order = {"test": 0, "negative_control": 1, "caller_control": 2, "titration": 3}
    rows = [r for r in rows if r[1] != "titration"]
    rows.sort(key=lambda r: (order.get(r[1], 9), r[0]))

    colours = {"test": C_TEST, "negative_control": C_DECOY, "caller_control": C_CALLER}
    labels = {"test": "clinical isolate", "negative_control": "read-level control",
              "caller_control": "caller-level control"}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.1, 2.9))

    # -- panel a: determinants called, per sample -----------------------------
    y = np.arange(len(rows))
    for i, (sid, role, n, _d) in enumerate(rows):
        col = colours.get(role, C_GREY)
        ax1.plot([0, n], [i, i], color=col, lw=0.8, zorder=1, alpha=0.6)
        ax1.plot([n], [i], "o", color=col, ms=5.5, zorder=2)
        # A zero-valued point still needs a visible mark at the baseline.
        ax1.annotate(f"{n}", (n, i), xytext=(5, 0), textcoords="offset points",
                     va="center", fontsize=MID, color=col)
    ax1.set_yticks(y)
    ax1.set_yticklabels([r[0].replace("_", " ") for r in rows])
    ax1.invert_yaxis()
    ax1.set_xlabel("resistance determinants called")
    ax1.set_title("Both controls return zero calls;\nevery isolate returns a full determinant set")
    ax1.margins(x=0.16, y=0.10)
    ax1.spines["left"].set_visible(False)
    ax1.tick_params(axis="y", length=0)

    seen = []
    for role in ("test", "negative_control", "caller_control"):
        if any(r[1] == role for r in rows):
            seen.append(plt.Line2D([], [], marker="o", ls="", ms=5.5,
                                   color=colours[role], label=labels[role]))
    ax1.legend(handles=seen, loc="lower right", handletextpad=0.4)
    panel_letter(ax1, "a")

    # -- panel b: assembly quality vs. the interpretability floor -------------
    stats = load_tsv("assembly_metrics.tsv")
    st = [s for s in stats if s["role"] in ("test", "negative_control")]
    if st:
        for s in st:
            col = colours.get(s["role"], C_GREY)
            d, n50 = float(s["mean_depth"]), int(s["n50"]) / 1e6
            ax2.plot([d], [n50], "o", color=col, ms=6,
                     mec="white", mew=0.5, zorder=3)
            ax2.annotate(s["sample_id"].replace("_", " "), (d, n50),
                         xytext=(4, 4), textcoords="offset points",
                         fontsize=SMALL, color=col)
        ax2.axvline(20, color=C_GREY, lw=0.7, ls="--", zorder=1)
        ax2.annotate("depth floor (20x):\nabsence not reported below this",
                     (20, 0.06), xytext=(6, 0), textcoords="offset points",
                     fontsize=SMALL, color=C_GREY, va="bottom")
        ax2.set_xlabel("realised depth (x), measured by remapping")
        ax2.set_ylabel("assembly N50 (Mb)")
        ax2.set_title("Isolates assemble to near-complete chromosomes;\nthe read-level control does not assemble")
        ax2.margins(0.14)
        ax2.text(1.0, -0.22, "higher = better", transform=ax2.transAxes,
                 ha="right", va="top", fontsize=SMALL, color=C_GREY)
    panel_letter(ax2, "b")

    out = os.path.join(FIGDIR, "fig1_controls_and_quality.png")
    fig.savefig(out)
    verify(fig, "fig1")
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
            and keep_roles.get(r["sample_id"]) == "test"]
    if not rows:
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
    cmap = mpl.cm.get_cmap("viridis").copy()
    cmap.set_bad("#f2f2f2")
    im = ax.imshow(M, cmap=cmap, vmin=90, vmax=100, aspect="auto")

    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels([s.replace("_", " ") for s in samples], rotation=30, ha="right")
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=SMALL, style="italic")
    ax.set_title("Shared core resistance genes plus lineage-specific determinants")
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
    ax1.annotate("complete recovery", (depths[0], 1.0), xytext=(2, -11),
                 textcoords="offset points", fontsize=SMALL, color=C_GREY)
    ax1.set_ylim(-0.04, 1.12)
    ax1.set_ylabel("determinants recovered\n(fraction of full-depth set)")
    ax1.set_xlabel("target depth (x)")
    ax1.set_title("Recovery of the full-depth determinant set\nfalls off below the depth floor")
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
    ax2.set_title("Contiguity collapses at low depth,\nwhich is what drives the missed calls")
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
