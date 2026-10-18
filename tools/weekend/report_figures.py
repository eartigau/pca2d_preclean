#!/usr/bin/env python
"""Les figures du rapport du week-end du 2026-10-17 (dix naines M NIRPS tirees
au hasard, et les sept cibles faibles).

    python tools/weekend/report_figures.py

Lit outputs/reports/weekend/weekend.json (analyse_weekend.py) et les rdb LBL
de l'arbre d'irrisor, ecrit les PDF dans outputs/reports/weekend/figures/.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
W = os.path.join(REPO, "outputs", "reports", "weekend")
OUT = os.path.join(W, "figures")
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "outputs", "reports", "nuit17"))
import analyse as A                                    # noqa: E402
from pca2d.lblscan import velocity_stats               # noqa: E402

A.TREE = "/Volumes/irrisor/pca2d_preclean/lbl/lblrdb"
NAMES = {"GJ1": "GJ 1", "GL699": "Barnard", "YZ_CMI": "YZ CMi", "GL406": "Wolf 359",
         "GJ643": "GJ 643", "GJ3737": "GJ 3737", "GJ3618": "GJ 3618",
         "GJ3622": "GJ 3622", "TOI5524": "TOI-5524", "TOI210": "TOI-210",
         "GJ3693": "GJ 3693", "LP_938M71": "LP 938-71", "TOI2407": "TOI-2407",
         "TOI756": "TOI-756", "SIPS_J1259M4336": "SIPS J1259", "TOI206": "TOI-206",
         "TOI782": "TOI-782"}
GROUP = {"random": ("#26a1ff", "o", "drawn at random (>100 nights)"),
         "faint": ("#e0564a", "s", "faint (EXTSN060 < 40, >50 nights)")}
FLOAT = matplotlib.ticker.FormatStrFormatter("%g")
NUL = matplotlib.ticker.NullFormatter()


def name(r):
    return NAMES.get(r["object"], r["object"])


def suffix(r):
    return "_W" if r.get("group", "random") == "random" else "_F"


def save(fig, stem):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, stem + ".pdf")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)


def plain_log(ax, xticks=None, yticks=None):
    for axis, ticks in ((ax.xaxis, xticks), (ax.yaxis, yticks)):
        if ticks is not None:
            axis.set_ticks(ticks)
        axis.set_major_formatter(FLOAT)
        axis.set_minor_formatter(NUL)


def label_points(ax, x, y, labels):
    order = np.argsort(x)
    for k, i in enumerate(order):
        ax.annotate(labels[i], (x[i], y[i]), fontsize=7, color="0.25",
                    xytext=(4, 6 if k % 2 == 0 else -10), textcoords="offset points")


def fig_predictor(t):
    """Le gain contre l'exces de dispersion du livre, et contre le SNR."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), sharey=True)
    spec = [("excess_rms", "excess scatter of the delivered series\n"
             r"$\sqrt{\mathrm{rms}^2-\sigma^2}/\sigma$", [1, 2, 5, 10, 20]),
            ("excess_robust", "same, robust sigma in place of the rms", [1, 2, 5, 10, 20]),
            ("snr", "median SNR per exposure (fit band)", [20, 50, 100, 200])]
    for ax, (key, xlabel, ticks) in zip(axes, spec):
        pts = [r for r in t if key in r and "nominal" in r]
        x = np.array([max(r[key], 0.5) for r in pts])
        g = np.array([r["nominal"]["gain"] for r in pts])
        for grp, (col, mk, lab) in GROUP.items():
            sel = np.array([r.get("group", "random") == grp for r in pts])
            if sel.any():
                lo = np.array([r["nominal"]["gain_16_84"][0] for r in pts])[sel]
                hi = np.array([r["nominal"]["gain_16_84"][1] for r in pts])[sel]
                ax.errorbar(x[sel], g[sel], yerr=[g[sel] - lo, hi - g[sel]], fmt=mk,
                            color=col, ms=6, capsize=2, lw=1, label=lab)
        label_points(ax, x, g, [name(r) for r in pts])
        rho, p = spearmanr(x, g)
        ax.set_title(r"Spearman $\rho = %+.2f$ ($p = %.3f$, $N = %d$)" % (rho, p, len(pts)),
                     fontsize=9.5)
        ax.axhline(1.0, color="0.35", lw=0.9)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_ylim(0.8, 3.0)
        plain_log(ax, ticks, [0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0])
        ax.set_xlabel(xlabel)
        ax.grid(alpha=0.25, which="both")
    axes[0].set_ylabel("gain in RV rms, delivered / precleaned (1-3)")
    axes[0].legend(fontsize=7.5, loc="upper left")
    fig.suptitle("What predicts the gain: the scatter the delivered series carries beyond its"
                 " photon noise, not the SNR", fontsize=11)
    fig.tight_layout()
    save(fig, "f1_predictor")


def fig_removed(t):
    """Ce que la correction retire en quadrature, en m/s, contre le rms livre."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, kind, title in ((axes[0], "rms", "rms"), (axes[1], "robust", "robust sigma")):
        for grp, (col, mk, lab) in GROUP.items():
            pts = [r for r in t if "nominal" in r and r.get("group", "random") == grp]
            if not pts:
                continue
            ref = np.array([r["nominal"]["ref"]["rms" if kind == "rms" else "robust"] for r in pts])
            rem = np.array([r["nominal"]["removed_%s" % kind] for r in pts])
            ax.plot(ref, rem, mk, color=col, ms=6, label=lab)
            label_points(ax, ref, rem, [name(r) for r in pts])
        xs = np.logspace(0, 2.5, 50)
        ax.plot(xs, xs, ls="--", color="0.6", lw=0.8, label="everything removed")
        ax.axhline(0, color="0.35", lw=0.9)
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=2)
        ax.set_xlabel("delivered %s (m/s)" % title)
        ax.set_ylabel(r"removed in quadrature, $\pm\sqrt{|a^2-b^2|}$ (m/s)")
        ax.set_title("%s: what the correction removes (positive) or adds (negative)" % title,
                     fontsize=9.5)
        ax.grid(alpha=0.25, which="both")
    axes[0].legend(fontsize=7.5)
    fig.tight_layout()
    save(fig, "f2_removed")


def fig_errors(t):
    """Le cout en barres d'erreur, avec et sans les coupures du residu."""
    fig, ax = plt.subplots(figsize=(7, 4.4))
    for key, col, mk, lab in (("nominal", "#e0564a", "o", "with residual cuts (nominal)"),
                              ("cutsoff", "#26a1ff", "s", "without residual cuts")):
        pts = [r for r in t if key in r]
        x = np.array([r["snr"] for r in pts])
        y = 100 * np.array([r[key]["median_error"] / r[key]["ref"]["median_error"] - 1 for r in pts])
        ax.plot(x, y, mk, color=col, ms=6, label=lab)
        if key == "nominal":
            label_points(ax, x, y, [name(r) for r in pts])
    ax.set_xscale("log")
    plain_log(ax, [20, 50, 100, 200])
    ax.axhline(0, color="0.35", lw=0.9)
    ax.set_xlabel("median SNR per exposure")
    ax.set_ylabel("median LBL error, corrected / delivered - 1 (%)")
    ax.set_title("The cost in error bars: systematic, and twice as large with the cuts",
                 fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    save(fig, "f3_errors")


def paired_cuts(t):
    """rms(avec coupures) / rms(sans), sur les memes poses, bootstrap apparie."""
    rng = np.random.default_rng(3)
    out = []
    for r in t:
        T = r["object"] + suffix(r)
        a = A.load(T + "_PCA2D_1-3", None)
        b = A.load(T + "_PCA2D_1-3_n17_cutsoff", None)
        if a is None or b is None:
            continue
        a, b = A.on_common(a, b)
        ratio = np.std(a["v"]) / np.std(b["v"])
        boot = [np.std(a["v"][k]) / np.std(b["v"][k])
                for k in (rng.integers(0, a["v"].size, a["v"].size) for _ in range(2000))]
        lo, hi = np.percentile(boot, [16, 84])
        out.append(dict(object=r["object"], name=name(r), group=r.get("group", "random"),
                        ratio=float(ratio), lo=float(lo), hi=float(hi), n=int(a["v"].size),
                        err_ratio=float(np.median(a["e"]) / np.median(b["e"]))))
    return out


def fig_cuts(cuts):
    if not cuts:
        return
    cuts = sorted(cuts, key=lambda c: c["ratio"])
    fig, ax = plt.subplots(figsize=(7.5, 0.45 * len(cuts) + 1.5))
    y = np.arange(len(cuts))
    for i, c in enumerate(cuts):
        col = "#147846" if c["hi"] < 1 else ("#b03a2e" if c["lo"] > 1 else "0.45")
        ax.errorbar(c["ratio"], i, xerr=[[c["ratio"] - c["lo"]], [c["hi"] - c["ratio"]]],
                    fmt=GROUP[c["group"]][1], color=col, ms=6, capsize=2)
    ax.set_yticks(y)
    ax.set_yticklabels([c["name"] for c in cuts], fontsize=8)
    ax.axvline(1.0, color="0.35", lw=0.9)
    g = np.exp(np.mean(np.log([c["ratio"] for c in cuts])))
    ax.axvline(g, color="0.6", ls="--", lw=0.8)
    ax.set_xlabel("rms with the residual cuts / rms without, same exposures")
    ax.set_title("Residual cuts, measured directly: green helps, red hurts\n"
                 "(geometric mean %.3f, dashed; error bars x%.2f in median with the cuts)"
                 % (g, np.median([c["err_ratio"] for c in cuts])), fontsize=9.5)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    save(fig, "f4_cuts")


def fig_core_tails(t):
    fig, ax = plt.subplots(figsize=(6.5, 5))
    pts = [r for r in t if "nominal" in r]
    x = np.array([r["nominal"]["removed_robust"] for r in pts])
    y = np.array([r["nominal"]["removed_rms"] for r in pts])
    for grp, (col, mk, lab) in GROUP.items():
        sel = np.array([r.get("group", "random") == grp for r in pts])
        if sel.any():
            ax.plot(x[sel], y[sel], mk, color=col, ms=6, label=lab)
    label_points(ax, x, y, [name(r) for r in pts])
    ax.axhline(0, color="0.35", lw=0.8)
    ax.axvline(0, color="0.35", lw=0.8)
    ax.set_xscale("symlog", linthresh=2)
    ax.set_yscale("symlog", linthresh=2)
    ax.set_xlabel("removed from the robust sigma (core), m/s")
    ax.set_ylabel("removed from the rms (tails included), m/s")
    ax.set_title("Upper left: tails trimmed, core widened", fontsize=10)
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    save(fig, "f5_core_tails")


def fig_compilation(t):
    """Les series par nuit de chaque cible: livree, nominale, sans coupures."""
    pts = [r for r in t if "nominal" in r]
    n = len(pts)
    cols = 3
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(14, 2.6 * rows), squeeze=False)
    for ax, r in zip(axes.flat, pts):
        T = r["object"] + suffix(r)
        series = [("delivered", T, T, "black"),
                  ("1-3", T + "_PCA2D_1-3", None, "#e0564a"),
                  ("1-3 no cuts", T + "_PCA2D_1-3_n17_cutsoff", None, "#26a1ff")]
        # LES MEMES POSES POUR LES TROIS SERIES: la livree couvre toutes les
        # poses, les corrigees seulement celles du fit; sans cette intersection
        # le rms par nuit du livre de Barnard sortait a 3.11 au lieu de 2.57
        loaded = [(lab, A.load(sci, tmpl), col) for lab, sci, tmpl, col in series]
        loaded = [x for x in loaded if x[1] is not None]
        common = loaded[0][1]
        for _, run, _ in loaded[1:]:
            common, _ = A.on_common(common, run)
        step = 4 * max(r["nominal"]["ref"]["robust"], 1.0)
        for k, (lab, run, col) in enumerate(loaded):
            run, _ = A.on_common(run, common)
            nt = A.nightly_means(run["t"], run["v"], run["e"])
            ax.errorbar(nt[0], nt[1] - np.median(nt[1]) - k * step, yerr=nt[2], fmt=".",
                        ms=2.5, lw=0.5, color=col,
                        label="%s %.2f" % (lab, velocity_stats(run["t"], run["v"], run["e"])["nightly_rms"]))
        ax.set_title("%s (H %.2f, SNR %.0f)" % (name(r), r["H"], r["snr"]), fontsize=8.5)
        ax.tick_params(labelsize=6.5)
        ax.legend(fontsize=5.8, loc="upper right", title="nightly rms, m/s", title_fontsize=5.8)
        ax.grid(alpha=0.2)
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel("rjd", fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel("nightly RV - median, offset (m/s)", fontsize=7.5)
    fig.suptitle("Every target: delivered (black), precleaned 1-3 (red), precleaned without the"
                 " residual cuts (blue), nightly means", fontsize=11)
    fig.tight_layout()
    save(fig, "f6_compilation")


def main():
    t = json.load(open(os.path.join(W, "weekend.json")))
    fig_predictor(t)
    fig_removed(t)
    fig_errors(t)
    cuts = paired_cuts(t)
    json.dump(cuts, open(os.path.join(W, "cuts_paired.json"), "w"), indent=1)
    fig_cuts(cuts)
    fig_core_tails(t)
    fig_compilation(t)


if __name__ == "__main__":
    main()
