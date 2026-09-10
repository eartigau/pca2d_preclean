#!/usr/bin/env python
"""LBL velocities of runs that differ only in their component counts, in one PDF.

    python -m pca2d.lblscan --object TOI2120 \
        --tags 1-3v 2-3v 3-3v 2-2v 2-4v 2-5v 2-6v 2-7v \
        --suffix _PCA2D_{tag}_P2 --star-default 2 --observer-default 3 \
        --out outputs/TOI2120/lbl_scan.pdf

A scan changes one count at a time. The runs at the default observer count are
the star scan, the runs at the default star count the observer scan, and the
run at both defaults is in each. Every run is compared on the exposures all of
them have, and against the delivered spectra on those same exposures.

    page 1  the numbers: rms, robust sigma, nightly rms and median error of
            each run, and the two things that say whether its fit can be
            trusted: how far its velocity term moved the star (rms of
            vrad_fit) and how closely a star coefficient follows the
            barycentric velocity
    page 2  those numbers against the component count, one scan per row
    page 3  the star scan, one panel of velocities per run
    page 4  the observer scan
    page 5  every star coefficient against the observing conditions
    page 6  LBL's STRPCA amplitudes against the fit's own coefficients
"""

from __future__ import annotations

import argparse
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import spearmanr

from .logger import log

#: delivered and corrected: slots 1 and 2 of the validated colour-blind safe
#: categorical palette the report's flux traces use (figures.sequence)
BEFORE, AFTER = "#2a78d6", "#eb6834"
#: a run's tag: star components, observer components, v for the velocity term
TAG = re.compile(r"^(\d+)-(\d+)(v?)$")
#: the numbers every run is summarised by, in the order they are shown
METRICS = (("rms", "rms"), ("robust", "robust sigma"),
           ("nightly_rms", "nightly rms"), ("median_error", "median error"))


# ------------------------------------------------------------- the velocities
def rdb_rows(path):
    """(rjd, vrad, svrad, table) of an LBL rdb, the finite rows only."""
    t = Table.read(path, format="ascii.rdb")
    rjd = np.asarray(t["rjd"], float)
    v = np.asarray(t["vrad"], float)
    e = np.asarray(t["svrad"], float)
    ok = np.isfinite(rjd) & np.isfinite(v) & np.isfinite(e) & (e > 0)
    return rjd[ok], v[ok], e[ok], t[ok]


def read_rdb(path):
    """(rjd, vrad, svrad) of an LBL rdb, the finite rows only."""
    return rdb_rows(path)[:3]


def nightly(rjd, v, e):
    """Weighted nightly means as (rjd, velocity, error), a night being a day of
    rjd + 0.5."""
    night = np.floor(rjd + 0.5).astype(int)
    out = []
    for n in np.unique(night):
        s = night == n
        w = 1 / e[s] ** 2
        out.append((rjd[s].mean(), np.sum(w * v[s]) / w.sum(), 1 / np.sqrt(w.sum())))
    return np.array(out).T


def velocity_stats(t, v, e):
    """rms, robust sigma (1.4826 MAD), median error and the rms of nightly
    weighted means, all about the median."""
    d = v - np.median(v)
    nt = nightly(t, v, e)
    return {"n": d.size, "rms": float(np.std(d)),
            "robust": float(1.4826 * np.median(np.abs(d - np.median(d)))),
            "median_error": float(np.median(e)),
            "nightly_rms": float(np.std(nt[1] - np.median(nt[1]))),
            "nights": int(nt.shape[1])}


# -------------------------------------------------------------------- the runs
def parse_tag(tag):
    """(star components, observer components) of a run's tag, 2-3v -> (2, 3)."""
    m = TAG.match(str(tag))
    if not m:
        raise SystemExit("%r is not a run tag like 2-3v" % tag)
    return int(m.group(1)), int(m.group(2))


def scans(tags, star_default=None, observer_default=None):
    """(star scan, observer scan): the tags at the default observer count, by
    star count, and the tags at the default star count, by observer count.
    A default not given is the count most of the tags share."""
    counts = [parse_tag(t) for t in tags]
    if observer_default is None:
        observer_default = max({c[1] for c in counts}, key=[c[1] for c in counts].count)
    if star_default is None:
        star_default = max({c[0] for c in counts}, key=[c[0] for c in counts].count)
    star = sorted((t for t, c in zip(tags, counts) if c[1] == observer_default),
                  key=lambda t: parse_tag(t)[0])
    observer = sorted((t for t, c in zip(tags, counts) if c[0] == star_default),
                      key=lambda t: parse_tag(t)[1])
    return star, observer


def load_run(obj, tag, suffix, outputs, lbl_dir):
    """What one run left: its LBL velocities, its fit and its coefficients.

    None when LBL has not measured it yet.
    """
    name = obj + suffix.format(tag=tag)
    rdb = os.path.join(lbl_dir, "lblrdb", "lbl_%s_%s.rdb" % (name, name))
    if not os.path.exists(rdb):
        log("no LBL velocities for %s yet (%s): left out" % (tag, rdb), "warn")
        return None
    t, v, e, table = rdb_rows(rdb)
    n_star, n_earth = parse_tag(tag)
    run = {"tag": tag, "name": name, "n_star": n_star, "n_earth": n_earth,
           "t": t, "v": v, "e": e, "table": table}
    folder = os.path.join(outputs, obj, tag)
    if os.path.exists(os.path.join(folder, "fit.npz")):
        run["fit"] = np.load(os.path.join(folder, "fit.npz"), allow_pickle=True)
    components = os.path.join(folder, "twoframe_components.fits")
    if os.path.exists(components):
        c = fits.getdata(components, "COEFFS")
        if "vrad_fit" in c.columns.names:
            keep = ~np.asarray(c["rejected"], bool) & np.isfinite(np.asarray(c["vrad_fit"], float))
            run["vrad_fit"] = np.asarray(c["vrad_fit"], float)[keep]
            run["berv"] = np.asarray(c["berv"], float)[keep]
    return run


def on_common(runs, delivered):
    """Every run, and the delivered spectra, cut to the exposures all have."""
    key = lambda t: np.round(t, 6)
    common = key(delivered["t"])
    for run in runs:
        common = np.intersect1d(common, key(run["t"]))
    for run in [delivered] + runs:
        keep = np.isin(key(run["t"]), common)
        for name in ("t", "v", "e"):
            run[name] = run[name][keep]
        run["table"] = run["table"][keep]
        run["stats"] = velocity_stats(run["t"], run["v"], run["e"])
    return len(common)


def star_conditions(run):
    """(labels, names, rho): Spearman rho of every star coefficient of a run's
    fit against the observing conditions the fit recorded."""
    from .plotting import rank_correlations
    fit = run.get("fit")
    if fit is None:
        return None
    keep = ~np.asarray(fit["rejected"], bool)
    a = np.asarray(fit["a"], float)[keep]
    comps = [("a%d" % (k + 1), a[:, k]) for k in range(a.shape[1])]
    rho = np.asarray(rank_correlations(comps, np.asarray(fit["anc_values"])[:, keep]))
    return [str(x) for x in fit["anc_labels"]], [n for n, _ in comps], rho


def trust(run):
    """(rms of vrad_fit in m/s, its rho with BERV, the star coefficient that
    follows BERV most closely and its rho)."""
    out = [np.nan, np.nan, "", np.nan]
    if "vrad_fit" in run:
        out[0] = float(np.std(run["vrad_fit"]))
        out[1] = float(spearmanr(run["vrad_fit"], run["berv"])[0])
    sc = star_conditions(run)
    if sc is not None and "BERV" in sc[0]:
        col = sc[2][:, sc[0].index("BERV")]
        k = int(np.argmax(np.abs(col)))
        out[2], out[3] = sc[1][k], float(col[k])
    return out


def strpca_pairs(run):
    """[(key, the fit's a_k less its mean, LBL's STRPCAk, its error)] of a run,
    exposure by exposure, matched by file name."""
    fit, table = run.get("fit"), run["table"]
    keys = sorted((c for c in table.colnames if re.fullmatch(r"STRPCA\d+", c)),
                  key=lambda c: int(c[6:]))
    if fit is None or not keys:
        return []
    rejected = np.asarray(fit["rejected"], bool)
    first = {}
    for r, name in enumerate(str(x) for x in fit["filename"]):
        first.setdefault(name.split("t.fits")[0], r)
    stems = [os.path.basename(str(f)).split("t_")[0] for f in table["FILENAME"]]
    hit = np.array([s in first and not rejected[first[s]] for s in stems])
    rows = np.array([first[s] for s, h in zip(stems, hit) if h], dtype=int)
    a = np.asarray(fit["a"], float)
    out = []
    for key in keys:
        k = int(key[6:]) - 1
        x = a[rows, k] - a[~rejected, k].mean()
        y = np.asarray(table[key], float)[hit]
        e = np.asarray(table["s" + key], float)[hit]
        ok = np.isfinite(x) & np.isfinite(y)
        out.append((key, x[ok], y[ok], e[ok]))
    return out


# ------------------------------------------------------------------- the pages
def summary_page(delivered, runs, n_common, title):
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.04, 0.95, title, fontsize=13, weight="bold")
    fig.text(0.04, 0.915, "%d exposures that every run has; velocities in m/s, each"
             " about its median; robust sigma = 1.4826 MAD" % n_common, fontsize=9)
    head = ("%-10s %5s %5s %8s %8s %9s %9s   %12s %9s   %s"
            % ("run", "star", "obs", "rms", "robust", "nightly", "median",
               "vrad_fit rms", "vs BERV", "star coef. most like BERV"))
    lines = [head, "-" * len(head)]
    s = delivered["stats"]
    lines.append("%-10s %5s %5s %8.2f %8.2f %9.2f %9.2f" % (
        "delivered", "", "", s["rms"], s["robust"], s["nightly_rms"], s["median_error"]))
    for run in runs:
        s = run["stats"]
        vr, vb, which, rho = trust(run)
        lines.append("%-10s %5d %5d %8.2f %8.2f %9.2f %9.2f   %12.1f %+9.2f   %s"
                     % (run["tag"], run["n_star"], run["n_earth"], s["rms"], s["robust"],
                        s["nightly_rms"], s["median_error"], vr, vb,
                        "%s %+.2f" % (which, rho) if which else ""))
    fig.text(0.04, 0.86, "\n".join(lines), family="monospace", fontsize=8.6, va="top")
    fig.text(0.04, 0.18,
             "vrad_fit is the star shift each fit's velocity term took out of the observer"
             " block's reach, written to every corrected file as PCASTR_V. A fit whose\n"
             "velocity term moves the star by far more than LBL measures, in step with the"
             " barycentric velocity, or whose star coefficient follows the barycentric\n"
             "velocity, the seeing or the sun, has let the star and the sky trade places, and"
             " its correction leaves the difference in the flux as velocity.",
             fontsize=8.5, va="top")
    return fig


def metrics_page(delivered, star, observer):
    fig, axes = plt.subplots(2, len(METRICS), figsize=(11, 6.4), squeeze=False)
    for row, (runs, count, label) in enumerate(((star, "n_star", "star components"),
                                                (observer, "n_earth", "observer components"))):
        for col, (key, name) in enumerate(METRICS):
            ax = axes[row, col]
            if runs:
                x = [r[count] for r in runs]
                y = [r["stats"][key] for r in runs]
                ax.plot(x, y, "o-", color=AFTER, lw=1.2, ms=5, label="corrected")
                ax.set_xticks(x)
            ax.axhline(delivered["stats"][key], color=BEFORE, ls="--", lw=1.0,
                       label="delivered")
            ax.set_xlabel(label, fontsize=8.5)
            ax.set_title(name + " (m/s)", fontsize=9)
            ax.tick_params(labelsize=8)
            ax.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=7.5, frameon=False)
    fig.suptitle("the numbers against the component count; top: the star scan,"
                 " bottom: the observer scan", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def series_page(delivered, runs, which):
    panels = [delivered] + runs
    fig, axes = plt.subplots(len(panels), 1, figsize=(11, 1.2 + 1.75 * len(panels)),
                             sharex=True, sharey=True, squeeze=False)
    axes = axes[:, 0]
    for k, (ax, run) in enumerate(zip(axes, panels)):
        colour = BEFORE if k == 0 else AFTER
        d = run["v"] - np.median(run["v"])
        ax.errorbar(run["t"], d, yerr=run["e"], fmt="o", ms=2.2, lw=0.5,
                    color=colour, ecolor=colour, alpha=0.55, capsize=0)
        nt = nightly(run["t"], run["v"], run["e"])
        ax.plot(nt[0], nt[1] - np.median(run["v"]), "o", ms=3.2, mfc="white",
                mec="0.15", mew=0.8, label="nightly mean")
        ax.axhline(0, color="0.6", lw=0.6)
        s = run["stats"]
        ax.set_title("%s: rms %.2f, robust sigma %.2f, nightly rms %.2f, median error"
                     " %.2f m/s" % (run.get("tag", "delivered spectra"), s["rms"], s["robust"],
                                    s["nightly_rms"], s["median_error"]), fontsize=8.5)
        ax.set_ylabel("m/s", fontsize=8)
        ax.tick_params(labelsize=7.5)
        ax.grid(alpha=0.15)
    lim = np.percentile(np.abs(np.concatenate([r["v"] - np.median(r["v"]) for r in panels])), 99.5)
    axes[0].set_ylim(-1.15 * lim, 1.15 * lim)
    axes[0].legend(fontsize=7.5, frameon=False, loc="upper right")
    axes[-1].set_xlabel("rjd (BJD - 2400000)", fontsize=9)
    fig.suptitle("the %s scan: LBL velocities of each run, and of the delivered spectra"
                 " on the same exposures" % which, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig


def conditions_page(runs):
    rows, names, labels = [], [], None
    for run in runs:
        sc = star_conditions(run)
        if sc is None:
            continue
        labels = sc[0]
        for name, r in zip(sc[1], sc[2]):
            rows.append(r)
            names.append("%s  %s" % (run["tag"], name))
    if not rows:
        return None
    rho = np.array(rows)
    fig, ax = plt.subplots(figsize=(11, 1.6 + 0.36 * len(rows)))
    im = ax.imshow(rho, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    for i in range(rho.shape[0]):
        for j in range(rho.shape[1]):
            if abs(rho[i, j]) >= 0.3:
                ax.text(j, i, "%+.2f" % rho[i, j], ha="center", va="center", fontsize=7.5,
                        color="white" if abs(rho[i, j]) > 0.6 else "black")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8, rotation=30, ha="right")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.01, label="Spearman rho")
    ax.set_title("every star coefficient against the observing conditions, |rho| printed"
                 " above 0.30: a star coefficient has no business following these",
                 fontsize=9.5)
    fig.tight_layout()
    return fig


def strpca_page(runs, aliased):
    pairs = [(run["tag"], p) for run in runs for p in strpca_pairs(run)]
    if not pairs:
        return None
    ncol = min(4, len(pairs))
    nrow = int(np.ceil(len(pairs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(11, 0.8 + 2.7 * nrow), squeeze=False)
    for ax in axes.ravel()[len(pairs):]:
        ax.axis("off")
    for ax, (tag, (key, x, y, e)) in zip(axes.ravel(), pairs):
        r = float(np.corrcoef(x, y)[0, 1]) if x.size > 2 else np.nan
        slope = float(np.sum(x * y) / np.sum(x * x)) if x.size else np.nan
        ax.errorbar(x, y, yerr=e, fmt="o", ms=2, lw=0.4, color=AFTER, ecolor=AFTER,
                    alpha=0.5, capsize=0)
        lim = 1.1 * np.percentile(np.abs(np.r_[x, y]), 99.5) if x.size else 1
        ax.plot([-lim, lim], [-lim, lim], color="0.45", lw=0.8)
        ax.set(xlim=(-lim, lim), ylim=(-lim, lim))
        late = aliased and int(key[6:]) >= 3
        ax.set_title("%s %s: r %.3f, slope %.2f%s" % (tag, key, r, slope,
                                                       "\nLBL divides this residual twice"
                                                       if late else ""), fontsize=8)
        ax.set_xlabel("the fit's a%s" % key[6:], fontsize=7.5)
        ax.set_ylabel("LBL's " + key, fontsize=7.5)
        ax.tick_params(labelsize=7)
    fig.suptitle("LBL's STRPCA amplitudes against the fit's own coefficients, exposure by"
                 " exposure", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--object", required=True)
    p.add_argument("--tags", nargs="+", required=True, help="the runs, e.g. 1-3v 2-3v")
    p.add_argument("--suffix", default="_PCA2D_{tag}",
                   help="the corrected objects' LBL names, as lbl.suffix spells them")
    p.add_argument("--star-default", type=int, default=None)
    p.add_argument("--observer-default", type=int, default=None)
    p.add_argument("--outputs", default="outputs")
    p.add_argument("--lbl-dir", default="lbl")
    p.add_argument("--title", default=None)
    p.add_argument("--out", required=True, help="the PDF to write")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    tags = list(dict.fromkeys(args.tags))
    star_tags, observer_tags = scans(tags, args.star_default, args.observer_default)
    delivered_rdb = os.path.join(args.lbl_dir, "lblrdb", "lbl_%s_%s.rdb"
                                 % (args.object, args.object))
    if not os.path.exists(delivered_rdb):
        raise SystemExit("no LBL velocities for the delivered spectra: %s" % delivered_rdb)
    t, v, e, table = rdb_rows(delivered_rdb)
    delivered = {"t": t, "v": v, "e": e, "table": table}
    loaded = {tag: load_run(args.object, tag, args.suffix, args.outputs, args.lbl_dir)
              for tag in tags}
    runs = [loaded[tag] for tag in tags if loaded[tag] is not None]
    if not runs:
        raise SystemExit("none of the runs has LBL velocities yet")
    n_common = on_common(runs, delivered)
    star = [loaded[t] for t in star_tags if loaded[t] is not None]
    observer = [loaded[t] for t in observer_tags if loaded[t] is not None]
    from .lbltemplate import resproj_divides_in_place
    title = args.title or "%s: LBL velocities across component counts" % args.object
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with PdfPages(args.out) as pdf:
        for fig in (summary_page(delivered, runs, n_common, title),
                    metrics_page(delivered, star, observer),
                    series_page(delivered, star, "star") if star else None,
                    series_page(delivered, observer, "observer") if observer else None,
                    conditions_page(runs),
                    strpca_page(runs, resproj_divides_in_place())):
            if fig is not None:
                pdf.savefig(fig)
                plt.close(fig)
    log("wrote %s: %d runs on %d common exposures" % (args.out, len(runs), n_common), "info")
    for run in [delivered] + runs:
        s = run["stats"]
        log("%-10s rms %6.2f  robust sigma %6.2f  nightly rms %6.2f  median error %5.2f"
            % (run.get("tag", "delivered"), s["rms"], s["robust"], s["nightly_rms"],
               s["median_error"]), "value")


if __name__ == "__main__":
    main()
