#!/usr/bin/env python
"""A run's report as a document: LaTeX, a table of contents, sections, numbers.

    python -m pca2d.texreport --run <run folder>       # write it again, now
    python -m pca2d.texreport --run <run folder> --lbl-dir <LBL tree>

The report a run left was the figures it drew, bound one after another with a
page of monospace text in front: a compilation of matplotlib outputs, which is
what it looked like (2026-09-16). This writes a report instead, from
templates/report.tex:

    abstract          what was run and what it did to the velocities, in words
    lists             the contents, and every figure by a short name
    1 Summary         the gains and losses of every star, one table each,
                      every row saying which way it went
    2 Velocities      per star: the time series (points, no lines), the same
                      points against BERV, d2v as an activity indicator, what
                      the correction changed, and the periodograms
    3 The run         every setting, the command, the code, the fit's numbers
    4 The correction  the figures of the figures stage, each with a caption
    5 Not drawn       what did not build, and why
    A Parameters      the resolved configuration, verbatim

It is written twice in a run's life: at the figures stage, when there are no
velocities yet, and after LBL, when there are. Both times from what is on disk:
the figures the figures stage keeps in <run>/report/figures, listed in
<run>/report/manifest.json, and LBL's rdb files. A run from before the report
existed has neither folder nor manifest; its bound PDF is cut along its own
bookmarks instead, so `--run` works on it too.

Nothing here is fatal to a run. Without pdflatex, or when it fails, the bound
PDF the figures stage wrote stays the report and the run says why.
"""

from __future__ import annotations

import argparse
import datetime
import glob
import json
import os
import re
import shutil
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

from . import bervbias
from .lblscan import AFTER, BEFORE, amplitude_at, nightly, rdb_rows, velocity_stats
from .logger import log

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE = os.path.join(HERE, "templates", "report.tex")
MANIFEST = "manifest.json"
#: the bound PDF's own text sections. They are written again here, as text,
#: so they are not taken as figures when an old bundle is cut up
TEXT_SECTIONS = ("Front matter", "Parameters", "What did not build",
                 "The velocities LBL measured")
#: where a TeX installation puts pdflatex when the window's PATH does not
#: know it: a program started from the Dock does not read the shell's profile
TEX_HOMES = ("/Library/TeX/texbin", "/usr/texbin", "/opt/homebrew/bin",
             "/usr/local/bin", "/usr/bin")
#: ink, for everything that is not data (dataviz: text never wears a series
#: colour); the two series colours are lblscan's, validated as a pair
INK, MUTED, GRID = "#1f1f1f", "#5f6368", "#e3e3e3"
#: a figure is drawn at the width it is printed at, so its text prints at the
#: size it was written at: A4 less two 2.2 cm margins
WIDTH = 6.6

#: what each figure of the figures stage shows, said under it. Keyed by the
#: bundle's own section titles
CAPTIONS = {
    "Every quantity of the correlations, against time":
        "Every quantity the correlations are measured against, over the"
        " campaign: the conditions the observer components are compared with.",
    "The sequence, step by step":
        "The whole chain in one window: the spectra as read, the fitted"
        " observer model, what was divided out, and what is left.",
    "Variance":
        "The weighted variance each component removes, star and observer"
        " blocks apart.",
    "Coefficients against time":
        "Every component's coefficient per exposure, over the campaign.",
    "Correlations":
        "Every coefficient against every measured quantity, the strongest"
        " correlations first.",
    "Coefficient periodogram":
        "Periodogram of every coefficient. A component varying at a known"
        " planet's period would subtract that planet: it is marked when the"
        " configuration gives the period.",
    "Weight spectrum":
        "The weights the fit gave each sample, which is where the model is"
        " trusted and where it is not.",
    "OH residual":
        "The airglow model against what the fit left, where SPIRou records it.",
}


# ================================================================= LaTeX ===
_TEX = {"\\": r"\textbackslash{}", "{": r"\{", "}": r"\}", "$": r"\$",
        "&": r"\&", "#": r"\#", "^": r"\^{}", "_": r"\_", "%": r"\%",
        "~": r"\textasciitilde{}", "<": r"\textless{}", ">": r"\textgreater{}",
        "|": r"\textbar{}"}
#: characters inputenc does not take in text mode, and what they become. The
#: dashes become hyphens on purpose: no em dash in anything written here
_UNICODE = {"\u00b2": r"\textsuperscript{2}", "\u00b3": r"\textsuperscript{3}",
            "\u03c3": r"$\sigma$", "\u00d7": r"$\times$", "\u00b1": r"$\pm$",
            "\u2265": r"$\geq$", "\u2264": r"$\leq$", "\u2013": "-",
            "\u2014": "-", "\u2212": "-", "\u2026": r"\ldots{}",
            "\u00b0": r"$^\circ$", "\u2192": r"$\rightarrow$",
            "\u0394": r"$\Delta$", "\u00b7": r"$\cdot$", "\u2019": "'",
            "\u2018": "`", "\u201c": "``", "\u201d": "''", "\u00a0": "~"}


def tex(text) -> str:
    """`text` as LaTeX prints it, whatever it holds.

    Object names carry underscores (SMETHELLS_20), run names plus signs and
    paths every character there is, and one unescaped one stops the document.
    A character nothing here knows becomes a question mark rather than an
    error: a report with one odd glyph is still a report.
    """
    out = []
    for ch in str(text):
        if ch in _TEX:
            out.append(_TEX[ch])
        elif ch in _UNICODE:
            out.append(_UNICODE[ch])
        elif ord(ch) < 256:
            out.append(ch)
        else:
            out.append("?")
    return "".join(out)


def tex_break(text) -> str:
    """tex(), with a place to break after every _ / + - . in it.

    A run's suffix (_PCA2D_{tag}_GL406_3b7955) is one word to LaTeX, and one
    word wider than its column runs over the next one.
    """
    out = tex(text)
    for mark in (r"\_", "/", "+", "-", "."):
        out = out.replace(mark, mark + r"\allowbreak{}")
    return out


def find_pdflatex():
    """The pdflatex this machine has, or None."""
    found = shutil.which("pdflatex")
    if found:
        return found
    for home in TEX_HOMES:
        path = os.path.join(home, "pdflatex")
        if os.access(path, os.X_OK):
            return path
    for path in sorted(glob.glob("/usr/local/texlive/*/bin/*/pdflatex"),
                       reverse=True):
        if os.access(path, os.X_OK):
            return path
    return None


def slug(text) -> str:
    """A file name LaTeX reads without complaint."""
    return re.sub(r"[^A-Za-z0-9]+", "-", str(text)).strip("-").lower() or "x"


def number(value, digits=2, unit=""):
    """A number as the tables print it, or a dash-free 'n/a'."""
    if value is None or not np.isfinite(value):
        return "n/a"
    text = "%.*f" % (digits, value)
    return text + (" " + unit if unit else "")


def wrap_lines(text, width=96):
    """Hard-wrapped, since verbatim never breaks a line itself."""
    out = []
    for line in str(text).splitlines():
        while len(line) > width:
            out.append(line[:width])
            line = "    " + line[width:]
        out.append(line)
    return "\n".join(out) + "\n"


# ============================================================ the figures ===
def pdf_pages(path) -> int:
    from pypdf import PdfReader
    try:
        return len(PdfReader(path).pages)
    except Exception:                                          # noqa: BLE001
        return 0


def read_manifest(folder):
    path = os.path.join(folder, MANIFEST)
    if not os.path.exists(path):
        return None
    with open(path) as handle:
        return json.load(handle)


def write_manifest(folder, manifest):
    os.makedirs(folder, exist_ok=True)
    tmp = os.path.join(folder, MANIFEST + ".tmp")
    with open(tmp, "w") as handle:
        json.dump(manifest, handle, indent=1)
    os.replace(tmp, os.path.join(folder, MANIFEST))


def keep_figures(order, folder):
    """Copy the figures stage's figures where the report can find them again.

    `order` is the bundle's [(section title, path)]. Returns the manifest's
    figure list, in that order, text sections left out.
    """
    figures = os.path.join(folder, "figures")
    os.makedirs(figures, exist_ok=True)
    kept = []
    for i, (title, path) in enumerate(order):
        if title in TEXT_SECTIONS or not os.path.exists(path):
            continue
        name = "%02d-%s.pdf" % (i, slug(title))
        shutil.copyfile(path, os.path.join(figures, name))
        kept.append({"title": title, "file": "figures/" + name,
                     "pages": pdf_pages(path)})
    return kept


def split_bundle(bundle, folder):
    """Cut a bound PDF from before the report along its own bookmarks.

    Returns the manifest's figure list. The bookmarks were written with page
    indices, so each section runs from its own index to the next one's.
    """
    from pypdf import PdfReader, PdfWriter

    reader = PdfReader(bundle)
    refs = {page.indirect_reference.idnum: i
            for i, page in enumerate(reader.pages)}
    marks = []
    for item in reader.outline:
        if isinstance(item, list):
            continue
        page = item.get("/Page")
        index = (int(page) if isinstance(page, (int, float))
                 else refs.get(getattr(page, "idnum", None)))
        if index is None:
            try:
                index = reader.get_destination_page_number(item)
            except Exception:                                   # noqa: BLE001
                index = None
        if index is not None:
            marks.append((int(index), str(item.title)))
    marks.sort()
    figures = os.path.join(folder, "figures")
    os.makedirs(figures, exist_ok=True)
    kept = []
    for n, (start, title) in enumerate(marks):
        end = marks[n + 1][0] if n + 1 < len(marks) else len(reader.pages)
        if title in TEXT_SECTIONS or end <= start:
            continue
        writer = PdfWriter()
        for page in reader.pages[start:end]:
            writer.add_page(page)
        name = "%02d-%s.pdf" % (n, slug(title))
        with open(os.path.join(figures, name), "wb") as handle:
            writer.write(handle)
        kept.append({"title": title, "file": "figures/" + name,
                     "pages": end - start})
    return kept


# ========================================================= the velocities ===
def rdb_file(tree, name):
    return os.path.join(tree, "lblrdb", "lbl_%s_%s.rdb" % (name, name))


def column(table, *names):
    """A column of an rdb as floats, whatever its case, or None."""
    lower = {c.strip().lower(): c for c in table.colnames}
    for name in names:
        if name.lower() in lower:
            values = np.asarray(table[lower[name.lower()]], float)
            return values if np.isfinite(values).any() else None
    return None


def load(tree, name, label):
    """One LBL object as the figures take it, or None when it has no rdb."""
    path = rdb_file(tree, name)
    if not os.path.exists(path):
        return None
    t, v, e, table = rdb_rows(path)
    if t.size < 3:
        return None
    return {"label": label, "name": name, "path": path, "t": t, "v": v,
            "e": e, "berv": column(table, "BERV"),
            "d2v": column(table, "d2v"), "sd2v": column(table, "sd2v"),
            "mtime": os.path.getmtime(path)}


def common(before, after):
    """Both series cut to the exposures they share, in time order."""
    key = lambda t: np.round(t, 6)                            # noqa: E731
    shared = np.intersect1d(key(before["t"]), key(after["t"]))
    out = []
    for run in (before, after):
        keep = np.isin(key(run["t"]), shared)
        order = np.argsort(run["t"][keep])
        cut = dict(run)
        for name in ("t", "v", "e", "berv", "d2v", "sd2v"):
            if cut.get(name) is not None:
                cut[name] = np.asarray(run[name])[keep][order]
        out.append(cut)
    return out[0], out[1]


def robust_sigma(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 2:
        return np.nan
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def inliers(x, nsig=5.0):
    """Where `x` is within nsig robust sigmas of its median."""
    x = np.asarray(x, float)
    ok = np.isfinite(x)
    if ok.sum() < 3:
        return ok
    sig = robust_sigma(x[ok])
    if not np.isfinite(sig) or sig == 0:
        return ok
    return ok & (np.abs(x - np.median(x[ok])) < nsig * sig)


def pearson(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 4 or np.std(x[ok]) == 0 or np.std(y[ok]) == 0:
        return np.nan
    return float(np.corrcoef(x[ok], y[ok])[0, 1])


def binned(x, y, width=2.0, least=3):
    """(centres, medians, errors of the medians) of y in bins of x."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < least:
        return np.array([]), np.array([]), np.array([])
    edges = np.arange(np.floor(x[ok].min() / width) * width,
                      x[ok].max() + width, width)
    centres, medians, errors = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        inside = ok & (x >= lo) & (x < hi)
        if inside.sum() >= least:
            centres.append(0.5 * (lo + hi))
            medians.append(np.median(y[inside]))
            errors.append(1.2533 * robust_sigma(y[inside])
                          / np.sqrt(inside.sum()))
    return np.array(centres), np.array(medians), np.array(errors)


def periodogram(t, y, e, pmin=1.1, samples=6000):
    """(periods, power, best period, its power, false-alarm probability)."""
    from astropy.timeseries import LombScargle

    t, y, e = (np.asarray(a, float) for a in (t, y, e))
    ok = np.isfinite(t) & np.isfinite(y) & np.isfinite(e) & (e > 0)
    baseline = float(np.ptp(t[ok])) if ok.sum() > 3 else 0.0
    if baseline <= 2 * pmin or ok.sum() < 8:
        return None
    fmin, fmax = 1.0 / baseline, 1.0 / pmin
    freq = np.linspace(fmin, fmax, samples)
    model = LombScargle(t[ok], y[ok], e[ok])
    power = model.power(freq)
    best = int(np.nanargmax(power))
    try:
        fap = float(model.false_alarm_probability(
            power[best], minimum_frequency=fmin, maximum_frequency=fmax))
    except Exception:                                          # noqa: BLE001
        fap = np.nan
    return 1.0 / freq, power, float(1.0 / freq[best]), float(power[best]), fap


def bias_row(fitted):
    """The table's view of one BERV-bias fit."""
    if fitted is None:
        return {}
    peak, lo, hi = fitted["peak"]
    return {"bias": fitted, "bias_peak": float(peak),
            "bias_err": float(0.5 * (hi - lo)),
            "bias_significance": fitted["significance"],
            "bias_detected": fitted["detected"],
            "bias_upper": fitted["upper"],
            "bias_width": float(fitted["sigma"][0]),
            "bias_width_err": float(0.5 * (fitted["sigma"][2]
                                           - fitted["sigma"][1])),
            "bias_amp_sigma_r": fitted["amp_sigma_r"],
            "jitter": float(fitted["jitter"][0])}


def star_numbers(before, after, planets=(), seed=0):
    """Everything the summary table says about one star, both ways."""
    out = {"n": int(before["t"].size),
           "nights": int(velocity_stats(after["t"], after["v"],
                                        after["e"])["nights"])}
    for key, run in (("before", before), ("after", after)):
        s = velocity_stats(run["t"], run["v"], run["e"])
        row = {"rms": s["rms"], "robust": s["robust"],
               "nightly_rms": s["nightly_rms"],
               "median_error": s["median_error"]}
        if run.get("berv") is not None:
            _, medians, _ = binned(run["berv"],
                                   run["v"] - np.median(run["v"]))
            row.update(berv_binned=float(np.std(medians))
                       if medians.size >= 3 else np.nan)
            # the bias's own shape, by MCMC; a straight line has no meaning
            row.update(bias_row(bervbias.fit(run["berv"], run["v"], run["e"],
                                             seed=seed)))
        if run.get("d2v") is not None:
            good = inliers(run["d2v"])
            row.update(d2v_sigma=robust_sigma(run["d2v"]),
                       d2v_error=float(np.nanmedian(run["sd2v"]))
                       if run.get("sd2v") is not None else np.nan,
                       d2v_r=pearson(run["d2v"][good], run["v"][good]))
        found = periodogram(run["t"], run["v"], run["e"])
        if found:
            row.update(peak_period=found[2], peak_power=found[3],
                       peak_fap=found[4])
        if run.get("d2v") is not None:
            good = inliers(run["d2v"])
            err = run.get("sd2v")
            spot = periodogram(run["t"][good], run["d2v"][good],
                               (err if err is not None
                                else np.ones_like(run["d2v"]))[good])
            if spot:
                row.update(d2v_peak_period=spot[2], d2v_peak_fap=spot[4])
        row["planets"] = [amplitude_at(run["t"], run["v"], run["e"], p)
                          for p in planets]
        out[key] = row
    change = after["v"] - before["v"]
    out["change_rms"] = float(np.std(change - np.median(change)))
    rb, ra = out["before"]["rms"], out["after"]["rms"]
    out["removed"] = float(np.sqrt(rb ** 2 - ra ** 2)) if rb > ra else np.nan
    out["planet_periods"] = [float(p) for p in planets]
    out["baseline"] = float(np.ptp(after["t"]))
    return out


def same_period(p1, p2, baseline):
    """Whether two periods are one peak: within a resolution element."""
    if not (p1 and p2 and baseline and np.isfinite(p1) and np.isfinite(p2)):
        return False
    return abs(1.0 / p1 - 1.0 / p2) < 1.0 / baseline


def verdict(before, after, better="lower", tolerance=0.02):
    """'gain', 'loss' or 'same' for one row of the table."""
    if not (np.isfinite(before) and np.isfinite(after)):
        return "same"
    if better == "lower":
        if after < before * (1 - tolerance):
            return "gain"
        if after > before * (1 + tolerance):
            return "loss"
        return "same"
    # closer to zero is better: a correlation, a slope
    if abs(after) < abs(before) - tolerance:
        return "gain"
    if abs(after) > abs(before) + tolerance:
        return "loss"
    return "same"


# --------------------------------------------------------------- drawing --
def _style(ax):
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=7.5, labelcolor=INK)


def _limits(*arrays, pad=1.25, q=99.0):
    values = np.concatenate([np.abs(np.asarray(a, float)) for a in arrays])
    values = values[np.isfinite(values)]
    top = float(np.percentile(values, q)) * pad if values.size else 1.0
    return (-top, top) if top > 0 else (-1.0, 1.0)


def _points(ax, x, y, e, colour, label, size=3.2, alpha=0.55, zorder=2):
    """Exposures as points: a marker and its bar, never a line between."""
    ax.errorbar(x, y, yerr=e, fmt="o", ls="none", ms=size, mfc=colour,
                mec="white", mew=0.4, ecolor=colour, elinewidth=0.5,
                alpha=alpha, capsize=0, label=label, zorder=zorder)


def _save(fig, path):
    fig.savefig(path)
    plt.close(fig)
    return path


def figure_time(before, after, path):
    fig, axes = plt.subplots(2, 1, figsize=(WIDTH, 5.4), sharex=True,
                             sharey=True)
    lim = _limits(before["v"] - np.median(before["v"]),
                  after["v"] - np.median(after["v"]))
    for ax, run, colour in zip(axes, (before, after), (BEFORE, AFTER)):
        middle = np.median(run["v"])
        nt = nightly(run["t"], run["v"], run["e"])
        # one exposure a night makes the nightly means the exposures again,
        # drawn over them: one layer then, not two that hide each other
        if nt.shape[1] < 0.8 * run["t"].size:
            _points(ax, run["t"], run["v"] - middle, run["e"], colour,
                    "exposures")
            label = "nightly means"
        else:
            label = "exposures (about one a night)"
        ax.errorbar(nt[0], nt[1] - middle, yerr=nt[2], fmt="o", ls="none",
                    ms=5.0, mfc=colour, mec=INK, mew=0.5, ecolor=INK,
                    elinewidth=0.6, capsize=0, label=label, zorder=3)
        s = velocity_stats(run["t"], run["v"], run["e"])
        ax.set_title("%s: rms %.2f m/s, nightly rms %.2f m/s, median error"
                     " %.2f m/s" % (run["label"], s["rms"], s["nightly_rms"],
                                    s["median_error"]),
                     fontsize=8.5, color=INK, loc="left")
        ax.axhline(0, color=MUTED, lw=0.6)
        ax.set_ylabel("velocity - median (m/s)", fontsize=8, color=INK)
        ax.set_ylim(*lim)
        ax.legend(fontsize=7, frameon=False, loc="upper right", ncol=2)
        _style(ax)
    axes[-1].set_xlabel("RJD (BJD - 2400000)", fontsize=8, color=INK)
    fig.tight_layout()
    return _save(fig, path)


def _bias_band(ax, fitted, grid, colour, label=None, offset=0.0):
    """The fitted bias: the band its 16th-84th percentiles span over `grid`,
    and its median when the bias is detected.

    Not the median of one that is not: the median of curves of every width
    the prior allows is a spike nobody fitted, and it reads as a feature.
    """
    lo, mid, hi = bervbias.envelope(fitted, grid)
    ax.fill_between(grid, lo + offset, hi + offset, color=colour,
                    alpha=0.22 if fitted["detected"] else 0.14, lw=0,
                    zorder=3, label=None if fitted["detected"] else label)
    if fitted["detected"]:
        ax.plot(grid, mid + offset, color=colour, lw=1.6, zorder=4,
                label=label)


def figure_berv(before, after, path):
    """The velocities against BERV with the fitted bias and its 1 sigma
    envelope, each series on its own, then both envelopes on one axis."""
    fb = (before.get("fits") or {})
    if before.get("berv") is None or after.get("berv") is None \
            or fb.get("before") is None or fb.get("after") is None:
        return None
    fits = {"before": fb["before"], "after": fb["after"]}
    fig = plt.figure(figsize=(WIDTH, 6.4))
    grid_spec = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.9])
    top = [fig.add_subplot(grid_spec[0, 0])]
    top.append(fig.add_subplot(grid_spec[0, 1], sharex=top[0], sharey=top[0]))
    both = fig.add_subplot(grid_spec[1, :], sharex=top[0])
    span = np.concatenate([before["berv"], after["berv"]])
    grid = np.linspace(np.nanmin(span), np.nanmax(span), 300)
    residuals = []
    for ax, run, key, colour in zip(top, (before, after), ("before", "after"),
                                    (BEFORE, AFTER)):
        fitted = fits[key]
        y = run["v"] - fitted["c"][0]
        residuals.append(y)
        _points(ax, run["berv"], y, run["e"], colour, "exposures", alpha=0.4)
        centres, medians, errors = binned(run["berv"], y)
        if centres.size:
            ax.errorbar(centres, medians, yerr=errors, fmt="s", ls="none",
                        ms=4.5, mfc="white", mec=INK, mew=0.8, ecolor=INK,
                        elinewidth=0.7, capsize=0, zorder=5,
                        label="2 km/s bins, median")
        _bias_band(ax, fitted, grid, INK, label="fitted bias, 1$\\sigma$")
        ax.set_title("%s\n%s" % (run["label"], bervbias.summary(fitted)),
                     fontsize=7, color=INK, loc="left")
        ax.axhline(0, color=MUTED, lw=0.6)
        ax.set_xlabel("BERV (km/s)", fontsize=8, color=INK)
        _style(ax)
    lim = _limits(*residuals)
    top[0].set_ylim(*lim)
    top[0].set_ylabel("velocity - fitted offset (m/s)", fontsize=8, color=INK)
    top[0].legend(fontsize=6, frameon=False, loc="upper left")
    for run, key, colour in ((before, "before", BEFORE),
                             (after, "after", AFTER)):
        fitted = fits[key]
        centres, medians, errors = binned(run["berv"],
                                          run["v"] - fitted["c"][0])
        if centres.size:
            both.errorbar(centres, medians, yerr=errors, fmt="s", ls="none",
                          ms=4.5, mfc="white", mec=colour, mew=0.9,
                          ecolor=colour, elinewidth=0.7, capsize=0, zorder=5)
        _bias_band(both, fitted, grid, colour,
                   label="%s: %s" % (run["label"], bervbias.summary(fitted)))
    both.axhline(0, color=MUTED, lw=0.6)
    both.set_xlabel("BERV (km/s)", fontsize=8, color=INK)
    both.set_ylabel("fitted bias (m/s)", fontsize=8, color=INK)
    both.set_title("both fits, 1$\\sigma$ envelopes, over the binned medians",
                   fontsize=7.5, color=INK, loc="left")
    both.legend(fontsize=6.5, frameon=False, loc="lower left",
                bbox_to_anchor=(0.0, 1.06), ncol=1, borderaxespad=0.0)
    _style(both)
    fig.tight_layout()
    return _save(fig, path)


def _corner(fig, cell, fitted, colour, title):
    """amp against sigma for one fit: the joint posterior and both
    marginals, with 1 and 2 sigma contours."""
    inner = cell.subgridspec(2, 2, width_ratios=[1.0, 0.35],
                             height_ratios=[0.35, 1.0], wspace=0.05,
                             hspace=0.05)
    joint = fig.add_subplot(inner[1, 0])
    top = fig.add_subplot(inner[0, 0], sharex=joint)
    side = fig.add_subplot(inner[1, 1], sharey=joint)
    amp = fitted["samples"][:, 0]
    lsig = np.log10(fitted["samples"][:, 1])
    alo, ahi = np.percentile(amp, [0.5, 99.5])
    pad = 0.05 * (ahi - alo) if ahi > alo else 1.0
    arange = (alo - pad, ahi + pad)
    srange = (np.log10(bervbias.SIGMA_MIN), np.log10(bervbias.SIGMA_MAX))
    counts, xe, ye = np.histogram2d(lsig, amp, bins=45, range=[srange, arange])
    ordered = np.sort(counts.ravel())[::-1]
    cumulative = np.cumsum(ordered) / ordered.sum()
    # the densities enclosing 39% and 86% of the mass: 1 and 2 sigma in 2D
    levels = sorted({float(ordered[np.searchsorted(cumulative, q)])
                     for q in (0.865, 0.393)})
    joint.imshow(counts.T, origin="lower", aspect="auto", cmap="Greys",
                 extent=[xe[0], xe[-1], ye[0], ye[-1]], alpha=0.55)
    if len(levels) >= 1 and levels[-1] > 0:
        joint.contour(0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]),
                      counts.T, levels=levels, colors=[colour],
                      linewidths=[0.9, 1.4][:len(levels)])
    joint.axhline(0, color=MUTED, lw=0.6)
    joint.set_xlabel("sigma (km/s)", fontsize=7.5, color=INK)
    joint.set_ylabel("amp ((m/s)/(km/s))", fontsize=7.5, color=INK)
    ticks = [t for t in (1, 2, 5, 10, 20, 50)
             if srange[0] <= np.log10(t) <= srange[1]]
    joint.set_xticks(np.log10(ticks))
    joint.set_xticklabels(["%g" % t for t in ticks])
    joint.set_xlim(*srange)
    joint.set_ylim(*arange)
    top.hist(lsig, bins=45, range=srange, color=colour, alpha=0.7)
    side.hist(amp, bins=45, range=arange, color=colour, alpha=0.7,
              orientation="horizontal")
    for ax in (top, side):
        ax.set_yticks([]) if ax is top else ax.set_xticks([])
        plt.setp(ax.get_xticklabels() if ax is top else ax.get_yticklabels(),
                 visible=False)
    top.set_title("%s\ncorrelation of amp and ln sigma: %.2f"
                  % (title, fitted["amp_sigma_r"]),
                  fontsize=7, color=INK, loc="left")
    for ax in (joint, top, side):
        _style(ax)


def figure_corner(before, after, path):
    fb = (before.get("fits") or {})
    if fb.get("before") is None or fb.get("after") is None:
        return None
    fig = plt.figure(figsize=(WIDTH, 3.6))
    cells = fig.add_gridspec(1, 2, wspace=0.35)
    _corner(fig, cells[0], fb["before"], BEFORE, before["label"])
    _corner(fig, cells[1], fb["after"], AFTER, after["label"])
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.14, top=0.84)
    return _save(fig, path)


def figure_d2v(before, after, path):
    if before.get("d2v") is None or after.get("d2v") is None:
        return None
    scale = 1e3                              # (m/s)^2 -> 10^3 (m/s)^2
    fig, axes = plt.subplots(2, 2, figsize=(WIDTH, 5.6))
    both = np.concatenate([before["d2v"], after["d2v"]]) / scale
    good = inliers(both)
    lo, hi = np.percentile(both[good], [0.5, 99.5]) if good.sum() > 4 \
        else (np.nanmin(both), np.nanmax(both))
    pad = 0.08 * (hi - lo) if hi > lo else 1.0
    dlim = (lo - pad, hi + pad)
    vlim = _limits(before["v"] - np.median(before["v"]),
                   after["v"] - np.median(after["v"]))
    for col, (run, colour) in enumerate(((before, BEFORE), (after, AFTER))):
        d2v = run["d2v"] / scale
        err = run["sd2v"] / scale if run.get("sd2v") is not None else None
        outside = int(np.sum((d2v < dlim[0]) | (d2v > dlim[1])))
        top = axes[0, col]
        _points(top, run["t"], d2v, err, colour, "exposures")
        top.set_ylim(*dlim)
        top.set_title("%s d2v\nrobust sigma %.0f, median error %.0f%s"
                      % (run["label"], robust_sigma(d2v),
                         np.nanmedian(err) if err is not None else np.nan,
                         "; %d off the scale" % outside if outside else ""),
                      fontsize=7.5, color=INK, loc="left")
        top.set_xlabel("RJD (BJD - 2400000)", fontsize=7.5, color=INK)
        _style(top)
        bottom = axes[1, col]
        keep = inliers(run["d2v"])
        _points(bottom, d2v[keep], (run["v"] - np.median(run["v"]))[keep],
                run["e"][keep], colour, "exposures")
        bottom.set_xlim(*dlim)
        bottom.set_ylim(*vlim)
        bottom.set_title("%s velocity against d2v\nr = %.2f"
                         % (run["label"],
                            pearson(run["d2v"][keep], run["v"][keep])),
                         fontsize=7.5, color=INK, loc="left")
        bottom.set_xlabel("d2v ($10^3$ m$^2$/s$^2$)", fontsize=7.5, color=INK)
        _style(bottom)
    axes[0, 0].set_ylabel("d2v ($10^3$ m$^2$/s$^2$)", fontsize=7.5, color=INK)
    axes[1, 0].set_ylabel("velocity - median (m/s)", fontsize=7.5, color=INK)
    fig.tight_layout()
    return _save(fig, path)


def figure_change(before, after, path):
    change = after["v"] - before["v"]
    change = change - np.median(change)
    panels = 2 if before.get("berv") is not None else 1
    fig, axes = plt.subplots(1, panels, figsize=(WIDTH, 3.0), sharey=True,
                             squeeze=False)
    lim = _limits(change)
    ax = axes[0, 0]
    _points(ax, before["t"], change, None, INK, "exposures", alpha=0.6)
    ax.set_xlabel("RJD (BJD - 2400000)", fontsize=8, color=INK)
    ax.set_ylabel("corrected - delivered (m/s)", fontsize=8, color=INK)
    ax.set_title("what the correction moved: rms %.2f m/s"
                 % np.std(change), fontsize=8, color=INK, loc="left")
    ax.set_ylim(*lim)
    ax.axhline(0, color=MUTED, lw=0.6)
    _style(ax)
    if panels == 2:
        ax = axes[0, 1]
        _points(ax, before["berv"], change, None, INK, "exposures", alpha=0.6)
        centres, medians, errors = binned(before["berv"], change)
        if centres.size:
            ax.errorbar(centres, medians, yerr=errors, fmt="s", ls="none",
                        ms=5, mfc="white", mec=INK, mew=0.9, ecolor=INK,
                        elinewidth=0.8, capsize=0, zorder=4)
        ax.set_xlabel("BERV (km/s)", fontsize=8, color=INK)
        ax.set_title("against BERV, 2 km/s bins", fontsize=8, color=INK,
                     loc="left")
        ax.axhline(0, color=MUTED, lw=0.6)
        _style(ax)
    fig.tight_layout()
    return _save(fig, path)


def figure_periodograms(before, after, planets, path):
    rows = [("velocity", "v", "e")]
    if before.get("d2v") is not None and after.get("d2v") is not None:
        rows.append(("d2v", "d2v", "sd2v"))
    fig, axes = plt.subplots(len(rows), 1, figsize=(WIDTH, 2.5 * len(rows) + 0.4),
                             sharex=True, squeeze=False)
    drawn = False
    for ax, (what, value, error) in zip(axes[:, 0], rows):
        for run, colour in ((before, BEFORE), (after, AFTER)):
            y, e = run[value], run.get(error)
            if e is None:
                e = np.ones_like(y)
            keep = inliers(y) if value == "d2v" else np.isfinite(y)
            found = periodogram(run["t"][keep], y[keep], e[keep])
            if not found:
                continue
            drawn = True
            periods, power, best, _, fap = found
            ax.plot(periods, power, color=colour, lw=0.9,
                    label="%s: peak %.3g d, FAP %.2g" % (run["label"], best, fap))
        top = ax.get_ylim()[1] if drawn else 1.0
        for letter, period in zip("bcdefgh", planets):
            ax.axvline(period, color=INK, lw=0.7, zorder=0)
            ax.text(period, top, " " + letter, fontsize=7, color=INK,
                    va="top", ha="left")
        for period in (365.25, 182.6, 29.53):
            ax.axvline(period, color=MUTED, lw=0.5, zorder=0, alpha=0.6)
        ax.set_xscale("log")
        ax.set_ylabel("%s power" % what, fontsize=8, color=INK)
        ax.legend(fontsize=6.5, frameon=False, loc="lower left",
                  bbox_to_anchor=(0.0, 1.0), ncol=2, borderaxespad=0.2)
        _style(ax)
    if not drawn:
        plt.close(fig)
        return None
    axes[-1, 0].set_xlabel("period (d); grey lines: 1 yr, 1/2 yr, 1 lunar"
                           " month%s" % ("; black: the known planets"
                                         if planets else ""),
                           fontsize=8, color=INK)
    fig.tight_layout()
    return _save(fig, path)


# ============================================================ the context ===
def find_tree(config, outdir, given=None):
    """The LBL tree this run's rdb files are in, or None.

    lbl.directory in full since 2026-09-16. Before, it was `lbl`, relative to
    wherever the run was started, which for the window is the configuration's
    folder: that one is read off the recorded command, then the working
    directory, then the output root. The first that holds an lblrdb folder
    wins.
    """
    if given:
        return os.path.abspath(given)
    chosen = str((config.get("lbl") or {}).get("directory") or "lbl")
    if os.path.isabs(chosen):
        return chosen
    starts = []
    command = (config.get("provenance") or {}).get("command") or ""
    found = re.search(r"--config\s+(\S+)", command)
    if found:
        starts.append(os.path.dirname(os.path.abspath(found.group(1))))
    starts.append(os.getcwd())
    root = (config.get("output") or {}).get("directory")
    if root:
        starts.append(os.path.abspath(root))
    for start in starts:
        tree = os.path.join(start, chosen)
        if os.path.isdir(os.path.join(tree, "lblrdb")):
            return tree
    return os.path.join(starts[0], chosen)


def stars_of(outdir, config):
    """The stars a run corrected: one, or the members of a joint run."""
    parent = os.path.dirname(os.path.normpath(outdir))
    if os.path.basename(os.path.dirname(parent)) == "joint":
        return os.path.basename(parent).split("+")
    return [str(config["input"]["object"])]


def planets_of(outdir, config, star, joint=False):
    """The published periods of one star, from its OWN configuration.

    A joint run's configuration is its first member's, target included, so
    it answers for that member only: read for the others, it gave GJ1
    Proxima's four planets. A member's own is its cube_config beside the run,
    else its `objects:` block in the configuration file the run was given.
    Nothing, rather than another star's periods, when neither says.
    """
    member = os.path.join(outdir, "cube_config_%s.yaml" % star)
    if os.path.exists(member):
        with open(member) as handle:
            target = (yaml.safe_load(handle) or {}).get("target") or {}
        return [float(p) for p in target.get("planets") or []]
    if not joint or str(star) == str(stars_first(config)):
        target = config.get("target") or {}
        return [float(p) for p in target.get("planets") or []]
    command = (config.get("provenance") or {}).get("command") or ""
    found = re.search(r"--config\s+(\S+)", command)
    if found and os.path.exists(found.group(1)):
        with open(found.group(1)) as handle:
            layered = yaml.safe_load(handle) or {}
        block = ((layered.get("objects") or {}).get(star) or {})
        target = block.get("target") or {}
        return [float(p) for p in target.get("planets") or []]
    return []


def stars_first(config):
    """The member a joint configuration was copied from."""
    command = (config.get("provenance") or {}).get("command") or ""
    found = re.search(r"--objects\s+(\S+)", command)
    if found:
        return found.group(1).split(",")[0]
    return str(config["input"].get("object") or "").split("+")[0]


def newest_input(outdir, star, joint):
    """When what LBL should have measured last changed: the fit, or the
    corrected spectra, whichever is newer. None when neither is there.

    Velocities older than that are from a previous fit of the same run name,
    which is what the report says of them at the figures stage of a run
    done again, before its LBL stage has caught up.
    """
    folder = os.path.join(outdir, "corrected", star) if joint \
        else os.path.join(outdir, "corrected")
    times = [os.path.getmtime(p) for p in glob.glob(os.path.join(folder, "*.fits"))]
    fit = os.path.join(outdir, "fit.npz")
    if os.path.exists(fit):
        times.append(os.path.getmtime(fit))
    return max(times) if times else None


# ============================================================== the stars ===
def first_header(config, star):
    """What one of this star's spectra says of it, as a dictionary, or {}.

    The primary header, and over it the first extension that carries APERO's
    PP_ keys: the object as APERO's database has it, each value with its
    source. PP_RV is in m/s whatever its comment says: NIRPS files label it
    [km/s] and hold GJ 1 at 25534 and Proxima at -22400, which are SIMBAD's
    25.534 and -22.4 km/s. It is given in km/s as PP_RV_KMS.
    """
    from astropy.io import fits

    from .config import spectra_dir

    root = (config.get("input") or {}).get("directory") or ""
    for folder in (os.path.join(root, star), spectra_dir(config)):
        files = sorted(glob.glob(os.path.join(folder, "*.fits")))
        files = [f for f in files if not os.path.basename(f).startswith(".")]
        if not files:
            continue
        try:
            with fits.open(files[0]) as hdulist:
                out = dict(hdulist[0].header)
                for hdu in hdulist[1:]:
                    if "PP_OBJN" in hdu.header:
                        out.update(dict(hdu.header))
                        rv = hdu.header.get("PP_RV")
                        if isinstance(rv, (int, float)):
                            out["PP_RV_KMS"] = rv / 1000.0
                        break
            return out
        except Exception:                                      # noqa: BLE001
            return {}
    return {}


def star_facts(star, config, folder):
    """SIMBAD's view of one star, kept beside the report for next time."""
    from . import simbad

    header = first_header(config, star)
    kept = os.path.join(folder, "simbad_%s.json" % slug(star))
    facts = simbad.star(star, header)
    if facts.get("ok"):
        with open(kept, "w") as handle:
            json.dump(facts, handle, indent=1)
    elif os.path.exists(kept):
        with open(kept) as handle:
            facts = json.load(handle)
        facts["stale"] = True
    facts["header"] = {k: header.get(k) for k in header
                       if (k.startswith("PP_") or k in (
                           "OBJECT", "DRSOBJN", "ESO OCS TARG SPTYPE",
                           "ESO OCS TARG JMAG", "OBJTEMP", "OBJMAG"))
                       and not isinstance(header.get(k), bytes)}
    return facts


def _sexa(value, hours):
    if value is None:
        return "n/a"
    value = value / 15.0 if hours else value
    sign = "-" if value < 0 else ("" if hours else "+")
    value = abs(value)
    d = int(value)
    m = int((value - d) * 60)
    sec = (value - d - m / 60.0) * 3600
    return "%s%02d %02d %0*.*f" % (sign, d, m, 5 if hours else 4,
                                   2 if hours else 1, sec)


def distance_text(facts):
    plx, err = facts.get("plx_value"), facts.get("plx_err")
    if not plx or plx <= 0:
        return "n/a"
    text = "%.3f mas" % plx
    if err:
        text += " $\\pm$ %.3f" % err
    return text + ", %.2f pc" % (1000.0 / plx)


def _pp(header, key, fmt="%s", unit=""):
    """An APERO value with its source, or ''."""
    value = header.get(key)
    if value in (None, "", "None"):
        return ""
    text = (fmt % value) if not isinstance(value, str) else value
    source = header.get(key + "S")
    return tex(text + (" " + unit if unit else "")) + (
        " \\muted{(%s)}" % tex(source)
        if source not in (None, "", "None") else "")


def star_table(star, facts, planets):
    """One star, SIMBAD's description beside the one APERO worked with."""
    header = facts.get("header") or {}
    rows = []

    def add(what, simbad_text, apero_text=""):
        rows.append((what, simbad_text, apero_text))

    apero_name = header.get("PP_OBJN") or header.get("DRSOBJN")
    add("APERO name", "", "\\textbf{%s}%s" % (
        tex(apero_name), " \\muted{(PP\\_OBJN)}" if header.get("PP_OBJN")
        else " \\muted{(DRSOBJN)}") if apero_name else "n/a")
    ok = facts.get("ok")
    if ok:
        add("name", "%s \\muted{(%s)}" % (
            tex(facts["main_id"].replace("NAME ", "")), tex(facts.get("otype"))),
            _pp(header, "PP_OBJNS"))
        found = "from %s" % tex(facts.get("resolved_from"))
        if facts.get("offset") is not None:
            found += (", %.1f arcsec from %s position at %.1f"
                      % (facts["offset"], facts.get("offset_from", "the"),
                         facts.get("offset_epoch", 2000)))
        add("found", found)
        if facts.get("ids"):
            add("identifiers", tex_break(", ".join(facts["ids"])))
    else:
        add("SIMBAD", "\\watch{not reached: %s}"
            % tex(facts.get("error") or "no answer"))
        add("name in the headers", tex(header.get("OBJECT", "n/a")))
    ra_apero = header.get("PP_RA")
    add("position (deg)",
        "%.6f, %+.6f \\muted{(ICRS, J2000)}" % (facts["ra"], facts["dec"])
        if ok and facts.get("ra") is not None else "",
        ("%.6f, %+.6f" % (ra_apero, header.get("PP_DEC"))
         + (" \\muted{(%s, JD %.1f)}" % (tex(header.get("PP_RAS")),
                                           header.get("PP_EPOCH"))
            if header.get("PP_EPOCH") else ""))
        if isinstance(ra_apero, (int, float)) else "")
    if ok:
        sptype = tex(facts.get("sp_type") or "n/a")
    else:
        sptype = tex(header.get("ESO OCS TARG SPTYPE", ""))
    add("spectral type", sptype, _pp(header, "PP_SPT"))
    add("parallax, distance", distance_text(facts) if ok else "",
        _pp(header, "PP_PLX", "%.3f", "mas") + (
            ", %.2f pc" % (1000.0 / header["PP_PLX"])
            if isinstance(header.get("PP_PLX"), (int, float))
            and header["PP_PLX"] > 0 else ""))
    add("proper motion (mas/yr)",
        "%+.2f, %+.2f" % (facts["pmra"], facts.get("pmdec") or 0.0)
        if ok and facts.get("pmra") is not None else "",
        ("%+.2f, %+.2f" % (header["PP_PMRA"], header.get("PP_PMDE", 0.0))
         + (" \\muted{(%s)}" % tex(header.get("PP_PMRAS"))
            if header.get("PP_PMRAS") else ""))
        if isinstance(header.get("PP_PMRA"), (int, float)) else "")
    rv = facts.get("rvz_radvel") if ok else None
    rv_apero = header.get("PP_RV_KMS")
    add("systemic velocity (km/s)",
        ("%.3f" % rv + (" $\\pm$ %.3f" % facts["rvz_err"]
                        if facts.get("rvz_err") else ""))
        if rv is not None else ("not in SIMBAD" if ok else ""),
        ("%.3f" % rv_apero + (" \\muted{(%s)}" % tex(header.get("PP_RVS"))
                              if header.get("PP_RVS") not in (None, "None")
                              else " \\muted{(no source: a placeholder)}"
                              if rv_apero == 0 else ""))
        if isinstance(rv_apero, (int, float)) else "")
    teff_simbad = ""
    measures = facts.get("teff") or [] if ok else []
    if measures:
        latest = measures[0]
        teff_simbad = "%.0f K" % latest["teff"]
        if latest.get("log_g") is not None:
            teff_simbad += ", log g %.2f" % latest["log_g"]
        if latest.get("fe_h") is not None:
            teff_simbad += ", [Fe/H] %+.2f" % latest["fe_h"]
        teff_simbad += " \\muted{(%s)}" % tex(latest["bibcode"])
        if len(measures) > 1:
            values = [m["teff"] for m in measures]
            teff_simbad += ("; %d values, %.0f to %.0f K"
                            % (len(values), min(values), max(values)))
    elif header.get("OBJTEMP"):
        teff_simbad = "%s K \\muted{(OBJTEMP)}" % tex(header["OBJTEMP"])
    add("Teff, latest", teff_simbad, _pp(header, "PP_TEFF", "%.0f", "K"))
    mags = (facts.get("mags") or {}) if ok else {}
    shown = ["%s %.2f" % (band, mags[band]) for band in
             ("B", "V", "G", "J", "H", "K") if mags.get(band) is not None]
    add("magnitudes", ", ".join(shown) or (
        "J %s \\muted{(headers)}" % tex(header["ESO OCS TARG JMAG"])
        if header.get("ESO OCS TARG JMAG") else "n/a"))
    add("known planets", ", ".join("%s %.6g d" % (letter, p) for letter, p in
                                   zip("bcdefgh", planets))
        or "none given in the configuration")
    add("read on", tex(facts.get("retrieved", "")) + (
        " \\watch{(kept from an earlier report: SIMBAD did not answer this"
        " time)}" if facts.get("stale") else ""),
        tex(str(header.get("PP_DDATE", ""))[:10]) + (
            " \\muted{(APERO's database)}" if header.get("PP_DDATE") else ""))
    body = "\n".join("%s & %s & %s \\\\" % (tex(what), a, b)
                     for what, a, b in rows)
    ragged = ">{\\raggedright\\arraybackslash}"
    columns = "".join(ragged + "p{%s\\linewidth}" % w
                      for w in ("0.2", "0.43", "0.3"))
    return ("\\subsection{%s}\n\\begin{longtable}{%s}\n\\toprule\n"
            " & SIMBAD & APERO (PP\\_ keys) \\\\\n\\midrule\n"
            "%s\n\\bottomrule\n\\end{longtable}\n"
            % (tex(star), columns, body))


def star_label(facts):
    """'Smethells 20, M1Ve, 50.7 pc' for the abstract, or ''."""
    if not facts.get("ok"):
        return ""
    parts = [facts["main_id"].replace("NAME ", "")]
    if facts.get("sp_type"):
        parts.append(facts["sp_type"])
    if facts.get("plx_value"):
        parts.append("%.1f pc" % (1000.0 / facts["plx_value"]))
    return ", ".join(parts)


# ============================================================== sections ===
def figure_block(path, caption, short=None, label=None, width=r"\linewidth",
                 page=None):
    """One figure. `short` is what the list of figures says of it: a caption
    is a paragraph, and a list of paragraphs is not a list anybody reads."""
    options = "width=%s,height=0.82\\textheight,keepaspectratio" % width
    if page is not None:
        options += ",page=%d" % page
    return ("\\begin{figure}\n\\centering\n\\includegraphics[%s]{%s}\n"
            "\\caption%s{%s}%s\n\\end{figure}\n"
            % (options, path, "[%s]" % short if short else "", caption,
               "\\label{%s}" % label if label else ""))


def mark(word):
    return {"gain": "\\gain{gain}", "loss": "\\loss{loss}",
            "same": "\\muted{same}", "watch": "\\watch{changed}"}[word]


def summary_table(star, numbers):
    b, a = numbers["before"], numbers["after"]
    rows = []

    def row(name, key, unit="m/s", digits=2, better="lower", tolerance=0.02):
        vb, va = b.get(key, np.nan), a.get(key, np.nan)
        if not (np.isfinite(vb) or np.isfinite(va)):
            return
        change = ("%.2f$\\times$" % (vb / va) if better == "lower"
                  and np.isfinite(vb) and np.isfinite(va) and va > 0 else "")
        rows.append("%s & %s & %s & %s & %s \\\\"
                    % (name, number(vb, digits), number(va, digits), change,
                       mark(verdict(vb, va, better, tolerance))))

    row("rms (m/s)", "rms")
    row("robust sigma (m/s)", "robust")
    row("nightly rms (m/s)", "nightly_rms")
    row("median error (m/s)", "median_error")
    if "bias_peak" in b:
        def said(side):
            if side["bias_detected"]:
                return "%.1f $\\pm$ %.1f (%.1f$\\sigma$)" % (
                    side["bias_peak"], side["bias_err"],
                    side["bias_significance"])
            return "$<$ %.1f" % side["bias_upper"]
        rows.append("BERV bias at its peak (m/s) & %s & %s & & %s \\\\"
                    % (said(b), said(a), mark(bias_verdict(b, a))))
        rows.append("its width sigma (km/s) & %s & %s & & \\\\" % tuple(
            ("%.1f $\\pm$ %.1f" % (side["bias_width"], side["bias_width_err"])
             if side["bias_detected"] else "\\muted{unconstrained}")
            for side in (b, a)))
        rows.append("amp-sigma correlation (posterior) & %s & %s & & \\\\"
                    % (number(b["bias_amp_sigma_r"]),
                       number(a["bias_amp_sigma_r"])))
        row("jitter beyond LBL's errors (m/s)", "jitter")
    if "berv_binned" in b:
        row("scatter of BERV-binned medians (m/s)", "berv_binned")
    if "d2v_sigma" in b:
        vb, va = b["d2v_sigma"] / 1e3, a["d2v_sigma"] / 1e3
        word = "watch" if verdict(vb, va, tolerance=0.1) != "same" else "same"
        rows.append("d2v robust sigma ($10^3$ m$^2$/s$^2$) & %s & %s & & %s \\\\"
                    % (number(vb, 1), number(va, 1), mark(word)))
        rows.append("velocity-d2v correlation r & %s & %s & & \\muted{activity}"
                    " \\\\" % (number(b["d2v_r"], 3), number(a["d2v_r"], 3)))
    if "peak_period" in b:
        rows.append("highest velocity peak (d) & %s & %s & & \\\\"
                    % (number(b["peak_period"], 3), number(a["peak_period"], 3)))
        rows.append("its false-alarm probability & %s & %s & & \\\\"
                    % (tex("%.2g" % b["peak_fap"]), tex("%.2g" % a["peak_fap"])))
    for letter, period, kb, ka in zip("bcdefgh", numbers["planet_periods"],
                                      b["planets"], a["planets"]):
        moved = abs(ka[0] - kb[0]) > 2 * np.hypot(ka[1], kb[1])
        rows.append("K at %s d (planet %s, m/s) & %s $\\pm$ %s & %s $\\pm$ %s"
                    " & & %s \\\\"
                    % (number(period, 4), letter, number(kb[0]), number(kb[1]),
                       number(ka[0]), number(ka[1]),
                       mark("watch" if moved else "same")))
    rows.append("\\midrule")
    rows.append("rms of corrected - delivered (m/s) & \\multicolumn{2}{c}{%s}"
                " & & \\\\" % number(numbers["change_rms"]))
    rows.append("noise removed in quadrature (m/s) & \\multicolumn{2}{c}{%s}"
                " & & \\\\" % (number(numbers["removed"])
                               if np.isfinite(numbers["removed"])
                               else "none: the rms grew"))
    return ("\\begin{longtable}{p{0.42\\linewidth}rrrl}\n"
            "\\caption{%s: the %d exposures both series share, over %d"
            " nights.}\\\\\n\\toprule\n & delivered & corrected & gain & \\\\\n"
            "\\midrule\n\\endhead\n%s\n\\bottomrule\n\\end{longtable}\n"
            % (tex(star), numbers["n"], numbers["nights"], "\n".join(rows)))


def bias_verdict(b, a):
    """gain, loss or same for the fitted BERV bias.

    Measured against the errors, not a fixed tolerance: a bias that was
    detected and is gone is a gain, one that appeared is a loss, and two
    detections are compared by their difference in units of its error.
    """
    if not ("bias_detected" in b and "bias_detected" in a):
        return "same"
    if b["bias_detected"] and not a["bias_detected"]:
        return "gain"
    if a["bias_detected"] and not b["bias_detected"]:
        return "loss"
    if not b["bias_detected"]:
        return "same"
    change = abs(b["bias_peak"]) - abs(a["bias_peak"])
    error = np.hypot(b["bias_err"], a["bias_err"])
    return "gain" if change > error else "loss" if -change > error else "same"


def bias_words(side):
    if side.get("bias_detected"):
        return "%.0f $\\pm$ %.0f m/s" % (side["bias_peak"], side["bias_err"])
    if "bias_upper" in side:
        return "none detected ($<$ %.0f m/s)" % side["bias_upper"]
    return "not fitted"


def verdict_lines(numbers):
    """(better, worse, moved) as lists of what, for the summary's bullets."""
    b, a = numbers["before"], numbers["after"]
    better, worse, moved = [], [], []
    for key, what, how in (("rms", "rms", "lower"),
                           ("robust", "robust sigma", "lower"),
                           ("nightly_rms", "nightly rms", "lower"),
                           ("median_error", "LBL's error bars", "lower"),
                           ("jitter", "the jitter beyond LBL's errors",
                            "lower"),
                           ("berv_binned", "the structure against BERV",
                            "lower")):
        if key not in b:
            continue
        word = verdict(b[key], a[key], how, 0.05 if how == "zero" else 0.02)
        (better if word == "gain" else worse if word == "loss"
         else []).append(what)
    word = bias_verdict(b, a)
    (better if word == "gain" else worse if word == "loss"
     else []).append("the fitted BERV bias")
    if "d2v_sigma" in b and verdict(b["d2v_sigma"], a["d2v_sigma"],
                                    tolerance=0.1) != "same":
        moved.append("d2v's scatter (%.0f to %.0f, in $10^3$ m$^2$/s$^2$)"
                     % (b["d2v_sigma"] / 1e3, a["d2v_sigma"] / 1e3))
    for letter, period, kb, ka in zip("bcdefgh", numbers["planet_periods"],
                                      b["planets"], a["planets"]):
        if abs(ka[0] - kb[0]) > 2 * np.hypot(ka[1], kb[1]):
            moved.append("the amplitude of planet %s (%.4g d)"
                         % (letter, period))
    return better, worse, moved


def activity_note(numbers):
    """A sentence when the velocity's strongest period is also d2v's."""
    a = numbers["after"]
    if same_period(a.get("peak_period"), a.get("d2v_peak_period"),
                   numbers.get("baseline")):
        return ("The corrected velocity's strongest peak, at %.3g d, is also"
                " d2v's: activity (the rotation or its alias) rather than a"
                " planet." % a["peak_period"])
    return ""


def star_headline(star, numbers):
    """The abstract's sentence for one star: short, and a planet only when
    the correction changed its amplitude."""
    b, a = numbers["before"], numbers["after"]
    ratio = b["rms"] / a["rms"] if a["rms"] > 0 else np.nan
    words = ("%s: rms %.2f to %.2f m/s (%s), nightly rms %.2f to %.2f m/s"
             % (tex(star), b["rms"], a["rms"],
                ("%.2f times lower" % ratio) if ratio >= 1
                else ("%.2f times higher" % (1 / ratio)),
                b["nightly_rms"], a["nightly_rms"]))
    if "bias_peak" in b:
        words += ", BERV bias %s to %s" % (bias_words(b), bias_words(a))
    words += "."
    _better, _worse, moved = verdict_lines(numbers)
    if moved:
        words += " \\watch{Changed where nothing should: %s.}" % "; ".join(moved)
    return words


def star_sentence(star, numbers):
    b, a = numbers["before"], numbers["after"]
    ratio = b["rms"] / a["rms"] if a["rms"] > 0 else np.nan
    words = ("%s: on the %d exposures both series share, the rms goes from"
             " %.2f to %.2f m/s (%s) and the nightly rms from %.2f to %.2f m/s."
             % (tex(star), numbers["n"], b["rms"], a["rms"],
                ("%.2f times lower" % ratio) if ratio >= 1
                else ("%.2f times higher" % (1 / ratio)),
                b["nightly_rms"], a["nightly_rms"]))
    if "bias_peak" in b:
        words += (" The bias that follows BERV, fitted, goes from %s to %s."
                  % (bias_words(b), bias_words(a)))
    for letter, period, kb, ka in zip("bcdefgh", numbers["planet_periods"],
                                      b["planets"], a["planets"]):
        words += (" At the %.4g d period of planet %s, K goes from %.2f to"
                  " %.2f m/s." % (period, letter, kb[0], ka[0]))
    return words


def velocity_section(star, before, after, numbers, folder, stale):
    base = "rv-" + slug(star)
    figures = os.path.join(folder, "figures")
    parts = ["\\subsection{%s}\n" % tex(star)]
    if stale:
        parts.append("\\watch{These velocities are older than this run's"
                     " fit or its corrected spectra: they are from an earlier"
                     " run under the same name, and LBL has not measured this"
                     " one yet.}\n\n")
    parts.append(star_sentence(star, numbers) + " " + activity_note(numbers)
                 + "\n")
    name = tex(star)
    drawn = [
        (figure_time(before, after, os.path.join(figures, base + "-time.pdf")),
         "%s: the velocities over the campaign" % name,
         "The velocities over the campaign, delivered above and corrected"
         " below, on the same scale. %s; no line joins them."
         % ("Small points are exposures, large ones the weighted nightly"
            " means" if numbers["nights"] < 0.8 * numbers["n"]
            else "One point per exposure, about one a night")),
        (figure_berv(before, after, os.path.join(figures, base + "-berv.pdf")),
         "%s: the BERV bias, fitted" % name,
         "The velocities against the barycentric velocity, with the bias a"
         " telluric line blended with the stellar lines produces, fitted by"
         " MCMC: $v = c + a\\,B\\,e^{-B^2/2\\sigma^2}$, with a jitter added"
         " to LBL's error bars. The line is the posterior median, the band its"
         " 1$\\sigma$ envelope (the band alone for a bias that is not"
         " detected); squares are medians in 2 km/s bins. Above,"
         " each series less its fitted offset; below, both envelopes on one"
         " axis. The peak is the bias at $B = \\pm\\sigma$,"
         " $a\\sigma e^{-1/2}$; below 3$\\sigma$ from zero, only an upper"
         " limit on it is quoted."),
        (figure_corner(before, after, os.path.join(figures, base + "-corner.pdf")),
         "%s: the covariance of amp and sigma" % name,
         "The joint posterior of the bias's amplitude $a$ and width"
         " $\\sigma$, delivered and corrected, with its 1 and 2$\\sigma$"
         " contours and both marginals. The two trade against each other,"
         " since a narrower bias needs a larger amplitude to reach the same"
         " points; a posterior filling the width's prior is a bias that is"
         " not there."),
        (figure_d2v(before, after, os.path.join(figures, base + "-d2v.pdf")),
         "%s: d2v, the activity indicator" % name,
         "d2v, LBL's second-derivative term, which follows the line width and"
         " is an activity indicator. Above, over the campaign; below, the"
         " velocity against it. A correction of the observer frame has no"
         " business changing it, and a velocity that correlates with it is"
         " activity rather than noise."),
        (figure_change(before, after, os.path.join(figures, base + "-change.pdf")),
         "%s: what the correction moved" % name,
         "What the correction changed, exposure by exposure: the corrected"
         " velocity less the delivered one, over time and against BERV."),
        (figure_periodograms(before, after, numbers["planet_periods"],
                             os.path.join(figures, base + "-periods.pdf")),
         "%s: periodograms of the velocity and of d2v" % name,
         "Lomb-Scargle periodograms of the velocity and of d2v, delivered and"
         " corrected. A known planet should keep its peak; a peak at a year or"
         " its harmonics is the Earth's."),
    ]
    for path, short, caption in drawn:
        if path:
            parts.append(figure_block("figures/" + os.path.basename(path),
                                      caption, short=short))
    parts.append("\\clearpage\n")
    return "".join(parts)


def run_section(config, manifest, folder):
    from .config import WINDOW_SETTINGS, setting_value

    parts = ["\\section{The run}\n"]
    rows = (manifest or {}).get("run") or []
    if rows:
        parts.append("\\begin{longtable}{p{0.34\\linewidth}p{0.6\\linewidth}}\n"
                     "\\toprule\n")
        for key, value in rows:
            parts.append("%s & %s \\\\\n" % (tex(key), tex_break(value)))
        parts.append("\\bottomrule\n\\end{longtable}\n")
    parts.append("\\subsection{Every setting the window can set}\n"
                 "\\begin{longtable}{p{0.28\\linewidth}p{0.2\\linewidth}"
                 "p{0.44\\linewidth}}\n\\toprule\nsetting & value & what it"
                 " does \\\\\n\\midrule\n\\endhead\n")
    for path, what in WINDOW_SETTINGS:
        value = setting_value(config, path)
        if isinstance(value, (list, tuple)):
            value = ", ".join(str(v) for v in value)
        elif isinstance(value, bool) or value is None:
            # as the configuration spells them, not as Python prints them
            value = {True: "true", False: "false", None: "null"}[value]
        parts.append("\\texttt{%s} & %s & %s \\\\\n"
                     % (tex_break(path), tex_break(value), tex(what)))
    parts.append("\\bottomrule\n\\end{longtable}\n")
    provenance = config.get("provenance") or {}
    command = provenance.get("command") or "not recorded (before 2026-09-15)"
    with open(os.path.join(folder, "command.txt"), "w") as handle:
        handle.write(wrap_lines(command, 90))
    parts.append("\\subsection{The command and the code}\n"
                 "The command this run was given:\n{\\footnotesize"
                 "\\verbatiminput{command.txt}}\n"
                 "Code: pca2d %s%s.\n"
                 % (tex((provenance.get("commit") or "unknown")[:12]),
                    tex(" with local changes") if provenance.get("dirty")
                    else ""))
    return "".join(parts)


def run_rows(outdir, cube=None):
    """The run's summary rows, for a run whose report folder has none.

    The figures stage records them since 2026-09-16; before, they were only
    printed on the bound PDF's first page. They are worked out again from the
    run's configuration and its fit, less the cube, which nothing recorded.
    """
    import types

    from .config import load_config
    from .figures.bundle import summary_rows

    path = os.path.join(outdir, "resolved_config.yaml")
    fit_path = os.path.join(outdir, "fit.npz")
    try:
        config = load_config(path)
        fit = np.load(fit_path) if os.path.exists(fit_path) else None
        rows = summary_rows(config, types.SimpleNamespace(
            outdir=outdir, config=path, cube=cube or "not recorded"), fit)
    except Exception as exc:                                   # noqa: BLE001
        log("  the run's summary could not be worked out again (%s)" % exc,
            "warn")
        return []
    if cube is None:
        rows = [(k, v) for k, v in rows
                if k not in ("cube", "this config's cube key")]
    # a block with no component has no shares to print
    rows = [(k, v) for k, v in rows if str(v).strip() not in ("[]", "")]
    return [[str(k), str(v)] for k, v in rows]


def correction_section(manifest, windows=()):
    figures = (manifest or {}).get("figures") or []
    if not figures:
        return ("\\section{The correction}\nThe figures stage has not drawn"
                " anything for this run yet.\n")
    parts = ["\\section{The correction}\n"]
    for entry in figures:
        title, path, pages = entry["title"], entry["file"], entry.get("pages", 1)
        parts.append("\\subsection{%s}\n" % tex(title))
        caption = tex(CAPTIONS.get(title, title))
        if pages <= 1:
            parts.append(figure_block(path, caption, short=tex(title)))
        else:
            # the sequence draws one page per window, in the config's order
            named = (list(windows) if title.startswith("The sequence")
                     and len(windows) == pages else [None] * pages)
            for page, window in zip(range(1, pages + 1), named):
                where = (" Window %s (centre:width, nm)." % tex(window)
                         if window else "")
                parts.append(figure_block(
                    path, "%s%s (%d of %d)" % (caption, where, page, pages),
                    short=("%s, window %s" % (tex(title), tex(window))
                           if window else "%s, %d of %d"
                           % (tex(title), page, pages)),
                    page=page))
        parts.append("\\clearpage\n")
    return "".join(parts)


def failure_section(manifest):
    failures = (manifest or {}).get("failures") or []
    if not failures:
        return ""
    parts = ["\\section{What did not build}\n\\begin{itemize}\n"]
    for item in failures:
        parts.append("\\item %s\\\\ \\muted{%s}\n"
                     % (("Not drawn: " if item.get("skipped") else "Failed: ")
                        + tex(item.get("what", "")),
                        tex(item.get("why", ""))[:600]))
    parts.append("\\end{itemize}\n")
    return "".join(parts)


def parameters_appendix(config, folder):
    body = yaml.safe_dump({k: v for k, v in config.items()},
                          sort_keys=False, default_flow_style=False, width=90)
    with open(os.path.join(folder, "parameters.yaml"), "w") as handle:
        handle.write(wrap_lines(body, 96))
    return ("\\appendix\n\\section{The parameters, as resolved}\n"
            "Every key of the configuration this run used, after the defaults,"
            " the instrument, the object and the command line were merged.\n"
            "{\\scriptsize\\verbatiminput{parameters.yaml}}\n")


# =============================================================== the build ===
def compile_tex(folder, name, pdflatex):
    """Run pdflatex twice (the contents need the second pass). True if done."""
    for _ in range(2):
        try:
            done = subprocess.run(
                [pdflatex, "-interaction=nonstopmode", "-halt-on-error",
                 name + ".tex"], cwd=folder, capture_output=True, text=True,
                timeout=600)
        except (OSError, subprocess.SubprocessError) as exc:
            log("pdflatex could not run: %s" % exc, "warn")
            return False
        if done.returncode != 0:
            lines = [line for line in (done.stdout or "").splitlines()
                     if line.startswith("!") or line.startswith("l.")]
            log("pdflatex stopped on %s.tex: %s"
                % (os.path.join(folder, name), " | ".join(lines[:4])
                   or (done.stdout or "")[-300:]), "warn")
            return False
    return os.path.exists(os.path.join(folder, name + ".pdf"))


def render(outdir, config=None, lbl_dir=None, out=None, stars=None):
    """Write the run's LaTeX report and compile it. Returns the PDF, or None.

    `outdir` is the run's folder. The report is written in <outdir>/report
    and copied to the run's own PDF name, <object>_<tag>.pdf, where the window
    and the run's last line look for it.
    """
    outdir = os.path.abspath(outdir)
    pdflatex = find_pdflatex()
    if pdflatex is None:
        log("no pdflatex on this machine, so the report stays the bound PDF of"
            " the figures. A TeX installation (MacTeX, or BasicTeX) gives the"
            " full report", "warn")
        return None
    if config is None:
        with open(os.path.join(outdir, "resolved_config.yaml")) as handle:
            config = yaml.safe_load(handle)
    folder = os.path.join(outdir, "report")
    os.makedirs(os.path.join(folder, "figures"), exist_ok=True)
    tag = os.path.basename(outdir)
    name = "%s_%s" % (config["input"].get("object") or "run", tag)
    out = out or os.path.join(outdir, name + ".pdf")

    manifest = read_manifest(folder)
    if manifest is None and os.path.exists(out):
        log("no report folder in %s: cutting the bound PDF along its"
            " bookmarks to write one" % outdir, "info")
        # the report is about to take its name: the old one is kept whole,
        # so that writing a run's report again can always be undone
        kept = os.path.join(folder, "bound_before_the_report.pdf")
        if not os.path.exists(kept):
            shutil.copyfile(out, kept)
            log("  the bound PDF as it was is kept: %s" % kept, "info")
        manifest = {"figures": split_bundle(out, folder), "failures": [],
                    "run": run_rows(outdir)}
        write_manifest(folder, manifest)

    stars = stars or stars_of(outdir, config)
    # <root>[/_name]/<object>/<tag>, or <root>[/_name]/joint/<A+B>/<tag>
    above = os.path.dirname(os.path.dirname(outdir))
    joint = os.path.basename(above) == "joint"
    run_name = os.path.basename(os.path.dirname(above) if joint else above)
    tree = find_tree(config, outdir, lbl_dir)
    from .lbl import object_names
    results, sections, missing = [], [], []
    for star in stars:
        before_name, after_name = object_names(config, star, tag)
        before = load(tree, before_name, "delivered")
        after = load(tree, after_name, "corrected")
        if before is None or after is None:
            missing.append((star, rdb_file(tree, before_name if before is None
                                           else after_name)))
            continue
        before, after = common(before, after)
        if before["t"].size < 4:
            missing.append((star, "fewer than four exposures in common"))
            continue
        planets = planets_of(outdir, config, star, joint)
        import zlib
        numbers = star_numbers(before, after, planets,
                               seed=zlib.crc32(star.encode()) & 0xffffffff)
        before["fits"] = {"before": numbers["before"].get("bias"),
                          "after": numbers["after"].get("bias")}
        newest = newest_input(outdir, star, joint)
        stale = newest is not None and after["mtime"] < newest
        results.append((star, numbers, stale))
        sections.append(velocity_section(star, before, after, numbers, folder,
                                         stale))

    facts = {star: star_facts(star, config, folder) for star in stars}
    body = ["\\section{The star%s}\n" % ("s" if len(stars) > 1 else "")]
    body.append("As SIMBAD describes %s, found from the names in the"
                " spectra's headers and the folder's own.\n"
                % ("them" if len(stars) > 1 else "it"))
    for star in stars:
        body.append(star_table(star, facts[star],
                               planets_of(outdir, config, star, joint)))
    body.append("\\section{Summary}\n")
    if results:
        body.append("One table per star: the delivered and corrected velocities"
                    " LBL measured, on the exposures both have. The last"
                    " column says which way each number went: \\gain{gain} or"
                    " \\loss{loss} where lower is better, \\watch{changed}"
                    " where the correction should have changed nothing.\n")
        for star, numbers, _stale in results:
            better, worse, moved = verdict_lines(numbers)
            items = []
            if better:
                items.append("\\item \\gain{Better}: %s." % ", ".join(better))
            if worse:
                items.append("\\item \\loss{Worse}: %s." % ", ".join(worse))
            if moved:
                items.append("\\item \\watch{Changed where nothing should}:"
                             " %s." % "; ".join(moved))
            note = activity_note(numbers)
            if note:
                items.append("\\item %s" % note)
            body.append("\\subsection*{%s}\n" % tex(star))
            if items:
                body.append("\\begin{itemize}\n%s\n\\end{itemize}\n"
                            % "\n".join(items))
            body.append(summary_table(star, numbers))
    for star, why in missing:
        body.append("\\muted{%s: no velocities to compare yet (%s).}\n\n"
                    % (tex(star), tex(why)))
    if sections:
        body.append("\\clearpage\n\\section{Velocities}\n")
        body.extend(sections)
    body.append(run_section(config, manifest, folder))
    body.append("\\clearpage\n")
    body.append(correction_section(
        manifest, (manifest or {}).get("windows")
        or (config.get("output") or {}).get("windows") or ()))
    body.append(failure_section(manifest))
    body.append(parameters_appendix(config, folder))

    tw = config.get("twoframe") or {}
    what = ("%d star and %d observer components, %s"
            % (tw.get("n_star", 0), tw.get("n_earth", 0),
               "one observer basis for %s" % ", ".join(stars)
               if len(stars) > 1 else "on %s alone" % stars[0]))
    named = ["%s (%s)" % (star, star_label(facts[star]))
             if star_label(facts[star]) else star for star in stars]
    abstract = [tex("Two-frame PCA correction of %s: %s."
                    % (" and ".join(named), what))]
    abstract.extend(star_headline(star, numbers) for star, numbers, _ in results)
    if not results:
        abstract.append("LBL has not measured this run yet, so the report"
                        " shows the correction and not its effect on the"
                        " velocities.")
    if any(stale for _, _, stale in results):
        abstract.append("\\watch{Some velocities predate the corrected spectra;"
                        " see their section.}")

    with open(TEMPLATE) as handle:
        document = handle.read()
    fill = {
        "TITLE": tex(" + ".join(stars)),
        "SUBTITLE": tex("pca2d-preclean run %s%s" % (
            tag, ", %s" % run_name.lstrip("_") if run_name.startswith("_")
            else "")),
        "DATE": tex(datetime.date.today().isoformat()),
        "RUNHEAD": tex("%s, %s" % (" + ".join(stars), tag)),
        "ABSTRACT": "\n\n".join(abstract),
        "BODY": "".join(body),
    }
    for key, value in fill.items():
        document = document.replace("<<%s>>" % key, value)
    with open(os.path.join(folder, name + ".tex"), "w") as handle:
        handle.write(document)
    if not compile_tex(folder, name, pdflatex):
        return None
    shutil.copyfile(os.path.join(folder, name + ".pdf"), out)
    log("report written: %s (%d star%s with velocities)"
        % (out, len(results), "" if len(results) == 1 else "s"), "value")
    return out


def parse_args(argv=None):
    p = argparse.ArgumentParser(prog="python -m pca2d.texreport",
                                description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", required=True,
                   help="the run's folder: <output root>/<object>/<tag>")
    p.add_argument("--lbl-dir", default=None,
                   help="LBL's tree, when the run's configuration does not say"
                        " where it is")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    done = render(args.run, lbl_dir=args.lbl_dir)
    return 0 if done else 1


if __name__ == "__main__":
    raise SystemExit(main())
