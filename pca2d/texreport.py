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


def weighted_line(x, y, e):
    """(slope, its error, intercept) of y against x, weighted by 1/e^2.

    The error is scaled up by the reduced chi2 when that is above one: LBL's
    error bars are known to be optimistic against the scatter of a night.
    """
    x, y, e = (np.asarray(a, float) for a in (x, y, e))
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (e > 0)
    if ok.sum() < 4 or np.ptp(x[ok]) == 0:
        return np.nan, np.nan, np.nan
    w = 1.0 / e[ok] ** 2
    A = np.column_stack([np.ones(ok.sum()), x[ok]])
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    p = cov @ ((A * w[:, None]).T @ y[ok])
    chi2r = float(np.sum(w * (y[ok] - A @ p) ** 2) / max(ok.sum() - 2, 1))
    return float(p[1]), float(np.sqrt(cov[1, 1] * max(chi2r, 1.0))), float(p[0])


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


def star_numbers(before, after, planets=()):
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
            slope, err, _ = weighted_line(run["berv"], run["v"], run["e"])
            _, medians, _ = binned(run["berv"],
                                   run["v"] - np.median(run["v"]))
            row.update(berv_r=pearson(run["berv"], run["v"]),
                       berv_slope=slope, berv_slope_err=err,
                       berv_binned=float(np.std(medians))
                       if medians.size >= 3 else np.nan)
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


def figure_berv(before, after, path):
    if before.get("berv") is None or after.get("berv") is None:
        return None
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 3.4), sharex=True,
                             sharey=True)
    lim = _limits(before["v"] - np.median(before["v"]),
                  after["v"] - np.median(after["v"]))
    for ax, run, colour in zip(axes, (before, after), (BEFORE, AFTER)):
        y = run["v"] - np.median(run["v"])
        _points(ax, run["berv"], y, run["e"], colour, "exposures", alpha=0.45)
        centres, medians, errors = binned(run["berv"], y)
        if centres.size:
            ax.errorbar(centres, medians, yerr=errors, fmt="s", ls="none",
                        ms=5, mfc="white", mec=INK, mew=0.9, ecolor=INK,
                        elinewidth=0.8, capsize=0, zorder=4,
                        label="2 km/s bins, median")
        slope, err, intercept = weighted_line(run["berv"], y, run["e"])
        if np.isfinite(slope):
            xs = np.array([np.nanmin(run["berv"]), np.nanmax(run["berv"])])
            ax.plot(xs, intercept + slope * xs, color=INK, lw=1.0, zorder=5,
                    label="weighted line")
        ax.set_title("%s\nslope %.2f $\\pm$ %.2f (m/s)/(km/s), r = %.2f"
                     % (run["label"], slope, err,
                        pearson(run["berv"], run["v"])),
                     fontsize=7.5, color=INK, loc="left")
        ax.axhline(0, color=MUTED, lw=0.6)
        ax.set_xlabel("BERV (km/s)", fontsize=8, color=INK)
        ax.set_ylim(*lim)
        _style(ax)
    axes[0].set_ylabel("velocity - median (m/s)", fontsize=8, color=INK)
    axes[0].legend(fontsize=6.5, frameon=False, loc="upper left")
    fig.tight_layout()
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
        ax.set_title("against BERV: r = %.2f" % pearson(before["berv"], change),
                     fontsize=8, color=INK, loc="left")
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
    if "berv_r" in b:
        row("velocity-BERV correlation r", "berv_r", digits=3, better="zero",
            tolerance=0.05)
        row("velocity-BERV slope (m/s per km/s)", "berv_slope", digits=3,
            better="zero", tolerance=0.02)
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


def verdict_lines(numbers):
    """(better, worse, moved) as lists of what, for the summary's bullets."""
    b, a = numbers["before"], numbers["after"]
    better, worse, moved = [], [], []
    for key, what, how in (("rms", "rms", "lower"),
                           ("robust", "robust sigma", "lower"),
                           ("nightly_rms", "nightly rms", "lower"),
                           ("median_error", "LBL's error bars", "lower"),
                           ("berv_r", "the correlation with BERV", "zero"),
                           ("berv_binned", "the structure against BERV",
                            "lower")):
        if key not in b:
            continue
        word = verdict(b[key], a[key], how, 0.05 if how == "zero" else 0.02)
        (better if word == "gain" else worse if word == "loss"
         else []).append(what)
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
    if "berv_r" in b:
        words += ", correlation with BERV r = %.2f to %.2f" % (b["berv_r"],
                                                               a["berv_r"])
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
    if "berv_r" in b:
        words += (" The correlation of the velocity with BERV goes from"
                  " r = %.2f to r = %.2f." % (b["berv_r"], a["berv_r"]))
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
         "%s: the velocities against BERV" % name,
         "The same velocities against the barycentric velocity. What the"
         " observer frame leaves in the velocities follows BERV, so a slope or"
         " a structure in the binned medians is telluric or instrumental, and"
         " the correction should remove it."),
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
        numbers = star_numbers(before, after, planets)
        newest = newest_input(outdir, star, joint)
        stale = newest is not None and after["mtime"] < newest
        results.append((star, numbers, stale))
        sections.append(velocity_section(star, before, after, numbers, folder,
                                         stale))

    body = ["\\section{Summary}\n"]
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
    abstract = [tex("Two-frame PCA correction of %s: %s."
                    % (" and ".join(stars), what))]
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
