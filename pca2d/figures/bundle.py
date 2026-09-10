#!/usr/bin/env python
"""Everything a run produced, as ONE multipage PDF.

    python -m pca2d.figures.bundle --config outputs/TOI2120/1-3v/resolved_config.yaml \
        --cube cache/cube_tfits_<key> --outdir outputs/TOI2120/1-3v

Scattering forty PDFs across three directories makes them hard to look at in
order, hard to send, and easy to compare across runs by accident: two files
with the same name in two directories are one careless click apart. A run now
leaves one document.

It opens with what was actually run, not what a config file says today: the
resolved parameters, after the package defaults are merged in, as plain
selectable text you can copy back into a YAML. Then the figures, grouped, with
PDF bookmarks so the sections are one click away.

The individual figures are built into a temporary directory and removed once
they are bound in. Pass --keep to leave them.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.backends.backend_pdf import PdfPages
from pypdf import PdfReader, PdfWriter

HERE = os.path.dirname(os.path.abspath(__file__))
# No sys.path edit here. There was one, pointing at the pca2d/ directory rather
# than at the repository, so it never made `import pca2d` work (the install
# does that); what it did do was put pca2d/ at the head of the path of the whole
# pipeline, since the CLI imports this module in-process. From then on the name
# `lbl` meant pca2d/lbl.py and not the LBL package, and a full run died at its
# last stage with "attempted relative import with no known parent package".

from pca2d.config import cache_key, load_config, spectra_dir  # noqa: E402
from pca2d.logger import log                              # noqa: E402

#: default windows, the ones this campaign looks at
WINDOWS = ["1200.3:2", "1220:4", "1267:2", "1593.6:2", "1669.5:5", "1700.5:3",
           "2196.2:2", "2450:5"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True)
    p.add_argument("--cube", required=True)
    p.add_argument("--outdir", required=True,
                   help="the fit's output directory; fit.npz is read from it"
                        " and the bundle is written into it")
    p.add_argument("--windows", nargs="+", default=WINDOWS)
    p.add_argument("--source-dir", default=None)
    p.add_argument("--out", default=None,
                   help="default: <outdir>/<object>_<tag>.pdf")
    p.add_argument("--keep", action="store_true",
                   help="leave the individual figures on disk as well")
    return p.parse_args(argv)


def wrap(text, width=96):
    """Break long lines so nothing runs off the page."""
    out = []
    for line in text.splitlines():
        while len(line) > width:
            cut = line.rfind(" ", 0, width)
            cut = cut if cut > width // 2 else width
            out.append(line[:cut])
            line = "  " + line[cut:].lstrip()
        out.append(line)
    return out


def text_pages(pdf, title, body, per_page=62, size=7.0):
    """`body` as selectable monospace text, paginated."""
    lines = wrap(body)
    for start in range(0, len(lines), per_page):
        chunk = lines[start:start + per_page]
        fig = plt.figure(figsize=(8.5, 11))
        fig.text(0.06, 0.965, title if start == 0 else title + " (continued)",
                 fontsize=11, weight="bold", va="top")
        fig.text(0.06, 0.93, "\n".join(chunk), fontsize=size, va="top",
                 family="monospace", linespacing=1.35)
        pdf.savefig(fig)
        plt.close(fig)


def object_name(config, args):
    """What to call this run.

    `input.object` when it is set, and otherwise the config file's own stem,
    which is how these are named. NOT the source directory's basename: that is
    `tfiles` for every target here, so a bundle would come out called tfiles.
    """
    named = config["input"].get("object")
    if named:
        return str(named)
    # outputs/<object>/<tag> is where this is being written, so the parent of
    # the output directory names the object. Tried before the config file's
    # stem, which is "config" whenever the generic config.yaml is used
    # directly, and before the source directory's basename, which is "tfiles"
    # for every target.
    parent = os.path.basename(os.path.dirname(os.path.normpath(args.outdir)))
    if parent and parent not in ("", ".", "outputs"):
        return parent
    stem = os.path.splitext(os.path.basename(args.config))[0]
    if stem and stem != "config":
        return stem
    return os.path.basename(os.path.normpath(spectra_dir(config)))


def summary(config, args, fit):
    """The few numbers a reader wants before anything else."""
    tw = config["twoframe"]
    rows = []
    rows.append(("object", object_name(config, args)))
    rows.append(("source", spectra_dir(config)))
    key = cache_key(config)
    here = os.path.basename(os.path.normpath(args.cube))
    rows.append(("cube", here))
    rows.append(("this config's cube key", key + (
        "" if key in here else
        "   <-- DOES NOT MATCH THE CUBE ABOVE. The cache key hashes the"
        " configuration, not the data, so a config edited after the cube was"
        " built points somewhere else. Every figure below came from the cube"
        " named above.")))
    rows.append(("domain", "%.1f to %.1f nm at dv = %.2f km/s"
                 % (config["domain"]["wave_min"], config["domain"]["wave_max"],
                    config["domain"]["dv"])))
    rows.append(("high-pass", "savgol, %d samples, order %d, %s"
                 % (config["highpass"]["window"], config["highpass"]["polyorder"],
                    config["highpass"]["mode"])))
    rows.append(("components", "%d star + %d observer"
                 % (tw["n_star"], tw["n_earth"])))
    rows.append(("nightly stacking", str(config["input"].get("nightly_stack"))))
    if fit is not None:
        n = int(np.asarray(fit["a"]).shape[0])
        rows.append(("rows fitted", "%d" % n))
        chi2n = float(fit["chi2_null"]) if "chi2_null" in fit.files else np.nan
        chi2b = float(fit["chi2_best"]) if "chi2_best" in fit.files else np.nan
        if np.isfinite(chi2n) and chi2n > 0:
            rows.append(("weighted variance removed", "%.2f%%"
                         % (100 * (1 - chi2b / chi2n))))
        ps = np.asarray(fit["power_star"]) if "power_star" in fit.files else None
        pe = np.asarray(fit["power_earth"]) if "power_earth" in fit.files else None
        if ps is not None and np.isfinite(chi2n) and chi2n > 0:
            rows.append(("star block shares",
                         np.array2string(100 * ps / chi2n, precision=2)))
            rows.append(("observer block shares",
                         np.array2string(100 * pe / chi2n, precision=2)))
        if "rejected" in fit.files:
            rows.append(("rejected by the MAD cut",
                         "%d" % int(np.asarray(fit["rejected"]).sum())))
    width = max(len(k) for k, _ in rows)
    return "\n".join("%-*s   %s" % (width, k, v) for k, v in rows)


def run(cmd, failures):
    """Run one diagnostic, and say so if it fails rather than dying.

    The second argument is the list of failures to append to. It was called
    `log` until it shadowed the logger of that name, which every call in here
    then went through, so a failing figure raised TypeError instead of being
    written to the bundle's last page.
    """
    log("   " + " ".join(os.path.basename(c) if c.endswith(".py") else c
                         for c in cmd[:4]) + " ...")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        failures.append((" ".join(cmd),
                         (r.stderr or r.stdout or "").strip()[-400:]))
        log("      failed, see the bundle's last page")
    return r.returncode == 0


def main(argv=None):
    args = parse_args(argv)
    config = load_config(args.config)
    fit_path = os.path.join(args.outdir, "fit.npz")
    fit = np.load(fit_path) if os.path.exists(fit_path) else None
    obj = object_name(config, args)
    tag = os.path.basename(os.path.normpath(args.outdir))
    out = args.out or os.path.join(args.outdir, "%s_%s.pdf" % (obj, tag))
    tmp = tempfile.mkdtemp(prefix="bundle_")
    failures = []

    # ---- the figures ---------------------------------------------------
    py = sys.executable
    d = lambda name: os.path.join(HERE, name)                 # noqa: E731
    log("  building the figures")
    run([py, d("sequence.py"), "--cube", args.cube, "--fit", fit_path,
         "--windows", *args.windows, "--source-dir",
         args.source_dir or spectra_dir(config),
         "--out", os.path.join(tmp, "sequence.pdf")], failures)
    run([py, d("sample_before_after.py"), "--cube", args.cube, "--fit", fit_path,
         "--windows", *args.windows, "--out",
         os.path.join(tmp, "before_after.pdf")], failures)
    run([py, d("weight_spectrum.py"), "--cube", args.cube, "--fit", fit_path,
         "--out", os.path.join(tmp, "weights.pdf")], failures)
    src = args.source_dir or spectra_dir(config)
    run([py, d("oh_residual.py"), "--cube", args.cube, "--fit", fit_path,
         "--source-dir", src, "--out", os.path.join(tmp, "oh.pdf")], failures)
    run([py, d("coeff_periodogram.py"), "--fit", fit_path,
         "--out", os.path.join(tmp, "periodogram.pdf")], failures)
    # Only scripts that take --windows are used, so the cube is loaded once per
    # script and not once per window. The per-window ones cost a full cube read
    # each: on 321 exposures that is 4.1 GB and a minute and a half, and with
    # three of them over eight windows the bundle spent an hour and a half
    # reading the same file twenty-four times. plot_snippet, star_frame_envelope
    # and compare_template are what went; sequence.py shows the whole chain for
    # every window in one read, and the stacks show the rest.
    per_window = []

    # ---- the front matter ----------------------------------------------
    front = os.path.join(tmp, "_front.pdf")
    with PdfPages(front) as pdf:
        fig = plt.figure(figsize=(8.5, 11))
        fig.text(0.5, 0.62, obj, fontsize=34, ha="center", weight="bold")
        fig.text(0.5, 0.565, "two-frame weighted PCA, %s" % tag, fontsize=13,
                 ha="center", color="0.3")
        fig.text(0.5, 0.52, "%d star + %d observer components"
                 % (config["twoframe"]["n_star"], config["twoframe"]["n_earth"]),
                 fontsize=11, ha="center", color="0.3")
        pdf.savefig(fig)
        plt.close(fig)
        text_pages(pdf, "What was run", summary(config, args, fit), size=8.5)
        text_pages(pdf, "Parameters, as resolved (copy this back into a YAML)",
                   yaml.safe_dump(config, sort_keys=False, default_flow_style=False,
                                  width=94))

    # ---- bind ------------------------------------------------------------
    order = [("Front matter", front)]
    # the whole chain first, one page per window: it is the overview every
    # other figure is a detail of
    order.append(("The sequence, step by step", os.path.join(tmp, "sequence.pdf")))
    order.append(("Variance", os.path.join(args.outdir, "variance.pdf")))
    # not the basis vectors: nine curves over half a million samples is a page
    # nobody can read, and what each component DOES is on the two pages that
    # follow, the coefficients against time and against everything measurable
    order.append(("Coefficients against time",
                  os.path.join(args.outdir, "coefficients_vs_time.pdf")))
    order.append(("Correlations", os.path.join(args.outdir, "correlations.pdf")))
    order.append(("Coefficient periodogram", os.path.join(tmp, "periodogram.pdf")))
    order.append(("Weight spectrum", os.path.join(tmp, "weights.pdf")))
    order.append(("OH residual", os.path.join(tmp, "oh.pdf")))
    order.append(("One exposure, before and after",
                  os.path.join(tmp, "before_after.pdf")))
    order.extend(per_window)

    writer = PdfWriter()
    page = 0
    bound = 0
    for title, path in order:
        if not os.path.exists(path):
            continue
        try:
            reader = PdfReader(path)
        except Exception as exc:                              # noqa: BLE001
            failures.append((path, str(exc)))
            continue
        writer.add_outline_item(title, page)
        for p in reader.pages:
            writer.add_page(p)
            page += 1
        bound += 1
    if failures:
        note = os.path.join(tmp, "_failures.pdf")
        with PdfPages(note) as pdf:
            text_pages(pdf, "What did not build",
                       "\n\n".join("%s\n    %s" % (c, e) for c, e in failures))
        reader = PdfReader(note)
        writer.add_outline_item("What did not build", page)
        for p in reader.pages:
            writer.add_page(p)
            page += 1

    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "wb") as fh:
        writer.write(fh)

    # A run leaves ONE document. The figures the fit itself wrote into the
    # output directory are now pages of it, and leaving them beside it is how
    # a directory fills with forty files that disagree about which run they
    # came from. The CSVs stay: they are data, not pages.
    if not args.keep:
        removed = 0
        for _title, path in order:
            if (path.startswith(os.path.abspath(args.outdir))
                    or os.path.dirname(os.path.abspath(path))
                    == os.path.abspath(args.outdir)) and path != out:
                try:
                    os.remove(path)
                    removed += 1
                except OSError:
                    pass
        if removed:
            log("  folded %d loose figures into the bundle and removed them"
                  % removed)
    if args.keep:
        keep = os.path.join(args.outdir, "figures")
        shutil.rmtree(keep, ignore_errors=True)
        shutil.copytree(tmp, keep)
        log("  individual figures kept in %s" % keep)
    shutil.rmtree(tmp, ignore_errors=True)
    log("wrote %s: %d pages from %d figures%s"
          % (out, page, bound,
             ", %d failed" % len(failures) if failures else ""))
    return None


if __name__ == "__main__":
    main()
