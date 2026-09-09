#!/usr/bin/env python
"""One exposure, before and after the correction, in each window of interest.

Drawn in the STAR'S rest frame, like every other figure in the set, so the same
wavelength means the same stellar feature from one page to the next.

    python diagnostics/sample_before_after.py --cube <cube> --fit <fit.npz> \
        --windows 1200.3:2 1593.6:2 1669.5:5 --out sample.pdf

The snippet figures show all the spectra as an image and the envelope figure
shows their spread. Neither shows what the correction did to ONE spectrum, which
is the thing an eye reads fastest: the flux as it arrived, the flux as it leaves,
and the model that was taken out between them.

One page per window. The exposure drawn is the one whose correction has the
median amplitude, so the page is representative rather than flattering: picking
the largest correction would show the method at its most impressive and tell you
nothing about a typical night.
"""

from __future__ import annotations

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

from pca2d.plotting import live_mask, sample_source_file, window_parity
from pca2d.twoframe import (LanczosShifter, load_cube, row_parity,
                                 star_model)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cube", required=True)
    p.add_argument("--fit", required=True)
    p.add_argument("--windows", nargs="+", required=True,
                   help="centre:width pairs in nm, e.g. 1200.3:2 1669.5:5")
    p.add_argument("--row", type=int, default=None,
                   help="cube row to draw; default is the median correction")
    p.add_argument("--frame", choices=["star", "observer"], default="star",
                   help="which rest frame to draw in; the star's by default")
    p.add_argument("--out", required=True)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    grid, data, w, meta = load_cube(args.cube)
    n, m = data.shape
    parity = row_parity(meta, n)
    # for windows an order overlap covers twice, see window_parity
    sample_path = sample_source_file(args.cube, meta)
    fit = np.load(args.fit)
    delta = -np.asarray(fit["berv"], dtype=float) / float(fit["dv"])
    shifter = LanczosShifter(m, a=8, max_shift=int(np.ceil(np.abs(delta).max())) + 2)
    model = (star_model(fit["P"], fit["a"], shifter, delta, n, m)
             + fit["b"] @ fit["Q"])


    # Into the star's frame, one array at a time so only one extra copy of an
    # (N, M) array is alive at once. The exposure drawn is then on the same
    # wavelength axis as every other figure in the set.
    if args.frame == "star":
        w = shifter.rows(w, -delta)
        w = np.where(live_mask(w), w, 0.0)   # see plotting.live_mask
        data = shifter.rows(data, -delta)
        model = shifter.rows(model, -delta)
        data[w <= 0] = 0.0

    live = live_mask(w)
    with np.errstate(invalid="ignore"):
        amplitude = np.array([np.nanstd(np.where(live[i], model[i], np.nan))
                              for i in range(n)])
    names = [str(v) for v in fit["filename"]] if "filename" in fit.files else []
    print("  median correction rms over rows: %.4f" % np.nanmedian(amplitude))

    pages = []
    for spec in args.windows:
        centre, _, width = spec.partition(":")
        pages.append((float(centre), float(width or 2.0)))

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with PdfPages(args.out) as pdf:
        for centre, width in pages:
            lo, hi = centre - 0.5 * width, centre + 0.5 * width
            # the row is chosen PER WINDOW: a cube row holds one order parity,
            # and the two parities cover different halves of the domain, so a
            # row picked once for the whole figure is empty in half the windows
            band = (grid >= lo) & (grid <= hi)
            covers = (live[:, band].sum(axis=1) > 0.5 * band.sum())
            # In an order overlap both parities cover the window, one near its
            # order's middle and one at an edge. Draw the middle one, so the
            # page is not a picture of the worse of two measurements.
            own, off, n_cov = (window_parity(centre, sample_path) if sample_path
                               else (None, float("nan"), 0))
            if own is not None and n_cov > 1 and (covers & (parity == own)).any():
                covers = covers & (parity == own)
                print("  %.1f-%.1f nm: two orders reach it, drawing parity %d"
                      " (%.2f of a half-width from the order centre)"
                      % (lo, hi, own, off))
            if args.row is not None:
                row = args.row
            elif covers.any():
                candidates = np.where(covers)[0]
                order = candidates[np.argsort(amplitude[candidates])]
                row = int(order[order.size // 2])
            else:
                print("  %.1f-%.1f nm: no row covers this window, skipped"
                      % (lo, hi))
                continue
            label = names[row] if row < len(names) else "row %d" % row
            win = band & live[row]
            if win.sum() < 10:
                print("  %.1f-%.1f nm: only %d live columns, skipped"
                      % (lo, hi, win.sum()))
                continue
            x = grid[win]
            before = data[row][win]
            corr = model[row][win]
            after = before - corr

            fig, axes = plt.subplots(2, 1, figsize=(9.0, 5.4), sharex=True,
                                     gridspec_kw={"height_ratios": [2.0, 1.0]})
            ax = axes[0]
            ax.plot(x, before, lw=0.9, color="0.35", label="before")
            ax.plot(x, after, lw=0.9, color="#1f4e9c", label="after")
            ax.axhline(0, color="0.85", lw=0.6)
            ax.set_ylabel(r"$\ln(f/\mathrm{savgol}\,f)$", fontsize=9)
            ax.legend(fontsize=8, loc="upper right", framealpha=0.9)
            ax.set_title("%s   %.2f-%.2f nm\nscatter %.4f before, %.4f after"
                         "  (%+.1f%%)"
                         % (label, lo, hi, np.std(before), np.std(after),
                            100 * (np.std(after) / max(np.std(before), 1e-30) - 1)),
                         fontsize=9)

            ax = axes[1]
            ax.plot(x, corr, lw=0.9, color="#b3261e")
            ax.axhline(0, color="0.85", lw=0.6)
            ax.set_ylabel("removed", fontsize=9)
            ax.set_xlabel("wavelength (nm), %s frame" % args.frame, fontsize=9)
            ax.text(0.005, 0.06, "rms %.4f" % np.std(corr), fontsize=7.5,
                    transform=ax.transAxes, color="0.35")

            for a in axes:
                a.tick_params(labelsize=8)
                a.grid(alpha=0.15)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
            print("  %.2f-%.2f nm: scatter %.4f -> %.4f (%+.1f%%), removed rms %.4f"
                  % (lo, hi, np.std(before), np.std(after),
                     100 * (np.std(after) / max(np.std(before), 1e-30) - 1),
                     np.std(corr)))
    print("wrote %s" % args.out)


if __name__ == "__main__":
    main()
