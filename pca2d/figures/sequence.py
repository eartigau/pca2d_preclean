#!/usr/bin/env python
"""Every step from the file to the corrected spectrum, on one page per window.

    python diagnostics/sequence.py --cube <cube> --fit <fit.npz> \
        --windows 1267:2 1669.5:5 --out sequence.pdf

The other figures each show one stage. This shows the chain, in the order the
pipeline applies it, on the same rows and the same wavelength axis, so that
what each stage removed is the difference between two panels the eye can put
side by side.

    1  the high pass: ln(f) minus a Savitzky-Golay of ln(f), and nothing else
    2  the reconstruction, M star + N observer components
    3  minus the OBSERVER block: what a corrected file holds
    4  minus everything: what nobody explained

Panel 3 is the product. The correction removes the observer block and nothing
else, because with no stellar template the first star component IS the star and
dividing it out would flatten the spectrum LBL is about to measure.

Everything is in the STAR'S rest frame with the rows ordered by barycentric
velocity, so a stellar feature is vertical and anything anchored to the Earth
slants.
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pca2d.plotting import (live_mask, nan_cmap, raw_log_flux_block,  # noqa: E402
                                 sample_source_file, window_parity)
from pca2d.twoframe import (LanczosShifter, fit_means, load_cube,  # noqa: E402
                                 mean_rows, row_parity, star_model)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cube", required=True)
    p.add_argument("--fit", required=True)
    p.add_argument("--windows", nargs="+", required=True,
                   help="centre:width pairs in nm")
    p.add_argument("--source-dir", default=None)
    p.add_argument("--n-overplot", type=int, default=5)
    p.add_argument("--out", required=True)
    return p.parse_args(argv)


def _filled(block, ok):
    """`block` with its gaps at each row's own median, ready to be shifted.

    A shift operator has to be handed finite numbers, and which finite number
    matters: in ln flux zero means a flux of one, which is fifty times the
    continuum here, and the Lanczos kernel spreads it over sixteen samples on
    either side. The row's own median is the value that disturbs nothing.
    """
    block = np.asarray(block, dtype=float)
    with np.errstate(invalid="ignore"):
        fill = np.nanmedian(np.where(ok, block, np.nan), axis=-1, keepdims=True)
    fill = np.where(np.isfinite(fill), fill, 0.0)
    return np.where(ok, block, np.broadcast_to(fill, block.shape))


def panel(ax, image, x, scale, title, cmap=None):
    im = ax.imshow(image, aspect="auto", cmap=cmap or nan_cmap("RdBu_r"),
                   vmin=-scale, vmax=scale, origin="upper",
                   extent=[x[0], x[-1], image.shape[0] - 0.5, -0.5])
    ax.set_title(title, fontsize=8.5)
    ax.set_ylabel("ordered by BERV", fontsize=8)
    ax.tick_params(labelsize=7)
    return im


def main(argv=None):
    args = parse_args(argv)
    grid, data, w0, meta = load_cube(args.cube)
    n, m = data.shape
    parity = row_parity(meta, n)
    berv = np.asarray(meta["berv"], dtype=float)
    fit = np.load(args.fit)
    dv = float(fit["dv"])
    delta = -np.asarray(fit["berv"], dtype=float) / dv
    shifter = LanczosShifter(m, a=8, max_shift=int(np.ceil(np.abs(delta).max())) + 2)

    # what the fit actually took out, read from the fit rather than recomputed
    means, group = fit_means(fit, meta, n, m)
    offset = mean_rows(means, group, n)
    n_star = int(fit["P"].shape[0])
    n_earth = int(fit["Q"].shape[0])
    star = star_model(fit["P"], fit["a"], shifter, delta, n, m)
    earth = fit["b"] @ fit["Q"]

    # The flux as delivered and the fit's own input are both gone from the
    # page. The first answered a question once, that the log of it lands on the
    # panel below and the filter eats nothing; the second differs from that
    # panel by an instrumental offset the eye cannot see.
    steps = [
        ("given", data - offset,
         "1. the high pass: $\\ln(f)$ minus a Savitzky-Golay of $\\ln(f)$"),
        ("model", star + earth,
         "2. the reconstruction, %d star + %d observer" % (n_star, n_earth)),
        ("corrected", data - offset - earth,
         "3. minus the OBSERVER block: what a corrected file holds"),
        ("resid", data - offset - star - earth,
         "4. minus everything: what nobody explained"),
    ]

    home = {}
    alive = live_mask(shifter.rows(w0, -delta))
    for name, arr, _ in steps:
        z = shifter.rows(arr, -delta)
        z[~alive] = np.nan
        home[name] = z

    names = [os.path.basename(str(v)) for v in meta["filename"]]
    sample_path = sample_source_file(args.cube, meta, args.source_dir)

    specs = []
    for spec in args.windows:
        c, _, w = spec.partition(":")
        specs.append((float(c), float(w or 2.0)))
    margin = int(np.ceil(np.abs(delta).max())) + 16
    blocks = {}
    for c, w in specs:
        j = np.where((grid >= c - 0.5 * w) & (grid <= c + 0.5 * w))[0]
        if j.size >= 10:
            blocks[c] = (max(0, j[0] - margin), min(m, j[-1] + 1 + margin))
    if blocks:
        cols = np.unique(np.concatenate([np.arange(a, b)
                                         for a, b in blocks.values()]))
        raw = raw_log_flux_block(args.cube, args.source_dir, names, parity,
                                 grid, cols)
    else:
        cols, raw = np.zeros(0, dtype=int), np.zeros((n, 0))

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with PdfPages(args.out) as pdf:
        for centre, width in specs:
            lo, hi = centre - 0.5 * width, centre + 0.5 * width
            win = (grid >= lo) & (grid <= hi)
            if win.sum() < 10:
                print("  %.1f-%.1f nm: too few columns, skipped" % (lo, hi))
                continue
            x = grid[win]
            cover = np.isfinite(home["given"][:, win]).sum(axis=1) / win.sum()
            keep = cover > 0.5 * max(cover.max(), 1e-9)
            own, off, ncov = (window_parity(centre, sample_path) if sample_path
                              else (None, float("nan"), 0))
            # Where two orders reach this window, keep the one whose middle it
            # sits nearest: the other measures the same wavelengths at an order
            # edge, where the blaze has fallen away. Done silently. It is not a
            # caveat, it is simply the right order to draw.
            if own is not None and ncov > 1 and (keep & (parity == own)).sum() >= 6:
                keep &= parity == own
                print("  %.1f-%.1f nm: two orders reach it, drawing the one"
                      " %.2f of a half-width from its centre" % (lo, hi, off))
            rows = np.where(keep)[0][np.argsort(berv[keep])]
            if rows.size < 6:
                print("  %.1f-%.1f nm: too few rows, skipped" % (lo, hi))
                continue
            jw = np.where(win)[0]
            block = {k: home[k][np.ix_(rows, jw)] for k in home}
            finite = block["given"][np.isfinite(block["given"])]
            scale = float(np.percentile(np.abs(finite), 98)) if finite.size else 1.0

            n_img = len(steps)
            fig, axes = plt.subplots(n_img + 1, 1,
                                     figsize=(9.4, 1.85 * n_img + 3.0),
                                     sharex=True,
                                     gridspec_kw={"height_ratios":
                                                  [1.0] * n_img + [1.7]})
            for r, (name, _, title) in enumerate(steps):
                extra = ("" if name in ("given", "model")
                         else "   scatter %.4f" % np.nanstd(block[name]))
                im = panel(axes[r], block[name], x, scale, title + extra)

            # the same rows in flux, before and after
            ax = axes[-1]
            pick = min(args.n_overplot, rows.size)
            spread = rows[np.linspace(0, rows.size - 1, pick).astype(int)]
            drawn = 0
            if centre in blocks:
                a0, b0 = blocks[centre]
                sub = LanczosShifter(b0 - a0, a=8, max_shift=margin)
                for cube_row in spread:
                    blk = raw[cube_row, np.searchsorted(cols, np.arange(a0, b0))]
                    if not np.isfinite(blk).any():
                        continue
                    ok = np.isfinite(blk)
                    lnf = sub.rows(_filled(blk[None, :], ok[None, :]),
                                   -delta[[cube_row]])[0]
                    live = sub.rows(ok.astype(float)[None, :],
                                    -delta[[cube_row]])[0] > 0.99
                    lnf = np.where(live, lnf, np.nan)
                    atm = sub.rows(earth[[cube_row], a0:b0], -delta[[cube_row]])[0]
                    f_in = np.exp(lnf[jw - a0])
                    f_out = np.exp((lnf - atm)[jw - a0])
                    norm = np.nanmedian(f_in)
                    if not np.isfinite(norm) or norm <= 0:
                        continue
                    ax.plot(x, f_in / norm, lw=0.8, color="#d62728", alpha=0.5,
                            label="before" if drawn == 0 else None)
                    ax.plot(x, f_out / norm, lw=0.8, color="#2ca02c", alpha=0.5,
                            label="observer block removed" if drawn == 0 else None)
                    drawn += 1
            if drawn:
                ax.set_ylabel("flux / its own median", fontsize=8.5)
                ax.legend(fontsize=8, loc="lower right", framealpha=0.9)
                ax.set_title("%d of those rows in FLUX, spread over BERV %+.1f to"
                             " %+.1f km/s" % (drawn, berv[spread].min(),
                                              berv[spread].max()), fontsize=8.5)
                ax.grid(alpha=0.15)
                ax.tick_params(labelsize=8)
            else:
                ax.text(0.5, 0.5, "could not re-read the source spectra",
                        ha="center", va="center", fontsize=8,
                        transform=ax.transAxes)
            axes[-1].set_xlabel("wavelength (nm), star frame", fontsize=9)
            # attached to EVERY axes, the trace panel included. A colourbar
            # steals width from the axes it is given, so leaving one out makes
            # it wider than the rest and the wavelength axes stop lining up
            # down the page, which is the whole point of the figure.
            fig.colorbar(im, ax=list(axes), fraction=0.022, pad=0.015,
                         label=r"$\ln(f/\mathrm{savgol}\,f)$")
            fig.suptitle("%.2f-%.2f nm: every step, in order, in the STAR'S REST"
                         " FRAME\nvertical structure belongs to the star;"
                         " anything slanted does not" % (lo, hi), fontsize=10)
            pdf.savefig(fig)
            plt.close(fig)
            print("  %.1f-%.1f nm: given %.4f | corrected %.4f (%+.0f%%) |"
                  " residual %.4f"
                  % (lo, hi, np.nanstd(block["given"]),
                     np.nanstd(block["corrected"]),
                     100 * (np.nanstd(block["corrected"]) /
                            max(np.nanstd(block["given"]), 1e-30) - 1),
                     np.nanstd(block["resid"])))
    print("wrote %s" % args.out)
    return None


if __name__ == "__main__":
    main()
