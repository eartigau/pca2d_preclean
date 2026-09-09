#!/usr/bin/env python
"""What the template's median actually sees, before it takes the median.

    python diagnostics/star_frame_stack.py --cube <cube> --fit <fit.npz> \
        --centre 1200.3 --width 2.0 --out stack.pdf

Every high-passed spectrum shifted into the star's rest frame and shown as an
image, with the rows ORDERED BY BARYCENTRIC VELOCITY. That ordering is the whole
point of the figure: in this frame a stellar feature sits at the same column in
every row and draws a vertical line, while anything fixed in the observer's
frame, an OH airglow residual or a telluric, moves by -BERV/dv columns from row
to row and draws a DIAGONAL. One glance separates the two.

Below the image, three things the median has to reconcile: the plain median over
all spectra, the hierarchical median that gives each BERV bin one vote, and the
individual bin medians. Where the bin medians disagree with each other, the
template is being asked to average things that are not the same, and the
disagreement is the diagnosis.

The bins are LBL's: 3000 m/s wide, at least 3 exposures each, see
lbl/recipes/lbl_template.py.
"""

from __future__ import annotations

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

from pca2d.plotting import (live_mask, nan_cmap, sample_source_file,
                                 window_parity)
from pca2d.twoframe import (LanczosShifter, load_cube, row_parity,
                                 star_model)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cube", required=True)
    p.add_argument("--fit", required=True)
    p.add_argument("--windows", nargs="+", required=True,
                   help="centre:width pairs in nm, e.g. 1200.3:2 1669.5:5;"
                        " one page per window, same convention as"
                        " sample_before_after.py")
    p.add_argument("--berv-bin", type=float, default=3000.0, help="m/s")
    p.add_argument("--min-entries", type=int, default=3)
    p.add_argument("--out", required=True,
                   help="a .pdf gives one multi-page file; a DIRECTORY gives one"
                        " stack_<lo>-<hi>nm.pdf per window, which is what you"
                        " want when the windows are looked at one at a time")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    grid, data, w, meta = load_cube(args.cube)
    n, m = data.shape
    parity = row_parity(meta, n)
    # for windows an order overlap covers twice, see window_parity
    sample_path = sample_source_file(args.cube, meta)
    fit = np.load(args.fit)
    dv = float(fit["dv"])
    berv = np.asarray(fit["berv"], dtype=float)
    delta = -berv / dv
    shifter = LanczosShifter(m, a=8, max_shift=int(np.ceil(np.abs(delta).max())) + 2)

    # the two blocks, so the same stack can be shown before and after them. The
    # template and the parity means are deliberately NOT removed: they are the
    # stellar spectrum, and taking them out would leave two images of noise with
    # nothing to line up by eye.
    blocks = (star_model(fit["P"], fit["a"], shifter, delta, n, m)
              + fit["b"] @ fit["Q"])
    home = shifter.rows(data, -delta)
    home_after = shifter.rows(data - blocks, -delta)
    home_w = shifter.rows(w, -delta)
    alive = live_mask(home_w)
    home[~alive] = np.nan
    home_after[~alive] = np.nan

    pages = []
    for spec in args.windows:
        centre, _, width = spec.partition(":")
        pages.append((float(centre), float(width or 2.0)))

    per_file = args.out.endswith(os.sep) or os.path.isdir(args.out) or not \
        args.out.lower().endswith(".pdf")
    if per_file:
        os.makedirs(args.out, exist_ok=True)
        for centre, width in pages:
            lo, hi = centre - 0.5 * width, centre + 0.5 * width
            path = os.path.join(args.out, "stack_%.0f-%.0fnm.pdf" % (lo, hi))
            with PdfPages(path) as pdf:
                one_page(pdf, grid, home, home_after, home_w, berv, centre,
                         width, args, parity=parity, sample_path=sample_path)
            print("     -> %s" % path)
    else:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with PdfPages(args.out) as pdf:
            for centre, width in pages:
                one_page(pdf, grid, home, home_after, home_w, berv, centre,
                         width, args, parity=parity, sample_path=sample_path)
        print("wrote %s" % args.out)


def one_page(pdf, grid, home, home_after, home_w, berv, centre, width, args,
             parity=None, sample_path=None):
    n = home.shape[0]
    lo, hi = centre - 0.5 * width, centre + 0.5 * width
    win = (grid >= lo) & (grid <= hi)
    if win.sum() < 10:
        print("  %.1f-%.1f nm: only %d columns, skipped" % (lo, hi, win.sum()))
        return
    x = grid[win]
    block = home[:, win]
    # Keep a row only if it covers a decent share of the window RELATIVE TO THE
    # BEST rows, not in absolute terms. The two order parities cover different
    # halves of the domain, and in a window one parity owns, the other's rows
    # are 20 to 25 per cent covered: above an absolute 20 per cent floor, so
    # they survive it, and then fill three quarters of the image with yellow
    # while carrying almost nothing. Judged against the best row they are
    # correctly dropped.
    cover = np.isfinite(block).sum(axis=1) / max(win.sum(), 1)
    keep = cover > 0.5 * np.max(cover) if cover.size and np.max(cover) > 0 else cover > 1
    # Where consecutive orders overlap, the coverage rule keeps both parities
    # because both genuinely cover the window: one near its order's middle and
    # one at an order edge, which is the same wavelengths measured worse. Keep
    # the middle one only.
    own, off, n_cov = (window_parity(centre, sample_path) if sample_path
                       else (None, float("nan"), 0))
    reason = "they cover under half of what the best row does"
    if (own is not None and n_cov > 1 and parity is not None
            and (keep & (parity == own)).sum() >= 6):
        keep &= parity == own
        reason = ("two orders reach this window and this one sits %.2f of a"
                  " half-width from its order centre" % off)
        print("  %.1f-%.1f nm: two orders reach it, drawing parity %d only"
              % (lo, hi, own))
    if keep.sum() < 6:
        print("  %.1f-%.1f nm: only %d rows cover it, skipped" % (lo, hi, keep.sum()))
        return
    after = home_after[:, win]
    block, after, b_berv = block[keep], after[keep], berv[keep]
    order = np.argsort(b_berv)
    block, after, b_berv = block[order], after[order], b_berv[order]


    labels = np.floor((b_berv - b_berv.min()) / (args.berv_bin / 1000.0)).astype(int)
    bins, bin_meds = [], []
    for value in np.unique(labels):
        rows = labels == value
        if rows.sum() < args.min_entries:
            continue
        bins.append((value, int(rows.sum()), float(np.median(b_berv[rows]))))
        with np.errstate(invalid="ignore"):
            bin_meds.append(np.nanmedian(block[rows], axis=0))
    bin_meds = np.vstack(bin_meds) if bin_meds else np.zeros((0, x.size))


    with np.errstate(invalid="ignore"):
        plain = np.nanmedian(block, axis=0)
    hier = (np.nanmedian(bin_meds, axis=0) if bin_meds.shape[0] >= 2 else plain)

    spread = np.nanstd(bin_meds, axis=0) if bin_meds.shape[0] >= 2 else np.zeros_like(x)
    scale = np.nanpercentile(np.abs(block[np.isfinite(block)]), 98) if np.isfinite(block).any() else 1.0

    # Stacked, not side by side: the two images then share the wavelength axis,
    # so a column that changes can be followed straight down the page instead of
    # being hunted for in the other panel.
    fig = plt.figure(figsize=(9.6, 12.4))
    gs = fig.add_gridspec(4, 1, height_ratios=[2.0, 2.0, 1.0, 0.9], hspace=0.10)
    ax_b = fig.add_subplot(gs[0])
    ax_a = fig.add_subplot(gs[1], sharex=ax_b, sharey=ax_b)
    axes = [ax_b, fig.add_subplot(gs[2], sharex=ax_b),
            fig.add_subplot(gs[3], sharex=ax_b)]

    with np.errstate(invalid="ignore"):
        sd_b = float(np.nanstd(block))
        sd_a = float(np.nanstd(after))
    for ax, img, name in (
            (ax_b, block, "before,  scatter %.4f" % sd_b),
            (ax_a, after, "after the two blocks,  %.4f  (%+.1f%%)"
             % (sd_a, 100 * (sd_a / max(sd_b, 1e-30) - 1)))):
        im = ax.imshow(img, aspect="auto", cmap=nan_cmap("RdBu_r"),
                       vmin=-scale, vmax=scale, origin="upper",
                       extent=[x[0], x[-1], img.shape[0] - 0.5, -0.5])
        ax.set_title(name, fontsize=9)
        ax.tick_params(labelsize=8)
        ax.set_xlim(x[0], x[-1])
    ax_b.tick_params(labelbottom=False)
    ax_a.tick_params(labelbottom=False)
    ax_a.set_ylabel("spectrum, ordered by BERV", fontsize=9)
    fig.colorbar(im, ax=[ax_b, ax_a], fraction=0.022, pad=0.015,
                 label=r"$\ln(f/\mathrm{savgol}\,f)$")
    ax = ax_b
    ax.set_ylabel("spectrum, ordered by BERV", fontsize=9)
    dropped = int((~keep).sum())
    if dropped:
        ax.text(0.988, 0.965, "%d rows of the other parity dropped\n(%s)"
                % (dropped, reason), fontsize=6.5, color="0.35", ha="right",
                va="top",
                transform=ax.transAxes,
                bbox=dict(fc="white", ec="0.85", lw=0.5, pad=2))
    # the slope an observer-frame feature would follow, drawn from the middle
    span = (b_berv.max() - b_berv.min())
    if span > 0:
        mid = 0.5 * (x[0] + x[-1])
        lam = mid * (1.0 + (b_berv - b_berv[block.shape[0] // 2]) / 299792.458)
        ax.plot(lam, np.arange(block.shape[0]), lw=1.0, color="#1f8f4e", ls="--",
                alpha=0.85)
        ax.text(0.012, 0.04, "green dashed: the track an observer-frame feature"
                             " follows here", fontsize=7.5, color="#1f8f4e",
                transform=ax.transAxes)
    fig.suptitle("%.2f-%.2f nm in the star's rest frame, before any median\n"
                 "vertical structure is the star; diagonal structure is not"
                 % (lo, hi), fontsize=10)

    ax = axes[1]
    for row, (_, count, centre) in zip(bin_meds, bins):
        ax.plot(x, row, lw=0.7, alpha=0.65,
                label="BERV %+.1f (%d)" % (centre, count))
    ax.plot(x, plain, lw=1.5, color="k", label="plain median")
    ax.plot(x, hier, lw=1.5, color="#b3261e", ls="--", label="hierarchical")
    ax.set_ylabel("median per BERV bin", fontsize=9)
    ax.legend(fontsize=6.5, ncol=3, loc="upper right", framealpha=0.9)

    ax = axes[2]
    ax.plot(x, spread, lw=1.0, color="#1f4e9c")
    ax.axhline(np.nanmedian(spread), color="0.6", lw=0.8, ls=":")
    ax.set_ylabel("scatter between\nbin medians", fontsize=9)
    ax.set_xlabel("wavelength (nm), star frame", fontsize=9)
    ax.text(0.005, 0.85, "high here = the BERV bins disagree, so the star block"
                         " is being asked\nto describe things that are not the"
                         " same", fontsize=7.5,
            transform=ax.transAxes, color="0.35", va="top")

    axes[1].tick_params(labelbottom=False)
    for a in axes[1:]:
        a.set_xlim(x[0], x[-1])
        a.tick_params(labelsize=8)
        a.grid(alpha=0.15)

    pdf.savefig(fig)
    plt.close(fig)
    print("  %.1f-%.1f nm: %d rows, %d bins, bin scatter %.5f median %.5f max,"
          " plain vs hierarchical rms %.5f"
          % (lo, hi, block.shape[0], len(bins), np.nanmedian(spread),
             np.nanmax(spread), np.sqrt(np.nanmean((plain - hier) ** 2))))


if __name__ == "__main__":
    main()
