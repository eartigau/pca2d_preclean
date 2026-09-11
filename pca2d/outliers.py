"""Samples the model leaves unexplained, found and set to NaN.

What panel 5 of the report's sequence figure shows, the data less everything
the fit models (the star block, the observer block and the parity offset), is
noise where the model is right and something else where it is not: a cosmic
ray the MAD cut let through, a hot pixel, a sky residual. Against the noise
around it, measured where it is, a sample far outside that noise is not a
measurement of the star, and a corrected file is better without it.

    running_stats      the local median and robust sigma (1.4826 MAD) of each
                       row, in boxes as wide as the high pass (highpass.window)
    clip               NaN wherever a value is more than nsig of them away
    residual_outliers  the same on the fit's own residual, per exposure, as
                       the (even, odd) grid mask the correct stage applies

The corrected files take it through reconstruct's NaN mask, beside the samples
the fit gave no weight, when correct.nsig_cut is set in config.yaml; null
leaves every sample in.
"""

from __future__ import annotations

import os
import warnings

import numpy as np

from .grids import pixel_shift


def running_stats(values, window, step=None, min_valid=0.5, origin=0):
    """(local median, robust sigma) of `values` along its last axis.

    The median and 1.4826 times the median absolute deviation about it, in
    boxes of `window` samples centred every `step` (window // 4 by default) and
    interpolated linearly in between. NaN are ignored; a box with fewer than
    `min_valid` of its samples finite gives no estimate, and its neighbours'
    are carried across it. Robust on purpose: the statistic has to describe
    the noise around an outlier, not be pulled by it. `origin` is the index of
    the first sample in a longer array: the box centres sit at multiples of
    `step` of that array, so a slice of it gets the same boxes as the whole.
    """
    x = np.atleast_2d(np.asarray(values, dtype=float))
    n, m = x.shape
    window = int(window) | 1
    half = window // 2
    step = max(1, int(step) if step else window // 4)
    centres = np.unique(np.r_[0, np.arange((-int(origin)) % step, m, step), m - 1])
    med = np.full((n, m), np.nan)
    sig = np.full((n, m), np.nan)
    index = np.arange(m)
    pad = np.full(half, np.nan)
    for r in range(n):
        boxes = np.lib.stride_tricks.sliding_window_view(
            np.concatenate([pad, x[r], pad]), window)[centres]
        count = np.isfinite(boxes).sum(axis=1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            mc = np.nanmedian(boxes, axis=1)
            sc = 1.4826 * np.nanmedian(np.abs(boxes - mc[:, None]), axis=1)
        ok = (count >= min_valid * window) & np.isfinite(sc) & (sc > 0)
        if ok.sum() < 2:
            continue
        med[r] = np.interp(index, centres[ok], mc[ok])
        sig[r] = np.interp(index, centres[ok], sc[ok])
    shape = np.shape(values)
    return med.reshape(shape), sig.reshape(shape)


def clip(values, window, nsig, **kwargs):
    """(values with NaN where they sit beyond nsig robust sigma of the local
    median, the mask of what was set to NaN)."""
    values = np.asarray(values, dtype=float)
    med, sig = running_stats(values, window, **kwargs)
    with np.errstate(invalid="ignore"):
        cut = np.abs(values - med) > float(nsig) * sig
    cut &= np.isfinite(values) & np.isfinite(sig)
    return np.where(cut, np.nan, values), cut


def residual_outliers(cube, fit, nsig, window, block=40000, chunk=16):
    """The samples of each exposure its residual puts beyond nsig: file -> (2, grid).

    The residual is panel 5's: the cube less the star block carried into each
    row's frame, the star-frame means if the fit has them, the observer block
    and the parity offset, NaN where the fit gave no weight. It is clipped row
    by row with running_stats over `window` samples. The cube is read `block`
    columns at a time, each with a margin covering the largest shift, the
    Lanczos kernel and half the window, so every column comes out as from the
    whole cube; `block=None` reads it whole. Row 0 of each mask is the even
    orders, row 1 the odd ones, as reconstruct.fit_weights_mask.
    """
    from .twoframe import (LanczosShifter, carried_means, cube_grid, fit_means,
                           fit_templates, load_cube, row_parity, star_model)
    m = np.asarray(cube_grid(cube)).size
    _, _, w1, meta = load_cube(cube, dtype=np.float32, columns=np.arange(0, 1))
    n = w1.shape[0]
    if len(fit["berv"]) != n:
        raise SystemExit("the fit has %d rows and the cube %d: it was not made on"
                         " this cube" % (len(fit["berv"]), n))
    parity = row_parity(meta, n) % 2
    names = [os.path.basename(str(v)) for v in meta["filename"]]
    delta = -pixel_shift(np.asarray(fit["berv"], dtype=float), float(fit["dv"]))
    reach = int(np.ceil(np.abs(delta).max())) + 2
    window = int(window) | 1
    margin = reach + 16 + window
    P, a = np.asarray(fit["P"], dtype=float), np.asarray(fit["a"], dtype=float)
    Q, b = np.asarray(fit["Q"], dtype=float), np.asarray(fit["b"], dtype=float)
    means, group = fit_means(fit, meta, n, m)
    T, tgroup = fit_templates(fit, meta, n, m)
    out = {name: np.zeros((2, m), dtype=bool) for name in names}
    step = m if block is None else int(block)
    for c0 in range(0, m, step):
        c1 = min(c0 + step, m)
        a0, b0 = max(0, c0 - margin), min(m, c1 + margin)
        _, data, w, _ = load_cube(cube, dtype=np.float32, columns=np.arange(a0, b0))
        shifter = LanczosShifter(b0 - a0, a=8, max_shift=reach)
        Tf = shifter.prepare(T[:, a0:b0]) if np.any(T[:, a0:b0]) else None
        inner = slice(c0 - a0, c1 - a0)
        for start in range(0, n, chunk):
            stop = min(start + chunk, n)
            rows = slice(start, stop)
            model = star_model(P[:, a0:b0], a[rows], shifter, delta[rows],
                               stop - start, b0 - a0)
            if Tf is not None:
                model += carried_means(Tf, tgroup, shifter, delta, start, stop)
            model += b[rows] @ Q[:, a0:b0]
            model += means[:, a0:b0][group[rows]]
            resid = np.where(w[rows] > 0, data[rows] - model, np.nan)
            _, cut = clip(resid, window, nsig, origin=a0)
            for i, r in enumerate(range(start, stop)):
                out[names[r]][int(parity[r])][c0:c1] = cut[i, inner]
        del data, w
    return out
