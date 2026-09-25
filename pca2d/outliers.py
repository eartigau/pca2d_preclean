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

ONE SAMPLE, AND THEN THE WHOLE COLUMN. Dividing the residual by that local
sigma turns it into a z, and a z is comparable from one exposure to the next:
so the same residual that clips a sample of one exposure also says, column by
column in the OBSERVER's frame, how the exposures agree there. A column where
too many of them are beyond the clip AND whose survivors still carry more
variance than noise is a column the model cannot describe in any exposure: a
detector defect, a telluric line the correction does not reach, a sky residual
that stands still. It goes for every exposure rather than for the few that
happened to be caught (asked for on 2026-09-25).

    frac    the fraction of exposures clipped in that column, > correct.column_frac
    chi2    the reduced chi2 of the SURVIVORS, > correct.column_chi2

Both, not either. What each threshold means, measured rather than assumed
(scratchpad/null_chi2.py, 4096 columns of 321 rows of pure Gaussian noise):

    the survivors' reduced chi2 is 0.973, not 1, because the clip takes the
    tails with it (expected_chi2 below is the closed form). Its median over
    columns was 0.979 and its highest 1.271, so 1.5 is well clear of noise.
    The clipped fraction reached 2.2% at most, so the 10% of column_frac
    cannot be met by noise alone: it is what keeps a column with a handful of
    genuinely bad exposures, whose survivors are clean, from being thrown away
    whole. Together they caught 39 of 40 columns given twice the noise of
    their neighbours, 40 of 40 at three times, and not one of the 4056 good
    ones. A column at 1.5 times the noise is NOT caught: its clipped fraction
    is 4.8%, under the 10%.
"""

from __future__ import annotations

import math
import os
import warnings

import numpy as np

from .grids import pixel_shift


def expected_chi2(nsig):
    """E[z^2] of what survives a clip at `nsig`, for unit Gaussian noise.

    The number the column test's threshold is measured against: clipping
    removes the tails, so the survivors of a 3 sigma clip carry 0.973 and not
    1. A threshold below this would throw away pure noise.
    """
    a = float(nsig)
    phi = math.exp(-0.5 * a * a) / math.sqrt(2.0 * math.pi)
    inside = math.erf(a / math.sqrt(2.0))
    return (inside - 2.0 * a * phi) / inside


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


def residual_outliers(cube, fit, nsig, window, block=40000, chunk=16,
                      column_frac=None, column_chi2=None, column_min_rows=10):
    """(file -> (2, grid) mask of what to NaN, a report on the columns).

    The residual is panel 5's: the cube less the star block carried into each
    row's frame, the star-frame means if the fit has them, the observer block
    and the parity offset, NaN where the fit gave no weight. It is clipped row
    by row with running_stats over `window` samples. The cube is read `block`
    columns at a time, each with a margin covering the largest shift, the
    Lanczos kernel and half the window, so every column comes out as from the
    whole cube; `block=None` reads it whole. Row 0 of each mask is the even
    orders, row 1 the odd ones, as reconstruct.fit_weights_mask.

    With `column_frac` and `column_chi2`, every exposure also loses the columns
    the exposures disagree on as a body: more than `column_frac` of them
    clipped there and the survivors' reduced chi2 still above `column_chi2`,
    over at least `column_min_rows` exposures of that parity. The report holds
    that mask, the two statistics per column and what they cost.
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
    # per order parity and per OBSERVER column, over every exposure: how many
    # were measured there, how many the clip took, and the z^2 of the rest
    wanted = bool(column_frac) and bool(column_chi2)
    seen = np.zeros((2, m), dtype=np.int32)
    taken = np.zeros((2, m), dtype=np.int32)
    kept = np.zeros((2, m), dtype=np.int32)
    sum_z2 = np.zeros((2, m), dtype=float)
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
            # the local median and sigma once, for the clip and for the z the
            # column test needs: dividing by that sigma is what makes one
            # exposure's residual comparable with another's
            med, sig = running_stats(resid, window, origin=a0)
            with np.errstate(invalid="ignore"):
                z = (resid - med) / sig
            measured = np.isfinite(z)
            cut = measured & (np.abs(z) > float(nsig))
            for i, r in enumerate(range(start, stop)):
                out[names[r]][int(parity[r])][c0:c1] = cut[i, inner]
            if wanted:
                rows_parity = parity[start:stop]
                for p in (0, 1):
                    which = np.flatnonzero(rows_parity == p)
                    if not which.size:
                        continue
                    ok = measured[which][:, inner]
                    gone = cut[which][:, inner]
                    survivor = np.where(ok & ~gone, z[which][:, inner], np.nan)
                    seen[p, c0:c1] += ok.sum(axis=0)
                    taken[p, c0:c1] += gone.sum(axis=0)
                    kept[p, c0:c1] += np.isfinite(survivor).sum(axis=0)
                    sum_z2[p, c0:c1] += np.nansum(survivor ** 2, axis=0)
        del data, w

    frac = np.divide(taken, seen, out=np.zeros((2, m)), where=seen > 0)
    chi2 = np.divide(sum_z2, kept, out=np.zeros((2, m)), where=kept > 0)
    columns = np.zeros((2, m), dtype=bool)
    if wanted:
        columns = ((seen >= int(column_min_rows)) & (frac > float(column_frac))
                   & (chi2 > float(column_chi2)))
        added = 0
        for name in names:
            added += int((columns & ~out[name]).sum())
            out[name] |= columns
        report = dict(columns=columns, frac=frac, chi2=chi2, seen=seen,
                      n_columns=int(columns.sum()), added=added,
                      thresholds=(float(column_frac), float(column_chi2)))
    else:
        report = dict(columns=columns, frac=frac, chi2=chi2, seen=seen,
                      n_columns=0, added=0, thresholds=None)
    return out, report
