#!/usr/bin/env python
"""Rebuild or correct an exposure from a two-frame fit, on its own wavelengths.

The `correct` stage of `pca2d-preclean`. Reads `twoframe_components.fits` and
puts the model back where the data came from: order by order, on the t.fits
file's native wavelengths, with each order taking the means of *its own
parity*.

    python -m pca2d.reconstruct --fits outputs/TOI2120/1-3v/twoframe_components.fits \
        --file data/TOI2120/2811170t.fits --out model.fits --plot model.pdf
    python -m pca2d.reconstruct --correct --all --cube cache/cube_tfits_<key> \
        --fits outputs/TOI2120/1-3v/twoframe_components.fits \
        --source-dir data/TOI2120 --corrected-dir outputs/TOI2120/1-3v/corrected

WHAT THE NUMBERS ARE. The cube was high-passed in the log, so every quantity
here is `ln f - savgol(ln f)`: zero in the continuum, negative in a line. The
Savitzky-Golay continuum was never stored, by anyone, so this cannot be turned
back into a flux. It is the right quantity to subtract from a fresh spectrum
that has been through the same high-pass, and the wrong one to multiply
anything by.

THE MODEL, per exposure n and order parity p:

    y_n = S_n (T_p + P^T a_n)  +  mu[p]  +  Q^T b_n

`S_n` carries the star block from the stellar rest frame into the observer
frame, a translation of -pixel_shift(BERV) samples because the grid is
log-uniform. `T_p` is the star-frame mean: zero by default, a median template
with the fit's --template, one per parity with --mean iterate. The observer
block and its mean `mu[p]` do not move.

WHAT A CORRECTED FILE HOLDS is panel 3 of the sequence figure, exactly: the
flux with the observer block and its per-parity mean divided out, and NaN
wherever the fit gave the sample no weight (fit_weights_mask). The star block
stays, T_p included. tests/test_panel3_is_the_correction.py holds them together.

THE PARITY WRAPPING is the part that is easy to get wrong. A t.fits row set is
split into even and odd orders, and `mu_earth` was fitted once per parity, so
order k must take `mean_even` if k is even and `mean_odd` if it is odd. Using
one mean for both puts a static even-minus-odd offset into every reconstructed
order, which is exactly the artefact the two-frame fit was fixed to avoid. The
star and Earth *components* have no parity: they are single vectors over the
whole grid, and only the means are per parity.
"""

from __future__ import annotations

import argparse
import os

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import CubicSpline

from .grids import pixel_shift
from .logger import log
from . import tfits as sptf
from . import twoframe as _bcd
from .progress import bar as _bar

#: the same treatment twoframe.load_cube gives the cube, so a refitted
#: coefficient is comparable with the fit's own rather than to its own scale
LN_CLIP_LOW = -0.5
RAMP_ZERO = 0.5


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fits", required=True,
                   help="the twoframe_components.fits the fit stage wrote")
    p.add_argument("--file", default=None,
                   help="the t.fits to act on; matched to COEFFS by basename")
    p.add_argument("--all", action="store_true",
                   help="every exposure in the fit, not just --file")
    p.add_argument("--out", default=None, help="write the model to this FITS")
    p.add_argument("--plot", default=None, help="write a check figure here")
    p.add_argument("--kernel-halfwidth", type=int, default=8,
                   help="Lanczos support, must match the fit (default 8)")

    p.add_argument("--correct", action="store_true",
                   help="divide the components out of the flux and write a new"
                        " t.fits named <stem>t_M-N.fits, M star components and N"
                        " Earth components removed. The template and the"
                        " observer means stay in: the filename names components,"
                        " so components are what comes out")
    p.add_argument("--n-star", type=int, default=None,
                   help="M, star components to remove (default: all in the fit)")
    p.add_argument("--n-earth", type=int, default=None,
                   help="N, Earth components to remove (default: all in the fit)")
    p.add_argument("--max-sky-ratio", type=float, default=None,
                   help="set to NaN any sample where the subtracted sky"
                        " emission exceeded this many times the stellar flux."
                        " Must match quality.max_sky_ratio of the fit, or the"
                        " corrected files carry samples the model never saw")
    p.add_argument("--corrected-dir", default="corrected",
                   help="where the t_M-N.fits files go")
    p.add_argument("--source-dir", default=None,
                   help="where to find the input t.fits when using --all or"
                        " --by-night")
    p.add_argument("--refit", action="store_true",
                   help="fit each exposure's OWN coefficients against the saved"
                        " basis instead of reusing its night's. The basis is the"
                        " campaign's and stays fixed; only the amplitudes are"
                        " re-estimated. This is what --by-night needs to be"
                        " honest: the sky changes within a night, so the"
                        " airglow amplitude of a 3-exposure visit is not one"
                        " number. Needs --config, the same YAML the cube was"
                        " built from, because the exposure has to be resampled"
                        " and high-passed exactly as the fit saw it")
    p.add_argument("--config", default=None,
                   help="the cube's YAML, required by --refit")
    p.add_argument("--by-night", action="store_true",
                   help="correct EVERY t.fits in --source-dir, not only the one"
                        " file per night that carries coefficients. With a"
                        " nightly-stacked cube the fit sees one row per night"
                        " but the night was built from several exposures; this"
                        " applies that night's coefficients to each of them,"
                        " using each exposure's own BERV to place the star"
                        " basis. Nights absent from the fit are skipped")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--cube", default=None,
                   help="the cube the fit was made from. With it, every sample"
                        " the fit gave no weight to is NaN in the corrected"
                        " file, exactly the samples panel 3 of the sequence"
                        " figure hides")
    args = p.parse_args(argv)
    if not args.file and not args.all and not args.by_night:
        p.error("give --file, or --all for every exposure in the fit, or"
                " --by-night for every exposure in --source-dir")
    return args


def load_model(path):
    """The basis, the means and the coefficient table, as plain arrays."""
    with fits.open(path) as h:
        head = h[0].header
        basis = h["BASIS"].data
        coeffs = Table(h["COEFFS"].data)
        n_star, n_earth = int(head["NSTAR"]), int(head["NEARTH"])
        model = {
            "grid": np.asarray(basis["wavelength"], dtype=np.float64),
            "template": np.asarray(basis["template"], dtype=np.float64),
            "P": np.array([basis["star_pc%d" % (k + 1)] for k in range(n_star)],
                          dtype=np.float64),
            "Q": np.array([basis["earth_pc%d" % (j + 1)] for j in range(n_earth)],
                          dtype=np.float64),
            "dv": float(head["DV"]),
            "n_star": n_star, "n_earth": n_earth,
            "content": head.get("CONTENT", "?"),
            "coeffs": coeffs,
        }
        names = [c for c in basis.columns.names if c.startswith("mean_")]
        # ordered even, odd -- the same order parity % 2 indexes
        model["means"] = {nm.split("_", 1)[1]: np.asarray(basis[nm], dtype=np.float64)
                          for nm in names}
        # the star-frame mean of each parity, from a --mean iterate fit; None
        # from an older one, which has the single `template` instead
        names = [c for c in basis.columns.names if c.startswith("template_")]
        model["templates"] = ({nm.split("_", 1)[1]:
                               np.asarray(basis[nm], dtype=np.float64)
                               for nm in names} or None)
    return model


def template_for(model, parity):
    """The star-frame mean for an order of this parity (0 even, 1 odd)."""
    per = model.get("templates")
    if per:
        names = list(per.keys())
        return per[names[parity % len(names)]]
    return model["template"]


def exposure_row(coeffs, filename):
    """The COEFFS row for this file, matched on basename."""
    want = os.path.basename(filename)
    names = [os.path.basename(str(v)) for v in coeffs["filename"]]
    if want not in names:
        stem = want.replace("t.fits", "")
        near = [n for n in names if stem[:20] in n][:3]
        raise SystemExit("%s is not in the fit%s"
                         % (want, ("; nearest: %s" % ", ".join(near)) if near else ""))
    row = coeffs[names.index(want)]
    if bool(row["rejected"]):
        log("WARNING: this exposure was rejected by the MAD cut. Its"
              " coefficients are zeros and the model below is the template"
              " plus the mean, nothing else.")
    return row


def model_on_grid(model, row):
    """The full model on the magic grid, split into what moves and what does not.

    Returned separately because only the star half is shifted, and keeping them
    apart is what lets the caller add the right parity's means afterwards:
    the star-frame one (template_for) before the carry, the observer one after.
    """
    n_star, n_earth = model["n_star"], model["n_earth"]
    a = np.array([row["a%d" % (k + 1)] for k in range(n_star)])
    b = np.array([row["b%d" % (j + 1)] for j in range(n_earth)])
    star_rest = a @ model["P"]           # the parity's template goes on per order
    earth = b @ model["Q"]
    return star_rest, earth, a, b


def carry_star(star_rest, berv, dv, halfwidth=8):
    """S_n applied to the star-frame vector: translate by -pixel_shift(BERV).

    Exactly the operator the fit used, taken from twoframe_bcd so the two can
    never drift apart. The sign is the one NOTES 11.11 settles: a stellar
    feature sits at lambda_bary = lambda_obs D(BERV), D the relativistic
    Doppler factor, so star -> observer is minus grids.pixel_shift(BERV, dv).
    """
    delta = np.array([-float(pixel_shift(float(berv), dv))])
    shifter = _bcd.SHIFTERS["lanczos"](
        star_rest.size, a=halfwidth,
        max_shift=int(np.ceil(abs(delta[0]))) + 2)
    prepared = shifter.prepare(star_rest[None, :])
    return shifter.carry(prepared, delta)[0, 0]


def reconstruct(model, row, halfwidth=8, path=None):
    """The model resampled onto every order of the file's own wavelength grid.

    Returns (wave, recon, parity): each (n_orders, n_pixels), parity per order.
    Samples the basis does not cover come back NaN rather than zero: zero is a
    meaningful value in a high-passed spectrum, "continuum", and handing back a
    flat continuum where the model has nothing to say would be a lie.
    """
    payload = sptf.read_tfits(path)
    wave = payload["wave"]
    n_orders, n_pixels = wave.shape
    grid = model["grid"]

    star_rest, earth, _, _ = model_on_grid(model, row)
    star_obs = [carry_star(template_for(model, p) + star_rest, row["berv"],
                           model["dv"], halfwidth) for p in (0, 1)]

    names = list(model["means"].keys())          # even, odd
    recon = np.full((n_orders, n_pixels), np.nan)
    parity = np.arange(n_orders) % 2
    for order in range(n_orders):
        mean = model["means"][names[parity[order]] if len(names) > 1 else names[0]]
        total = star_obs[parity[order]] + earth + mean
        # the basis has no support where the mean is exactly zero, and that is
        # per parity, so the mask has to be rebuilt for each order's parity
        support = mean != 0.0
        w = wave[order]
        ok = np.isfinite(w)
        if not ok.any():
            continue
        lo = max(np.searchsorted(grid, np.nanmin(w)) - 4, 0)
        hi = min(np.searchsorted(grid, np.nanmax(w)) + 4, grid.size)
        if hi - lo < 4:
            continue
        seg_g, seg_v, seg_s = grid[lo:hi], total[lo:hi], support[lo:hi]
        spline = CubicSpline(seg_g, np.where(seg_s, seg_v, 0.0), extrapolate=False)
        values = spline(w)
        # drop any sample whose neighbourhood on the grid was not supported
        nearest = np.clip(np.searchsorted(seg_g, w) - 1, 0, seg_s.size - 2)
        covered = seg_s[nearest] & seg_s[nearest + 1]
        recon[order] = np.where(ok & covered, values, np.nan)
    return wave, recon, parity


def correction_on_grid(model, row, n_star=None, n_earth=None, halfwidth=8):
    """The part of the model a correction removes, on the magic grid.

    The components only: M star and N Earth, the pair `t_M-N.fits` is named
    after. Never the template, the star's mean spectrum, whose removal would
    not be a correction but a different product. The observer-frame mean the
    fit took out of each parity is not here either, only because which one
    applies depends on the order: order_correction adds it, order by order,
    together with the observer components.

    M and N may be smaller than the fit's, which is the point of naming the file
    after them: t_5-5, t_0-5 and t_5-0 are three different corrections of the
    same exposure and the filename is what tells them apart.
    """
    k = model["n_star"] if n_star is None else int(n_star)
    j = model["n_earth"] if n_earth is None else int(n_earth)
    if not 0 <= k <= model["n_star"] or not 0 <= j <= model["n_earth"]:
        raise SystemExit("the fit has %d star and %d Earth components; asked for %d and %d"
                         % (model["n_star"], model["n_earth"], k, j))
    star_rest = np.zeros(model["grid"].size)
    if k:
        a = np.array([row["a%d" % (i + 1)] for i in range(k)])
        star_rest = a @ model["P"][:k]
    earth = np.zeros(model["grid"].size)
    if j:
        b = np.array([row["b%d" % (i + 1)] for i in range(j)])
        earth = b @ model["Q"][:j]
    star_obs = (carry_star(star_rest, row["berv"], model["dv"], halfwidth)
                if k else star_rest)
    return star_obs + earth, k, j


def _field(row, name):
    """row[name] as a float, None when the row has no such field.

    Not `name in row`: on a FITS or astropy table row that asks whether the
    name is one of the row's VALUES, which it never is, and the PCASTR_V card
    it used to guard was never written to a corrected file.
    """
    try:
        return float(row[name])
    except (KeyError, ValueError, IndexError, TypeError):
        return None


def coefficient_cards(model, row, k, j):
    """The fit's coefficients for one exposure, as (keyword, value, comment).

    Every component the fit HAS, not only the ones divided out. With
    correct.n_star = 0 the star coefficients are exactly what stays in the
    corrected flux, and they are the thing a velocity is worth correlating
    against afterwards; nothing downstream carries them otherwise, since an rdb
    knows only what LBL measured.

    Every keyword here is PCASTR or PCAOBS and then either a two-digit index or
    _N and _D, so one prefix names the frame and the suffix names the number.
    The PCA2xxxx keywords beside them are the ones that belong to the run as a
    whole rather than to one frame: PCA2BERV, PCA2NPIX, PCA2REJ.

    Two digits and not three. A FITS keyword is eight characters: PCASTR001
    would be nine, and astropy would write it as a HIERARCH card, which is not
    what a file whose whole point is to stay an ordinary t.fits should carry.
    Ninety-nine components is already twice what the validator calls sane.
    """
    cards = []
    # The velocity the fit took out of the observer block's reach. Reported and
    # NOT applied: the corrected flux is the delivered flux with the observer
    # components divided out, and this number is what would have been written
    # into it had the term not been fitted.
    vrad = _field(row, "vrad_fit")
    if vrad is not None and np.isfinite(vrad):
        cards.append(("PCASTR_V", vrad, "m/s, star shift the fit absorbed"))
    for prefix, letter, count, removed, frame in (
            ("PCASTR", "a", model["n_star"], k, "star"),
            ("PCAOBS", "b", model["n_earth"], j, "observer")):
        if count > 99:
            log("%d %s components is more than the %s01..99 keywords can name;"
                " the rest are not written to the header"
                % (count, frame, prefix), "warn")
        listed = min(int(count), 99)
        # Two counts, and they are different numbers: _N is how many amplitude
        # cards follow, so a reader loops without guessing, and _D is how many
        # of those components were actually divided out of the flux. With
        # correct.n_star = 0 the first is 2 and the second is 0.
        cards.append(("%s_N" % prefix, listed,
                      "%s-frame comps listed, %s01..%02d"
                      % (frame, prefix, listed)))
        cards.append(("%s_D" % prefix, int(removed),
                      "%s-frame comps divided out of the flux" % frame))
        for i in range(listed):
            key = "%s%02d" % (prefix, i + 1)
            value = float(row["%s%d" % (letter, i + 1)])
            cards.append((key, value, "%s-frame comp %d amplitude, %s"
                          % (frame, i + 1, "divided out" if i < removed
                             else "left in the flux")))
    return cards


def order_correction(model, correction, n_earth, order, wave):
    """What correct_file divides out of one order, on that order's own pixels.

    Returns (values, live): the model in ln f at each pixel, and where it may
    be applied, or None when the order misses the grid. `correction` is
    correction_on_grid's, the components. To it goes the mean the fit
    subtracted from this order's parity before solving for them, whenever any
    observer component is divided out: panel 3 of the sequence figure, "what a
    corrected file holds", is the data minus that mean minus the observer
    block, and a file that kept the mean kept a static observer-frame pattern
    the panel shows removed (2026-09-10, the stripes at 1267 nm). It belongs
    to the observer block, so a file that divides none of it out, t_k-0,
    keeps it.
    """
    names = list(model["means"].keys())
    grid = model["grid"]
    mean = model["means"][names[order % 2] if len(names) > 1 else names[0]]
    # the basis has no support where the mean is exactly zero, per parity
    support = mean != 0.0
    ok = np.isfinite(wave)
    if not ok.any():
        return None
    lo = max(np.searchsorted(grid, np.nanmin(wave)) - 4, 0)
    hi = min(np.searchsorted(grid, np.nanmax(wave)) + 4, grid.size)
    if hi - lo < 4:
        return None
    total = correction + mean if n_earth else correction
    seg_g, seg_s = grid[lo:hi], support[lo:hi]
    spline = CubicSpline(seg_g, np.where(seg_s, total[lo:hi], 0.0),
                         extrapolate=False)
    values = spline(wave)
    nearest = np.clip(np.searchsorted(seg_g, wave) - 1, 0, seg_s.size - 2)
    live = ok & seg_s[nearest] & seg_s[nearest + 1] & np.isfinite(values)
    return values, live


def fit_weights_mask(cube):
    """The grid samples the fit weighted, per exposure: file name -> (2, grid).

    Exactly what panel 3 of the sequence figure shows and nothing else: the
    weights the fit read from the cube (twoframe.load_cube, with the figure's
    own defaults) and, as the figure does, not a sample stranded between gaps
    (plotting.live_mask). Row 0 is the even orders, row 1 the odd ones.
    """
    from .plotting import live_mask
    _, _, w, meta = _bcd.load_cube(cube, dtype=np.float32)
    parity = _bcd.row_parity(meta, w.shape[0])
    names = [os.path.basename(str(v)) for v in meta["filename"]]
    out = {}
    for i, name in enumerate(names):
        out.setdefault(name, np.zeros((2, w.shape[1]), dtype=bool))
        out[name][int(parity[i]) % 2] = live_mask(w[i:i + 1])[0]
    return out


def weighted_on_pixels(alive, wave, grid):
    """Which of an order's pixels fall on grid samples the fit weighted.

    A pixel counts as weighted when both grid samples around it were, the rule
    order_correction uses for the basis support. Pixels outside the grid are
    outside the fit altogether and come back True, since there is no panel 3
    there to match.
    """
    inside = np.isfinite(wave) & (wave >= grid[0]) & (wave <= grid[-1])
    near = np.clip(np.searchsorted(grid, np.where(inside, wave, grid[0])) - 1,
                   0, grid.size - 2)
    return ~inside | (alive[near] & alive[near + 1])


def correct_file(model, row, path, outdir, n_star=None, n_earth=None,
                 halfwidth=8, overwrite=False, max_sky=None, alive=None):
    """Write a t.fits with the model divided out. Returns the new path.

    The model lives in `ln f - savgol(ln f)`, so removing it from the flux is a
    division, not a subtraction:

        f_corrected = f * exp(-model)

    with the model the components and, with any observer component, the parity
    mean the fit removed first (order_correction). That leaves the
    Savitzky-Golay continuum exactly where it was. That matters:
    the continuum was never part of the fit and never stored, so anything that
    claimed to reconstruct it would be inventing it.

    Every other extension is copied untouched, including Recon and the blaze, so
    the result is still a t.fits and still reducible by anything that reads one.
    Samples outside the basis support keep their original flux rather than being
    blanked: the correction has nothing to say there, which is not the same as
    the flux being unknown.
    """
    correction, k, j = correction_on_grid(model, row, n_star, n_earth, halfwidth)

    with sptf.robust_open(path) as hdulist:
        out = fits.HDUList([h.copy() for h in hdulist])
    # the same instrument table the reader uses; a SPIRou file must have its
    # AB extensions corrected, not the single-fibre A that also exists there
    e_flux, e_wave, _, _ = sptf.extensions_for(out)
    flux = np.asarray(out[e_flux].data, dtype=np.float64)
    wave = np.asarray(out[e_wave].data, dtype=np.float64)
    n_orders = flux.shape[0]
    touched = 0
    blanked = 0
    for order in range(n_orders):
        got = order_correction(model, correction, j, order, wave[order])
        if got is None:
            continue
        values, live = got
        flux[order] = np.where(live, flux[order] * np.exp(-values), flux[order])
        if alive is not None:
            # NaN wherever the fit gave the sample no weight: exactly the
            # samples panel 3 of the sequence figure hides. The model there was
            # constrained by nothing, and a flux corrected by it is not one.
            bad = ~weighted_on_pixels(alive[order % 2], wave[order], model["grid"])
            blanked += int((bad & np.isfinite(flux[order])).sum())
            flux[order][bad] = np.nan
        touched += int(live.sum())

    # Propagate the sky mask: samples the airglow drowned were excluded from
    # the fit, so leaving them in the corrected file would hand LBL data the
    # model never saw and could not correct. NaN, not zero: a zero is a
    # measurement of no flux and would be believed.
    masked = 0
    if max_sky:
        e_sky = sptf.sky_extension_for(out)
        if e_sky:
            sky = np.asarray(out[e_sky].data, dtype=np.float64)
            if sky.shape == flux.shape:
                with np.errstate(invalid="ignore", divide="ignore"):
                    ratio = np.abs(sky) / np.where(flux > 0, flux, np.nan)
                drown = np.isfinite(ratio) & (ratio > float(max_sky))
                flux[drown] = np.nan
                masked = int(drown.sum())

    out[e_flux].data = flux
    head = out[0].header
    if max_sky:
        head["PCA2SKYR"] = (float(max_sky), "sky/flux above this was set to NaN")
        head["PCA2SKYN"] = (masked, "samples removed as sky-dominated")
    head["PCA2REF"] = (True, "two-frame PCA correction applied")
    head["PCA2BERV"] = (float(row["berv"]), "km/s used to carry the star basis")
    head["PCA2NPIX"] = (touched, "samples corrected")
    head["PCA2MEAN"] = (bool(j), "observer-frame parity mean divided out")
    if alive is not None:
        head["PCA2WNAN"] = (blanked, "samples the fit gave no weight, set to NaN")
    head["PCA2REJ"] = (bool(row["rejected"]), "exposure was MAD-rejected")
    for key, value, comment in coefficient_cards(model, row, k, j):
        head[key] = (value, comment)
    head.add_history("two-frame PCA: %d star-frame + %d observer-frame"
                     " components%s divided out"
                     % (k, j, " and the parity mean" if j else ""))
    head.add_history("f_corrected = f * exp(-model), model in ln f - savgol(ln f)")

    stem = os.path.basename(path)
    stem = stem[:-len("t.fits")] if stem.endswith("t.fits") else stem.rsplit(".", 1)[0]
    new = os.path.join(outdir, "%st_%d-%d.fits" % (stem, k, j))
    os.makedirs(outdir, exist_ok=True)
    if os.path.exists(new) and not overwrite:
        raise SystemExit("%s exists; pass --overwrite" % new)
    out.writeto(new, overwrite=True)
    out.close()
    return new, touched, flux.size


def write_fits(path, wave, recon, parity, row, model, source):
    primary = fits.PrimaryHDU()
    primary.header["SOURCE"] = (os.path.basename(source), "t.fits rebuilt")
    primary.header["CONTENT"] = (model["content"], "ln f - savgol(ln f)")
    primary.header["BERV"] = (float(row["berv"]), "km/s, used for the star shift")
    primary.header["BJD"] = (float(row["bjd"]), "")
    primary.header["REJECTED"] = (bool(row["rejected"]), "MAD cut")
    hdus = [primary,
            fits.ImageHDU(recon, name="MODEL"),
            fits.ImageHDU(wave, name="WAVE"),
            fits.ImageHDU(parity.astype(np.int16), name="PARITY")]
    fits.HDUList(hdus).writeto(path, overwrite=True)
    log("wrote %s" % path)


def plot(path, wave, recon, row, source):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    finite = np.isfinite(recon)
    if not finite.any():
        log("nothing to plot: the model covers none of this file")
        return
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25})
    fig, axes = plt.subplots(2, 1, figsize=(12, 6))
    axes[0].plot(wave[finite], recon[finite], ".", ms=0.6, alpha=0.4, color="tab:blue")
    axes[0].set_xlabel("wavelength (nm)")
    axes[0].set_ylabel("model, ln f - savgol(ln f)")
    axes[0].set_title("%s   BERV %.3f km/s" % (os.path.basename(source), row["berv"]))
    # one order of each parity, side by side, to show the wrapping is right
    orders = [o for o in range(recon.shape[0]) if np.isfinite(recon[o]).sum() > 100]
    for o in orders[len(orders) // 2: len(orders) // 2 + 2]:
        axes[1].plot(wave[o], recon[o], lw=0.7,
                     label="order %d (%s)" % (o, "even" if o % 2 == 0 else "odd"))
    axes[1].set_xlabel("wavelength (nm)")
    axes[1].set_ylabel("model")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    log("wrote %s" % path)


def night_index(coeffs):
    """Map a night number to its COEFFS row.

    The night is floor(BJD - 0.5), the same convention the cube builder uses, so
    an exposure taken after midnight belongs to the night it started in.
    """
    nights = {}
    for i, bjd in enumerate(np.asarray(coeffs["bjd"], dtype=float)):
        nights[int(np.floor(bjd - 0.5))] = i
    return nights


def rows_by_night(model, args):
    """One (path, row) pair per t.fits in --source-dir, matched on its night.

    The row's coefficients are the night's; the BERV is the exposure's own,
    read from its header, because the star basis has to be placed where that
    exposure actually saw the star and exposures within a night differ by a
    fraction of a km/s.
    """
    import glob
    coeffs = model["coeffs"]
    index = night_index(coeffs)
    out, missing = [], 0
    for path in sorted(glob.glob(os.path.join(args.source_dir, "*t.fits"))):
        try:
            payload = sptf.read_tfits(path)
        except Exception:                                     # noqa: BLE001
            missing += 1
            continue
        meta = payload["meta"]
        night = int(np.floor(float(meta["bjd"]) - 0.5))
        if night not in index:
            missing += 1
            continue
        row = Table(coeffs)[index[night]]
        row["berv"] = float(meta["berv"])      # this exposure, not the night
        out.append((path, row))
    return out, missing


def refit_row(model, path, config, shifter, base_row):
    """This exposure's own coefficients, against the campaign's fixed basis.

    The alternative, which this replaces, was to hand every exposure of a night
    the night's coefficients and change only its BERV. That is defensible for
    the star, which does not vary over three hours, and wrong for everything
    the observer block describes: OH and O2 airglow rise and fall through a
    night by factors, and the water column moves with the airmass. A visit of
    four exposures was being corrected with one sky amplitude.

    The basis P and Q is NOT refitted, only the amplitudes a and b. Refitting
    the basis per exposure would be a different method, and a single spectrum
    does not constrain five components of anything.

    Returns a copy of `base_row` with fresh a1..ak, b1..bj and this exposure's
    BERV, or None if the file cannot be resampled onto the grid.
    """
    grid = model["grid"]
    payload = sptf.read_tfits(path)
    values, sigma, trans, good, n_written = sptf.resample_exposure(
        payload, grid, config)
    if not n_written or not good.any():
        return None

    # The SAME weights the fit used, not merely 1/sigma^2. twoframe.load_cube
    # also zeroes anything below ln_clip_low, ramps by the telluric
    # transmission and blanks the grid edges, and the transmission ramp in
    # particular matters here: it is near zero exactly where the atmosphere is
    # opaque, which is where the observer block does most of its work. Leaving
    # it out shifted the refitted amplitude by tens of per cent against the
    # fit's own, which is a method difference masquerading as a measurement.
    with np.errstate(divide="ignore", invalid="ignore"):
        w = np.where(good & np.isfinite(sigma) & (sigma > 0),
                     1.0 / np.where(sigma > 0, sigma, 1.0) ** 2, 0.0)
    w[~np.isfinite(values)] = 0.0
    w[values < LN_CLIP_LOW] = 0.0
    if trans is not None:
        w *= np.clip((trans - RAMP_ZERO) / (1.0 - RAMP_ZERO), 0.0, 1.0)
    w[:, :_bcd.EDGE] = 0.0
    w[:, -_bcd.EDGE:] = 0.0
    data = np.where(w > 0, values, 0.0)

    # take out exactly what the fit took out before it solved for coefficients:
    # the mean it stored per order parity, and the template carried into this
    # exposure's frame
    names = list(model["means"].keys())
    for parity in (0, 1):
        data[parity] -= model["means"][names[parity % len(names)]]
    berv = float(payload["meta"]["berv"])
    delta = np.full(2, -float(pixel_shift(berv, model["dv"])))
    for parity in (0, 1):
        template = template_for(model, parity)
        if template is not None and np.any(template):
            data[parity] -= _bcd.carry_template(template, shifter,
                                                delta[parity:parity + 1])[0]
    data[w <= 0] = 0.0

    tie = np.zeros(2, dtype=int)          # the two parities are one exposure
    # No velocity column here on purpose: this path re-solves ONE exposure
    # against a basis the fit already fixed, and the shift it would find has
    # nowhere to go. What the fit found is in the COEFFS row, and that is what
    # the header reports.
    a, b, _, _ = _bcd.joint_coeffs(data, w, model["P"], model["Q"], shifter,
                                   delta, exposure=tie)
    if base_row is None:
        return None
    # a plain dict, not a Table Row: a Row is a view into its table, so writing
    # coefficients into it would edit the fit's own coefficient table and the
    # next exposure would inherit them
    try:
        row = {name: base_row[name] for name in base_row.colnames}
    except AttributeError:
        row = dict(base_row)
    for k in range(model["n_star"]):
        row["a%d" % (k + 1)] = float(a[0, k])
    for j in range(model["n_earth"]):
        row["b%d" % (j + 1)] = float(b[0, j])
    row["berv"] = berv
    return row


def correct_many(model, args):
    """--correct over one file or over every exposure in the fit."""
    if (args.all or args.by_night) and not args.source_dir:
        raise SystemExit("--source-dir is required with --all or --by-night:"
                         " where the delivered t.fits are")
    if args.by_night:
        pairs, missing = rows_by_night(model, args)
        log("  %d exposures matched to a fitted night, %d skipped"
              % (len(pairs), missing))
    elif args.all:
        names = [os.path.basename(str(v)) for v in model["coeffs"]["filename"]]
        pairs = [(os.path.join(args.source_dir, n), None) for n in names]
    else:
        pairs = [(args.file, None)]
    covered = []
    shifter = config = None
    if args.refit:
        if not args.config:
            raise SystemExit("--refit needs --config, the YAML the cube was"
                             " built from: the exposure has to be resampled and"
                             " high-passed exactly as the fit saw it")
        from .config import load_config
        config = load_config(args.config)
        m = model["grid"].size
        top = int(np.ceil(np.abs(np.asarray(model["coeffs"]["berv"],
                                            dtype=float)).max()
                          / model["dv"])) + 2
        shifter = _bcd.LanczosShifter(m, a=8, max_shift=top)
        log("  refitting coefficients per exposure against the fixed basis")
    elif args.by_night:
        log("  WARNING: every exposure of a night gets that night's"
              " coefficients. The sky is not constant over a night; --refit"
              " solves for each exposure's own amplitudes")

    alive_by_file = None
    if getattr(args, "cube", None):
        log("  reading the fit's weights from %s: every sample the fit gave no"
            " weight to will be NaN, as panel 3 of the sequence figure hides it"
            % args.cube)
        alive_by_file = fit_weights_mask(args.cube)
    else:
        log("  no --cube: samples the fit gave no weight to keep their flux,"
            " which panel 3 of the sequence figure does not show", "warn")
    written = skipped = refitted = 0
    progress = _bar(pairs, desc="correcting", unit="file")
    for path, preset in progress:
        if not os.path.exists(path):
            log("  missing, skipped: %s" % path)
            skipped += 1
            continue
        row = preset if preset is not None else exposure_row(model["coeffs"], path)
        if args.refit:
            fresh = refit_row(model, path, config, shifter, row)
            if fresh is None:
                log("  could not resample, skipped: %s" % os.path.basename(path))
                skipped += 1
                continue
            row = fresh
            refitted += 1
        alive = (alive_by_file.get(os.path.basename(path))
                 if alive_by_file is not None else None)
        new, touched, total = correct_file(
            model, row, path, args.corrected_dir, args.n_star, args.n_earth,
            args.kernel_halfwidth, args.overwrite, max_sky=args.max_sky_ratio,
            alive=alive)
        written += 1
        # onto the bar, not onto its own line: three hundred of these scroll
        # the narration off the screen and say nothing a total cannot
        covered.append(100.0 * touched / max(total, 1))
        if hasattr(progress, "set_postfix_str"):
            progress.set_postfix_str("%s  %.0f%% covered"
                                     % (os.path.basename(new), covered[-1]))
    if covered:
        import numpy as _np
        log("  samples corrected per file: median %.1f%%, worst %.1f%%"
              % (float(_np.median(covered)), float(_np.min(covered))))
    log("wrote %d file%s to %s%s%s"
          % (written, "" if written == 1 else "s", args.corrected_dir,
             ", %d missing" % skipped if skipped else "",
             ", %d refitted individually" % refitted if refitted else ""))
    return written


def main(argv=None):
    args = parse_args(argv)
    model = load_model(args.fits)
    log("%s: %d star + %d earth components, dv = %.4f km/s, content = %s"
          % (os.path.basename(args.fits), model["n_star"], model["n_earth"],
             model["dv"], model["content"]))
    if args.correct:
        k = model["n_star"] if args.n_star is None else args.n_star
        j = model["n_earth"] if args.n_earth is None else args.n_earth
        log("correcting: removing %d star + %d Earth components -> t_%d-%d.fits"
              % (k, j, k, j))
        # not `return correct_many(...)`: this is a console_script entry
        # point, so the count of files written would become the exit status.
        # 316 corrected spectra exited 316.
        correct_many(model, args)
        return None
    row = exposure_row(model["coeffs"], args.file)
    log("  BJD %.5f   BERV %+.4f km/s   star shift %+.2f samples"
          % (row["bjd"], row["berv"], -row["berv"] / model["dv"]))
    wave, recon, parity = reconstruct(model, row, args.kernel_halfwidth, args.file)
    good = np.isfinite(recon)
    log("  rebuilt %d of %d samples (%.1f%%) over %d orders"
          % (good.sum(), recon.size, 100 * good.mean(), recon.shape[0]))
    log("  model rms %.5f in ln flux" % np.nanstd(recon[good]))
    if args.out:
        write_fits(args.out, wave, recon, parity, row, model, args.file)
    if args.plot:
        plot(args.plot, wave, recon, row, args.file)
    return None          # exit status, not a value; see above


if __name__ == "__main__":
    main()
