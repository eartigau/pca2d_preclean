"""Shared plotting conventions, so every figure says the same thing the same way.

One rule so far, and it earns its own module because getting it wrong is
invisible: **a missing sample is yellow.**

Matplotlib's default is to paint a NaN in the axes background colour, which on a
white figure is white and on a diverging colour map reads as "zero", the value
that means "continuum" in a high-passed spectrum. A gap and a flat continuum
then look identical, and a row of a cube that carries no data at all becomes a
pale stripe indistinguishable from a genuinely featureless spectrum. Yellow sits
outside every diverging map used here, so a gap can never be mistaken for a
measurement.
"""

from __future__ import annotations

import matplotlib
from .logger import log
import numpy as np

NAN_COLOUR = "#ffd400"      # banana


def nan_cmap(name="RdBu_r", under=None, over=None):
    """A colour map that paints NaN yellow.

    Returns a copy, never the registered map: mutating a global colour map
    changes every later figure in the process, including ones drawn by code
    that never asked for this convention.
    """
    cmap = matplotlib.colormaps[name].copy()
    cmap.set_bad(NAN_COLOUR)
    if under is not None:
        cmap.set_under(under)
    if over is not None:
        cmap.set_over(over)
    return cmap


#: how far a gap may be and still condemn a valid sample, in samples. The same
#: number the cube is built with; see config quality.isolated_window.
ISOLATED_WINDOW = 3


def live_mask(weights, window=ISOLATED_WINDOW):
    """Which samples of a shifted array are worth drawing.

    `weights > 0` is not enough. Carrying an array into another frame smears
    the mask by the shift, and what comes out has valid samples stranded inside
    regions that are otherwise gone; on screen they are speckles of colour in a
    band of yellow, and they read as measurements. This drops any sample with a
    gap within `window` on both sides, the same rule the cube is built with.
    """
    from .preprocess import drop_isolated

    return drop_isolated(np.asarray(weights) > 0, window)


def drop_empty_rows(*arrays, reference=None):
    """Remove rows that are entirely missing, rather than drawing them blank.

    A spectrum that contributes nothing to the window on screen should not
    occupy a line of the image. Left in, it draws a flat stripe that a reader
    reasonably takes for a measurement of something; the honest picture has one
    row per spectrum that is actually there. Returns the trimmed arrays plus the
    boolean mask of what was kept, so a caller can report how many went.

    `reference` chooses which array decides; by default the first one.
    """
    import numpy as np

    ref = arrays[0] if reference is None else reference
    keep = np.isfinite(np.asarray(ref)).any(axis=1)
    return [np.asarray(a)[keep] for a in arrays] + [keep]


def snap_to_order_centre(wavelength_nm, path):
    """Move a requested wavelength to the centre of the echelle order nearest it.

    A 2 nm window placed at an arbitrary wavelength can land on an order edge,
    where the blaze has fallen and two orders overlap, which is the least
    representative place in the spectrum to look. Snapping to an order centre
    puts the window where the throughput peaks and where exactly one order
    contributes.

    Returns (snapped_nm, order_index, requested_nm). If the file cannot be read
    the request is returned unchanged rather than raising: this is a convenience,
    not a requirement.
    """
    import numpy as np

    try:
        from .io import robust_open
        from .tfits import extensions_for
        with robust_open(path) as hdulist:
            # NOT a hardcoded "WaveA": on a SPIRou file that is the single-fibre
            # wavelength map, and the order centres would be read off the wrong
            # extension without any error
            _, e_wave, _, _ = extensions_for(hdulist)
            wave = np.asarray(hdulist[e_wave].data, dtype=float)
    except Exception:                                         # noqa: BLE001
        return float(wavelength_nm), None, float(wavelength_nm)

    centres = np.array([np.nanmedian(w[np.isfinite(w) & (w > 0)])
                        if np.isfinite(w).any() else np.nan for w in wave])
    good = np.isfinite(centres)
    if not good.any():
        return float(wavelength_nm), None, float(wavelength_nm)
    idx = np.where(good)[0][np.argmin(np.abs(centres[good] - wavelength_nm))]
    return float(centres[idx]), int(idx), float(wavelength_nm)


def order_bounds(path):
    """(lo, hi) in nm for every order of a file's wavelength extension.

    One row per order, in the order the extension stores them, which is the
    same counting the cube uses: its builder assigns `parity = order % 2`, so
    the index returned here maps straight onto a cube row's parity. An order
    with no finite wavelength comes back as NaN rather than raising.
    """
    import numpy as np

    from .io import robust_open
    from .tfits import extensions_for

    with robust_open(path) as hdulist:
        # NOT a hardcoded extension name, for the reason given above
        _, e_wave, _, _ = extensions_for(hdulist)
        wave = np.asarray(hdulist[e_wave].data, dtype=float)
    wave = np.where(np.isfinite(wave) & (wave > 0), wave, np.nan)
    lo = np.full(wave.shape[0], np.nan)
    hi = np.full(wave.shape[0], np.nan)
    ok = np.isfinite(wave).any(axis=1)
    if ok.any():
        lo[ok] = np.nanmin(wave[ok], axis=1)
        hi[ok] = np.nanmax(wave[ok], axis=1)
    return lo, hi


def window_parity(wavelength_nm, path):
    """Which order parity owns a wavelength where two orders overlap.

    Consecutive echelle orders share their ends, so a window that falls in an
    overlap is measured twice: once near the middle of one order, where the
    blaze peaks and the extraction is at its best, and once at the very edge of
    its neighbour, where the throughput has fallen off and the same wavelengths
    are far noisier. Since the two parities are separate cube rows, both sets
    are drawn, and the figure then shows a measurement stacked above a poor
    copy of itself, which is worse than showing one.

    The order kept is the one whose middle the wavelength sits closest to,
    measured in units of that order's own half-width so that orders of
    different widths compare: 0 is dead centre, 1 is an edge.

    Returns (parity, offset, n_covering). `n_covering` is how many orders reach
    this wavelength at all, so 1 means there was nothing to choose and the
    caller can say so. (None, nan, 0) if the file cannot be read, in which case
    the caller should fall back to its coverage rule rather than fail: this
    sharpens a figure, it does not gate one.
    """
    try:
        lo, hi = order_bounds(path)
    except Exception:                                         # noqa: BLE001
        return None, float("nan"), 0
    return parity_from_bounds(wavelength_nm, lo, hi)


def parity_from_bounds(wavelength_nm, lo, hi):
    """The arithmetic of `window_parity`, with the file already read.

    Separated so it can be checked without a spectrum on disk: the order index
    of the winner is what decides the parity, and getting that off by one puts
    every overlapping window on the wrong side.
    """
    import numpy as np

    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    half = 0.5 * (hi - lo)
    with np.errstate(invalid="ignore"):
        offset = np.abs(wavelength_nm - 0.5 * (lo + hi)) / np.where(half > 0,
                                                                   half, np.nan)
        covers = np.isfinite(offset) & (offset <= 1.0)
    if not covers.any():
        # outside every order, or an unreadable solution: nothing to arbitrate
        return None, float("nan"), 0
    idx = np.where(covers)[0][np.argmin(offset[covers])]
    return int(idx % 2), float(offset[idx]), int(covers.sum())


def raw_log_flux_block(cube_dir, source_dir, names, parities, grid, cols):
    """ln(flux) BEFORE the high-pass, every cube row, on a slice of the grid.

    The panels a fit produces all live in `ln f - savgol(ln f)`, which has no
    flux units and no continuum. To show what the correction does to the
    spectrum as it sits in the file, the same exposures are resampled again
    with the high-pass switched off, which leaves ln(flux) on the same grid.

    Each exposure is opened once and both of its order parities taken from that
    read, and only the columns asked for are kept, so a page of windows costs
    one pass over the files rather than one pass per window.
    """
    import os

    import numpy as np

    from . import tfits as sptf
    from .config import load_config, spectra_dir

    config = load_config(os.path.join(cube_dir, "cube_config.yaml"))
    config["highpass"] = dict(config["highpass"], method="none")
    directory = source_dir or spectra_dir(config)
    out = np.full((len(names), cols.size), np.nan)
    by_name = {}
    for i, name in enumerate(names):
        by_name.setdefault(os.path.basename(str(name)), []).append(i)
    read = 0
    for name, rows in by_name.items():
        try:
            payload = sptf.read_tfits(os.path.join(directory, name))
            values, _sigma, _trans, good, _n = sptf.resample_exposure(
                payload, grid, config)
        except Exception:                                     # noqa: BLE001
            continue
        for i in rows:
            p = int(parities[i])
            out[i] = np.where(good[p][cols], values[p][cols], np.nan)
        read += 1
    log("  re-read %d of %d exposures without the high-pass, %d columns"
          % (read, len(by_name), cols.size))
    return out


def sample_source_file(cube_path, meta=None, source_dir=None):
    """One readable t.fits behind a cube, for anything that needs its geometry.

    The order boundaries are a property of the instrument, not of the night, so
    any exposure of the set answers the question and the first one that opens
    will do. Returns None if the config names no directory, the metadata holds
    no filenames, or none of them are on this machine, so that a figure drawn
    from a cube alone still draws.
    """
    import os

    directory = source_dir or source_directory(cube_path)
    if not directory or meta is None:
        return None
    try:
        names = [os.path.basename(str(v)) for v in meta["filename"]]
    except Exception:                                         # noqa: BLE001
        return None
    for name in names[:20]:
        path = os.path.join(directory, name)
        if os.path.exists(path):
            return path
    return None


# ---------------------------------------------------------------------------
# Correlations against ancillary quantities
# ---------------------------------------------------------------------------
# A component carries no label. That a1 varies on 83 days is a result; that it
# varies with the precipitable water column would say plainly that a "stellar"
# component is describing the atmosphere. This is the diagnostic that catches
# that, which is why the fit writes it itself rather than leaving it to a
# separate script somebody may forget to run.

#: cube metadata column -> label
ANCILLARY = [
    ("berv", "BERV"),
    ("airmass", "airmass"),
    ("snr_band", "band SNR"),
    ("seeing", "seeing"),
    ("humidity", "humidity"),
    ("exptime", "exp. time"),
    ("sun_elevation", "sun elev."),
    ("median_flux", "median flux"),
    ("mjdmid", "time"),
]
#: APERO's telluric-preclean water exponent, in the FluxA header, not the cube
WATER_KEY = "TLPEH2O"


def water_column(filenames, source_dir):
    """TLPEH2O per exposure, read from the t.fits themselves.

    NaN where the file is missing or the key absent, so a partial read shows as
    gaps in one row rather than silently dropping exposures and shifting every
    other correlation.
    """
    import os
    from astropy.io import fits

    out = np.full(len(filenames), np.nan)
    if not source_dir:
        return out
    for i, name in enumerate(filenames):
        path = os.path.join(source_dir, os.path.basename(str(name)))
        if not os.path.exists(path):
            continue
        try:
            with fits.open(path) as hdulist:
                for hdu in hdulist:
                    value = hdu.header.get(WATER_KEY)
                    if value is not None:
                        out[i] = float(value)
                        break
        except Exception:                                     # noqa: BLE001
            continue
    return out


def source_directory(cube_path):
    """Where the t.fits behind a cube live, per the config the build copied in.

    `pca2refs-cube` drops the config it used into the cache directory, so the
    fit can recover the source folder without being told twice. Returns None if
    the file is absent, in which case the water row is simply omitted rather
    than the whole matrix failing.
    """
    import os

    path = os.path.join(cube_path, "cube_config.yaml")
    if not os.path.isdir(cube_path) or not os.path.exists(path):
        return None
    try:
        import yaml
        with open(path) as fh:
            cfg = yaml.safe_load(fh) or {}
        return (cfg.get("input") or {}).get("directory")
    except Exception:                                         # noqa: BLE001
        return None


def ancillary_table(meta, filenames, source_dir=None):
    """(labels, values) for every ancillary quantity that varies."""
    import os

    labels, values = [], []
    if meta is not None:
        by_name = {}
        for row in meta:
            by_name.setdefault(os.path.basename(str(row["filename"])), row)
        for key, label in ANCILLARY:
            if key not in getattr(meta, "colnames", []):
                continue
            col = np.array([float(by_name[n][key]) if n in by_name else np.nan
                            for n in filenames])
            if np.isfinite(col).sum() > 10 and np.nanstd(col) > 0:
                labels.append(label)
                values.append(col)
    water = water_column(filenames, source_dir)
    if np.isfinite(water).sum() > 10:
        labels.append("water (%s)" % WATER_KEY)
        values.append(water)
    return labels, (np.array(values) if values else np.zeros((0, len(filenames))))


def rank_correlations(comps, values):
    """Spearman rho for every (component, ancillary) pair, NaN where undefined.

    Rank rather than Pearson: a coefficient has no reason to depend linearly on
    airmass, and one outlying night should not manufacture a correlation.
    """
    from scipy.stats import spearmanr

    rho = np.full((len(comps), values.shape[0]), np.nan)
    for i, (_, c) in enumerate(comps):
        for j in range(values.shape[0]):
            good = np.isfinite(c) & np.isfinite(values[j])
            if good.sum() > 10 and np.std(values[j][good]) > 0:
                rho[i, j] = spearmanr(c[good], values[j][good]).statistic
    return rho


def plot_correlations(rho, comps, labels, path, title="", flag=0.3):
    """The matrix, with the star and Earth blocks separated by a rule."""
    import matplotlib.pyplot as plt

    if rho.size == 0:
        log("  no ancillary quantity varies; skipping %s" % path)
        return
    plt.rcParams.update({"font.size": 9, "axes.grid": False})
    fig, ax = plt.subplots(figsize=(1.2 + 0.62 * len(labels), 1.4 + 0.34 * len(comps)))
    im = ax.imshow(rho, cmap=nan_cmap("RdBu_r"), vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(comps)))
    ax.set_yticklabels([n for n, _ in comps], fontsize=8)
    n_star = sum(1 for n, _ in comps if n.startswith("a"))
    if 0 < n_star < len(comps):
        ax.axhline(n_star - 0.5, color="k", lw=1.0)
    for i in range(rho.shape[0]):
        for j in range(rho.shape[1]):
            if np.isfinite(rho[i, j]) and abs(rho[i, j]) >= flag:
                ax.text(j, i, "%.2f" % rho[i, j], ha="center", va="center",
                        fontsize=7, color="white" if abs(rho[i, j]) > 0.6 else "black")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02,
                 label="Spearman rank correlation")
    ax.set_title(title or "coefficients against ancillary quantities", fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    log("  wrote %s" % path)


def write_correlation_table(rho, comps, labels, path):
    """The same numbers as plot_correlations, machine-readable."""
    import csv

    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["component", "quantity", "spearman_rho"])
        for i, (name, _) in enumerate(comps):
            for j, label in enumerate(labels):
                value = rho[i, j]
                writer.writerow([name, label,
                                 "" if not np.isfinite(value) else "%.4f" % value])
    log("  wrote %s" % path)


def report_correlations(rho, comps, labels, flag=0.7):
    """Print the pairs worth looking at, so a run says it without a PDF."""
    if rho.size == 0:
        return
    strong = [(abs(rho[i, j]), comps[i][0], labels[j], rho[i, j])
              for i in range(rho.shape[0]) for j in range(rho.shape[1])
              if np.isfinite(rho[i, j]) and abs(rho[i, j]) >= flag]
    if not strong:
        log("  no |rho| above %.2f: no component tracks an ancillary quantity"
              % flag)
        return
    log("  |rho| >= %.2f, i.e. a component that follows something known:" % flag)
    for _, name, label, value in sorted(strong, reverse=True):
        log("    %-4s %-22s %+.2f" % (name, label, value))
