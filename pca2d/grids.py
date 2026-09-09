"""The destination wavelength grid, and the nominal photometric bands.

Two things live here because both are properties of the *output* of the
pipeline rather than of any particular input file:

  * the "magic grid": a log-uniform wavelength grid with a constant velocity
    step, defined analytically by an anchor wavelength and that step. APERO's
    NIRPS s1d_v products are sampled on exactly this grid with wave0 = 965 nm
    and dv = 0.5 km/s, so building it from those two numbers reproduces the s1d
    sampling to 5e-13 nm and needs no file on disk. That matters for the t.fits
    path, where there is no s1d grid to recycle;

  * the nominal photometric bands. The PCA is fitted on those and only those:
    between the bands the atmosphere is opaque, the blaze has collapsed, or
    both, and columns like that contribute nothing but noise and telluric
    residuals to the eigenvectors.
"""

from __future__ import annotations

import numpy as np

from .logger import log

C_KMS = 299792.458

# MKO half-power points, in nm. K is included on purpose even though NIRPS
# stops at 1.95 um: the same configuration has to serve an instrument that
# does reach K, and a band that falls outside the grid is simply dropped.
DEFAULT_BANDS = {
    "Y": [970.0, 1070.0],
    "J": [1170.0, 1330.0],
    "H": [1490.0, 1780.0],
    "K": [2030.0, 2370.0],
}


def magic_grid(wave0: float, dv: float, wave_min: float, wave_max: float) -> np.ndarray:
    """Log-uniform grid  lambda_i = wave0 * exp(i * dv / c), cropped to the band.

    The anchor matters as much as the step: it is what makes two runs, two
    instruments' worth of reduction, or a t.fits and an s1d of the same
    exposure land on identical wavelengths, so that products can be compared
    sample by sample instead of being re-splined into each other.
    """
    wave0 = float(wave0)
    step = float(dv) / C_KMS
    if step <= 0:
        raise ValueError("domain.dv must be positive, got %r" % dv)
    first = int(np.ceil(np.log(float(wave_min) / wave0) / step))
    last = int(np.floor(np.log(float(wave_max) / wave0) / step))
    if last < first:
        raise ValueError(
            "empty grid: %.3f-%.3f nm at %.4f km/s from an anchor of %.3f nm"
            % (wave_min, wave_max, dv, wave0)
        )
    return wave0 * np.exp(step * np.arange(first, last + 1))


def grid_dv(grid: np.ndarray) -> float:
    """Velocity step of a log-uniform grid, in km/s."""
    return float(np.median(np.diff(np.log(np.asarray(grid, dtype=np.float64)))) * C_KMS)


def band_mask(grid: np.ndarray, bands: dict, names=None):
    """Boolean mask of the grid columns inside the requested bands.

    Returns (mask, used) where `used` lists the bands that actually intersect
    the grid. A band that misses the grid entirely -- K on NIRPS -- is reported
    and ignored rather than treated as an error, because the same config is
    meant to be reusable across instruments.
    """
    grid = np.asarray(grid, dtype=np.float64)
    mask = np.zeros(grid.size, dtype=bool)
    used, missing = [], []
    for name in (names if names else list(bands)):
        if name not in bands:
            log("band '%s' is requested but not defined; ignored" % name, "warn")
            continue
        low, high = [float(v) for v in bands[name]]
        inside = (grid >= low) & (grid <= high)
        if not np.any(inside):
            missing.append(name)
            continue
        mask |= inside
        used.append(name)
        log("band %s: %.1f-%.1f nm, %d grid columns" % (name, low, high, inside.sum()),
            "value")
    if missing:
        log("bands outside the grid, skipped: %s" % ", ".join(missing), "warn")
    if not np.any(mask):
        log("no band overlaps the grid; fitting on the whole domain instead", "warn")
        mask[:] = True
        used = ["<all>"]
    return mask, used
