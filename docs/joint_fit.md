# One observer basis for several stars

The atmosphere and the instrument belong to the night, not to the target. Two
stars observed in the same campaign share them; what they do not share is the
star itself, its BERV coverage and its systemic velocity. Fitting several
objects together therefore buys a cleaner observer basis: **what is common to
all of them cannot be a star, and what follows one star cannot be common.**

    pca2d-preclean --objects PROXIMA,GJ1,GJ3090 --n-star 0

Everything that run writes lands under `outputs/joint/PROXIMA+GJ1+GJ3090/0-3/`,
and each object is measured by LBL on its own, as `<object>_PCA2D_<M-N>_joint`.

## What is shared and what is not

| | shared | per object |
|---|---|---|
| observer components (the atmosphere, the instrument) | yes | |
| star spectrum, per order parity | | yes |
| exposure amplitudes b_n | | each row has its own |
| corrected spectra, LBL template, velocities | | yes |

The cube is the objects' own cubes, concatenated row for row on the one grid
they share, with three labels rewritten (`pca2d/joint.py`):

- **parity** becomes `2 * (the object's index) + the order parity`. The fit
  estimates one star spectrum per label, so each object gets its own, while
  `parity % 2` still gives the order parity, **which means the same thing for
  every object**: it is a property of the detector, not of the star.
- **exposure** is offset per object, so the two rows of an exposure stay tied
  to each other and no two objects share an id.
- **object** carries the name, which the metadata already had.

The correction of one object is told which pair of star spectra is its own
(`reconstruct --star-group`, `2 * index`); the observer components it divides
out are the shared ones.

## The nominal parameters

Those of `config.yaml`, as they stand after the TOI-2120 campaign, plus a
fixed star:

- **star:** `--n-star 0`, so the star is the per-parity template alone, taken
  out with a coefficient of exactly 1. A coefficient in front of a term in the
  log is an exponent on the flux, and an exponent on the star's mean spectrum
  describes no star. On TOI-4552 dropping the star component took the robust
  scatter from 14.56 to 11.20 m/s.
- **star spectra:** `mean: star`, one per parity per object, a BERV-binned
  median taken out once before any component.
- **observer components:** 3, on a cubic B-spline star basis (`star_basis:
  spline`), no velocity term.
- **correction:** every observer component divided out only where it is
  significant, all exposures together (`correct.shrink`), and one mask for
  every epoch (`correct.mask: common`).
- **LBL:** the corrected spectra are measured against the template LBL builds
  from them (`lbl.star_template: false`).
- **high pass:** Savitzky-Golay of 100 km/s; grid step 0.5 km/s.

## What the objects must have in common

The joint cube refuses anything else: one instrument, one domain, one grid
step, one high pass, one set of quality cuts. The run compares the objects'
cube keys with the object taken out of them and stops if they differ. Choose
targets of comparable SNR: a much brighter star would dominate the basis by
weight alone.

## What it does not do

- **The water column of the correlation matrix** is read from one object's
  spectra, so it is dropped in a joint run.
- **The velocity term** has not been tried jointly.
- **LBL** never sees the joint cube: each object is prepared, measured and
  compiled on its own, exactly as in a solo run.
