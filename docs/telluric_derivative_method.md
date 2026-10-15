# Telluric fit by the derivative method: hand-off note

Written on 2026-10-15 by the Claude session working on `papier_helium` (LHS 1140
He I referee report). Several of Étienne's messages meant for the pca2d session
landed in that one by mistake, and the work they triggered was done there on
WINERED data. This note passes both on. Nothing in `pca2d_preclean` was changed
apart from adding this file.

## Messages meant for this session (verbatim, in order)

1. « non, juste non, les telluriques sont au mauvais endroit! Fais une correlation
   entre tapas et le spectre de l'étoile non-corrigé, tu vas le voir »
2. « autour de 686-689 nm »
3. « en fait celle là est ok TAPAS dans le repere du SOL (BERV +23.7 km/s) »
4. « celle là est fausse TAPAS sur la grille barycentrique (ce que LBL fait) »
5. « le projet telluclean a une meilleure facon de converger avec la dérivée du
   spectre, reprendre cette méthode de convergence plus robuste, tester et
   implémenter si ça marche »
   Clarifications given afterwards: not the LBL `lbl_telluclean`, it is the
   SPIRou/NIRPS code of `telluric_fit` (formerly `tapas_tellu`), « c'est le code
   avec predict abso ».
6. « la masse d'air aussi, ça sera plus précis »
7. « ne pas laisser libre la vitesse »
8. « non pas global par nuit »
9. « la colonne d'eau change trop vite »
10. « les telluriques sont faites par fichier, j'insiste »

Messages 2 to 4 refer to the O2 B band and to plots this session never saw
(BERV +23.7 km/s matches none of the helium data).

## The method in `telluric_fit`

`telluric_fit/tellu_tools.py`, function `optimize_exponents` (called from
`predict_abso.py`):

- The model transmission T is TAPAS raised to per-molecule exponents
  (H2O free, the dry species CH4, CO2 and O2 bounded to ±10% of the airmass),
  convolved to the instrument resolution after the exponentiation.
- The objective is the weighted standard deviation of the pixel-to-pixel
  gradient of the corrected spectrum, `np.gradient(sp / T, axis=1)`. A wrong
  exponent leaves residuals at every telluric line, which inflate that scatter;
  stellar lines contribute the same amount whatever the exponents, so they do
  not bias the minimum.
- Weights:
  - smooth transmission weights from `get_transmission_weights`, with H2O
    parameters `depth_max` 0.6, `depth_saturated` 0.5, `reject_saturated` 0.8
    and a sigmoid width of 0.02 (`telluric_config_nirps.yaml`); these drop
    saturated lines and their wings without hard masks that make the
    minimizer oscillate;
  - a blaze weight (0 below 30% of the 90th percentile of the blaze);
  - an outlier roll-off on the gradient, `1 / (1 + exp(dev - 6))`, with sigma
    taken from the 16th-84th percentiles.
- `weighted_nanstd` uses reliability weights:
  `sum(w (x - m)^2) / (sum(w) - sum(w^2) / sum(w))`.
- The minimizer is Nelder-Mead, tolerance 5e-4, with H2O bounds (0.001, 20).
- In that code the velocity is not a fit parameter.

## What was tested (WINERED, order 163, 1079-1086 nm, 190 frames over 3 nights)

The test used the TAPAS La Silla model (`~/tapas/European Southern Obs (ESO) La
Silla Chile_2024-01-01_1000-1250nm/`), convolved to R = 68000. Only H2O matters
in that range (O2 absorbs 0.3% at most). Results:

- **Velocity and H2O exponent both free, per frame:**
  - converges in about 0.3 s per frame, with a clear minimum in velocity;
  - the velocity scatters by 0.15 to 0.20 km/s from frame to frame;
  - without bounds on the exponent, a few frames diverge (exponents of
    -1.6 and -157679), so the bounds are needed.
- **Comparison with the cross-correlation method:** the CCF of a
  high/low-airmass ratio against TAPAS gives velocities 0.3 to 0.66 km/s lower
  than the derivative method on all three nights.
- **One velocity per night, H2O exponent = water column x airmass:**
  - the nod A and nod B frames give the same velocity within 0.06 km/s;
  - the first and second halves of each night differ by 0.15 to 0.28 km/s;
  - the per-frame H2O exponent does not follow the airmass closely, which
    agrees with message 9: the water column changes too fast for a model
    shared over a night.
- **Scheme following messages 6 to 10, applied before being dropped:**
  - the velocity is not free;
  - the dry absorbers are fixed at the airmass of each file;
  - the H2O exponent is the only free parameter, fitted file by file
    (bounded, Brent);
  - each file is divided by its own T, and pixels whose transmission weight
    is below 0.5 are set to NaN, as in `telluric_fit`.

On WINERED the velocity is a real quantity: the zero point of the WARP
wavelength solution, which is +1.15, +1.20 and -4.52 km/s for the three
nights. For APERO-reduced NIRPS or SPIRou spectra the wavelength solution is
already in the ground frame, so a fixed velocity means TAPAS at rest in the
ground frame (message 3), not on the barycentric grid (message 4).

## Core code (stand-alone, numpy + scipy)

```python
import numpy as np
from scipy.optimize import minimize_scalar

C_KMS = 299792.458
H2O_WEIGHTS = dict(depth_max=0.6, depth_saturated=0.5, reject_saturated=0.8,
                   transition_sigma=0.02)


def doppler_ratio(v):
    """Relativistic wavelength ratio for a velocity v (km/s, positive receding)."""
    b = v / C_KMS
    return np.sqrt((1 + b) / (1 - b))


def make_model(wave_tapas, trans_h2o, trans_dry, resolution, step_kms=0.25):
    """Optical depths on a uniform ln(lambda) grid and the LSF kernel."""
    step = np.arctanh(step_kms / C_KMS)
    lnw = np.arange(np.log(wave_tapas[0]), np.log(wave_tapas[-1]), step)
    tau = -np.log(np.clip(np.interp(lnw, np.log(wave_tapas), trans_h2o), 1e-6, 1))
    tau_dry = -np.log(np.clip(np.interp(lnw, np.log(wave_tapas), trans_dry), 1e-6, 1))
    sig = np.arctanh(1 / resolution / 2.3548) / step
    x = np.arange(-int(5 * sig), int(5 * sig) + 1)
    k = np.exp(-0.5 * (x / sig) ** 2)
    return dict(wave=np.exp(lnw), tau=tau, tau_dry=tau_dry, kernel=k / k.sum())


def transmission(model, wave, alpha_h2o, v, airmass):
    """exp(-alpha tau_H2O - airmass tau_dry), convolved after the power,
    shifted by v."""
    conv = np.convolve(np.exp(-alpha_h2o * model['tau'] - airmass * model['tau_dry']),
                       model['kernel'], mode='same')
    return np.interp(wave / doppler_ratio(v), model['wave'], conv)


def transmission_weights(trans, depth_max, depth_saturated, reject_saturated,
                         transition_sigma):
    """Same as telluric_fit.tellu_tools.get_transmission_weights."""
    w_depth = 1 / (1 + np.exp(-(trans - depth_max) / transition_sigma))
    sat, above = trans < depth_saturated, trans >= reject_saturated
    contam = np.zeros(len(trans), bool)
    for idx in (range(len(trans)), range(len(trans) - 1, -1, -1)):
        inside = False
        for i in idx:
            if sat[i]:
                inside = True
            elif inside and above[i]:
                inside = False
            contam[i] |= inside
    return w_depth * (~contam)


def weighted_nanstd(d, w):
    ok = np.isfinite(d) & np.isfinite(w) & (w > 0)
    d, w = d[ok], w[ok]
    m = np.sum(w * d) / np.sum(w)
    return np.sqrt(np.sum(w * (d - m) ** 2) / (np.sum(w) - np.sum(w ** 2) / np.sum(w)))


def objective(model, wave, sp, alpha_h2o, v, airmass, outlier_sigma=6.0):
    trans = transmission(model, wave, alpha_h2o, v, airmass)
    grad = np.gradient(sp / trans)
    w = transmission_weights(trans, **H2O_WEIGHTS)
    valid = w > 0.5
    p16, p84 = np.nanpercentile(grad[valid], [16, 84])
    dev = np.abs(grad - np.nanmedian(grad[valid])) / ((p84 - p16) / 2)
    w = w / (1 + np.exp(np.clip(dev - outlier_sigma, -500, 500)))
    return weighted_nanstd(grad, w)


def correct_file(model, wave, sp, v_fixed, airmass, bounds=(1e-3, 5)):
    """Per-file telluric correction: H2O exponent the only free parameter."""
    good = np.isfinite(sp) & (sp > 0)
    alpha = minimize_scalar(
        lambda a: objective(model, wave[good], sp[good], a, v_fixed, airmass),
        bounds=bounds, method='bounded', options=dict(xatol=1e-4)).x
    trans = transmission(model, wave, alpha, v_fixed, airmass)
    corr = sp / trans
    corr[transmission_weights(trans, **H2O_WEIGHTS) < 0.5] = np.nan
    return corr, alpha
```

This is 1D, for a single order. For an e2ds cube, take the gradient along
axis 1 and add the blaze weight, as `optimize_exponents` does.
