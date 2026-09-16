"""The velocity bias that follows BERV, fitted by MCMC.

A telluric line blended with a stellar line pulls the velocity measured on it
by an amount that depends on how far apart the two are, which is the BERV.
Superposed, the blend is symmetric and nothing moves; far apart, nothing is
blended; in between, the pull is the derivative of a Gaussian:

    v(BERV) = c + amp * BERV * exp(-0.5 * (BERV / sigma)**2)

odd in BERV, zero at 0, extreme at BERV = +-sigma, where it is worth

    peak = amp * sigma * exp(-1/2)

(written `amp*np.exp(-0.5*BERV/sigma)*BERV` when it was asked for, on
2026-09-16: the square is what makes the bias die away on both sides; taken
as written it grows without bound for negative BERV). A straight line against
BERV, which the report drew before, has no such shape and no meaning here.

`amp` is signed, and the fit explores both signs: its prior is symmetric and
half the walkers start on each side. `amp` and `sigma` are what the fit is
for; `c`, the velocities' arbitrary
zero, and `jitter`, the scatter LBL's error bars do not account for, are
fitted beside them and marginalised. Without the jitter the posterior would be
as narrow as LBL's error bars are optimistic: 13 m/s against 44 m/s of scatter
on SMETHELLS_20. The sampler is Goodman and Weare's affine-invariant stretch
move, emcee's algorithm, in a few lines of numpy: emcee is not in the
pipeline's environment, and four parameters do not need it.

    amp       modified Jeffreys, p(amp) ~ 1 / (|amp| + a0), (m/s)/(km/s),
              within +-1e4 (asked for on 2026-09-16). Scale-invariant above
              the knee a0, flat below it, and signed, since a bias pulls
              either way; 1/|amp| itself has no normalisation at zero.
              The knee is the noise level (Gregory 2005, ApJ 631, 1198, for
              RV amplitudes): the amp whose peak, at the prior's median
              width sqrt(1 * 60) km/s, is the median error over sqrt(N)
    sigma     log-uniform over 1 to 60 km/s. Below 1 km/s is narrower than
              any blend of two lines can be at a resolution of 70 000 to
              80 000 (a 4 km/s FWHM): allowed, the few points within +-sigma
              let amp run to hundreds and the upper limit on a bias that is
              not there doubles. Above 60, with BERV within +-30, the
              Gaussian is a straight line the data cannot tell apart
    c         uniform, m/s
    jitter    log-uniform over 1e-3 to 1e4 m/s
"""

from __future__ import annotations

import numpy as np

SIGMA_MIN, SIGMA_MAX = 1.0, 60.0
JITTER_MIN, JITTER_MAX = 1e-3, 1e4
AMP_MAX = 1e4
NAMES = ("amp", "sigma", "c", "jitter")
#: standard deviations from zero before a bias is called a bias. Below it the
#: posterior of sigma is the prior's, and a peak and a width quoted from it
#: would describe the prior; an upper limit on the peak is what is said then
DETECTED = 3.0
UNITS = ("(m/s)/(km/s)", "km/s", "m/s", "m/s")


def shape(berv, amp, sigma):
    """The bias alone, without the offset: amp * B * exp(-B^2 / 2 sigma^2)."""
    berv = np.asarray(berv, float)
    return amp * berv * np.exp(-0.5 * (berv / sigma) ** 2)


def peak(amp, sigma):
    """The bias at BERV = sigma, its largest."""
    return amp * sigma * np.exp(-0.5)


def amp_knee(e, n):
    """The modified Jeffreys prior's knee, from the data's own noise."""
    width = np.sqrt(SIGMA_MIN * SIGMA_MAX)
    return float(np.median(e) / np.sqrt(max(n, 1)) / (width * np.exp(-0.5)))


def log_probability(theta, berv, v, e, knee=None):
    """ln posterior of (amp, ln sigma, c, ln jitter), one row per walker."""
    theta = np.atleast_2d(theta)
    amp, lsig, c, ljit = theta.T
    inside = ((np.abs(amp) < AMP_MAX)
              & (lsig > np.log(SIGMA_MIN)) & (lsig < np.log(SIGMA_MAX))
              & (ljit > np.log(JITTER_MIN)) & (ljit < np.log(JITTER_MAX)))
    out = np.full(theta.shape[0], -np.inf)
    if not inside.any():
        return out
    sigma = np.exp(lsig[inside])[:, None]
    model = c[inside][:, None] + amp[inside][:, None] * berv[None, :] \
        * np.exp(-0.5 * (berv[None, :] / sigma) ** 2)
    var = e[None, :] ** 2 + np.exp(2 * ljit[inside])[:, None]
    out[inside] = -0.5 * np.sum((v[None, :] - model) ** 2 / var + np.log(var),
                                axis=1)
    if knee is None:
        knee = amp_knee(e, e.size)
    out[inside] -= np.log(np.abs(amp[inside]) + knee)
    return out


def stretch(log_prob, start, steps, rng, a=2.0):
    """Goodman and Weare's stretch move. Returns (chain, acceptance).

    `start` is (walkers, ndim) with an even number of walkers; `log_prob`
    takes a (n, ndim) array. The walkers are updated in two halves, each
    against the other, which is what keeps the move valid when it is done
    for many walkers at once.
    """
    pos = np.array(start, float)
    walkers, ndim = pos.shape
    lp = log_prob(pos)
    chain = np.empty((steps, walkers, ndim))
    accepted = 0
    halves = (np.arange(walkers // 2), np.arange(walkers // 2, walkers))
    for step in range(steps):
        for moving, fixed in (halves, halves[::-1]):
            z = ((a - 1.0) * rng.random(moving.size) + 1.0) ** 2 / a
            partners = pos[fixed[rng.integers(0, fixed.size, moving.size)]]
            proposal = partners + z[:, None] * (pos[moving] - partners)
            lp_new = log_prob(proposal)
            keep = np.log(rng.random(moving.size)) < \
                (ndim - 1) * np.log(z) + lp_new - lp[moving]
            pos[moving[keep]] = proposal[keep]
            lp[moving[keep]] = lp_new[keep]
            accepted += int(keep.sum())
        chain[step] = pos
    return chain, accepted / float(steps * walkers)


def starting_point(berv, v, e):
    """The best (amp, ln sigma, c, ln jitter) on a grid of sigma.

    For each sigma the model is linear in (c, amp), so it is solved rather
    than searched; the jitter is the scatter the error bars leave over.
    """
    excess = max(float(np.var(v - np.median(v)) - np.median(e) ** 2), 1.0)
    w = 1.0 / (e ** 2 + excess)
    best = None
    for sigma in np.geomspace(SIGMA_MIN * 1.2, SIGMA_MAX / 1.2, 40):
        A = np.column_stack([np.ones_like(berv), shape(berv, 1.0, sigma)])
        coef, *_ = np.linalg.lstsq(A * np.sqrt(w)[:, None], v * np.sqrt(w),
                                   rcond=None)
        chi2 = float(np.sum(w * (v - A @ coef) ** 2))
        if best is None or chi2 < best[0]:
            best = (chi2, coef[1], sigma, coef[0])
    _, amp, sigma, c = best
    return np.array([amp, np.log(sigma), c, 0.5 * np.log(excess)])


def fit(berv, v, e, walkers=32, steps=2500, burn=1000, seed=0):
    """Posterior of the BERV bias of one velocity series.

    Returns a dictionary: `samples` (n, 4) in (amp, sigma, c, jitter),
    `acceptance`, and for each of amp, sigma, c, jitter and peak its median
    and 16th and 84th percentiles, plus the amp-sigma correlation of the
    posterior. None when there are too few points to fit four parameters.
    """
    berv, v, e = (np.asarray(x, float) for x in (berv, v, e))
    ok = np.isfinite(berv) & np.isfinite(v) & np.isfinite(e) & (e > 0)
    if ok.sum() < 12 or np.ptp(berv[ok]) < 2 * SIGMA_MIN:
        return None
    berv, v, e = berv[ok], v[ok], e[ok]
    rng = np.random.default_rng(seed)
    knee = amp_knee(e, e.size)
    centre = starting_point(berv, v, e)
    scale = np.array([max(abs(centre[0]) * 0.1, 1e-2), 0.1,
                      max(np.std(v) * 0.01, 1e-2), 0.1])
    start = centre + scale * rng.normal(size=(walkers, 4))
    # half the walkers start on the other sign: a bias pulls either way, and
    # an ensemble that starts on one side only explores the other if it
    # happens to wander there. On SMETHELLS_20 the mirrored half rejoins the
    # delivered bias's negative mode within 500 steps, and the corrected
    # posterior keeps both signs, P(amp > 0) = 0.15 either way
    start[walkers // 2:, 0] *= -1
    start[:, 1] = np.clip(start[:, 1], np.log(SIGMA_MIN) + 1e-3,
                          np.log(SIGMA_MAX) - 1e-3)
    chain, acceptance = stretch(
        lambda theta: log_probability(theta, berv, v, e, knee), start,
        steps, rng)
    flat = chain[burn:].reshape(-1, 4)
    samples = np.column_stack([flat[:, 0], np.exp(flat[:, 1]), flat[:, 2],
                               np.exp(flat[:, 3])])
    out = {"samples": samples, "acceptance": acceptance, "n": int(ok.sum()),
           "knee": knee}
    for i, name in enumerate(NAMES):
        out[name] = np.percentile(samples[:, i], [50, 16, 84])
    peaks = peak(samples[:, 0], samples[:, 1])
    out["peak"] = np.percentile(peaks, [50, 16, 84])
    out["peak_sigma"] = float(np.std(peaks))
    # how many standard deviations the peak is from no bias at all
    out["significance"] = (float(abs(np.median(peaks)) / np.std(peaks))
                           if np.std(peaks) > 0 else 0.0)
    out["amp_sigma_r"] = float(np.corrcoef(samples[:, 0],
                                           np.log(samples[:, 1]))[0, 1])
    out["p_positive"] = float(np.mean(samples[:, 0] > 0))
    # what can be said when nothing is seen: the bias is below this
    out["upper"] = float(np.percentile(np.abs(peaks), 95))
    out["detected"] = out["significance"] >= DETECTED
    return out


def envelope(result, grid, draws=400, seed=1):
    """(16th, 50th, 84th) percentiles of the bias curve over `grid`."""
    rng = np.random.default_rng(seed)
    samples = result["samples"]
    pick = samples[rng.integers(0, samples.shape[0], draws)]
    curves = pick[:, 0][:, None] * grid[None, :] \
        * np.exp(-0.5 * (grid[None, :] / pick[:, 1][:, None]) ** 2)
    return np.percentile(curves, [16, 50, 84], axis=0)


def summary(result):
    """'peak -70.6 +9.6/-10.0 m/s at 6.7 km/s (7.1 sigma)', or the limit."""
    if result is None:
        return "not fitted"
    if not result["detected"]:
        return ("none detected (%.1f sigma): |peak| < %.1f m/s at 95%%"
                % (result["significance"], result["upper"]))
    p, lo, hi = result["peak"]
    s = result["sigma"][0]
    return ("peak %.1f +%.1f/-%.1f m/s at %.1f km/s (%.1f sigma)"
            % (p, hi - p, p - lo, s, result["significance"]))
