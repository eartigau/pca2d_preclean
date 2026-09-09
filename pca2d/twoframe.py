#!/usr/bin/env python
"""Two-frame (star + Earth) weighted PCA -- block coordinate descent prototype.

The design is NOTES.md section 11. This is the stage-one build of 11.7: the
blocks are updated one at a time (star first), with the coefficients of both
blocks solved *jointly* per spectrum so the fit can decide which frame a feature
belongs to. It is a benchmark and a behaviour check, not production code: the
means are handled by a plain observer-frame subtraction rather than the
component-zero trick of 11.1, and nothing is written out but numbers.

    pca2refs-fit --iters 6
    python pca2d/twoframe.py --order earth_first     # the losing order

Working frame is the OBSERVER frame, so the weights (photon noise, telluric
ramp) never move. The star basis is carried into each spectrum's frame by a
Fourier phase ramp: on the log-uniform magic grid a Doppler shift is a pure
translation, so the operator is exact and unitary and its adjoint is simply the
opposite shift (NOTES.md 11.3).
"""

from __future__ import annotations

import argparse
import os
import time
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .logger import log
import numpy as np
from astropy.io import fits
from astropy.table import Table
from numpy.lib.stride_tricks import sliding_window_view
from scipy.fft import irfft, next_fast_len, rfft

C_KMS = 299792.458
PARITY_NAMES = ("even", "odd")
EDGE = 64          # samples zeroed at each end; the FFT wraps, the shifts are <=43


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=None,
                   help="YAML whose `twoframe:` section supplies every option"
                        " below. Anything given on the command line wins over"
                        " it, and anything in neither falls back to"
                        " pca2d/config.py:DEFAULTS")
    p.add_argument("--cube", default="cache/cube_06b3ec43de9f.npz",
                   help="observer-frame cube, log_sub high-pass (build with registration.frame: observer)")
    p.add_argument("--order", choices=["star_first", "earth_first"], default=None,
                   help="which block is updated first; star_first for a reason, see 11.7")
    p.add_argument("--iters", type=int, default=None)
    p.add_argument("-k", "--n-star", type=int, default=None,
                   help="M, components in the STELLAR rest frame")
    p.add_argument("-j", "--n-earth", type=int, default=None,
                   help="N, components in the OBSERVER frame")
    p.add_argument("--leakage", action="store_true", default=None,
                   help="report the cross-block leakage matrix at the end (slow)")
    p.add_argument("--shift", choices=["lanczos", "fft"], default=None,
                   help="shift operator; fft is kept only for comparison (NOTES 11.3)")
    p.add_argument("-a", "--kernel-halfwidth", type=int, default=None,
                   help="Lanczos support in samples")
    p.add_argument("--chunk", type=int, default=None,
                   help="spectra per carry block; default keeps the (chunk, K, M)"
                        " buffer near 128 MB, which is what stops the K=10 run"
                        " from swapping")
    p.add_argument("--no-gap-guard", dest="gap_guard", action="store_false", default=None,
                   help="keep samples whose tap window reaches into a hole in the"
                        " basis support (see gap_guard); only useful to measure"
                        " what the guard is buying")
    p.add_argument("--min-snr-frac", type=float, default=None,
                   help="drop spectra below this fraction of the median band SNR")
    p.add_argument("--clip", type=float, default=None,
                   help="soft-clip threshold in sigma; 0 disables")
    p.add_argument("--mean", choices=("offset", "full"), default="offset",
                   help="how much of the observer-frame mean to take out"
                        " before the fit. 'offset' (default) removes only the"
                        " even-minus-odd half-difference, the instrumental"
                        " term of NOTES 13.9, and leaves the part both"
                        " parities share in the data so the Earth block can"
                        " describe it and the correction can remove it."
                        " 'full' removes the whole mean per parity, which is"
                        " what every fit before 2026-09-09 did; the static"
                        " atmosphere is then outside the model and survives"
                        " the correction untouched")
    p.add_argument("--no-template", dest="template", action="store_false", default=None,
                   help="skip the star-frame median template. The star block"
                        " then spends its first component rebuilding the mean"
                        " spectrum instead of describing variability")
    p.add_argument("--outdir", default="_obsolete/outputs/PROXIMA/twoframe",
                   help="where to write coefficients.csv and coefficients_vs_time.pdf")
    p.add_argument("--dtype", choices=("float64", "float32"), default=None,
                   help="storage for the two (rows, pixels) arrays the fit"
                        " lives in. float32 halves them, and their two model"
                        " temporaries with them, which is what lets every"
                        " exposure be fitted rather than nightly means. The"
                        " weighted sums accumulate in float64 either way")
    p.add_argument("--no-tie-parities", dest="tie_parities", action="store_false", default=None,
                   help="solve a separate coefficient vector for each parity"
                        " row. The two rows of an exposure are one measurement,"
                        " so this is only meaningful as a diagnostic")
    # --max-mad is kept as an alias so every existing command, config and log
    # line still means exactly what it meant. The name was wrong, not the
    # number: the cut is in robust sigmas, 1.4826 x MAD, so the default of 10
    # has always been 14.83 plain MADs. Renaming without changing the value is
    # deliberate; changing the value would make every fit already on disk
    # incomparable with the next one.
    # LBL's own values, so the two steps of the chain bin the same way; see
    # lbl/recipes/lbl_template.py, "bin cube by BERV (to give equal weighting
    # to epochs)". 0 restores the plain median over all spectra.
    p.add_argument("--template-berv-bin", type=float, default=None,
                   help="BERV bin width in m/s for the hierarchical template;"
                        " 0 or unset gives the plain median")
    p.add_argument("--template-berv-min-entries", type=int, default=3,
                   help="exposures a BERV bin needs before it gets a vote")
    p.add_argument("--max-sigma", "--max-mad", dest="max_mad", type=float,
                   default=None,
                   help="drop a spectrum whose coefficient on ANY of the K+J"
                        " components sits more than this many robust sigmas"
                        " (1.4826 x MAD) from that component's median, then"
                        " refit; 0 disables the cut. 10 sigma = 14.83 MAD")
    p.add_argument("--max-mad-rounds", type=int, default=None,
                   help="how many times to reject and refit at most")
    p.add_argument("--replot", action="store_true",
                   help="regenerate the figures from <outdir>/fit.npz and exit,"
                        " without loading a cube or refitting")
    args = p.parse_args(argv)
    return _resolve(args)


# every option the YAML can supply; the CLI still wins over it
CONFIGURABLE = ("n_star", "n_earth", "iters", "order", "tie_parities", "max_mad",
                "template_berv_bin", "template_berv_min_entries",
                "max_mad_rounds", "clip", "min_snr_frac", "template", "shift",
                "kernel_halfwidth", "gap_guard", "leakage", "chunk", "dtype")


def _resolve(args):
    """Fill anything the command line left unset from the YAML, then DEFAULTS.

    Three layers, and the order matters: an explicit flag must beat the config
    file, or a one-off experiment silently inherits the nominal settings and the
    run is not what its command line says it is. Options default to None rather
    than to a value so "not given" is distinguishable from "given the default".
    """
    from .config import DEFAULTS, load_config

    section = dict(DEFAULTS["twoframe"])
    if args.config:
        section.update(load_config(args.config).get("twoframe", {}) or {})
        log("options from %s: %s" % (args.config, ", ".join(
            "%s=%s" % (k, section[k]) for k in ("n_star", "n_earth", "iters",
                                                "tie_parities", "max_mad"))))
    for name in CONFIGURABLE:
        if getattr(args, name, None) is None:
            setattr(args, name, section[name])
    return args


# --------------------------------------------------------------------------
# the shift operator
# --------------------------------------------------------------------------
class FourierShifter:
    """Translation by a fractional number of samples, via a Fourier phase ramp.

    KEPT FOR COMPARISON ONLY -- see NOTES.md 11.3, amendment of 2026-08-24, and
    `--shift fft`. A phase ramp is a valid translation only for a band-limited,
    gap-free signal. The spectra qualify (7.5 samples per FWHM); the weight maps
    do not, being a mask with a step edge at each of ~22 gaps per row. Shifting w
    this way drives 0.63% of samples negative, and the exact diag(S^T W S) under
    a global operator is flat straight through a gap rather than gapped.
    """

    banded = False

    def __init__(self, n_pixels, pad=256, **_):
        self.n = n_pixels
        self.n_pad = next_fast_len(n_pixels + pad)
        self.freq = np.fft.rfftfreq(self.n_pad)

    def rows(self, arr, pix, chunk=64):
        out = np.empty_like(arr)
        for start in range(0, arr.shape[0], chunk):
            stop = min(start + chunk, arr.shape[0])
            buf = np.zeros((stop - start, self.n_pad))
            buf[:, :self.n] = arr[start:stop]
            spec = rfft(buf, axis=1)
            spec *= np.exp(-2j * np.pi * self.freq[None, :] * pix[start:stop, None])
            out[start:stop] = irfft(spec, n=self.n_pad, axis=1)[:, :self.n]
        return out

    def prepare(self, basis):
        buf = np.zeros((basis.shape[0], self.n_pad))
        buf[:, :self.n] = basis
        return rfft(buf, axis=1)

    def carry(self, prepared, pix):
        ramp = np.exp(-2j * np.pi * self.freq[None, None, :] * pix[:, None, None])
        return irfft(prepared[None] * ramp, n=self.n_pad, axis=2)[:, :, :self.n]

    def adjoint(self, arr, pix):
        """S^T x. Only equal to shift(-pix) because the operator is (nearly) unitary."""
        return self.rows(arr, -pix)

    def diag_normal(self, w, pix):
        """diag(S^T W S), approximated by the shifted weight map. This is the
        step that is wrong; the clip exists only to hide the ringing."""
        return np.clip(self.rows(w, -pix), 0.0, None)


class LanczosShifter:
    """Translation by a compactly supported windowed-sinc kernel.

    out[m] = sum_t c_t y[m - base + t],  base = floor(p), t in [-a, a)

    so S[m, j] = c_{j - m + base} and the adjoint is a scatter with the *same*
    taps, exact by construction. Nothing here relies on unitarity, which is what
    11.3 originally wanted the FFT for and which the FFT does not actually
    deliver on real input at fractional shift anyway.

    Two properties the FFT does not have, both measured in NOTES.md 11.3:

      * a gap contaminates only +-a samples, exactly zero beyond, against a 61 px
        span with a 1/k tail for the FFT;
      * S^T W S is banded with bandwidth 2a, so its diagonal -- which is what the
        M-step uses -- is a good approximation rather than a crude one, and is
        computable exactly as sum_t c_t^2 w[j + base - t]. It is non-negative by
        construction, so no clipping is needed anywhere.

    Accuracy costs nothing at this sampling: a = 8 gives an rms error 1.2% of the
    per-pixel photon noise.

    Everything is done on a zero-padded copy so that no sample wraps around; the
    FFT version needs apodisation for the same reason and gets it from EDGE.
    """

    banded = True

    #: worker threads for carry(); set to 1 to serialise, None to autodetect
    threads = None

    def __init__(self, n_pixels, a=8, max_shift=64, threads=None):
        self.n = n_pixels
        self.a = int(a)
        self.pad = int(max_shift) + self.a + 2
        if threads is None:
            threads = min(8, (os.cpu_count() or 1))
        self.threads = int(threads)

    def _taps(self, p):
        """(base, c) with c[i] the weight on y[m - base + t], t = i - a."""
        base = int(np.floor(p))
        frac = p - base
        t = np.arange(-self.a, self.a)
        x = -(frac + t)
        c = np.sinc(x) * np.sinc(x / self.a)
        return base, c / c.sum()

    def _padded(self, arr):
        return np.pad(arr, [(0, 0)] * (arr.ndim - 1) + [(self.pad, self.pad)])

    def _apply(self, padded, base, c, sign):
        """sum_t c_t padded[..., off + t], with off placing the shift.

        sign=+1 is the forward operator, sign=-1 the adjoint (same taps, the
        offset reflected), which is exactly the transpose of the same band.

        Written as one contraction over a sliding-window view rather than as 2a
        accumulating passes. The taps are a contiguous, overlapping band, so the
        window view costs nothing to build and the sum becomes a single matmul
        over the tap axis. Measured 1.7x faster on the 421782-column grid, and
        bit-identical: the same 16 products in the same order.
        """
        off = self.pad - sign * base
        taps = c if sign > 0 else c[::-1]
        lo = off - self.a if sign > 0 else off - self.a + 1
        segment = padded[..., lo:lo + self.n + 2 * self.a - 1]
        window = sliding_window_view(segment, 2 * self.a, axis=-1)
        return window @ taps

    def rows(self, arr, pix):
        padded = self._padded(arr)
        out = np.empty((arr.shape[0], self.n))
        for i in range(arr.shape[0]):
            base, c = self._taps(pix[i])
            out[i] = self._apply(padded[i], base, c, +1)
        return out

    def adjoint(self, arr, pix):
        """S^T x, exactly: the same taps scattered rather than gathered."""
        padded = self._padded(arr)
        out = np.empty((arr.shape[0], self.n))
        for i in range(arr.shape[0]):
            base, c = self._taps(pix[i])
            out[i] = self._apply(padded[i], base, c, -1)
        return out

    def diag_normal(self, w, pix):
        """diag(S^T W S)[j] = sum_t c_t^2 w[j + base - t]. Exact, never negative."""
        padded = self._padded(w)
        out = np.empty((w.shape[0], self.n))
        for i in range(w.shape[0]):
            base, c = self._taps(pix[i])
            out[i] = self._apply(padded[i], base, c ** 2, -1)
        return out

    def prepare(self, basis):
        return self._padded(basis)

    def carry(self, prepared, pix, threads=None):
        """S_n applied to a whole basis, for a chunk of spectra.

        The rows are independent, and after the sliding-window rewrite of
        `_apply` the work is a matmul that releases the GIL, so a thread pool
        helps: measured 5.6x on ten cores for a 421782-column grid, with an
        identical result since each row writes its own slice. Before that
        rewrite the same pool bought only 2.5x, the loop then being limited by
        allocating sixteen temporaries per row rather than by arithmetic.
        """
        out = np.empty((pix.size, prepared.shape[0], self.n))
        n_threads = self.threads if threads is None else threads
        if n_threads and n_threads > 1 and pix.size > 1:
            from concurrent.futures import ThreadPoolExecutor

            def one(i):
                base, c = self._taps(pix[i])
                out[i] = self._apply(prepared, base, c, +1)

            with ThreadPoolExecutor(min(n_threads, pix.size)) as pool:
                list(pool.map(one, range(pix.size)))
            return out
        for i in range(pix.size):
            base, c = self._taps(pix[i])
            out[i] = self._apply(prepared, base, c, +1)
        return out


SHIFTERS = {"fft": FourierShifter, "lanczos": LanczosShifter}


# --------------------------------------------------------------------------
# the three steps
# --------------------------------------------------------------------------
def joint_coeffs(data, w, P, Q, shifter, delta, chunk=64, exposure=None):
    """Step 1 of 11.2: solve [a_n ; b_n] together, per spectrum.

    B_n = [ S_n P^T , Q^T ] is M x (K+J); the off-diagonal blocks of
    B_n^T W_n B_n are the cross-frame overlaps and are the whole point.

    With `exposure` given, the rows sharing an id are solved as ONE measurement.
    For a t.fits cube those rows are the even and the odd orders of a single
    exposure: one spectrum, of one star, at one instant, split across the
    detector. It has one coefficient vector, and this is how that is imposed.

    THE CONSTRAINT. Per row r the normal equations are

        (B_r^T W_r B_r) c = B_r^T W_r y_r,   B_r = [ S_r P^T | Q^T ]

    with W_r the diagonal weights of that row. Rather than solving that twice
    per exposure, the two systems are SUMMED and solved once,

        (sum_r B_r^T W_r B_r) c_n = sum_r B_r^T W_r y_r

    and the single c_n is written back to both rows.

    WHY THAT IS EXACT, not an approximation. B_r depends on the row only through
    S_r, and S_r depends only on the BERV. Both rows are the same exposure, so
    they carry the same BERV to machine precision and B_even = B_odd. The two
    systems therefore differ only through their weights, which say which columns
    each parity covers, and those supports are disjoint outside the order
    overlaps. Summing them is then exactly the normal equation of the single
    least-squares problem in which every sample of the exposure, from either
    parity, constrains one coefficient vector. The sum reassembles the exposure;
    it does not average two estimates of it.

    NOTE what is and is not per parity. The components are not: P and Q are
    single vectors over the whole grid, and there is no such thing as an even
    component or an odd one. Only the observer-frame mean is fitted per parity,
    because the even-minus-odd offset is a static property of the instrument,
    and only the weights know which columns a row covers.
    """
    n_spectra = data.shape[0]
    n_star, n_earth = P.shape[0], Q.shape[0]
    n_tot = n_star + n_earth
    coeffs = np.zeros((n_spectra, n_tot))
    cond = np.zeros(n_spectra)
    eye = np.eye(n_tot)
    Pf = shifter.prepare(P)
    B = np.empty((n_tot, data.shape[1]))
    Bw = np.empty((n_tot, data.shape[1]))
    tied = exposure is not None
    amats = np.zeros((n_spectra, n_tot, n_tot)) if tied else None
    bvecs = np.zeros((n_spectra, n_tot)) if tied else None
    for start in range(0, n_spectra, chunk):
        stop = min(start + chunk, n_spectra)
        SP = shifter.carry(Pf, delta[start:stop])
        for i in range(stop - start):
            n = start + i
            B[:n_star] = SP[i]
            B[n_star:] = Q
            # into a buffer allocated once: B * w[n] is (K+J, n_pixels), 34 MB
            # at these sizes, and allocating it per row is 11% of the loop
            np.multiply(B, w[n], out=Bw)
            amat = Bw @ B.T
            bvec = Bw @ data[n]
            if tied:
                amats[n], bvecs[n] = amat, bvec
                continue
            trace = np.trace(amat)
            if trace <= 0:
                continue
            cond[n] = np.linalg.cond(amat)
            coeffs[n] = np.linalg.solve(amat + 1e-12 * trace * eye, bvec)
    if tied:
        exposure = np.asarray(exposure)
        for value in np.unique(exposure):
            rows = np.where(exposure == value)[0]
            amat = amats[rows].sum(axis=0)
            trace = np.trace(amat)
            if trace <= 0:
                continue
            cond[rows] = np.linalg.cond(amat)
            # the ridge is 1e-12 of the trace: far below any real curvature, but
            # enough to keep the solve finite when a component is unconstrained
            # for this exposure, which happens when its support is fully masked
            coeffs[rows] = np.linalg.solve(amat + 1e-12 * trace * eye,
                                           bvecs[rows].sum(axis=0))
    return coeffs[:, :n_star], coeffs[:, n_star:], cond


def gap_guard(w, delta, a, verbose=True):
    """Zero the weight wherever the carried star basis would reach into a hole.

    This is what it takes to run the two-frame model over a domain that is
    mostly holes. The M-step side is already safe: `diag_normal` reports an
    honest near-zero star-frame weight inside a gap and `mstep` floors it away,
    so the basis is simply not defined there. The *forward* side is not. The
    model at sample m is sum_t c_t P[m - base + t], t in [-a, a), so a sample
    within a of a hole edge is built partly from basis samples that were never
    constrained and are therefore zero. The model comes out too small, the fit
    sees a residual it cannot explain in the star frame, and the Earth block
    obligingly absorbs a rim of star-frame structure around every gap.

    The exact fix is to stop using those samples: a sample is admissible only if
    its whole tap window lands on constrained basis samples. Note this acts on
    the weights, never on the operator, so the forward map, its adjoint and
    `diag_normal` stay the same object -- which matters, because an operator
    that is not its own adjoint is what NOTES.md 11.9 spent a day diagnosing.

    Apodising the taps at the hole edge, or renormalising them over the live
    subset, would retain those samples instead of dropping them. It is not worth
    it here: on the 965-1950 nm NIRPS domain the holes are a few dozen wide
    blocks, so the rim is a fraction of a percent of the data, and the price
    would be a shift operator whose adjoint has to be rederived.
    """
    live = w.sum(axis=0) > 0
    n_pixels = live.size
    window = np.lib.stride_tricks.sliding_window_view(
        np.concatenate([np.zeros(a, bool), live, np.zeros(a, bool)]), 2 * a)
    eroded = window.all(axis=1)[:n_pixels]      # eroded[j] = all(live[j-a:j+a])

    before = int(np.count_nonzero(w > 0))
    base = np.floor(delta).astype(int)
    index = np.arange(n_pixels)
    for value in np.unique(base):
        rows = base == value
        source = index - value                  # sample m needs eroded[m - base]
        ok = (source >= 0) & (source < n_pixels)
        admissible = np.zeros(n_pixels, dtype=bool)
        admissible[ok] = eroded[source[ok]]
        w[rows] *= admissible
    after = int(np.count_nonzero(w > 0))

    if verbose:
        edges = int(np.count_nonzero(np.diff(live.astype(np.int8)) != 0))
        log("gap guard: %d holes in the basis support (%.1f%% of the grid),"
              " %d distinct integer shifts" % (edges // 2 + 1,
                                               100 * np.mean(~live),
                                               np.unique(base).size))
        log("  dropped %d of %d weighted samples (%.3f%%) within %d of a hole"
              % (before - after, before, 100 * (before - after) / max(before, 1), a))
    return w


def star_model(P, a, shifter, delta, n_spectra, n_pixels, chunk=64):
    """sum_k a_nk (S_n P_k), never materialising S_n P^T for all n at once."""
    out = np.zeros((n_spectra, n_pixels))
    Pf = shifter.prepare(P)
    for start in range(0, n_spectra, chunk):
        stop = min(start + chunk, n_spectra)
        out[start:stop] = np.einsum(
            "nk,nkm->nm", a[start:stop], shifter.carry(Pf, delta[start:stop])
        )
    return out


def mstep(residual, weights, coeffs, basis, floor_frac=1e-2):
    """One EMPCA eigenvector sweep, in whatever frame the residual lives in.

    `denominator > 0` is not a sufficient guard. With an honest banded operator
    the star-frame weight really does fall to ~1e-15 of nominal inside a telluric
    gap (that is the whole point of NOTES 11.3), and dividing by it puts a huge
    spike in the eigenvector: it drove cond(A) from 6 to 5e10 the first time this
    ran, and a 1e-8 floor still left max|v|/rms at 196 against 17 for the input.
    The FFT operator never tripped this only because its smearing filled the gaps
    back in with fictitious weight.

    floor_frac = 1e-2 costs essentially nothing here, because the denominator is
    strongly bimodal: constrained samples sit within a factor of ~1.3 of the
    maximum (median 6.4e4 against max 8.7e4) and only **0.11%** of samples fall
    below 1% of it. Measured counts below a given fraction of the max, out of
    38696 samples: 1e-8 -> 29, 1e-4 -> 36, 1e-2 -> 44, 1e-1 -> 46. There is no
    population in between to lose; the spikes come entirely from the handful of
    samples sitting just above whatever floor is chosen.

    The floor is a fraction of the **median** of the positive denominators, not
    of their maximum. Over the 1500-1600 nm band those two differ by a factor
    1.3 and the choice is immaterial, which is why the original spelling was
    fine. Over 965-1950 nm they do not: the per-order extraction SNR runs from
    108 to 294, so the weight alone spans a factor of seven before the telluric
    ramp and the per-column exposure count are applied, and a floor set at 1% of
    the global maximum would silently delete legitimately constrained samples at
    the blue end of Y. What the floor has to separate is real data from
    numerical residue, and residue sits at ~1e-15 of nominal; the typical
    constrained weight is the right yardstick for that, the largest one is not.
    """
    basis = basis.copy()
    # NOTE: `residual` is consumed in place. Every caller passes a freshly built
    # temporary, and copying it costs an N x M float64 array we cannot spare.
    for k in range(basis.shape[0]):
        ck = coeffs[:, k]
        numerator = (weights * residual).T @ ck
        denominator = weights.T @ (ck ** 2)
        positive = denominator[denominator > 0]
        scale = float(np.median(positive)) if positive.size else 0.0
        floor = floor_frac * scale
        with np.errstate(invalid="ignore", divide="ignore"):
            vec = np.where(denominator > floor,
                           numerator / np.maximum(denominator, floor), 0.0)
        norm = np.linalg.norm(vec)
        if norm > 1e-12:
            basis[k] = vec / norm
        residual -= np.outer(ck, basis[k])
    for i in range(basis.shape[0]):          # modified Gram-Schmidt
        for j in range(i):
            basis[i] -= np.dot(basis[i], basis[j]) * basis[j]
        basis[i] /= np.linalg.norm(basis[i])
    return basis


def update_star(data, w, P, Q, a, b, shifter, delta, chunk):
    """Step 3: deflate the Earth model, carry the weighted residual home.

    The normal equation wants sum_n c^2 S_n^T W_n S_n; we use its diagonal. With
    a banded operator that matrix has bandwidth 2a and the diagonal is a fair
    approximation to it; with the FFT it is dense and the diagonal is not (11.3).
    Either way the diagonal is now computed exactly by `diag_normal` rather than
    approximated by a shifted weight map, and needs no clipping.
    """
    # built in place: `data - b @ Q`, then `w * that`, then the adjoint. Spelled
    # the obvious way this holds four N x M float64 arrays at once.
    resid = b @ Q
    resid *= -1.0
    resid += data
    resid *= w
    resid_star = shifter.adjoint(resid, delta)              # S^T (W r)
    del resid
    w_star = shifter.diag_normal(w, delta)                  # diag(S^T W S)
    with np.errstate(invalid="ignore", divide="ignore"):
        resid_star /= np.where(w_star > 1e-12, w_star, 1.0)
    resid_star[w_star <= 1e-12] = 0.0
    return mstep(resid_star, w_star, a, P)


def update_earth(data, w, P, Q, a, b, shifter, delta, chunk):
    """Step 2: deflate the star model. No shifting of the weights at all."""
    model = star_model(P, a, shifter, delta, data.shape[0], data.shape[1], chunk)
    model *= -1.0
    model += data                       # in place: data - model
    return mstep(model, w, b, Q)


# --------------------------------------------------------------------------
def read_cube_files(path, dtype=np.float64):
    """Read a cached cube, in either of the two formats `cube.py` has written.

    Until 2026-08-26 the cache was a single compressed `.npz`. It is now a
    directory holding `grid.npy`, `data.npy`, `sigma.npy`, an optional
    `trans.npy` and `meta.fits`: a full-domain t.fits cube is several GB, and
    `savez_compressed` on that spends minutes in zlib for a payload that is
    float32 noise and does not compress. Both are read here so the older caches
    on disk stay usable.
    """
    if os.path.isdir(path):
        def _load(name, required=True):
            full = os.path.join(path, name)
            if os.path.exists(full):
                return np.load(full)
            if required:
                raise FileNotFoundError("%s has no %s" % (path, name))
            return None
        grid = _load("grid.npy")
        data = np.asarray(_load("data.npy"), dtype=dtype)
        sigma = np.asarray(_load("sigma.npy"), dtype=dtype)
        trans = _load("trans.npy", required=False)
        meta = Table.read(os.path.join(path, "meta.fits"))
    else:
        blob = np.load(path, allow_pickle=True)
        grid = blob["grid"]
        data = np.asarray(blob["data"], dtype=dtype)
        sigma = np.asarray(blob["sigma"], dtype=dtype)
        trans = blob["trans"] if "trans" in blob.files else None
        meta = Table(blob["meta"])
    if trans is not None:
        trans = np.asarray(trans, dtype=dtype)
    return grid, data, sigma, trans, meta


def row_parity(meta, n_rows):
    """Which half of the echellogram each cube row came from.

    A t.fits cube carries two rows per exposure, the even orders and the odd
    ones (NOTES.md 13.2). An s1d cube, and any cube built before 2026-08-26,
    has one row per exposure and no `parity` column; those get a single group,
    so everything downstream can be written once.
    """
    if "parity" in getattr(meta, "colnames", []):
        return np.asarray(meta["parity"], dtype=int)
    return np.zeros(n_rows, dtype=int)


def row_exposure(meta, n_rows):
    """Which original exposure each cube row came from.

    The two parity rows of a t.fits exposure share this id, which is what ties
    them together in the coefficient solve. Without an `exposure` column every
    row is its own exposure and the tie is a no-op.
    """
    if "exposure" in getattr(meta, "colnames", []):
        return np.asarray(meta["exposure"], dtype=int)
    return np.arange(n_rows, dtype=int)


def parity_means(data, w, parity):
    """Weighted mean spectrum per parity, one row per group, in the cube frame.

    Not a refinement: with a t.fits cube a single global mean leaves a static
    even-minus-odd offset in the residual, which is the order overlaps not being
    at the same resolution, and it is the largest coherent structure in the
    cube. In the main pipeline it took the whole of PC1 (NOTES.md 13.9). It is
    fixed in the instrument frame, so the observer-frame mean is exactly where
    it belongs and removing it per parity kills it outright.
    """
    groups = np.unique(parity)
    means = np.zeros((groups.size, data.shape[1]))
    for i, value in enumerate(groups):
        rows = parity == value
        total = w[rows].sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            means[i] = np.where(total > 0,
                                (w[rows] * data[rows]).sum(axis=0)
                                / np.where(total > 0, total, 1.0), 0.0)
    return means, groups


def fit_means(blob, meta, n_rows, n_pixels):
    """The mean spectrum a saved fit removed, as (means, group_per_row).

    `means[group]` is what was subtracted from each cube row. Since 2026-08-26
    a fit.npz stores one mean per order parity and the parity of every row;
    older files store a single `mean`, which is read back as one group so the
    caller needs no special case.
    """
    files = list(getattr(blob, "files", []))
    if "means" in files:
        means = np.atleast_2d(blob["means"])
        parity = np.asarray(blob["parity"], dtype=int) if "parity" in files else None
        if parity is None or parity.size != n_rows:
            # the fit and this cube disagree on which rows survived the SNR cut;
            # the parity of a row is a property of the cube, so take it from there
            parity = row_parity(meta, n_rows)
        groups = np.unique(parity)
        if groups.size != means.shape[0]:
            raise SystemExit(
                "the fit removed %d mean spectra but this cube has %d order"
                " parities: mismatched runs" % (means.shape[0], groups.size))
        return means, np.searchsorted(groups, parity)
    mean = blob["mean"] if "mean" in files else np.zeros(n_pixels)
    return np.atleast_2d(mean), np.zeros(n_rows, dtype=int)


def mean_rows(means, group, n_rows):
    """`means[group]` as a broadcast view when there is only one group.

    With a single mean this is a zero-copy view; fancy-indexing it instead
    would allocate a full (N, M) array, which on the 764 x 92426 cube of
    section 11 is 565 MB of pure waste.
    """
    if means.shape[0] == 1:
        return np.broadcast_to(means[0], (n_rows, means.shape[1]))
    return means[group]


def subtract_means(data, means, group):
    """data -= means[group], in place and without an (N, M) temporary."""
    if means.shape[0] == 1:
        data -= means[0][None, :]
        return data
    for i in range(means.shape[0]):
        rows = group == i
        data[rows] -= means[i][None, :]
    return data


def load_cube(path, ln_clip_low=-0.5, ramp_zero=0.5, min_snr_frac=0.5,
              dtype=np.float64):
    """The observer-frame cube plus the weights, condensed from cube.py.

    `min_snr_frac` drops spectra whose band SNR is below that fraction of the
    median. The pipeline's `quality.min_snr` is an *absolute* floor, which only
    catches outright extraction failures; this is the relative cut. A spectrum at
    a third of the usual SNR is not wrong, it is just nine times less
    informative, and the weights already know that -- but it still gets a full
    vote in the median template and contributes a row of mostly-noise
    coefficients, so it is cheaper to drop it than to carry it.
    """
    grid, data, sigma, trans, meta = read_cube_files(path, dtype)

    if min_snr_frac:
        snr = np.asarray(meta["snr_band"], dtype=float)
        threshold = float(min_snr_frac * np.nanmedian(snr))
        keep = np.isfinite(snr) & (snr >= threshold)
        if not keep.all():
            log("dropping %d / %d spectra with band SNR < %.0f%% of the median"
                  " (%.1f)" % ((~keep).sum(), keep.size, 100 * min_snr_frac,
                               threshold))
            data, sigma, meta = data[keep], sigma[keep], meta[keep]
            if trans is not None:
                trans = trans[keep]

    with np.errstate(invalid="ignore", divide="ignore"):
        w = np.where(np.isfinite(sigma) & (sigma > 0), 1.0 / sigma ** 2, 0.0)
    w[~np.isfinite(data)] = 0.0
    w[data < ln_clip_low] = 0.0
    if trans is not None:
        w *= np.clip((trans - ramp_zero) / (1.0 - ramp_zero), 0.0, 1.0)
    w[:, :EDGE] = 0.0
    w[:, -EDGE:] = 0.0
    data = np.where(w > 0, data, 0.0)
    return grid, data, w, meta


def leakage(data, w, P, Q, shifter, delta, chunk=64):
    """How much of each block can the *other* block reproduce, averaged over n."""
    n_spectra, n_pixels = data.shape
    n_star, n_earth = P.shape[0], Q.shape[0]
    to_earth = np.zeros(n_star)
    to_star = np.zeros(n_earth)
    Pf = shifter.prepare(P)
    for start in range(0, n_spectra, chunk):
        stop = min(start + chunk, n_spectra)
        SP = shifter.carry(Pf, delta[start:stop])
        for i in range(stop - start):
            wi = w[start + i]
            gram_q = (Q * wi) @ Q.T + 1e-10 * np.eye(n_earth)
            gram_p = (SP[i] * wi) @ SP[i].T + 1e-10 * np.eye(n_star)
            for k in range(n_star):
                vec = SP[i, k]
                norm = float(vec @ (wi * vec))
                if norm > 0:
                    rhs = (Q * wi) @ vec
                    to_earth[k] += float(np.linalg.solve(gram_q, rhs) @ rhs) / norm
            for j in range(n_earth):
                vec = Q[j]
                norm = float(vec @ (wi * vec))
                if norm > 0:
                    rhs = (SP[i] * wi) @ vec
                    to_star[j] += float(np.linalg.solve(gram_p, rhs) @ rhs) / norm
    return to_earth / n_spectra, to_star / n_spectra



# --------------------------------------------------------------------------
# post-processing and output
# --------------------------------------------------------------------------
def block_power(w, coeffs, carried, chunk=64):
    """Weighted power carried by each component of one block.

    sum_nm w_nm (c_nk V_nkm)^2, with V the block's vectors as they appear in the
    observer frame: fixed for the Earth block, shifted per spectrum for the
    star block. Same quantity `wpca.component_power` computes for the
    single-frame case, and the right thing to rank on for the same reason: the
    pixels a component acts on can carry very different weights from one
    component to the next.
    """
    n_spectra = w.shape[0]
    n_comp = coeffs.shape[1]
    power = np.zeros(n_comp)
    for start in range(0, n_spectra, chunk):
        stop = min(start + chunk, n_spectra)
        vecs = carried(start, stop)                      # (n, k, m)
        if vecs.ndim == 2:                               # fixed basis
            power += ((coeffs[start:stop] ** 2)
                      * (w[start:stop] @ (vecs ** 2).T)).sum(axis=0)
        else:
            power += ((coeffs[start:stop] ** 2)
                      * np.einsum("nm,nkm->nk", w[start:stop], vecs ** 2)).sum(axis=0)
    return power


def tidy_block(basis, coeffs, power):
    """Sort by decreasing weighted power and pin the sign of each vector.

    The sign of an eigenvector is arbitrary; fixing it (largest excursion
    positive) is what `wpca.sign_convention` does, and keeps runs comparable.
    """
    order = np.argsort(power)[::-1]
    basis, coeffs, power = basis[order].copy(), coeffs[:, order].copy(), power[order]
    for k in range(basis.shape[0]):
        if basis[k][np.argmax(np.abs(basis[k]))] < 0:
            basis[k] *= -1.0
            coeffs[:, k] *= -1.0
    return basis, coeffs, power


def nightly_bin(bjd, coeffs):
    night = np.floor(bjd - 0.5).astype(int)
    unique = np.unique(night)
    out_t = np.zeros(unique.size)
    out_c = np.zeros((unique.size, coeffs.shape[1]))
    for i, value in enumerate(unique):
        sel = night == value
        out_t[i] = np.mean(bjd[sel])
        out_c[i] = np.mean(coeffs[sel], axis=0)
    return out_t, out_c


def plot_coeffs(bjd, a, b, pa, pb, chi2_null, path, err_a=None, err_b=None):
    """Both blocks' coefficients against time, star on the left, Earth on the right.

    Side by side and on a shared time axis on purpose: the question this figure
    exists to answer is whether a feature in one block also appears in the other,
    which is what cross-frame leakage looks like from the coefficient side.
    """
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25})
    n_rows = max(a.shape[1], b.shape[1])
    t = bjd - 2400000.5
    fig, axes = plt.subplots(n_rows, 2, figsize=(13, 1.7 * n_rows + 1.2),
                             sharex=True, squeeze=False)
    blocks = [("star rest frame", a, pa, "tab:blue", err_a),
              ("observer (Earth) rest frame", b, pb, "tab:red", err_b)]
    for col, (title, coeffs, power, colour, err) in enumerate(blocks):
        bt, bc = nightly_bin(bjd, coeffs)
        # the nightly error is the per-exposure error added in quadrature and
        # divided by the count, i.e. sqrt(sum e^2)/n -- not e/sqrt(n), because
        # the exposures in a night do not all carry the same weight
        if err is not None:
            night = np.floor(bjd - 0.5).astype(int)
            uniq = np.unique(night)
            be = np.zeros((uniq.size, coeffs.shape[1]))
            for i, v in enumerate(uniq):
                sel = night == v
                be[i] = np.sqrt(np.nansum(err[sel] ** 2, axis=0)) / sel.sum()
        for row in range(n_rows):
            ax = axes[row][col]
            if row >= coeffs.shape[1]:
                ax.axis("off")
                continue
            ax.plot(t, coeffs[:, row], ".", ms=2.5, alpha=0.35, color=colour)
            if err is not None:
                ax.errorbar(bt - 2400000.5, bc[:, row], yerr=be[:, row],
                            fmt="none", ecolor="0.35", elinewidth=0.6,
                            capsize=0, alpha=0.8, zorder=2)
            ax.plot(bt - 2400000.5, bc[:, row], "o", ms=3.2, mfc="none",
                    mew=0.8, color="k", label="nightly mean", zorder=3)
            ax.axhline(0.0, color="k", lw=0.6, alpha=0.4)
            ax.set_ylabel("%s%d" % ("a" if col == 0 else "b", row + 1))
            ax.text(0.995, 0.95, "%.1f%% of weighted variance"
                    % (100 * power[row] / chi2_null),
                    transform=ax.transAxes, ha="right", va="top", fontsize=7,
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7",
                              lw=0.5, alpha=0.85))
            ax.margins(y=0.18)
            if row == 0:
                ax.set_title(title)
                ax.legend(fontsize=7, loc="lower left")
        axes[-1][col].set_xlabel("MJD (days)")
    fig.suptitle("Two-frame PCA: coefficients versus time\n"
                 "left = star-anchored block (fitted first), "
                 "right = Earth-anchored block, same fit", y=0.998, fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    log("  wrote %s" % path)


def plot_variance(power_star, power_earth, chi2_null, chi2_best, path):
    """Weighted variance carried by each component of both blocks.

    A single-frame PCA has one bar per component because it fits one basis.
    This model fits two, so the figure carries K + J bars: how
    the power divides between the star-anchored and the Earth-anchored block is
    the whole question the two-frame model exists to answer.

    Everything is quoted against `chi2_null`, the cube after the star-frame
    median template and the per-parity observer mean have been taken out. That
    is what the blocks actually fit and it is the only reference that means
    anything here: the raw pre-template variance is dominated by the stellar
    spectrum itself, which no component is trying to describe, so a share of it
    says nothing about the model. The headline number is the RMS of the residual
    against the RMS of that median-subtracted cube.

    DO NOT read these bars against the single-frame runs in `_obsolete/outputs/`.
    Neither half of the fraction is the same quantity. The numerator here is
    `block_power`, a component's own weighted power, which that pipeline used
    only to *rank* components; the bars it plotted were the incremental chi2
    removed when a component is added on top of the previous ones. The two agree
    only when the coefficients are weighted-uncorrelated, which they are not.
    The denominator differs too: there, the mean-subtracted cube; here, the cube
    after a star-frame template that already takes out 95% of the variance.
    """
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25})
    n_star, n_earth = power_star.size, power_earth.size
    power = np.concatenate([power_star, power_earth])
    share = 100 * power / chi2_null                  # of what the blocks fit
    # RMS of the residual over RMS of the median-subtracted cube
    rms_ratio = np.sqrt(max(chi2_best / chi2_null, 0.0)) if chi2_best else np.nan
    labels = (["a%d" % (k + 1) for k in range(n_star)]
              + ["b%d" % (j + 1) for j in range(n_earth)])
    colours = ["tab:blue"] * n_star + ["tab:red"] * n_earth
    x = np.arange(power.size)

    fig, ax = plt.subplots(figsize=(1.0 * power.size + 3.0, 4.2))
    ax.bar(x, share, color=colours, alpha=0.85)
    for xi, s in zip(x, share):
        ax.text(xi, s, "%.2f" % s, ha="center", va="bottom",
                fontsize=6.5, color="0.25")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("component  (a = star rest frame, b = observer rest frame)")
    ax.set_ylabel("component's own weighted power\n(% of the median-subtracted cube)")
    ax.margins(y=0.16)

    twin = ax.twinx()
    twin.plot(x, np.cumsum(share), "k.-", lw=1)
    twin.set_ylabel("cumulative (%)")
    twin.grid(False)

    handles = [plt.Rectangle((0, 0), 1, 1, color="tab:blue", alpha=0.85),
               plt.Rectangle((0, 0), 1, 1, color="tab:red", alpha=0.85)]
    ax.legend(handles, ["star block, %d components (%.2f%% together)"
                        % (n_star, share[:n_star].sum()),
                        "Earth block, %d components (%.2f%% together)"
                        % (n_earth, share[n_star:].sum())],
              fontsize=7.5, loc="upper right")
    fig.suptitle("Two-frame PCA: contribution to variance, %d + %d components"
                 % (n_star, n_earth), fontsize=10)
    ax.set_title("bars are each component's own weighted power, as %% of the"
                 " median-subtracted cube.\nThe %d together take the residual RMS"
                 " to %.4f of its starting value: a %.2f%% reduction."
                 % (power.size, rms_ratio, 100 * (1 - rms_ratio)),
                 fontsize=7.5, color="0.35")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    log("  wrote %s" % path)


def write_components_fits(path, grid, P, Q, template, means, power_star,
                          power_earth, chi2_null, chi2_best, table, dv,
                          tied=True, max_mad=0.0, frame="observer",
                          highpass="log_sub", weights=None, parity=None,
                          min_fraction=0.2):
    """Everything needed to rebuild a spectrum, in one file.

    Three extensions. BASIS carries the magic grid, the star-frame template,
    one observer-frame mean per parity and the K + J basis vectors, all on the
    same grid and all in the same units as the cube. COEFFS carries one row per
    exposure: the coefficients, their errors, the BERV needed to place the star
    basis, and the MAD-cut flag. VARIANCE carries the weighted power of each
    component, which is what "eigenvalue" means for a weighted PCA where the
    vectors are orthonormal in wavelength but their coefficients need not be.

    The quantity is the *high-passed* log flux, ln f - savgol(ln f), because
    that is what the cube holds and therefore all the fit ever saw. The
    Savitzky-Golay continuum is not stored anywhere, so nothing here can be
    turned back into a flux; see pca2d/reconstruct.py (pca2refs-apply).
    """
    n_star, n_earth = P.shape[0], Q.shape[0]
    primary = fits.PrimaryHDU()
    h = primary.header
    h["INFORMAT"] = ("tfits", "input format the cube was built from")
    h["CONTENT"] = ("hp-log-flux", "ln f - savgol(ln f); no continuum stored")
    h["NSTAR"] = (n_star, "star-frame components")
    h["NEARTH"] = (n_earth, "observer-frame components")
    h["NPARITY"] = (means.shape[0], "means stored, one per order parity")
    h["FRAME"] = (frame, "frame the cube was registered in")
    h["HIGHPASS"] = (highpass, "high-pass mode")
    h["WAVE0"] = (float(grid[0]), "nm, first sample of the magic grid")
    h["DV"] = (float(dv), "km/s per sample; a shift is a translation")
    h["TIED"] = (bool(tied), "parity rows share one coefficient set")
    h["MAXMAD"] = (float(max_mad), "MAD cut applied to the coefficients")
    h["CHI2NULL"] = (float(chi2_null), "weighted var of the median-subtracted cube")
    h["CHI2BEST"] = (float(chi2_best), "weighted var of the residual")
    ratio = float(np.sqrt(max(chi2_best / chi2_null, 0.0)))
    h["RMSRATIO"] = (ratio, "residual RMS / median-subtracted RMS")
    h["CONSTRAI"] = ("n_spectra/n_rows >= %.2f" % min_fraction,
                     "definition of the `constrained` flag")
    h["R2"] = (1.0 - chi2_best / chi2_null, "on the median-subtracted cube")

    cols = [fits.Column(name="wavelength", format="D", unit="nm", array=grid),
            fits.Column(name="template", format="D", array=template)]
    names = PARITY_NAMES if means.shape[0] == 2 else ["all"]
    for i, nm in enumerate(names[:means.shape[0]]):
        cols.append(fits.Column(name="mean_%s" % nm, format="D", array=means[i]))
    for k in range(n_star):
        cols.append(fits.Column(name="star_pc%d" % (k + 1), format="D", array=P[k]))
    for j in range(n_earth):
        cols.append(fits.Column(name="earth_pc%d" % (j + 1), format="D", array=Q[j]))

    # How well each column is determined, and by how much of the data. Without
    # this a reader cannot tell a column carried by 250 spectra from one carried
    # by three, and the vectors look equally trustworthy everywhere. Counted per
    # parity because half the domain is reached by the even orders only.
    if weights is not None:
        live = np.asarray(weights) > 0
        groups = ([0] if parity is None
                  else list(np.unique(np.asarray(parity))))
        total = 0
        for i, value in enumerate(groups):
            rows = (slice(None) if parity is None
                    else np.asarray(parity) == value)
            name = (PARITY_NAMES[i] if len(groups) == 2 and i < 2 else "all")
            count = live[rows].sum(axis=0).astype(np.int32)
            wsum = np.asarray(weights)[rows].sum(axis=0)
            cols.append(fits.Column(name="n_spectra_%s" % name, format="J",
                                    array=count))
            cols.append(fits.Column(name="weight_sum_%s" % name, format="D",
                                    array=wsum))
            total += count
        fraction = total / float(live.shape[0])
        cols.append(fits.Column(name="constrained", format="L",
                                array=fraction >= min_fraction))
    basis = fits.BinTableHDU.from_columns(cols, name="BASIS")

    # one row per exposure, not per cube row: with the parities tied the two
    # rows of an exposure carry byte-identical coefficients, errors and chi2,
    # so keeping both would only invite someone to average a number with itself
    tab = Table(table)
    if tied and "filename" in tab.colnames:
        _, first = np.unique(np.asarray(tab["filename"]), return_index=True)
        tab = tab[np.sort(first)]
    coeffs = fits.BinTableHDU(tab, name="COEFFS")

    power = np.concatenate([power_star, power_earth])
    var = Table()
    var["component"] = (["a%d" % (k + 1) for k in range(n_star)]
                        + ["b%d" % (j + 1) for j in range(n_earth)])
    var["block"] = ["star"] * n_star + ["earth"] * n_earth
    var["weighted_power"] = power
    var["percent_of_median_subtracted"] = 100 * power / chi2_null
    variance = fits.BinTableHDU(var, name="VARIANCE")

    fits.HDUList([primary, basis, coeffs, variance]).writeto(path, overwrite=True)
    log("  wrote %s" % path)


def plot_components(grid, P, Q, power_star, power_earth, chi2_null, path):
    """The basis vectors themselves, star block left, Earth block right.

    The coefficients say when something happened; only this says what. Reading
    the two columns against each other is how a component is identified: a star
    vector should carry stellar lines, an Earth vector should carry telluric
    residuals, and a vector that looks like the other block's is leakage.

    Samples no row covers are stored as exact zeros. They are masked to NaN
    before plotting so the line breaks at a gap instead of diving to zero and
    drawing 622 spurious spikes across the domain.
    """
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25})
    n_rows = max(P.shape[0], Q.shape[0])
    fig, axes = plt.subplots(n_rows, 2, figsize=(13, 1.7 * n_rows + 1.2),
                             sharex=True, squeeze=False)
    blocks = [("star rest frame", P, power_star, "tab:blue", "a"),
              ("observer (Earth) rest frame", Q, power_earth, "tab:red", "b")]
    for col, (title, basis, power, colour, letter) in enumerate(blocks):
        for row in range(n_rows):
            ax = axes[row][col]
            if row >= basis.shape[0]:
                ax.axis("off")
                continue
            v = np.where(basis[row] != 0.0, basis[row], np.nan)
            ax.plot(grid, v, lw=0.25, color=colour)
            ax.axhline(0.0, color="k", lw=0.6, alpha=0.4)
            ax.set_ylabel("%s%d" % (letter, row + 1))
            ax.text(0.995, 0.95, "%.2f%% of the median-subtracted cube"
                    % (100 * power[row] / chi2_null),
                    transform=ax.transAxes, ha="right", va="top", fontsize=7,
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7",
                              lw=0.5, alpha=0.85))
            if row == 0:
                ax.set_title(title)
        axes[-1][col].set_xlabel("wavelength (nm)")
    fig.suptitle("Two-frame PCA: the basis vectors", y=0.998, fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    log("  wrote %s" % path)


def write_variance_table(power_star, power_earth, chi2_null, chi2_best, path):
    """The same numbers as plot_variance, machine-readable."""
    n_star = power_star.size
    power = np.concatenate([power_star, power_earth])
    share = 100 * power / chi2_null
    table = Table()
    table["component"] = (["a%d" % (k + 1) for k in range(n_star)]
                          + ["b%d" % (j + 1) for j in range(power_earth.size)])
    table["block"] = ["star"] * n_star + ["earth"] * power_earth.size
    table["weighted_power"] = power
    table["percent_of_median_subtracted"] = share
    table["cumulative_percent"] = np.cumsum(share)
    table.meta["residual_rms_ratio"] = float(np.sqrt(max(chi2_best / chi2_null, 0.0)))
    table.write(path, format="csv", overwrite=True)
    log("  wrote %s" % path)



def _hierarchical_median(home, berv, bin_kms, min_entries):
    """Median within each BERV bin, then median across bins.

    A plain median over spectra is robust to a minority, and that is exactly
    the assumption that fails here. Anything fixed in the OBSERVER frame, an OH
    airglow residual above all, lands at the same STAR-frame column for every
    exposure sharing a barycentric velocity. On TOI-2120, 23 of 80 exposures sit
    within one resolution element of each other, so 29 per cent of the sample
    votes for the same wrong value and the median is pulled toward it.

    Giving each BERV bin one vote removes that. Measured on TOI-2120 against a
    template built without the dominant group, the excess deviation on
    OH-affected columns falls from 1.46 to 1.20 times the deviation on control
    columns, and the template's own noise rises by 5 per cent.

    This is LBL's approach, not ours: see lbl/recipes/lbl_template.py, "bin cube
    by BERV (to give equal weighting to epochs)", whose BERVBIN_SIZE of 3000 m/s
    and BERVBIN_MIN_ENTRIES of 3 are the defaults used here so that the two
    steps of the chain bin the same way.
    """
    labels = np.floor((berv - np.nanmin(berv)) / bin_kms).astype(int)
    per_bin = []
    for value in np.unique(labels):
        rows = labels == value
        if rows.sum() < min_entries:
            continue
        per_bin.append(np.nanmedian(home[rows], axis=0))
    if len(per_bin) < 2:                       # not enough bins to be worth it
        return np.nanmedian(home, axis=0)
    return np.nanmedian(np.vstack(per_bin), axis=0)

def star_frame_template(data, w, shifter, delta, min_spectra=20, berv=None,
                        berv_bin=None, berv_min_entries=3):
    """Component zero of the star block: the median spectrum in the star frame.

    NOTES.md 11.12: with only an observer-frame mean removed, the star block
    spends its first component (32% of the weighted variance, correlation -0.85)
    turning that mean into this template. Removing it up front hands that
    component back to the fit.

    Median rather than mean so flares and cosmics do not set the template, and
    samples a spectrum does not constrain are masked before the median rather
    than contributing a zero, which would pull the template toward zero exactly
    under the tellurics.
    """
    home = shifter.rows(data, -delta)          # observer -> star
    home_w = shifter.rows(w, -delta)
    home[home_w <= 0] = np.nan
    with warnings.catch_warnings():
        # a t.fits cube has columns no row covers at all (dead orders, and the
        # half of the domain the other parity owns); those are all-NaN by
        # construction and are zeroed by the min_spectra cut below
        warnings.simplefilter("ignore", RuntimeWarning)
        if berv_bin and berv is not None:
            template = _hierarchical_median(home, np.asarray(berv, dtype=float),
                                            float(berv_bin) / 1000.0,
                                            int(berv_min_entries))
        else:
            template = np.nanmedian(home, axis=0)
    count = np.sum(np.isfinite(home), axis=0)
    template = np.where(np.isfinite(template), template, 0.0)
    template[count < min_spectra] = 0.0
    return template


def carry_template(template, shifter, delta):
    """S_n T for every spectrum: the template in each observer frame."""
    return shifter.carry(shifter.prepare(template[None, :]), delta)[:, 0, :]



try:                                    # optional, and worth having
    import bottleneck as _bn
except ImportError:                     # pragma: no cover
    _bn = None


def _nanmedian_over_rows(z):
    """Median over spectra, per wavelength column, ignoring NaN.

    This is the single most expensive operation in an iteration after the carry:
    two calls per clip step over an (n_spectra, 421782) array, and numpy's
    nanmedian sorts every column. Bottleneck's does the same job five to six
    times faster, but only reduces over the LAST axis, so the array is
    transposed first. That copy costs 0.05 s against the 1.9 s it saves, and it
    is exact: measured difference from numpy, 0.

    Falls back to numpy when bottleneck is absent, so the package keeps working
    without it and merely runs slower.
    """
    if _bn is None:
        return np.nanmedian(z, axis=0)
    return _bn.nanmedian(np.ascontiguousarray(z.T), axis=1)


def clip_weights(w0, residual, clip=3.0, min_spectra=20):
    """Iterative soft down-weighting of local >clip-sigma outliers.

    Returns w0 scaled by a factor that is 1 inside the clip and falls as
    (clip/|z|)^2 outside it. That exponent is not arbitrary: it is what you get
    by inflating the variance of an outlying sample until it *is* a clip-sigma
    sample, so a 6-sigma point is kept with a quarter of its weight rather than
    thrown away. Nothing is ever hard-rejected, which matters because a flare is
    real data and a hard cut on residual would sculpt the very variability the
    PCA is here to find.

    "Local" is per wavelength column: z is the residual in units of the
    per-pixel photon sigma, rescaled by a robust MAD over spectra *at that
    wavelength*. So a region the model fits badly does not have its whole
    column clipped, only the spectra that are outliers within it.

    Always applied to the original weights, never compounded, so this is
    iterative re-weighting and not a ratchet that slowly deletes the data.
    """
    good = w0 > 0
    # z = residual * sqrt(w0), built in place: the obvious spelling allocates a
    # sigma array, a z array, a |z - centre| array and a factor array, four
    # N x M float64 temporaries. At M = 92426 that is 2.3 GB of garbage per call
    # and it is what drove the K = J = 10 run into swap.
    z = np.sqrt(w0, dtype=np.float64)
    z *= residual
    z[~good] = np.nan

    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        # columns no row constrains are all-NaN; the min_spectra test below is
        # what handles them, so the median's complaint is noise
        warnings.simplefilter("ignore", RuntimeWarning)
        centre = _nanmedian_over_rows(z)
        z -= centre[None, :]
        np.abs(z, out=z)
        scale = 1.4826 * _nanmedian_over_rows(z)
        enough = np.sum(good, axis=0) >= min_spectra
        fallback = np.nanmedian(scale[enough]) if np.any(enough) else 1.0
    scale = np.where(enough & np.isfinite(scale) & (scale > 0), scale, fallback)
    if not np.isfinite(fallback) or fallback <= 0:
        scale = np.ones_like(scale)

    # z already holds |z - centre|; finish the normalisation in place
    with np.errstate(invalid="ignore", divide="ignore"):
        z /= scale[None, :]
    hit = float(np.count_nonzero(z[good] > clip) / max(np.count_nonzero(good), 1))
    with np.errstate(invalid="ignore", divide="ignore"):
        np.maximum(z, clip, out=z)          # inside the clip -> exactly clip
        z **= -2.0
        z *= clip ** 2                      # -> 1 inside, (clip/|z|)^2 outside
    z[~np.isfinite(z)] = 1.0
    z *= w0
    return z, hit



def coefficient_errors(data, w, P, Q, shifter, delta, a, b, chunk=64,
                       exposure=None):
    """Per-spectrum 1-sigma errors on the coefficients, two flavours.

    The coefficient step is a weighted least squares, so under "the weights are
    inverse variances" the covariance of c_n is just the inverse of the normal
    matrix already formed in joint_coeffs:

        Cov(c_n) = (B_n^T W_n B_n)^-1 = A_n^-1

    That is `sigma_formal`. It is optimistic for two reasons, so a second column
    is reported alongside it:

    1. The weights are NOT purely inverse variances. The telluric ramp of 2.3 and
       the soft clip of 11.14 both multiply w by a distrust factor. The clip is
       self-consistent (it *is* a variance inflation, by construction), the ramp
       is not.
    2. The noise model of 2.2 is reconstructed rather than read -- APERO's eflux
       column is all zeros -- so its absolute normalisation is only as good as
       the header SNR it was tied to. Measured residual/sigma is 1.07, not 1.00.

    `sigma_scaled` multiplies by sqrt(chi2_red) per spectrum, the standard fix,
    which absorbs both. Use it.

    With `exposure` given the covariance is formed on the summed normal matrix
    of the group, and chi2_red on the group's pooled residual and pooled degrees
    of freedom, matching the tied solve in joint_coeffs. Inverting one row's
    normal matrix while the coefficient was fitted to two would overstate the
    error by roughly sqrt(2).

    **Neither includes the uncertainty in P and Q.** The bases are estimated from
    the same 764 spectra, and these errors are conditional on them being exactly
    right. Getting that honestly needs a bootstrap over spectra, refitting the
    bases each time. Expect these to be underestimates.
    """
    n_spectra, n_pixels = data.shape
    n_star, n_earth = P.shape[0], Q.shape[0]
    n_tot = n_star + n_earth
    eye = np.eye(n_tot)
    var = np.full((n_spectra, n_tot), np.nan)
    chi2_red = np.full(n_spectra, np.nan)
    coeffs = np.hstack([a, b])
    Pf = shifter.prepare(P)
    B = np.empty((n_tot, n_pixels))
    tied = exposure is not None
    amats = np.zeros((n_spectra, n_tot, n_tot)) if tied else None
    chi2_rows = np.zeros(n_spectra)
    dof_rows = np.zeros(n_spectra)
    for start in range(0, n_spectra, chunk):
        stop = min(start + chunk, n_spectra)
        SP = shifter.carry(Pf, delta[start:stop])
        for i in range(stop - start):
            n = start + i
            B[:n_star] = SP[i]
            B[n_star:] = Q
            wn = w[n]
            amat = (B * wn) @ B.T
            resid = data[n] - coeffs[n] @ B
            chi2_rows[n] = float(np.sum(wn * resid ** 2, dtype=np.float64))
            dof_rows[n] = int(np.count_nonzero(wn))
            if tied:
                amats[n] = amat
                continue
            trace = np.trace(amat)
            if trace <= 0:
                continue
            try:
                cov = np.linalg.inv(amat + 1e-12 * trace * eye)
            except np.linalg.LinAlgError:
                cov = np.linalg.pinv(amat)
            var[n] = np.clip(np.diag(cov), 0.0, None)
            dof = dof_rows[n] - n_tot
            if dof > 0:
                chi2_red[n] = chi2_rows[n] / dof
    if tied:
        exposure = np.asarray(exposure)
        for value in np.unique(exposure):
            rows = np.where(exposure == value)[0]
            amat = amats[rows].sum(axis=0)
            trace = np.trace(amat)
            if trace <= 0:
                continue
            try:
                cov = np.linalg.inv(amat + 1e-12 * trace * eye)
            except np.linalg.LinAlgError:
                cov = np.linalg.pinv(amat)
            var[rows] = np.clip(np.diag(cov), 0.0, None)
            dof = dof_rows[rows].sum() - n_tot
            if dof > 0:
                chi2_red[rows] = chi2_rows[rows].sum() / dof
    sigma_formal = np.sqrt(var)
    sigma_scaled = sigma_formal * np.sqrt(chi2_red)[:, None]
    return sigma_formal, sigma_scaled, chi2_red


def mad_outliers(coeffs, parity, threshold, already=None):
    """Rows whose coefficient on ANY component sits too far from the others.

    Robust z-score, |c - median| / (1.4826 MAD), component by component; a row
    is flagged if any single component exceeds the threshold.

    THE THRESHOLD IS IN SIGMAS, NOT IN MADS. The 1.4826 puts the MAD on the
    scale of a standard deviation for Gaussian data, so `threshold=10` rejects
    beyond 10 robust sigmas, which is 14.83 plain MADs. The option was called
    --max-mad for a while and that name was simply wrong; it survives as an
    alias for --max-sigma so that existing configs keep their exact behaviour.
    Nothing here reports a "MAD" that carries the factor: where a MAD is
    reported, in the RV compilation and the run summaries, it is
    median(|x - median x|) with no scaling.

    Computed within each parity group rather than over the pooled rows. The two
    parities are different populations here -- the star block's coefficients are
    very nearly pure even-minus-odd -- and pooling them makes each component's
    distribution bimodal, which inflates the MAD by roughly the separation of
    the two modes and hides every real outlier inside it. With an s1d cube there
    is one group and this reduces to the plain cut.
    """
    n_rows = coeffs.shape[0]
    flag = np.zeros(n_rows, dtype=bool)
    if threshold <= 0:
        return flag
    live = np.ones(n_rows, dtype=bool) if already is None else ~np.asarray(already)
    for value in np.unique(parity):
        rows = (parity == value) & live
        if rows.sum() < 5:                 # too few for a robust scale
            continue
        block = coeffs[rows]
        centre = np.median(block, axis=0)
        mad = 1.4826 * np.median(np.abs(block - centre), axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            z = np.abs(block - centre) / np.where(mad > 0, mad, np.inf)
        flag[np.where(rows)[0][z.max(axis=1) > threshold]] = True
    return flag


def correlate_and_plot(outdir, a, b, labels, values, title=""):
    """Write correlations.pdf and correlations.csv for one fit.

    This runs on every fit rather than living in a diagnostic script somebody
    has to remember. A component carries no label: that a1 varies on 83 days is
    a result only for as long as it is not simply tracking the water column, and
    the cheapest way to keep that honest is to draw the matrix every time.
    """
    from .plotting import (plot_correlations, rank_correlations,
                           report_correlations, write_correlation_table)

    comps = [("a%d" % (k + 1), a[:, k]) for k in range(a.shape[1])]
    comps += [("b%d" % (j + 1), b[:, j]) for j in range(b.shape[1])]
    if not labels:
        log("  no ancillary quantity available, skipping correlations.pdf")
        return
    rho = rank_correlations(comps, np.asarray(values, dtype=float))
    plot_correlations(rho, comps, list(labels),
                      os.path.join(outdir, "correlations.pdf"), title=title)
    write_correlation_table(rho, comps, list(labels),
                            os.path.join(outdir, "correlations.csv"))
    report_correlations(rho, comps, list(labels))


def replot(outdir, cube=None):
    """Regenerate every figure from a saved fit.npz, without touching the cube.

    fit.npz has always been advertised as "re-plot without refitting"; this is
    the switch that makes good on it. It buys more than convenience: the figures
    are written only after the last iteration, so until now a half-hour fit that
    produced a plot nobody liked had to be run again from the top to change a
    label. Everything plot_coeffs and plot_variance need is in the archive.
    """
    path = os.path.join(outdir, "fit.npz")
    fit = np.load(path)
    n_star = fit["a"].shape[1]
    log("re-plotting %s: %d star + %d Earth components, %d spectra, iterate %d"
          % (path, n_star, fit["b"].shape[1], fit["bjd"].size, fit["best_iter"]))
    sigma = fit["sigma_scaled"]
    chi2_null = float(fit["chi2_null"])
    chi2_best = float(fit["chi2_best"]) if "chi2_best" in fit.files else np.nan
    keep = ~fit["rejected"] if "rejected" in fit.files else np.ones(
        fit["bjd"].size, dtype=bool)
    if not keep.all():
        log("  %d of %d spectra were rejected by the MAD cut and are not"
              " plotted" % (int((~keep).sum()), keep.size))
    plot_coeffs(fit["bjd"][keep], fit["a"][keep], fit["b"][keep],
                fit["power_star"], fit["power_earth"], chi2_null,
                os.path.join(outdir, "coefficients_vs_time.pdf"),
                err_a=sigma[keep][:, :n_star], err_b=sigma[keep][:, n_star:])
    plot_variance(fit["power_star"], fit["power_earth"], chi2_null, chi2_best,
                  os.path.join(outdir, "variance.pdf"))
    write_variance_table(fit["power_star"], fit["power_earth"], chi2_null,
                         chi2_best, os.path.join(outdir, "variance.csv"))
    if "grid" in fit.files:
        plot_components(fit["grid"], fit["P"], fit["Q"], fit["power_star"],
                        fit["power_earth"], chi2_null,
                        os.path.join(outdir, "components.pdf"))
    else:
        log("  no grid in the archive, skipping components.pdf")
    # the ancillary quantities are stored in the archive, so the correlation
    # matrix redraws without the cube or the source spectra being present
    if "anc_values" in fit.files and fit["anc_values"].size:
        correlate_and_plot(outdir, fit["a"][keep], fit["b"][keep],
                           [str(x) for x in fit["anc_labels"]],
                           fit["anc_values"][:, keep],
                           title="%d exposures" % int(keep.sum()))
    elif cube and "filename" in fit.files:
        # an archive written before the ancillary table was stored: rebuild it
        # from the cube, matching on filename so a rejected row cannot shift the
        # columns against the coefficients
        from astropy.table import Table
        from .plotting import ancillary_table, source_directory
        meta_path = os.path.join(cube, "meta.fits")
        if not os.path.exists(meta_path):
            log("  no %s, skipping correlations.pdf" % meta_path)
            return
        names = [os.path.basename(str(v)) for v in fit["filename"]]
        labels, values = ancillary_table(Table.read(meta_path), names,
                                         source_dir=source_directory(cube))
        correlate_and_plot(outdir, fit["a"][keep], fit["b"][keep], labels,
                           values[:, keep] if len(labels) else values,
                           title="%d exposures" % int(keep.sum()))
    else:
        log("  no ancillary quantities in the archive and no --cube given,"
              " skipping correlations.pdf")


def main(argv=None):
    args = parse_args(argv)
    if args.replot:
        replot(args.outdir, args.cube)
        return None                  # exit status, not a value; see below
    # float32 halves the two (rows, pixels) arrays the fit lives in, and the
    # two temporaries the model costs with them. On 321 individual exposures of
    # TOI-2120 that is the difference between 11.9 GB and 6 GB, which is the
    # difference between fitting every file and coadding nights. Every weighted
    # sum still accumulates in float64: a float32 accumulator over 3.7e8 terms
    # loses the low bits of the very totals the fit is judged by.
    grid, data, w, meta = load_cube(args.cube, min_snr_frac=args.min_snr_frac,
                                    dtype=np.dtype(args.dtype))
    log("fit arrays in %s: %.1f GB for the data and the weights together"
          % (data.dtype.name, 2 * data.nbytes / 1e9))
    n_spectra, n_pixels = data.shape
    parity = row_parity(meta, n_spectra)
    n_parities = np.unique(parity).size
    if n_parities > 1:
        log("t.fits cube: %d rows = %d exposures x %d order parities, two"
              " measurements of the same photons at two places on the detector,"
              " carrying the same BERV." % (n_spectra, n_spectra // n_parities,
                                            n_parities))
    tie = row_exposure(meta, n_spectra) if args.tie_parities else None
    if tie is not None and n_parities > 1:
        log("tying the %d parity rows of each exposure to one set of"
              " coefficients (--no-tie-parities to undo)" % n_parities)
    dv = float(np.median(np.diff(np.log(grid))) * C_KMS)
    # delta is the STAR -> OBSERVER shift in samples, i.e. what S_n applies.
    #
    # APERO's s1d_v wavelengths are in the observer frame, and a stellar feature
    # sits at lambda_bary = lambda_obs * (1 + BERV/c). On a log-uniform grid that
    # is a translation of +BERV/dv samples, so observer -> star is +BERV/dv and
    # the inverse, star -> observer, is **minus** BERV/dv.
    #
    # This had the wrong sign until 2026-08-24 and the error was silent: the fit
    # still converged, to a "star" block anchored to a mirror frame moving at
    # -BERV, which corresponds to nothing physical. See NOTES.md 11.11. Verified
    # against the pipeline's own barycentric cube: rows(observer, +BERV/dv)
    # reproduces it to 1.4% of the data rms, rows(observer, -BERV/dv) does not.
    delta = -np.asarray(meta["berv"], dtype=np.float64) / dv
    shifter = SHIFTERS[args.shift](
        n_pixels, a=args.kernel_halfwidth,
        max_shift=int(np.ceil(np.abs(delta).max())) + 2,
    )
    if args.gap_guard and getattr(shifter, "banded", False):
        gap_guard(w, delta, args.kernel_halfwidth)
        data = np.where(w > 0, data, 0.0)
    elif args.gap_guard:
        log("gap guard skipped: it is defined for a compactly supported"
              " kernel, and %s is not one" % args.shift)
    n_max = max(args.n_star, args.n_earth)
    # The carried basis is (chunk, K, n_pixels) float64 and is the single largest
    # transient in the run: at K=5 and 421782 columns one spectrum costs 17 MB,
    # so the chunk is sized to keep that buffer near 128 MB rather than to any
    # round number. Bigger is not faster; the loop is bandwidth-bound, which is
    # also why threading carry() helps and why the tap loop was rewritten as a
    # single contraction over a sliding-window view.
    chunk = args.chunk or max(4, int(128e6 / (8 * n_max * n_pixels)))
    log("carry chunk %d spectra -> %.0f MB per (chunk, K, M) buffer"
          % (chunk, 8 * chunk * n_max * n_pixels / 1e6))
    log("cube %s: N=%d M=%d dv=%.3f km/s  shifts %.1f..%.1f pix"
          % (os.path.basename(args.cube), n_spectra, n_pixels, dv,
             delta.min(), delta.max()))
    log("shift operator: %s%s"
          % (args.shift, " (a=%d)" % args.kernel_halfwidth
             if args.shift == "lanczos" else ""))

    # Everything is scored against the raw high-passed cube so that the numbers
    # are comparable across runs and against the template baseline of 11.12.
    chi2_raw_rows = np.sum(w * data ** 2, axis=1, dtype=np.float64)
    chi2_raw = float(chi2_raw_rows.sum())
    w0 = w.copy()          # modified only when the MAD cut retires a spectrum

    if args.template:
        template = star_frame_template(
            data, w, shifter, delta,
            berv=np.asarray(meta["berv"], dtype=float),
            berv_bin=args.template_berv_bin,
            berv_min_entries=args.template_berv_min_entries)
        template_model = carry_template(template, shifter, delta)
        data = data - template_model
        data[w <= 0] = 0.0
        log("subtracted the star-frame median template: %.4f of the raw"
              " weighted variance left" % (float(np.sum(w * data ** 2, dtype=np.float64)) / chi2_raw))
    else:
        template = np.zeros(n_pixels)
        template_model = None

    # component zero of the *Earth* block, on whatever the template left behind.
    # One mean per parity, not one overall: a t.fits cube holds two rows per
    # exposure and they carry a static even-minus-odd offset, fixed in the
    # instrument frame, which a single mean would leave in the residual for a
    # component to spend itself on (NOTES.md 13.9). With an s1d cube there is
    # one group and this is the old code path exactly.
    means, parity_groups = parity_means(data, w, parity)
    mean = means[0] if means.shape[0] == 1 else np.average(
        means, axis=0, weights=[np.sum(parity == g) for g in parity_groups])
    group = np.searchsorted(parity_groups, parity)
    # What that mean holds splits in two, and only one half should come out
    # here. The even-minus-odd offset is instrumental and has to go, for the
    # reason above. The part the two parities SHARE is whatever is static in
    # the observer's frame, which is the atmosphere: the airglow emission that
    # sits at the same wavelength every night lands in it at full amplitude.
    # Subtracting that before the fit hides it from the components, and since
    # `reconstruct.correction_on_grid` removes the components and never the
    # means, it is then never removed from anything. Leaving it in costs one
    # observer component and lets the correction reach it.
    common = means.mean(axis=0)
    if args.mean == "offset":
        means = means - common
    subtract_means(data, means, group)
    data[w <= 0] = 0.0
    if args.mean == "offset":
        live_any = w.sum(axis=0) > 0
        log("left the shared part of the observer-frame mean IN the data,"
              " %.4g rms, for the Earth block to describe"
              % float(np.std(common[live_any])))
    if means.shape[0] > 1:
        both = np.ones(n_pixels, dtype=bool)
        for i, value in enumerate(parity_groups):
            both &= w[parity == value].sum(axis=0) > 0
        log("subtracted one mean per parity; the even/odd offset it removed is"
              " %.4g rms over the %d columns both cover"
              % (float(np.std((means[0] - means[1])[both])) if both.any() else np.nan,
                 int(both.sum())))
    chi2_null = float(np.sum(w * data ** 2, dtype=np.float64))
    log("after the observer-frame mean as well:      %.4f of the raw"
          " weighted variance left" % (chi2_null / chi2_raw))
    if args.clip > 0:
        log("soft clip at %.1f sigma, local per wavelength column" % args.clip)

    # Random orthonormal start. Orthonormal because the update steps
    # re-orthonormalise anyway and starting inside the constraint set avoids a
    # first iteration spent getting there; random because the fit converges in a
    # dozen iterations from anywhere, so a cleverer seed buys nothing and would
    # tie the result to whatever produced it. The seed is fixed, so a rerun with
    # the same cube and options reproduces the same basis.
    rng = np.random.default_rng(0)
    P = np.linalg.qr(rng.normal(size=(n_pixels, args.n_star)))[0].T.copy()
    Q = np.linalg.qr(rng.normal(size=(n_pixels, args.n_earth)))[0].T.copy()

    # This alternation is NOT guaranteed to decrease chi2, and on this dataset it
    # does not: it peaks around iteration 6 and then degrades while the condition
    # number runs away (see NOTES.md 11.9). Two approximations are responsible --
    # the diagonal stand-in for S_n^T W_n S_n in update_star, and the
    # re-orthonormalisation after each M-step. Until the star update is made an
    # exact minimiser, keep the best iterate rather than the last one, and stop
    # once it has clearly turned over.
    rejected = np.zeros(n_spectra, dtype=bool)
    for mad_round in range(int(args.max_mad_rounds) + 1):
        best = None
        worse = 0
        hit = 0.0
        a_prev = b_prev = None      # carried between iterations, see the clip step
        for iteration in range(args.iters):
            # Re-weight before fitting, using the model as it stands. Iteration 0 has
            # no model yet, so the residual is the data itself, which is the right
            # thing: it catches cosmics and flares against the template before the
            # components ever see them.
            if args.clip > 0:
                # The coefficients this needs are the ones the previous
                # iteration ended on: same P, same Q, same w, so solving again
                # returns the same numbers. Reusing them removes one of the four
                # coefficient solves per iteration, which is the single most
                # expensive call in the loop. Only iteration 0 has nothing to
                # reuse.
                if a_prev is None:
                    a_prev, b_prev, _ = joint_coeffs(data, w, P, Q, shifter,
                                                     delta, chunk=chunk,
                                                     exposure=tie)
                model_prev = star_model(P, a_prev, shifter, delta, n_spectra, n_pixels, chunk) \
                    + b_prev @ Q
                w, hit = clip_weights(w0, data - model_prev, clip=args.clip)

            # ---- E-step: coefficients for both blocks at fixed bases --------
            # Solved JOINTLY, never one block then the other. The off-diagonal
            # part of the normal matrix is the overlap between the carried star
            # basis and the fixed Earth basis, and it is exactly what lets the
            # fit attribute a feature to one frame rather than the other. Solve
            # the blocks separately and each one claims whatever it can reach.
            tick = time.time()
            a, b, cond = joint_coeffs(data, w, P, Q, shifter, delta, chunk=chunk,
                                 exposure=tie)
            t_coeff = time.time() - tick

            # ---- M-step: one basis, then re-solve, then the other -----------
            # The re-solve in the middle is not optional. update_earth needs the
            # coefficients that go with the basis update_star has just produced;
            # feeding it the stale ones makes the second update fit a residual
            # that no longer exists, and the alternation stops descending.
            tick = time.time()
            if args.order == "star_first":
                P = update_star(data, w, P, Q, a, b, shifter, delta, chunk)
                a, b, _ = joint_coeffs(data, w, P, Q, shifter, delta, chunk=chunk,
                                 exposure=tie)
                Q = update_earth(data, w, P, Q, a, b, shifter, delta, chunk)
            else:
                Q = update_earth(data, w, P, Q, a, b, shifter, delta, chunk)
                a, b, _ = joint_coeffs(data, w, P, Q, shifter, delta, chunk=chunk,
                                 exposure=tie)
                P = update_star(data, w, P, Q, a, b, shifter, delta, chunk)
            t_basis = time.time() - tick

            a, b, cond = joint_coeffs(data, w, P, Q, shifter, delta, chunk=chunk,
                                 exposure=tie)
            # these are what the next iteration's clip step would recompute
            a_prev, b_prev = a, b
            model = star_model(P, a, shifter, delta, n_spectra, n_pixels, chunk) + b @ Q
            # scored on w0: the clip weights change every iteration, so scoring on
            # them would be a moving target and the monotonicity check meaningless
            chi2 = float(np.sum(w0 * (data - model) ** 2, dtype=np.float64))
            flag = ""
            if best is None or chi2 < best[0]:
                best = (chi2, P.copy(), Q.copy(), iteration)
                worse = 0
            else:
                worse += 1
                flag = "   <- worse than iter %d" % best[3]
            # `data` here is already template- and mean-subtracted, so this chi2 is
            # the residual of the *full* model and can be quoted against chi2_raw.
            log("  iter %d  R2=%.6f  left/raw=%.4f  clipped=%.3f%%"
                  "  cond(A) med/max %.1f/%.1f  [%.1fs coeff, %.1fs bases]%s"
                  % (iteration, 1 - chi2 / chi2_null, chi2 / chi2_raw, 100 * hit,
                     np.median(cond), cond.max(), t_coeff, t_basis, flag))
            if worse >= 2:
                log("  chi2 has turned over; stopping and keeping iteration %d" % best[3])
                break

        chi2, P, Q, best_iter = best
        log("  best iterate: %d, R2 = %.6f, %.4f of the raw weighted variance"
              " left" % (best_iter, 1 - chi2 / chi2_null, chi2 / chi2_raw))

        # ---- MAD cut: retire a spectrum the components disagree about ------
        # After the fit, not before: an outlier is only definable against the
        # other spectra's coefficients, which is what the fit just produced.
        # The rejected rows keep their place in every array and simply lose
        # their weight, so nothing downstream has to be re-indexed and
        # joint_coeffs leaves their coefficients at zero of its own accord.
        if args.max_mad <= 0 or mad_round == int(args.max_mad_rounds):
            break
        a_now, b_now, _ = joint_coeffs(data, w, P, Q, shifter, delta, chunk=chunk,
                                 exposure=tie)
        fresh = mad_outliers(np.hstack([a_now, b_now]), parity, args.max_mad,
                             rejected)
        if not fresh.any():
            log("  MAD cut at %.1f sigma: nothing left to reject" % args.max_mad)
            break
        rejected |= fresh
        log("  MAD cut at %.1f sigma: rejecting %d spectra this round, %d of %d"
              " in total (%.1f%%); refitting"
              % (args.max_mad, int(fresh.sum()), int(rejected.sum()), n_spectra,
                 100 * rejected.mean()))
        w0[fresh] = 0.0
        w = w0.copy()
        data[fresh] = 0.0
        # the per-parity mean moved when those rows left; take out the
        # difference, and in offset mode only the part of that difference that
        # separates the parities, so the shared term stays in as before
        residual_means, _ = parity_means(data, w0, parity)
        if args.mean == "offset":
            residual_means = residual_means - residual_means.mean(axis=0)
        subtract_means(data, residual_means, group)
        data[w0 <= 0] = 0.0
        means = means + residual_means
        chi2_raw = float(chi2_raw_rows[~rejected].sum())
        chi2_null = float(np.sum(w0 * data ** 2, dtype=np.float64))

    chi2, P, Q, best_iter = best
    if rejected.any():
        log("  MAD cut kept %d of %d spectra (%d rejected at %.1f sigma)"
              % (n_spectra - int(rejected.sum()), n_spectra,
                 int(rejected.sum()), args.max_mad))

    if args.leakage:
        to_earth, to_star = leakage(data, w, P, Q, shifter, delta, chunk)
        log("  leakage star -> Earth block: %s"
            % np.array2string(np.round(to_earth, 4)), "value")
        log("  leakage Earth -> star block: %s"
            % np.array2string(np.round(to_star, 4)), "value")

    # ------------------------------------------------------- outputs -------
    a, b, _ = joint_coeffs(data, w, P, Q, shifter, delta, chunk=chunk,
                                 exposure=tie)
    Pf = shifter.prepare(P)
    power_star = block_power(w, a, lambda s0, s1: shifter.carry(Pf, delta[s0:s1]))
    power_earth = block_power(w, b, lambda s0, s1: Q)
    P, a, power_star = tidy_block(P, a, power_star)
    Q, b, power_earth = tidy_block(Q, b, power_earth)

    sig_formal, sig_scaled, chi2_red = coefficient_errors(
        data, w0, P, Q, shifter, delta, a, b, chunk, exposure=tie)
    n_star = a.shape[1]
    log("  median reduced chi2 per spectrum: %.3f   (1.0 would mean the noise"
          " model of 2.2 is exactly right)" % np.nanmedian(chi2_red))
    log("  median error scaling sqrt(chi2_red): %.3f"
          % np.nanmedian(np.sqrt(chi2_red)))
    with np.errstate(invalid="ignore", divide="ignore"):
        snr_a = np.nanmedian(np.abs(a) / sig_scaled[:, :n_star], axis=0)
        snr_b = np.nanmedian(np.abs(b) / sig_scaled[:, n_star:], axis=0)
    log("  median |coefficient| / sigma, star block : %s"
          % np.array2string(snr_a, precision=1))
    log("  median |coefficient| / sigma, Earth block: %s"
          % np.array2string(snr_b, precision=1))

    os.makedirs(args.outdir, exist_ok=True)
    bjd = np.asarray(meta["bjd"], dtype=float)
    table = Table()
    table["filename"] = meta["filename"]
    table["bjd"] = bjd
    table["berv"] = np.asarray(meta["berv"], dtype=float)
    table["airmass"] = np.asarray(meta["airmass"], dtype=float)
    table["snr_band"] = np.asarray(meta["snr_band"], dtype=float)
    table["chi2_red"] = chi2_red
    table["rejected"] = rejected      # MAD cut; their coefficients are zeros
    for k in range(a.shape[1]):
        table["a%d" % (k + 1)] = a[:, k]              # star-frame block
        table["ea%d" % (k + 1)] = sig_scaled[:, k]    # 1-sigma, chi2-scaled
    for j in range(b.shape[1]):
        table["b%d" % (j + 1)] = b[:, j]              # observer-frame block
        table["eb%d" % (j + 1)] = sig_scaled[:, n_star + j]
    csv_path = os.path.join(args.outdir, "coefficients.csv")
    table.write(csv_path, format="csv", overwrite=True)
    log("  wrote %s" % csv_path)

    log("  star  block weighted variance share: %s"
          % np.array2string(100 * power_star / chi2_null, precision=2))
    log("  Earth block weighted variance share: %s"
          % np.array2string(100 * power_earth / chi2_null, precision=2))
    # the ancillary quantities, gathered once and stored in the archive so that
    # --replot can redraw the correlation matrix on its own
    from .plotting import ancillary_table, source_directory
    anc_labels, anc_values = ancillary_table(
        meta, [os.path.basename(str(v)) for v in meta["filename"]],
        source_dir=source_directory(args.cube))

    npz_path = os.path.join(args.outdir, "fit.npz")
    np.savez_compressed(npz_path, bjd=bjd, a=a, b=b, P=P, Q=Q,
                        anc_labels=np.asarray(anc_labels, dtype="U32"),
                        anc_values=np.asarray(anc_values, dtype=float),
                        power_star=power_star, power_earth=power_earth,
                        chi2_null=chi2_null, chi2_raw=chi2_raw,
                        best_iter=best_iter, template=template, mean=mean,
                        means=means, parity=parity,
                        used_template=bool(args.template),
                        sigma_formal=sig_formal, sigma_scaled=sig_scaled,
                        chi2_red=chi2_red, rejected=rejected,
                        max_mad=float(args.max_mad), chi2_best=chi2, grid=grid,
                        tied=bool(args.tie_parities), dv=dv,
                        filename=np.asarray(meta["filename"], dtype="U64"),
                        berv=np.asarray(meta["berv"], dtype=float))
    log("  wrote %s (re-plot without refitting)" % npz_path)

    # the rejected rows carry zero weight, so their coefficients are zeros and
    # plotting them would draw a spurious line on the axis
    keep = ~rejected
    plot_coeffs(bjd[keep], a[keep], b[keep], power_star, power_earth, chi2_null,
                os.path.join(args.outdir, "coefficients_vs_time.pdf"),
                err_a=sig_scaled[keep][:, :n_star],
                err_b=sig_scaled[keep][:, n_star:])
    plot_variance(power_star, power_earth, chi2_null, chi2,
                  os.path.join(args.outdir, "variance.pdf"))
    write_variance_table(power_star, power_earth, chi2_null, chi2,
                         os.path.join(args.outdir, "variance.csv"))
    plot_components(grid, P, Q, power_star, power_earth, chi2_null,
                    os.path.join(args.outdir, "components.pdf"))
    correlate_and_plot(args.outdir, a[keep], b[keep], anc_labels,
                       anc_values[:, keep] if len(anc_labels) else anc_values,
                       title="%d exposures" % int(keep.sum()))
    write_components_fits(os.path.join(args.outdir, "twoframe_components.fits"),
                          grid, P, Q, template, means, power_star, power_earth,
                          chi2_null, chi2, table, dv, tied=bool(args.tie_parities),
                          max_mad=float(args.max_mad),
                          frame="observer", highpass="log_sub",
                          weights=w0, parity=parity)
    # Nothing, on purpose. This is a console_script entry point, so whatever it
    # returns is handed to sys.exit: returning the two bases printed a pair of
    # arrays and exited 1 on a successful fit, which is why nothing could be
    # chained after it and why a `set -e` script died at the first one.
    return None


if __name__ == "__main__":
    main()
