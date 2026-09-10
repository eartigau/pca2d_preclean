"""YAML configuration handling with defaults and a reproducibility hash."""

from __future__ import annotations

import copy
import hashlib
import json
import math
import os
import re

import yaml


class _ConfigLoader(yaml.SafeLoader):
    """yaml.SafeLoader without YAML 1.1's base-60 numbers.

    PyYAML follows YAML 1.1, where digits separated by a colon are a number in
    base 60: an unquoted `1267:2` is read as 1267 * 60 + 2 = 76022, and even
    `2450:0.5` as 147000.5. The windows in output.windows are written exactly
    that way, centre:width in nm, so three of the eight became numbers in the
    tens of thousands, landed outside the grid, and were skipped by the figure
    scripts without a word. Only the ones whose centre had a decimal point
    survived, because `1200.3:2` is not a base-60 pattern.

    Everything else is SafeLoader: the same ints, floats, booleans and nulls,
    so no other value in a config changes type.
    """


_ConfigLoader.yaml_implicit_resolvers = {
    first: [(tag, rx) for tag, rx in entries
            if tag not in ("tag:yaml.org,2002:int", "tag:yaml.org,2002:float")]
    for first, entries in yaml.SafeLoader.yaml_implicit_resolvers.items()}
# PyYAML 6's own two patterns, each minus its base-60 alternative
_ConfigLoader.add_implicit_resolver(
    "tag:yaml.org,2002:int",
    re.compile(r"""^(?:[-+]?0b[0-1_]+
                    |[-+]?0[0-7_]+
                    |[-+]?(?:0|[1-9][0-9_]*)
                    |[-+]?0x[0-9a-fA-F_]+)$""", re.X),
    list("-+0123456789"))
_ConfigLoader.add_implicit_resolver(
    "tag:yaml.org,2002:float",
    re.compile(r"""^(?:[-+]?(?:[0-9][0-9_]*)\.[0-9_]*(?:[eE][-+][0-9]+)?
                    |\.[0-9][0-9_]*(?:[eE][-+][0-9]+)?
                    |[-+]?\.(?:inf|Inf|INF)
                    |\.(?:nan|NaN|NAN))$""", re.X),
    list("-+0123456789."))


def read_yaml(stream):
    """A YAML document, read the way this package's configs mean it."""
    return yaml.load(stream, Loader=_ConfigLoader)


def check_windows(windows, domain):
    """Say out loud about any window no figure can be drawn in.

    The figure scripts skip such a window in a subprocess whose output is only
    shown when it fails, so a window that could not be drawn used to vanish from
    the PDF with nothing said. A number where a centre:width was meant is the
    signature of the base-60 reading above, from a config read some other way.
    """
    from .logger import log

    lo, hi = float(domain["wave_min"]), float(domain["wave_max"])
    for spec in windows or []:
        centre, _, width = str(spec).partition(":")
        try:
            c, w = float(centre), float(width or 2.0)
        except ValueError:
            log("output.windows: %r is not centre:width in nm, so it is not"
                " drawn" % (spec,), "warn")
            continue
        if not lo <= c <= hi:
            hint = ("; a number here is what YAML 1.1 makes of an unquoted"
                    " centre:width, 1267:2 being read as 76022"
                    if isinstance(spec, (int, float)) else "")
            log("output.windows: %s is centred at %.1f nm, outside the domain"
                " %.1f-%.1f nm, so it is not drawn%s" % (spec, c, lo, hi, hint),
                "warn")

# Every tunable lives here; config.yaml only needs to override what differs.
DEFAULTS = {
    "input": {
        "format": "tfits",           # 'tfits' = APERO order-by-order spectra
                                     #           (the nominal input)
                                     # 's1d'   = the resampled s1d_v products
        # The input ROOT, and not the folder of spectra: one run reads
        # <directory>/<object>/, so the same config serves every target on the
        # disk and a run needs nothing but a name. See spectra_dir().
        "directory": "data",
        "pattern": "*t.fits",
        "max_files": None,           # None = all; small int for quick tests
        "object": None,              # keep only this OBJECT (None = no filter)
        # "auto" (the default), True or False. Coadding a night is not a
        # modelling choice, it exists so the fit holds in memory, and fitting
        # the individual spectra is strictly more information. "auto" weighs
        # the fit's footprint against output.max_memory_gb, a config value and
        # not the machine's RAM, so that the decision stays a pure function of
        # the configuration the cube cache key hashes. See cube.resolve_stacking.
        "nightly_stack": "auto",     # coadd the exposures of each night before
                                     # the PCA. Every exposure is still
                                     # registered with its own BERV first, so
                                     # nothing is smeared; NIRPS takes 3-4
                                     # exposures a visit, so this makes the cube
                                     # (and the decomposition) 3-4 times smaller
        # --- 'tfits' only -------------------------------------------------
        "orders": None,              # None = every order; or an explicit list
        "min_order_finite_fraction": 0.05,   # skip an order with fewer usable
                                     # pixels than this. NIRPS loses orders
                                     # 44-45 and 71-74 entirely to water; that
                                     # is expected, not an error.
        "blaze_min_frac": 0.05,      # reject pixels where the blaze has fallen
                                     # below this fraction of its order peak
    },
    "domain": {
        "wave_min": 965.0,           # nm
        "wave_max": 1950.0,          # nm
        "grid_source": "magic",      # 'magic'  = build the log-uniform grid from
                                     #            wave0 + dv below (the only
                                     #            option that works for t.fits,
                                     #            which carry no common grid)
                                     # 'native' = recycle the grid stored in the
                                     #            first s1d_v file
        "wave0": 965.0,              # nm, anchor of the magic grid
        # Where dv comes from. False: the number below. True: the finest pixel
        # step in the first spectrum, times SMART_DV_FRACTION, and THE NUMBER
        # BELOW IS NOT READ AT ALL. A float here is that fraction. The measured
        # step is kept as domain.pixel_dv for the record.
        "smart_dv": False,
        # km/s per pixel (500 m/s). Read only when smart_dv is false, and may
        # be null when it is true, since nothing would look at it.
        "dv": 0.5,
        # Nominal photometric bands (MKO half-power points, nm). The PCA is
        # fitted on these and only these; see pca.fit_bands.
        "bands": {
            "Y": [970.0, 1070.0],
            "J": [1170.0, 1330.0],
            "H": [1490.0, 1780.0],
            "K": [2030.0, 2370.0],
        },
    },
    "registration": {
        "frame": "barycentric",      # 'barycentric' (stellar rest frame) or
                                     # 'observer' (telluric rest frame, no shift)
        "target_berv": 0.0,          # km/s; frame the spectra are registered to
        "spline_order": 3,           # cubic spline; the BERV shift is fractional
        "mask_threshold": 0.999,     # good-pixel fraction required after resampling
    },
    "highpass": {
        "method": "savgol",          # 'savgol' or 'none'
        "window": 31,                # pixels (odd); 15.5 km/s at dv = 0.5 km/s,
                                     # about 4 NIRPS resolution elements
        "polyorder": 2,
        "mode": "divide",            # 'divide'  -> y = ln(f / lowpass(f))
                                     # 'log_sub' -> y = ln(f) - lowpass(ln f)
        "frame": "observer",         # s1d only: apply before ('observer') or
                                     # after ('target') the BERV registration.
                                     # t.fits always filters on the destination
                                     # grid, per order -- the native sampling is
                                     # not uniform in velocity.
    },
    "quality": {
        "min_snr": 1.0,              # reject if the band SNR is below this
        "require_positive_flux": True,   # reject frames with median flux <= 0
        "max_nan_fraction": 0.5,     # reject if more of the band than this is NaN
        # Epoch window, as rjd = jd - 2400000. Both default to null, i.e. keep
        # everything: excluding part of a time series is a claim about the data,
        # not a hygiene default, and it must be made in the config where it is
        # visible. For NIRPS the first season sits before rjd 60200
        # (2023-09-12); on Proxima it carries most of the excess velocity
        # scatter. `quality` is hashed into the cube cache key, so setting
        # either of these rebuilds the cube rather than silently reusing one
        # built over a different window.
        # Discard a sample where the sky emission the pipeline subtracted was
        # more than this many times the stellar flux. null keeps everything.
        # Hashed into the cube cache key only when set.
        "max_sky_ratio": None,
        # Drop a sample that has a masked neighbour within this many samples on
        # BOTH sides: it survived the cuts but its neighbourhood did not, so its
        # high pass and any shift of it were built from samples that are not
        # there. 3 is the width that was asked for. null disables the check.
        # Hashed into the cube cache key only when set.
        "isolated_window": None,
        "min_rjd": None,
        "max_rjd": None,
    },
    "weights": {
        "mode": "photon",            # 'photon' or 'uniform'
        "snr_pixel_scale": 1.0,      # EXTSN is per e2ds pixel; rescale if wanted
        "ln_clip_low": -0.5,         # y below this -> zero weight (deep lines)
        "ln_clip_high": None,        # y above this -> zero weight (None = keep)
        "sigma_floor_frac": 0.0,     # add this (in ln units) in quadrature
        "telluric_mode": "recon",        # 'recon'     = the paired
                                         #               *_s1d_v_recon_A.fits
                                         # 'from_data' = airmass regression
                                         # 'file'      = external model
                                         # 'none'      = no telluric weighting
        "telluric_ramp_zero": 0.5,       # transmission giving zero weight
        "telluric_ramp_one": 1.0,        # transmission giving full weight
        "telluric_file": None,           # for mode 'file': (wavelength_nm, T)
        "telluric_power": 2.0,           # for modes 'from_data' / 'file'
        "telluric_min_transmission": None,
        # Width in pixels of the boxcar used to measure the noise from the
        # sample-to-sample scatter, or null to trust the photon model alone.
        # The larger of the two sigmas is used. Hashed into the cube cache key
        # only when set, so turning it on rebuilds the cube and leaving it off
        # changes no existing key.
        "empirical_noise_box": None,
        "min_good_fraction": 0.2,    # drop grid columns with fewer good spectra
    },
    "output": {
        # The cube's storage dtype is NOT here: it is cube.CUBE_DTYPE, hard
        # float32, for the reasons written beside it.
        # The output ROOT: every product of a run lands in
        # <directory>/<object>/<M>-<N>/, corrected spectra included, so that
        # nothing is ever written beside the input files.
        "directory": "outputs",
        "tag": "run",
        # Rebuildable intermediates, deliberately NOT under the output root:
        # a cube costs twenty minutes and does not depend on where the products
        # of a run are asked to go, so pointing --out-dir somewhere new must
        # not orphan it.
        "cache_directory": "cache",
        # Where a run's products are kept, when that is not the internal disk:
        # each run folder, and LBL's, is a link to its place under it, made
        # before anything is written (storage.py).
        # None here, so that a clone on another machine does not aim at a disk
        # it does not have; config.yaml names the one this campaign uses.
        "fits_directory": None,
        "use_cache": True,
        "max_memory_gb": 12.0,       # refuse to allocate a cube larger than this
                                     # rather than driving the machine to swap
    },
    # ---------------------------------------------------------------- fit ---
    # The two-frame fit. Two component counts, never one: M in the stellar rest
    # frame and N in the observer frame, and they need not be equal.
    "twoframe": {
        # 2 and 7 since 2026-09-09, measured on TOI-2120 with the static sky
        # left in the data for the fit to describe. Against 5+5: the observer
        # block explains 25.8% of the weighted variance instead of 16.5%, its
        # first component 17.4% instead of 9.7%, and the leakage from it into
        # the star block falls from 0.086 to 0.0116. Five star components were
        # enough rope for the star block to help describe an observer-frame
        # line, through the (feature, its derivative x shift) pair; two are not.
        "n_star": 2,
        "n_earth": 7,
        "iters": 16,
        "tie_parities": True,
        # Rejection threshold in ROBUST SIGMAS, 1.4826 x MAD, not in MADs:
        # 10 here is 14.83 plain MADs. `max_sigma` is the name; `max_mad` is
        # accepted as an alias because configs and logs on disk use it and the
        # value must keep meaning exactly what it meant.
        # Hierarchical template: median within each BERV bin, then median
        # across bins, so a clump of exposures sharing a barycentric velocity
        # cannot vote twice. LBL's values, so both steps of the chain bin the
        # same way. null restores the plain median.
        "template_berv_bin": None,        # m/s
        "template_berv_min_entries": 3,
        "max_mad": 10.0,
        "max_mad_rounds": 2,
        "clip": 3.0,
        "min_snr_frac": 0.5,
        # Storage for the two (rows, pixels) arrays the fit lives in. float32
        # halves them and their model temporaries, which is what makes fitting
        # every exposure possible instead of nightly means; the weighted sums
        # accumulate in float64 whatever this says.
        "dtype": "float64",
        # The static part of the model (twoframe --mean). "offset": the
        # one-shot per-parity observer-frame offset of the fit behind the
        # reports, which the correction divides out with the observer block.
        # "iterate", one mean per order parity in EACH frame re-estimated at
        # every sweep with the components centred, separates the frames on a
        # synthetic cube but does not converge yet on TOI2120 (2026-09-10: its
        # star-frame mean moved more at every sweep and chi2 turned over after
        # one), so it is not the default.
        "mean": "offset",
        # The one-shot star-frame median of the older modes; "iterate" makes
        # its own star-frame mean, one per parity, and ignores this.
        "template": False,
        "shift": "lanczos",
        "kernel_halfwidth": 8,
        "gap_guard": True,
        "leakage": True,
        "chunk": None,
        "order": "star_first",
        # One velocity per exposure, fitted beside the components as the
        # derivative of the reconstructed star carried to that exposure's BERV.
        # It exists to keep the star's own motion OUT of the observer block:
        # without it the block describes the shift, and dividing the block out
        # of the flux writes that velocity into the corrected spectrum. Measured
        # on TOI2120 it accounts for 41% of the variance of what the correction
        # did to LBL's velocities. It is fitted and never divided out, since
        # taking it out of the data would remove the signal being looked for.
        "velocity_term": True,
        # The shift is measured only on columns inside a photometric band whose
        # telluric transmission stays above this in 90% of the exposures. Zero
        # turns the transmission cut off and keeps the band cut. Fitting the
        # star's velocity where the flux is mostly atmosphere is fitting it to
        # the wrong thing: unmasked, the term came out 1.8 times larger than
        # the velocity LBL measures on the same spectra.
        "velocity_min_transmission": 0.95,
    },
    # ---------------------------------------------------------------- lbl ---
    # Handing both sets of spectra to LBL, the delivered ones and the corrected
    # ones, as two objects in one LBL tree. See pca2d/lbl.py.
    "lbl": {
        "prepare": True,             # write LBL's config and its run script
        "run": False,                # and run it. Hours, so it is asked for.
        "directory": "lbl",          # LBL's DATA_DIR, its own tree
        # The corrected object's name. {tag} becomes the run's <M>-<N>, and
        # it is in the default because without it two runs at different
        # component counts write their corrected spectra into ONE LBL science
        # folder, where LBL's glob takes the mixture and says nothing.
        "suffix": "_PCA2D_{tag}",
        "before": True,              # the delivered spectra as their own object
        "after": True,               # and the corrected ones as another
        # LBL picks the stellar model its mask comes from by this, and stops
        # without it. 'auto' reads it from the spectra, where APERO writes it
        # as OBJTEMP and PP_TEFF; a number here overrides the header.
        "teff": "auto",
        "template": None,            # OBJECT_COMPARISON; None = each its own
        # The corrected object's template is the fit's first star component
        # at its mean amplitude (lbltemplate.py), written by LBL's own writer
        # where LBL looks for it, so LBL's template step finds it and skips.
        # A template LBL made for that object is never replaced.
        "star_template": True,
        # With two or more star components, the ones past the first go to
        # LBL as RESPROJ tables STRPCA2..N, in the place of its DTEMP
        # gradients, and each exposure's projection on them is an rdb column.
        "strpca": True,
        "steps": ["template", "mask", "compute", "compile"],
        "link": "symlink",           # 'symlink' or 'copy' into LBL's tree
        "input_file": None,          # LBL's glob inside a science folder
        "instrument": None,          # LBL's name for this spectrograph, which
        "data_source": None,         # with the source selects its reader
    },
}


def _deep_update(base: dict, new: dict) -> dict:
    out = copy.deepcopy(base)
    for key, value in (new or {}).items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_update(out[key], value)
        else:
            out[key] = value
    return out


# What `domain.smart_dv: true` means: sample the finest pixel on the detector
# at this fraction of its own step. Below 1 the grid is finer than the data,
# which is the only thing resampling has to guarantee.
SMART_DV_FRACTION = 0.7


def measure_pixel_dv(directory: str, pattern: str = "*t.fits") -> float:
    """The smallest step between adjacent pixels, in km/s, from the first file.

    Read from the wavelength EXTENSION, order by order, and never from header
    polynomials. The minimum is taken over every order, since it is the finest
    pixel anywhere in the domain that decides how finely the common grid has to
    be sampled for no part of the spectrum to be smoothed by the resampling.

    One file, the first that opens, and not a median over many: the answer has
    to be a pure function of the configuration and the data, because it becomes
    `domain.dv` and `domain.dv` is hashed into the cube cache key.
    """
    import glob

    import numpy as np

    from .io import robust_open
    from .grids import velocity
    from .tfits import extensions_for

    paths = sorted(glob.glob(os.path.join(directory, pattern)))
    if not paths:
        raise SystemExit("no file matching %s in %s" % (pattern, directory))
    for path in paths[:20]:
        try:
            with robust_open(path) as hdulist:
                _, e_wave, _, _ = extensions_for(hdulist)
                wave = np.asarray(hdulist[e_wave].data, dtype=np.float64)
        except Exception:                                     # noqa: BLE001
            continue
        wave = np.atleast_2d(wave)
        with np.errstate(invalid="ignore", divide="ignore"):
            step = velocity(np.diff(np.log(wave), axis=-1))
        # a NaN, a repeated wavelength or a decreasing one is not a pixel step
        step = step[np.isfinite(step) & (step > 0)]
        if step.size:
            finest, typical = float(step.min()), float(np.median(step))
            if finest < 0.5 * typical:
                from .logger import log
                log("the finest pixel step in %s is %.4f km/s against a median"
                    " of %.4f: that looks like a glitch in the wavelength"
                    " solution rather than a pixel, and it is about to set the"
                    " grid for the whole run"
                    % (os.path.basename(path), finest, typical), "warn")
            return finest
    raise SystemExit("none of the first files in %s carries a usable wavelength"
                     " grid, so domain.smart_dv has nothing to measure"
                     % directory)


def smart_dv_from_step(step: float, value=True) -> float:
    """The grid step for a measured pixel step: a fraction of it, rounded down.

    Down to the nearest 10 m/s, and down rather than to the nearest so that the
    grid stays at or below the requested fraction of a pixel. Rounding at all
    is what keeps a wavelength solution that moved by a hair between two
    reductions from re-keying the cube and orphaning twenty minutes of work.
    """
    fraction = SMART_DV_FRACTION if value is True else float(value)
    if not 0 < fraction < 1:
        raise ValueError("domain.smart_dv must be true or a fraction in (0, 1),"
                         " got %r" % (value,))
    dv = math.floor(fraction * float(step) * 100.0) / 100.0
    if dv <= 0:
        raise ValueError("smart dv resolves to %g km/s from a pixel step of %g:"
                         " below 10 m/s there is nothing left to round to"
                         % (dv, step))
    return dv


def resolve_smart_dv(cfg: dict) -> dict:
    """`domain.dv` measured from the data instead of chosen.

    A magic grid is uniform in velocity and a spectrograph is not, so one dv
    has to serve the finest pixel in the domain; anything below that is paid
    for on every sample of every exposure and buys nothing. On SPIRou the step
    is around 2.27 km/s, which a fixed dv of 0.5 oversamples more than fourfold
    in every array the run allocates and every sweep the fit makes over them.

    `domain.dv` in the config is then DEAD: it is overwritten here, not
    combined with anything, and the run says so out loud rather than leaving a
    number in the file to be believed. It may also simply be null.

    Resolved here, at load, and written into `domain.dv` as a plain number, so
    that the cube cache key hashes the step that was actually used and the
    saved config says what it was. `smart_dv` may also be a number, the
    fraction to use in place of the default 0.7.
    """
    from .logger import log

    value = cfg["domain"].get("smart_dv")
    if not value:
        return cfg
    if cfg["domain"].get("grid_source", "magic") == "native":
        log("domain.smart_dv is ignored with grid_source 'native': the grid"
            " comes from the file and brings its own step", "warn")
        return cfg
    step = measure_pixel_dv(spectra_dir(cfg),
                            cfg["input"].get("pattern", "*t.fits"))
    dv = smart_dv_from_step(step, value)
    fraction = SMART_DV_FRACTION if value is True else float(value)
    was = cfg["domain"].get("dv")
    log("smart dv: finest pixel is %.4f km/s, so dv = %.2f km/s, %.0f%% of it"
        % (step, dv, 100 * fraction), "value")
    if isinstance(was, (int, float)):
        log("domain.dv: the %.2f km/s in the config is NOT used, smart_dv"
            " measured the step above instead" % float(was), "warn")
    cfg["domain"]["dv"] = float(dv)
    cfg["domain"]["pixel_dv"] = float(step)

    # Every window below is counted in SAMPLES, so a coarser grid widens all of
    # them in velocity without a line of the config changing. Said out loud
    # rather than rescaled: what these should cover is a modelling decision,
    # and one this function has no business making on its own.
    for key, section in (("window", "highpass"),
                         ("empirical_noise_box", "weights")):
        n = cfg[section].get(key)
        if n:
            log("  %s.%s is %d samples, so %.0f km/s at this step"
                % (section, key, n, n * dv), "warn")
    return cfg


def spectra_dir(config: dict) -> str:
    """The folder holding one object's spectra: input.directory/input.object.

    Resolved on every call rather than stored back into `input.directory`.
    Storing it would give that one key two meanings, the root in a config a
    person writes and the object's folder in a config a run saved, and loading
    a saved config would then append the object a second time.
    """
    inp = config.get("input") or {}
    root = inp.get("directory") or "data"
    obj = inp.get("object")
    if not obj:
        return root
    # a config written before the root/object split already ends in the object
    if os.path.basename(os.path.normpath(root)) == str(obj):
        return root
    return os.path.join(root, str(obj))


def _apply_run_overrides(cfg: dict, object_name, data_dir, out_dir) -> dict:
    """What the command line said, on top of whichever layer merged last."""
    if object_name:
        cfg["input"]["object"] = object_name
    if data_dir:
        cfg["input"]["directory"] = data_dir
    if out_dir:
        cfg["output"]["directory"] = out_dir
    return cfg


def detect_instrument(directory: str, pattern: str = "*t.fits") -> str:
    """Read INSTRUME from the first file that opens, and refuse to guess.

    The instrument decides which FITS extensions hold the flux, the wavelength
    solution and the blaze, and getting that wrong fails silently: reading
    FluxA from a SPIRou file raises nothing, it returns half the light with a
    different blaze and SNR. So it is read from the data and never offered as
    a command-line flag.
    """
    import glob

    from .io import robust_open

    paths = sorted(glob.glob(os.path.join(directory, pattern)))
    if not paths:
        raise SystemExit("no file matching %s in %s" % (pattern, directory))
    for path in paths[:20]:
        try:
            with robust_open(path) as hdulist:
                name = str(hdulist[0].header.get("INSTRUME", "")).strip().upper()
        except Exception:                                     # noqa: BLE001
            continue
        if name:
            return name
    raise SystemExit("none of the first files in %s declares INSTRUME"
                     % directory)


def load_config(path: str | None, object_name: str | None = None,
                data_dir: str | None = None, out_dir: str | None = None,
                instrument: str | None = None) -> dict:
    """The three layers of config.yaml, merged, with the instrument resolved.

    `general`, then `instruments.<INSTRUME>`, then `objects.<object_name>`,
    each winning over the last, all on top of the package DEFAULTS. The
    instrument is read from the files themselves unless one is named, and the
    extension table it carries is registered with `tfits` so every reader in
    the package agrees on which extension is which.

    A file written in the old flat form, with `input:` and `domain:` at the top
    level, still loads: it is treated as one more layer on top of `general`.
    """
    user = {}
    if path is not None:
        if not os.path.exists(path):
            raise FileNotFoundError("config file not found: %s" % path)
        with open(path, "r") as handle:
            user = read_yaml(handle) or {}

    layered = any(k in user for k in ("general", "instruments", "objects"))
    cfg = _deep_update(DEFAULTS, user.get("general", {}) if layered else user)
    if layered and user.get("general") is None:
        cfg = _deep_update(cfg, {k: v for k, v in user.items()
                                 if k not in ("general", "instruments", "objects")})

    cfg = _apply_run_overrides(cfg, object_name, data_dir, out_dir)

    table = (user.get("instruments") or {}) if layered else {}
    if table:
        from . import tfits as _tfits
        _tfits.register_instruments(
            {name: dict(block.get("extensions") or {})
             for name, block in table.items() if block.get("extensions")})
    # Everything below this point reads a spectrum, so say plainly here what a
    # missing folder means rather than letting a glob come back empty inside
    # the instrument detection or the dv measurement.
    if object_name and (table or cfg["domain"].get("smart_dv")):
        directory = spectra_dir(cfg)
        if not os.path.isdir(directory):
            raise SystemExit(
                "no directory %s. The object names a folder under the input"
                " root, which is input.directory (%s) unless --data-dir"
                " overrides it." % (directory, cfg["input"]["directory"]))

    # only look at the data when an object was named: loading the file to read
    # it, which the tests and the docs do, must not need a telescope
    if instrument is None and table and object_name:
        instrument = detect_instrument(spectra_dir(cfg),
                                       cfg["input"].get("pattern", "*t.fits"))
    if instrument:
        block = dict(table.get(instrument.upper()) or {})
        block.pop("extensions", None)
        if not block and table:
            raise SystemExit(
                "config.yaml has no `instruments: %s:` block. The instrument is"
                " read from INSTRUME in the data, so add one rather than"
                " renaming anything." % instrument)
        cfg = _deep_update(cfg, block)
        cfg["input"]["instrument"] = instrument.upper()

    if object_name and layered:
        cfg = _deep_update(cfg, (user.get("objects") or {}).get(object_name, {}))
    # last word to the command line: an object block may move the roots, a flag
    # passed to this run overrules it
    cfg = _apply_run_overrides(cfg, object_name, data_dir, out_dir)

    # dv from the data, once the instrument's domain and the object's folder
    # are both known. Only when an object was named, for the same reason the
    # instrument is only detected then: loading the file to read it must not
    # need the data to be on this disk. A run saves the resolved number, so
    # every stage after this one reads a plain dv.
    if object_name and cfg["domain"].get("smart_dv"):
        cfg = resolve_smart_dv(cfg)

    # A knob that used to exist. Removing it silently would mean a config
    # asking for float64 got float32 and no word about it.
    if cfg["output"].pop("cube_dtype", None) is not None:
        from .logger import log
        log("output.cube_dtype is no longer a setting and was dropped from this"
            " configuration: the cube is float32, in cube.CUBE_DTYPE", "warn")

    # max_sigma is the correct name for the cut; max_mad is what it was called
    # first. Both drive the same key so a config written either way behaves
    # identically, and a config carrying both is a contradiction worth refusing.
    tf = (user or {}).get("twoframe") or {}
    for new, old in (("max_sigma", "max_mad"),
                     ("max_sigma_rounds", "max_mad_rounds")):
        if new in tf:
            if old in tf and tf[old] != tf[new]:
                raise ValueError("config sets both twoframe.%s and twoframe.%s"
                                 " to different values" % (new, old))
            cfg["twoframe"][old] = tf[new]
    return cfg


def cache_key(config: dict) -> str:
    """Hash of the config entries that affect the registered data cube.

    Anything downstream of the cube (PCA settings, plots) is deliberately left
    out so that re-running with a different n_components reuses the cache.
    """
    relevant = {
        "input": config["input"],
        "quality": config["quality"],
        # domain.bands only chooses which columns the *fit* sees, so it must not
        # invalidate a cube that took twenty minutes to build
        "domain": {k: v for k, v in config["domain"].items() if k != "bands"},
        "registration": config["registration"],
        "highpass": config["highpass"],
        # weights that are baked into the cube rather than applied later
        "weights": {
            k: config["weights"][k]
            for k in ("mode", "snr_pixel_scale", "sigma_floor_frac",
                      "telluric_mode", "empirical_noise_box")
        },
    }
    # Entries left unset are stripped before hashing. Without this, adding a new
    # optional knob to `quality` changes the key of every existing config and
    # orphans cubes that took twenty minutes to build, purely because a default
    # of None appeared in the dictionary. A knob that is not set describes no
    # behaviour, so it must not describe a different cube either.
    # The epoch window is dropped from the hash when unset. Without this, adding
    # these two knobs would change the key of every existing config and orphan
    # cubes that took twenty minutes to build, purely because a default of None
    # appeared in the dictionary. Only these two are treated this way: the older
    # entries are hashed exactly as before, including their Nones, so keys
    # computed by earlier versions still match.
    relevant["quality"] = {k: v for k, v in relevant["quality"].items()
                           if not (k in ("min_rjd", "max_rjd", "max_sky_ratio",
                                         "isolated_window")
                                   and v is None)}
    relevant["weights"] = {k: v for k, v in relevant["weights"].items()
                           if not (k == "empirical_noise_box" and v is None)}
    # The Doppler shift that registers each spectrum became relativistic
    # (grids.doppler) on 2026-09-10; it was the first-order 1 + v/c, 0.003
    # pixel (1.5 m/s) off at a 30 km/s BERV. A cube registered the old way must
    # not be reused, so the formula is in the key. An observer-frame cube
    # applies no shift at all, and keeps its key and its telluric map.
    if relevant["registration"].get("frame") != "observer":
        relevant["doppler"] = "relativistic"
    blob = json.dumps(relevant, sort_keys=True, default=str).encode()
    return hashlib.sha1(blob).hexdigest()[:12]

