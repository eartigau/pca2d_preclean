"""YAML configuration handling with defaults and a reproducibility hash."""

from __future__ import annotations

import copy
import hashlib
import json
import os

import yaml

# Every tunable lives here; config.yaml only needs to override what differs.
DEFAULTS = {
    "input": {
        "format": "tfits",           # 'tfits' = APERO order-by-order spectra
                                     #           (the nominal input)
                                     # 's1d'   = the resampled s1d_v products
        "directory": "data/tfiles",
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
        "dv": 0.5,                   # km/s per pixel (500 m/s)
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
        # storage dtype of the cached cube. Lives here and not under the fit:
        # it is a property of the cube, and float64 doubles the memory for
        # nothing the fit can use.
        "cube_dtype": "float32",
        "directory": "outputs",
        "tag": "run",
        "cache_directory": "cache",
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
        # False since 2026-09-09. The template was a star-frame median
        # subtracted before the components ran, and it is the mirror of the
        # observer-frame mean that used to be subtracted before them too. That
        # one was removed for a reason that applies here word for word: what is
        # taken out before the fit is outside the model, and
        # reconstruct.correction_on_grid removes components and never a
        # template, so it was fitted by nothing and removed from nothing. Now
        # neither frame gets a free zeroth order and each block describes its
        # own static content. The only thing still subtracted is the
        # even-minus-odd instrumental offset, which belongs to neither.
        "template": False,
        "shift": "lanczos",
        "kernel_halfwidth": 8,
        "gap_guard": True,
        "leakage": True,
        "chunk": None,
        "order": "star_first",
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
                data_dir: str = "data", instrument: str | None = None) -> dict:
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
            user = yaml.safe_load(handle) or {}

    layered = any(k in user for k in ("general", "instruments", "objects"))
    cfg = _deep_update(DEFAULTS, user.get("general", {}) if layered else user)
    if layered and user.get("general") is None:
        cfg = _deep_update(cfg, {k: v for k, v in user.items()
                                 if k not in ("general", "instruments", "objects")})

    if object_name:
        cfg["input"]["object"] = object_name
        cfg["input"]["directory"] = os.path.join(data_dir, object_name)

    table = (user.get("instruments") or {}) if layered else {}
    if table:
        from . import tfits as _tfits
        _tfits.register_instruments(
            {name: dict(block.get("extensions") or {})
             for name, block in table.items() if block.get("extensions")})
    # only look at the data when an object was named: loading the file to read
    # it, which the tests and the docs do, must not need a telescope
    if instrument is None and table and object_name:
        instrument = detect_instrument(cfg["input"]["directory"],
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
        cfg["input"]["object"] = object_name
        cfg["input"]["directory"] = os.path.join(data_dir, object_name)

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
    blob = json.dumps(relevant, sort_keys=True, default=str).encode()
    return hashlib.sha1(blob).hexdigest()[:12]

