#!/usr/bin/env python
"""The figures of the site, from one run, as SVG next to the page.

    python docs/make_figures.py --run outputs/TOI2120/1-3v
    python docs/make_figures.py --run outputs/TOI2120/1-3v --only lbl \
        --lbl-objects "TOI2120_PCA2D_1-3v=1-3v, LBL's own template" \
                      "TOI2120_PCA2D_1-3v_PT=1-3v, the fit's star template"

--format pdf or png draws the same figures in another format, to look at them
somewhere that does not read SVG.

Every figure is drawn by the code the run's own report uses where there is
one (pca2d.figures.sequence.draw_window, twoframe.plot_coeffs,
plotting.plot_correlations), so the site cannot show something the pipeline
does not. The corrected-file figure puts panel 3 of the report beside the
corrected files themselves, through the same high pass, and beside what the
files held before 2026-09-10, when the parity mean was left in them. The
velocities come from the run's LBL stage: the delivered spectra and the
corrected ones, on the same exposures.
"""

from __future__ import annotations

import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.signal import savgol_filter

from pca2d.config import load_config, spectra_dir
from pca2d.figures.sequence import AFTER, BEFORE, draw_window, load_context, window_arrays
from pca2d.grids import doppler, pixel_shift
from pca2d.lblscan import nightly, read_rdb, velocity_stats  # noqa: F401
from pca2d.logger import log
from pca2d.plotting import nan_cmap, plot_correlations, rank_correlations
from pca2d.reconstruct import correction_on_grid, load_model, order_correction
from pca2d.twoframe import exposures_label, plot_coeffs

HERE = os.path.dirname(os.path.abspath(__file__))
C_KMS = 299792.458


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", default="outputs/TOI2120/1-3v")
    p.add_argument("--cube", default=None,
                   help="the cube the run was fitted on; the only cache/cube_* by default")
    p.add_argument("--lbl-dir", default="lbl")
    p.add_argument("--out", default=os.path.join(HERE, "figures"))
    p.add_argument("--only", nargs="+", default=None,
                   choices=("sequence", "corrected", "coefficients", "lbl", "strpca"),
                   help="draw only these figures")
    p.add_argument("--format", default="svg", choices=("svg", "pdf", "png"),
                   help="svg for the site; pdf or png to look at them elsewhere")
    p.add_argument("--lbl-objects", nargs="+", default=None, metavar="NAME[=LABEL]",
                   help="the LBL objects set against the delivered spectra, one"
                        " panel each; the run's own corrected object by default")
    return p.parse_args(argv)


#: what save() writes; main sets it from --format
FORMAT = "svg"


def save(fig, out, name):
    path = os.path.join(out, name + "." + FORMAT)
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    log("wrote %s (%.0f kB)" % (path, os.path.getsize(path) / 1e3), "value")


# ----------------------------------------------------------- corrected file
def highpass_at(wave, flux, x, berv, window=151, polyorder=2, dv=0.5):
    """The pipeline's high pass of one order (tfits + preprocess.highpass):
    onto a log grid of dv km/s in the observer frame, gaps bridged linearly in
    flux, ln f minus its Savitzky-Golay; read at star-frame wavelengths x."""
    good = np.isfinite(flux) & (flux > 0)
    u = np.arange(np.log(wave[good][0]), np.log(wave[good][-1]), dv / C_KMS)
    lnf = np.log(np.interp(u, np.log(wave[good]), flux[good]))
    lnf -= savgol_filter(lnf, window, polyorder)
    lnf[np.interp(u, np.log(wave), (~good).astype(float)) > 0] = np.nan
    return np.interp(np.log(x) - np.log(doppler(berv)), u, lnf,
                     left=np.nan, right=np.nan)


def corrected_file_figure(ctx, run, source, centre=1267.0, width=2.0):
    """Panel 3, the corrected files, and the files before the parity mean fix."""
    lo, hi = centre - width / 2, centre + width / 2
    model = load_model(os.path.join(run, "twoframe_components.fits"))
    coeffs = {os.path.basename(str(r["filename"])): r for r in model["coeffs"]}
    corrected = sorted(glob.glob(os.path.join(run, "corrected", "*.fits")))
    first = os.path.join(source, os.path.basename(corrected[0]).split("t_")[0] + "t.fits")
    w_all = fits.getdata(first, "WaveAB")
    order = int(np.argmin(np.abs(w_all[:, w_all.shape[1] // 2] - centre)))
    arr = window_arrays(ctx["cube"], ctx["fit"], ctx["means"], ctx["group"], ctx["grid"],
                        ctx["dv"], ctx["delta"], centre, width, templates=ctx["templates"])
    g = arr["grid"]
    win = (g >= lo) & (g <= hi)
    gw = g[win]
    rows = {n: r for r, n in enumerate(ctx["names"]) if ctx["parity"][r] == order % 2}
    k = model["n_earth"]
    # the high pass the depicted run used, not today's default: a run saved
    # before 2026-09-11 filtered over 151 samples, a later one over 100 km/s
    from pca2d.config import load_config
    hp = load_config(os.path.join(run, "resolved_config.yaml"))["highpass"]
    cut = dict(window=int(hp["window"]), polyorder=int(hp["polyorder"]), dv=ctx["dv"])
    views = {"panel3": [], "files": [], "before": []}
    bervs = []
    for path in corrected:
        raw = os.path.basename(path).split("t_")[0] + "t.fits"
        if raw not in rows:
            continue
        with fits.open(os.path.join(source, raw)) as h:
            f_raw = h["FluxAB"].data[order].astype(float)
            wave = h["WaveAB"].data[order].astype(float)
            berv = float(h["FluxAB"].header["BERV"])
        with fits.open(path) as h:
            f_cor = h["FluxAB"].data[order].astype(float)
        # the files as they were written before 2026-09-10: the observer
        # components divided out and the parity mean left in (n_earth=0 in
        # order_correction is the call that leaves it out of the sum)
        corr, _, _ = correction_on_grid(model, coeffs[raw], 0, k)
        values, live = order_correction(model, corr, 0, order, wave)
        f_old = np.where(live, f_raw * np.exp(-values), f_raw)
        views["panel3"].append(arr["home"]["corrected"][rows[raw]][win])
        views["files"].append(highpass_at(wave, f_cor, gw, berv))
        views["before"].append(highpass_at(wave, f_old, gw, berv))
        bervs.append(berv)
    o = np.argsort(bervs)
    views = {name: np.array(v)[o] for name, v in views.items()}
    views = {name: v - np.nanmedian(v, axis=1, keepdims=True) for name, v in views.items()}
    rms = lambda a: float(np.sqrt(np.nanmean(a ** 2)))
    stats = {"rows": len(o), "order": order,
             "files_vs_panel3": rms(views["files"] - views["panel3"]),
             "before_vs_panel3": rms(views["before"] - views["panel3"])}
    scale = np.nanpercentile(np.abs(views["panel3"]), 98)
    fig, axes = plt.subplots(3, 1, figsize=(9.4, 9.6), sharex=True)
    ext = [gw[0], gw[-1], len(o) - 0.5, -0.5]
    titles = [("panel3", "panel 3 of the report: the data minus the observer block and its parity mean"),
              ("files", "the corrected files, through the same high pass   (minus panel 3: rms %.4f)"
               % stats["files_vs_panel3"]),
              ("before", "the files before 2026-09-10, parity mean left in   (minus panel 3: rms %.4f)"
               % stats["before_vs_panel3"])]
    for ax, (name, title) in zip(axes, titles):
        im = ax.imshow(views[name], aspect="auto", cmap=nan_cmap("RdBu_r"),
                       vmin=-scale, vmax=scale, extent=ext)
        ax.set_title(title, fontsize=8.5)
        ax.set_ylabel("ordered by BERV", fontsize=8)
        ax.tick_params(labelsize=7)
    axes[-1].set_xlabel("wavelength (nm), star frame, order %d" % order, fontsize=9)
    fig.colorbar(im, ax=list(axes), fraction=0.022, pad=0.015,
                 label=r"$\ln(f/\mathrm{savgol}\,f)$")
    fig.suptitle("%.1f-%.1f nm: a corrected file holds panel 3; yellow is what the fit gave"
                 " no weight, NaN in the file" % (lo, hi), fontsize=10)
    return fig, stats


# ----------------------------------------------------------------- velocities
# read_rdb, nightly and velocity_stats live in pca2d.lblscan, the scan report,
# so the site and the report cannot count the same velocities two ways


def velocity_figure(series):
    """One panel per LBL object, on the exposures every one of them has.

    `series` is [(label, rdb)], the delivered spectra first. Returns the
    figure and [(label, stats)].
    """
    data = [read_rdb(path) for _, path in series]
    common = np.round(data[0][0], 6)
    for t, _, _ in data[1:]:
        common = np.intersect1d(common, np.round(t, 6))
    rows = []
    for (label, _), (t, v, e) in zip(series, data):
        keep = np.isin(np.round(t, 6), common)
        rows.append((label, t[keep], v[keep], e[keep],
                     velocity_stats(t[keep], v[keep], e[keep])))
    fig, axes = plt.subplots(len(rows), 1, figsize=(9.4, 1.0 + 2.3 * len(rows)),
                             sharex=True, sharey=True, squeeze=False)
    axes = axes[:, 0]
    for k, (ax, (label, t, v, e, s)) in enumerate(zip(axes, rows)):
        colour = BEFORE if k == 0 else AFTER
        ax.errorbar(t, v - np.median(v), yerr=e, fmt="o", ms=3, lw=0.6,
                    color=colour, ecolor=colour, alpha=0.75, capsize=0)
        ax.axhline(0, color="0.6", lw=0.6)
        ax.set_title("%s: rms %.2f m/s, robust sigma %.2f m/s, median error"
                     " %.2f m/s, nightly rms %.2f m/s"
                     % (label, s["rms"], s["robust"], s["median_error"],
                        s["nightly_rms"]), fontsize=8.5)
        ax.set_ylabel("velocity - median (m/s)", fontsize=8.5)
        ax.grid(alpha=0.15)
        ax.tick_params(labelsize=8)
    lim = np.percentile(np.abs(np.concatenate([v - np.median(v)
                                               for _, _, v, _, _ in rows])), 99.5)
    axes[0].set_ylim(-1.15 * lim, 1.15 * lim)
    axes[-1].set_xlabel("rjd (BJD - 2400000)", fontsize=9)
    fig.suptitle("LBL velocities of the same %d exposures" % len(common), fontsize=10)
    return fig, [(label, s) for label, _, _, _, s in rows]


def strpca_figure(rdb_path, fit):
    """LBL's STRPCA amplitudes against the fit's own coefficients.

    LBL projects each exposure's residual about the template on the tables,
    line by line; the fit measured the same thing as a_k, on the whole grid at
    once, and the template is the star at the mean a_k. The two agreeing,
    exposure by exposure, is the check that the tables say what they are meant
    to. Returns the figure and [(key, exposures, pearson r, slope)].
    """
    t = Table.read(rdb_path, format="ascii.rdb")
    keys = sorted((c for c in t.colnames if c.startswith("STRPCA")),
                  key=lambda c: int(c[len("STRPCA"):]))
    rejected = np.asarray(fit["rejected"], dtype=bool)
    first = {}
    for r, name in enumerate(str(x) for x in fit["filename"]):
        first.setdefault(name.split("t.fits")[0], r)
    stems = [os.path.basename(str(f)).split("t_")[0] for f in t["FILENAME"]]
    hit = np.array([s in first and not rejected[first[s]] for s in stems])
    rows = np.array([first[s] for s, h in zip(stems, hit) if h])
    rjd = np.asarray(t["rjd"], float)[hit]
    a = np.asarray(fit["a"], float)
    stats = []
    fig, axes = plt.subplots(len(keys), 2, figsize=(9.4, 0.9 + 3.0 * len(keys)),
                             squeeze=False, gridspec_kw={"width_ratios": [2.3, 1]})
    for (ax_t, ax_s), key in zip(axes, keys):
        k = int(key[len("STRPCA"):]) - 1
        x = a[rows, k] - a[~rejected, k].mean()
        y = np.asarray(t[key], float)[hit]
        e = np.asarray(t["s" + key], float)[hit]
        ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(e)
        r = float(np.corrcoef(x[ok], y[ok])[0, 1])
        slope = float(np.sum(x[ok] * y[ok]) / np.sum(x[ok] ** 2))
        stats.append((key, int(ok.sum()), r, slope))
        ax_t.errorbar(rjd[ok], y[ok], yerr=e[ok], fmt="o", ms=3, lw=0.6, color=AFTER,
                      ecolor=AFTER, alpha=0.7, capsize=0,
                      label="LBL: every line projected on %s" % key)
        ax_t.plot(rjd[ok], x[ok], "o", ms=2.4, color=BEFORE,
                  label="the fit: a%d less its mean" % (k + 1))
        ax_t.axhline(0, color="0.6", lw=0.6)
        ax_t.set_ylabel("amplitude along star component %d" % (k + 1), fontsize=8.5)
        ax_t.legend(fontsize=7.5, frameon=False, loc="upper left")
        ax_t.tick_params(labelsize=8)
        ax_t.grid(alpha=0.15)
        ax_s.errorbar(x[ok], y[ok], yerr=e[ok], fmt="o", ms=2.4, lw=0.5, color=AFTER,
                      ecolor=AFTER, alpha=0.55, capsize=0)
        lim = 1.1 * np.percentile(np.abs(np.r_[x[ok], y[ok]]), 99.5)
        ax_s.plot([-lim, lim], [-lim, lim], color="0.45", lw=0.8)
        ax_s.set(xlim=(-lim, lim), ylim=(-lim, lim))
        ax_s.set_aspect("equal")
        ax_s.set_xlabel("the fit's a%d" % (k + 1), fontsize=8.5)
        ax_s.set_ylabel("LBL's %s" % key, fontsize=8.5)
        ax_s.set_title("%d exposures: r = %.3f, slope %.3f" % (ok.sum(), r, slope),
                       fontsize=8.5)
        ax_s.tick_params(labelsize=8)
    axes[-1, 0].set_xlabel("rjd (BJD - 2400000)", fontsize=9)
    fig.suptitle("the star's variability, measured twice: by the fit on the whole"
                 " grid, by LBL line by line", fontsize=10)
    return fig, stats


def main(argv=None):
    global FORMAT
    args = parse_args(argv)
    FORMAT = args.format
    os.makedirs(args.out, exist_ok=True)
    cfg = load_config(os.path.join(args.run, "resolved_config.yaml"))
    source = spectra_dir(cfg)
    cube = args.cube or sorted(glob.glob(os.path.join(cfg["output"]["cache_directory"], "cube_*")))[0]
    ctx = load_context(cube, os.path.join(args.run, "fit.npz"), source)

    want = lambda name: args.only is None or name in args.only
    for centre, width, name in ((1669.5, 5.0, "sequence_1669nm"), (1267.0, 2.0, "sequence_1267nm")):
        fig = draw_window(ctx, centre, width) if want("sequence") else None
        if fig is not None:
            save(fig, args.out, name)

    if want("corrected"):
        fig, stats = corrected_file_figure(ctx, args.run, source)
        save(fig, args.out, "corrected_file_1267nm")
        log("corrected files against panel 3: rms %.4f now, %.4f with the parity"
            " mean left in" % (stats["files_vs_panel3"], stats["before_vs_panel3"]), "value")

    fit = ctx["fit"]
    keep = ~fit["rejected"]
    n_star = fit["a"].shape[1]
    sigma = fit["sigma_scaled"]
    if want("coefficients"):
        _coefficient_figures(fit, keep, n_star, sigma, args.out)

    obj, tag = os.path.normpath(args.run).split(os.sep)[-2:]
    rdb = lambda name: os.path.join(args.lbl_dir, "lblrdb", "lbl_%s_%s.rdb" % (name, name))
    series = [("delivered spectra", rdb(obj))]
    for item in args.lbl_objects or ["%s_PCA2D_%s=corrected, %s" % (obj, tag, tag)]:
        name, _, label = item.partition("=")
        series.append((label or name, rdb(name)))
    missing = [path for _, path in series if not os.path.exists(path)]
    if not want("lbl"):
        pass
    elif missing:
        log("no LBL velocities yet in %s: the velocity figure waits"
            % ", ".join(missing), "warn")
    else:
        fig, stats = velocity_figure(series)
        save(fig, args.out, "lbl_velocities")
        for label, s in stats:
            log("%-44s %d exposures: rms %.2f m/s, robust sigma %.2f, median error"
                " %.2f, nightly rms %.2f over %d nights"
                % (label, s["n"], s["rms"], s["robust"], s["median_error"],
                   s["nightly_rms"], s["nights"]), "value")

    own = rdb("%s_PCA2D_%s" % (obj, tag))
    if want("strpca") and os.path.exists(own) and any(
            c.startswith("STRPCA") for c in Table.read(own, format="ascii.rdb").colnames):
        fig, stats = strpca_figure(own, fit)
        save(fig, args.out, "strpca_%s" % tag)
        for key, n, r, slope in stats:
            log("%s against the fit's own coefficient, %d exposures: r = %.3f,"
                " slope %.3f" % (key, n, r, slope), "value")
    elif want("strpca") and args.only:
        log("no STRPCA columns in %s: the fit has one star component, or LBL has"
            " not measured it yet" % own, "warn")
    return None


def _coefficient_figures(fit, keep, n_star, sigma, out):
    """The run's own coefficient and correlation figures."""
    plot_coeffs(fit["bjd"][keep], fit["a"][keep], fit["b"][keep], fit["power_star"],
                fit["power_earth"], float(fit["chi2_null"]),
                os.path.join(out, "coefficients_vs_time." + FORMAT),
                err_a=sigma[keep][:, :n_star], err_b=sigma[keep][:, n_star:])
    plt.close("all")
    comps = [("a%d" % (i + 1), fit["a"][keep][:, i]) for i in range(n_star)]
    comps += [("b%d" % (j + 1), fit["b"][keep][:, j]) for j in range(fit["b"].shape[1])]
    labels = [str(x) for x in fit["anc_labels"]]
    rho = rank_correlations(comps, fit["anc_values"][:, keep])
    plot_correlations(rho, comps, labels, os.path.join(out, "correlations." + FORMAT),
                      title=exposures_label(fit["filename"], keep))
    plt.close("all")
    log("wrote coefficients_vs_time.%s and correlations.%s" % (FORMAT, FORMAT), "value")


if __name__ == "__main__":
    main()
