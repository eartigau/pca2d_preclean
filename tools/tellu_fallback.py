"""CE QUE LE REPLI FAIT VRAIMENT, et pourquoi il ne faut pas s'en contenter.

Quand le critere de la derivee n'a pas de quoi trancher, `fit_exponents`
renvoie `fallback` (0.35) au lieu d'une valeur mesuree: la courbe de cout est
plate, ou son minimum tombe sur un bord du balayage. C'est le bon
comportement, puisqu'il y a toujours de l'eau et qu'un zero serait un echec
deguise en mesure.

MAIS CE N'EST PAS NEUTRE. Le spectre est bel et bien divise par une
transmission a cet exposant-la, donc une correction REELLE est appliquee a
partir d'un chiffre qui ne vient pas de ces spectres. Cette figure le montre.

    page 1  la courbe de cout de la pose, celle qui a declenche le repli,
            a cote de celle d'une pose qui converge: l'une a un creux franc,
            l'autre monte des la plus petite valeur
    page 2  le spectre avant et apres, sur deux bandes, et le RAPPORT des
            deux en dessous, qui est exactement ce que le repli a change
    page 3  le meme rapport pour plusieurs exposants, pour voir de combien
            la correction se trompe si la vraie colonne d'eau n'est pas 0.35

Mesure du 2026-10-16 sur GJ 581 (M3.0 V): ses 120 poses recoivent toutes le
repli, et au fond de la bande B le nettoyage remonte le flux d'un facteur
2.68. Voir docs/telluric_derivative_test.md.

    python tools/tellu_fallback.py [--raw DIR] [--tc DIR] [--out PDF]
"""
import argparse
import glob
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.io import fits  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from scipy.interpolate import InterpolatedUnivariateSpline as IUS  # noqa: E402

sys.path.insert(0, "/Users/eartigau/lbl-rapide")

C = 299792.458
#: les bandes montrees, choisies pour que l'une soit profonde et l'autre non
ZOOMS = ((688.0, 690.5, "bande B, saturee"),
         (594.0, 595.5, "eau 594, moderee"))
#: les exposants compares page 3, autour du repli
ALPHAS = (0.05, 0.35, 1.0)


def splines(tapas):
    table = fits.getdata(tapas)
    wave = np.array(table["WAVELENGTH"], dtype=float)
    water = IUS(wave, np.array(table["ABSO_WATER"], dtype=float), k=1, ext=3)
    dry = IUS(wave, np.array(table["ABSO_OTHERS"], dtype=float), k=1, ext=3)
    return water, dry


def ground_frame(path):
    """Le spectre, sa grille au repere du SOL, et ce qu'il faut des en-tetes."""
    with fits.open(path) as hdul:
        flux = hdul[1].data.astype(float)
        wave = hdul[4].data.astype(float) / 10.0
        header = hdul[0].header
    berv = float(header.get("HIERARCH ESO QC BERV", 0.0))
    return dict(
        flux=flux,
        wave=wave / np.sqrt((1.0 + berv / C) / (1.0 - berv / C)),
        airmass=float(header.get("HIERARCH ESO TEL AIRM START", np.nan)),
        snr=float(header.get("HIERARCH ESO QC ORDER50 SNR", np.nan)),
        berv=berv,
        water=header.get("TLPEH2O"),
        name=os.path.basename(path))


def cost_curve(frame, deriv, spl_water, spl_dry, scan):
    """Le critere contre l'exposant, sur les bandes d'eau."""
    keep = deriv.band_mask(frame["wave"], deriv.WATER_BANDS_HARPS)
    w = frame["wave"][keep].ravel()
    f = frame["flux"][keep].ravel()
    good = np.isfinite(w) & np.isfinite(f) & (f > 0)
    w, f = w[good], f[good]
    order = np.argsort(w)
    w, f = w[order], f[order]
    model = deriv.model_grid(spl_water, spl_dry, w.min() * 0.999,
                             w.max() * 1.001, 115000.0)
    return np.array([deriv.objective(f, deriv.transmission(
        model, w, a, frame["airmass"])) for a in scan])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    base = "/Users/eartigau/tmp/pca2d_harps_v2/lbl_data"
    parser.add_argument("--raw", default=base + "/science/GJ581")
    parser.add_argument("--tc", default=base + "/science/GJ581_tc")
    parser.add_argument("--ok", default=base + "/science/HD10700",
                        help="une campagne qui converge, pour comparer")
    parser.add_argument("--tapas", default=base + "/other/tapas_lbl.fits")
    parser.add_argument("--out",
                        default="/Users/eartigau/tmp/diagnostic_harps/"
                                "tellu_fallback.pdf")
    args = parser.parse_args(argv)

    from lbl.science import tellu_derivative as deriv

    spl_water, spl_dry = splines(args.tapas)
    cleaned = sorted(glob.glob(os.path.join(args.tc, "*.fits")))
    if not cleaned:
        raise SystemExit("aucun spectre nettoye dans %s" % args.tc)
    tc_path = cleaned[0]
    raw_path = os.path.join(args.raw, os.path.basename(tc_path))
    if not os.path.exists(raw_path):
        raise SystemExit("pas de spectre brut pour %s" % tc_path)
    frame = ground_frame(raw_path)
    after = fits.getdata(tc_path, 1).astype(float)
    good = sorted(glob.glob(os.path.join(args.ok, "*.fits")))
    other = ground_frame(good[0]) if good else None

    scan = np.array([0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.6,
                     1.0, 1.6, 2.5])
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with PdfPages(args.out) as pdf:
        _page_cost(pdf, frame, other, deriv, spl_water, spl_dry, scan)
        _page_effect(pdf, frame, after)
        _page_alphas(pdf, frame, deriv, spl_water, spl_dry)
    print("%s" % args.out)


def _page_cost(pdf, frame, other, deriv, spl_water, spl_dry, scan):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, f, title in ((axes[0], frame, "la pose qui recoit le REPLI"),
                         (axes[1], other, "une pose qui CONVERGE")):
        if f is None:
            continue
        curve = cost_curve(f, deriv, spl_water, spl_dry, scan)
        ax.plot(scan, curve / curve.min(), "o-", ms=4)
        ax.axvline(0.35, color="C3", ls="--", lw=1.2, label="le repli, 0.35")
        best = scan[int(np.argmin(curve))]
        ax.axvline(best, color="C2", ls=":", lw=1.2,
                   label="le minimum, %.3f" % best)
        span = (curve.max() - curve.min()) / curve.min()
        ax.set_xscale("log")
        ax.set_xlabel("exposant d'eau")
        ax.set_ylabel("cout / son minimum")
        ax.set_title("%s\n%s   SNR %.0f   variation %.1f%%"
                     % (title, f["name"][7:26], f["snr"], 100 * span),
                     fontsize=9)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
    fig.suptitle("LE CRITERE A-T-IL DE QUOI TRANCHER? A gauche la courbe monte"
                 " des la plus petite valeur:\nil n'y a pas de minimum, donc"
                 " rien a mesurer, et fit_exponents rend le repli.",
                 fontsize=10)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def _page_effect(pdf, frame, after):
    fig, axes = plt.subplots(2, len(ZOOMS), figsize=(11, 6), squeeze=False,
                             gridspec_kw=dict(height_ratios=[2, 1]))
    for col, (low, high, what) in enumerate(ZOOMS):
        band = ((frame["wave"] > low) & (frame["wave"] < high)
                & np.isfinite(frame["flux"]) & (frame["flux"] > 0)
                & np.isfinite(after))
        if band.sum() < 30:
            continue
        w = frame["wave"][band]
        before = frame["flux"][band]
        now = after[band]
        order = np.argsort(w)
        w, before, now = w[order], before[order], now[order]
        scale = np.nanpercentile(before, 95)
        top, bottom = axes[0, col], axes[1, col]
        top.plot(w, before / scale, "k-", lw=0.9, label="avant")
        top.plot(w, now / scale, "-", lw=0.9, color="C3", label="apres")
        top.set_title("%s" % what, fontsize=9)
        top.legend(fontsize=7)
        bottom.plot(w, now / np.clip(before, 1e-9, None), "-", lw=0.9,
                    color="C3")
        bottom.axhline(1.0, color="k", lw=0.8, ls=":")
        bottom.set_xlabel("longueur d'onde (nm), repere du sol")
        for ax in (top, bottom):
            ax.grid(alpha=0.3)
        top.set_ylabel("flux normalise")
        bottom.set_ylabel("apres / avant")
    fig.suptitle("CE QUE LE REPLI A FAIT AU SPECTRE. Le rapport du bas est la"
                 " correction appliquee,\nissue d'un exposant qui ne vient pas"
                 " de ce spectre.", fontsize=10)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def _page_alphas(pdf, frame, deriv, spl_water, spl_dry):
    """De combien la correction se trompe si la vraie colonne n'est pas 0.35."""
    fig, axes = plt.subplots(1, len(ZOOMS), figsize=(11, 4.2), squeeze=False)
    for col, (low, high, what) in enumerate(ZOOMS):
        band = ((frame["wave"] > low) & (frame["wave"] < high)
                & np.isfinite(frame["flux"]) & (frame["flux"] > 0))
        if band.sum() < 30:
            continue
        w = np.sort(frame["wave"][band])
        model = deriv.model_grid(spl_water, spl_dry, w.min() * 0.999,
                                 w.max() * 1.001, 115000.0)
        ax = axes[0, col]
        for alpha in ALPHAS:
            trans = deriv.transmission(model, w, alpha, frame["airmass"])
            ax.plot(w, 1.0 / np.clip(trans, 1e-3, None), lw=1.1,
                    label="alpha = %.2f" % alpha,
                    ls="--" if abs(alpha - 0.35) > 1e-9 else "-")
        ax.set_title(what, fontsize=9)
        ax.set_xlabel("longueur d'onde (nm)")
        ax.set_ylabel("facteur applique au flux")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
    fig.suptitle("DE COMBIEN LE REPLI SE TROMPE. Le trait plein est le repli;"
                 " les pointilles sont ce que\nla correction aurait fait si la"
                 " vraie colonne d'eau avait ete trois fois plus faible ou"
                 " plus forte.", fontsize=10)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


if __name__ == "__main__":
    main()
