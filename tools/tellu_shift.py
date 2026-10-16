"""LE MODELE TELLURIQUE DOIT ETRE A VITESSE ZERO, et cette figure le verifie.

L'atmosphere ne bouge pas. Evalue dans le repere du sol, le modele TAPAS doit
donc coller aux raies telluriques du spectre SANS qu'on ait a le decaler, et
la vitesse n'est pas, et ne doit pas etre, un parametre libre de l'ajustement
(lbl/science/tellu_derivative.py n'en a aucun).

Si un decalage apparait quand meme, c'est que la grille de longueur d'onde du
spectre n'est pas dans le repere qu'on croit, et c'est un defaut des DONNEES ou
de leur lecture, pas une chose a absorber dans l'ajustement. Cette figure est
la pour le debugger.

CE QU'ELLE TRACE, par pose, le critere de la derivee etant evalue a exposant
fixe pour une serie de decalages du modele:

  page 1  le decalage prefere de chaque pose contre la date, le BERV, la masse
          d'air et le SNR, avec les modes du spectrographe en couleurs
          differentes. Une population a zero et une population ailleurs, c'est
          deux reperes, donc deux instruments
  page 2  la courbe de cout contre le decalage, une courbe par pose, pour voir
          si le minimum est franc ou s'il n'y a rien a trouver
  page 3  pour les poses qui s'ecartent le plus, le spectre et le modele
          superposes dans une bande d'eau, au decalage nul et au decalage
          prefere, c'est-a-dire ce que l'oeil doit trancher

Mesure du 2026-10-15 sur les 100 poses HARPS de HD 10700: 92 poses preferent
0.0 km/s, et les 8 poses du mode EGGS preferent +16.0 km/s. Voir
quality.drs_modes dans config.yaml.

   python tools/tellu_shift.py [dossier] [sortie.pdf]
"""
import glob
import os
import sys

import numpy as np
from astropy.io import fits
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)

C = 299792.458
#: le decalage balaye, en km/s. Large exprès: un vrai defaut de repere vaut
#: plusieurs dizaines de km/s, pas quelques centaines de m/s
SHIFTS = np.arange(-30.0, 30.1, 1.0)
#: l'exposant auquel le critere est evalue. Il n'a pas a etre juste: on cherche
#: OU le modele colle, pas combien il y a d'eau, et un exposant plausible
#: suffit pour que les raies soient a la bonne place
ALPHA = 0.3
#: la profondeur d'eau minimale d'un pixel retenu. Sans elle, les raies
#: stellaires, qui ne bougent pas avec le modele, diluent le signal
MIN_DEPTH = 0.02


def _spl(tapas_path):
    """Les deux splines de TAPAS, eau et sec, en nanometres."""
    table = fits.getdata(tapas_path)
    wave = np.array(table["WAVELENGTH"], dtype=float)
    water = IUS(wave, np.array(table["ABSO_WATER"], dtype=float), k=1, ext=3)
    dry = IUS(wave, np.array(table["ABSO_OTHERS"], dtype=float), k=1, ext=3)
    return water, dry


def _frame(path, spl_water, bands):
    """Un spectre, sa grille ramenee au repere du SOL, et ses en-tetes.

    Le S2D de l'ESO porte WAVEDATA_VAC_BARY, en angstroms: une grille deja
    ramenee au barycentre. L'atmosphere etant au sol, on defait cette
    correction, relativistiquement.
    """
    with fits.open(path) as hdul:
        flux = hdul[1].data.astype(float)
        wave = hdul[4].data.astype(float) / 10.0        # angstrom -> nm
        hdr = hdul[0].header
    berv = float(hdr.get("HIERARCH ESO QC BERV", 0.0))
    ground = wave / np.sqrt((1.0 + berv / C) / (1.0 - berv / C))
    keep = np.zeros(ground.shape, dtype=bool)
    for low, high in bands:
        keep |= (ground > low) & (ground < high)
    keep &= (1.0 - spl_water(ground)) > MIN_DEPTH
    w = ground[keep].ravel()
    f = flux[keep].ravel()
    good = np.isfinite(w) & np.isfinite(f) & (f > 0)
    w, f = w[good], f[good]
    order = np.argsort(w)
    return dict(
        wave=w[order], flux=f[order], berv=berv,
        airmass=float(hdr.get("HIERARCH ESO TEL AIRM START", np.nan)),
        snr=float(hdr.get("HIERARCH ESO QC ORDER50 SNR", np.nan)),
        mjd=float(hdr.get("MJD-OBS", np.nan)),
        mode=str(hdr.get("HIERARCH ESO INS MODE", "?")).strip(),
        name=os.path.basename(path))


def cost_curve(frame, deriv, model):
    """Le critere de la derivee contre le decalage, pour une pose."""
    out = np.full(SHIFTS.size, np.nan)
    for i, v in enumerate(SHIFTS):
        ratio = np.sqrt((1.0 + v / C) / (1.0 - v / C))
        trans = deriv.transmission(model, frame["wave"] * ratio, ALPHA,
                                   frame["airmass"])
        out[i] = deriv.objective(frame["flux"], trans)
    return out


def main(directory=None, out=None):
    directory = directory or os.path.join(
        HERE, "..", "..", "tmp", "pca2d_harps_v2", "lbl_data", "science",
        "HD10700")
    directory = os.path.abspath(directory)
    out = out or os.path.join(HERE, "outputs", "tellu_shift.pdf")
    os.makedirs(os.path.dirname(out), exist_ok=True)

    sys.path.insert(0, "/Users/eartigau/lbl-rapide")
    from lbl.science import tellu_derivative as deriv

    tapas = os.path.join(os.path.dirname(os.path.dirname(directory)),
                         "other", "tapas_lbl.fits")
    spl_water, spl_dry = _spl(tapas)
    bands = deriv.WATER_BANDS_HARPS

    files = sorted(glob.glob(os.path.join(directory, "*.fits")))
    if not files:
        raise SystemExit("aucun spectre dans %s" % directory)
    print("%d poses dans %s" % (len(files), directory), flush=True)

    frames, curves = [], []
    for i, path in enumerate(files, 1):
        frame = _frame(path, spl_water, bands)
        if frame["wave"].size < 500:
            print("  %3d/%d %s: trop peu de pixels d'eau"
                  % (i, len(files), frame["name"][:40]), flush=True)
            continue
        model = deriv.model_grid(spl_water, spl_dry,
                                 frame["wave"].min() * 0.996,
                                 frame["wave"].max() * 1.004, 115000.0)
        curve = cost_curve(frame, deriv, model)
        frames.append(frame)
        curves.append(curve)
        best = SHIFTS[int(np.nanargmin(curve))]
        print("  %3d/%d %s  mode %-5s BERV %+7.2f -> decalage %+6.1f km/s"
              % (i, len(files), frame["name"][7:26], frame["mode"],
                 frame["berv"], best), flush=True)

    curves = np.array(curves)
    best = SHIFTS[np.nanargmin(curves, axis=1)]
    contrast = np.nanmax(curves, axis=1) / np.nanmin(curves, axis=1)
    modes = np.array([f["mode"] for f in frames])
    _draw(out, frames, curves, best, contrast, modes, deriv, spl_water, spl_dry)
    print("\n%s" % out)
    for mode in sorted(set(modes)):
        sel = modes == mode
        print("  mode %-6s n=%3d  decalage median %+6.1f km/s  (de %+.1f a %+.1f)"
              % (mode, int(sel.sum()), float(np.median(best[sel])),
                 best[sel].min(), best[sel].max()))
    zero = np.abs(best) <= 1.0
    print("  a vitesse nulle (a 1 km/s pres): %d poses sur %d"
          % (int(zero.sum()), len(best)))


def _draw(out, frames, curves, best, contrast, modes, deriv,
          spl_water, spl_dry):
    colours = {m: c for m, c in zip(sorted(set(modes)),
                                    ("tab:blue", "tab:red", "tab:green",
                                     "tab:orange"))}
    with PdfPages(out) as pdf:
        # page 1: le decalage prefere, contre tout ce qui pourrait l'expliquer
        fig, axes = plt.subplots(2, 2, figsize=(11, 8))
        xs = [("date (MJD)", np.array([f["mjd"] for f in frames])),
              ("BERV (km/s)", np.array([f["berv"] for f in frames])),
              ("masse d'air", np.array([f["airmass"] for f in frames])),
              ("SNR ordre 50", np.array([f["snr"] for f in frames]))]
        for ax, (label, x) in zip(axes.ravel(), xs):
            for mode in sorted(set(modes)):
                sel = modes == mode
                ax.plot(x[sel], best[sel], "o", ms=5, alpha=0.8,
                        color=colours[mode], label="mode %s" % mode)
            ax.axhline(0.0, color="k", lw=1.0, ls="--")
            ax.set_xlabel(label)
            ax.set_ylabel("decalage prefere (km/s)")
            ax.grid(alpha=0.3)
        axes[0, 0].legend(fontsize=8)
        fig.suptitle("LE MODELE TELLURIQUE DOIT ETRE A ZERO. Une population"
                     " ailleurs, c'est un autre repere.\n"
                     "%d poses, critere de la derivee a exposant %.2f fixe"
                     % (len(frames), ALPHA))
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # page 2: les courbes de cout, pour voir si le minimum est franc
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        for curve, mode in zip(curves, modes):
            norm = (curve - np.nanmin(curve)) / (np.nanmax(curve)
                                                 - np.nanmin(curve))
            axes[0].plot(SHIFTS, norm, lw=0.8, alpha=0.5,
                         color=colours[mode])
        axes[0].axvline(0.0, color="k", lw=1.0, ls="--")
        axes[0].set_xlabel("decalage du modele (km/s)")
        axes[0].set_ylabel("cout, normalise par pose")
        axes[0].set_title("une courbe par pose")
        axes[0].grid(alpha=0.3)
        for mode in sorted(set(modes)):
            sel = modes == mode
            axes[1].plot(best[sel], contrast[sel], "o", ms=5, alpha=0.8,
                         color=colours[mode], label="mode %s" % mode)
        axes[1].axvline(0.0, color="k", lw=1.0, ls="--")
        axes[1].set_xlabel("decalage prefere (km/s)")
        axes[1].set_ylabel("contraste du cout (max/min)")
        axes[1].set_title("un minimum franc, ou rien a trouver")
        axes[1].legend(fontsize=8)
        axes[1].grid(alpha=0.3)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # page 3: l'oeil tranche. Les poses les plus decalees, spectre et
        # modele a zero et au decalage prefere
        worst = np.argsort(-np.abs(best))[:4]
        fig, axes = plt.subplots(len(worst), 1, figsize=(11, 2.6 * len(worst)),
                                 squeeze=False)
        for ax, idx in zip(axes[:, 0], worst):
            frame = frames[idx]
            model = deriv.model_grid(spl_water, spl_dry,
                                     frame["wave"].min() * 0.996,
                                     frame["wave"].max() * 1.004, 115000.0)
            band = (frame["wave"] > 687.0) & (frame["wave"] < 690.0)
            if band.sum() < 50:
                band = np.ones(frame["wave"].shape, dtype=bool)
            w = frame["wave"][band]
            f = frame["flux"][band]
            f = f / np.nanpercentile(f, 95)
            ax.plot(w, f, "k-", lw=0.8, label="spectre (repere du sol)")
            for v, style, colour in ((0.0, "-", "tab:blue"),
                                     (best[idx], "--", "tab:red")):
                ratio = np.sqrt((1.0 + v / C) / (1.0 - v / C))
                trans = deriv.transmission(model, w * ratio, ALPHA,
                                           frame["airmass"])
                ax.plot(w, trans, style, lw=1.2, color=colour, alpha=0.8,
                        label="modele a %+.1f km/s" % v)
            ax.set_title("%s  mode %s  BERV %+.2f km/s"
                         % (frame["name"][7:26], frame["mode"], frame["berv"]),
                         fontsize=9)
            ax.set_ylabel("flux normalise")
            ax.legend(fontsize=7, ncol=3)
            ax.grid(alpha=0.3)
        axes[-1, 0].set_xlabel("longueur d'onde (nm), repere du sol")
        fig.suptitle("LES QUATRE POSES LES PLUS DECALEES: ou le modele"
                     " tombe-t-il vraiment?")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


if __name__ == "__main__":
    main(*sys.argv[1:3])
