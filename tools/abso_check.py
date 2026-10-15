"""La correction tellurique, avec TAPAS dans le BON repere.

Trois courbes par fenetre:
  le spectre brut, normalise par son continu
  le spectre nettoye par lbl_telluclean
  le modele TAPAS evalue a lambda / D(BERV), c'est-a-dire dans le repere du
  sol, ou l'atmosphere se trouve

Et dessous, ce que le nettoyage a divise contre ce que le modele demande.
"""
import glob
import sys

import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, "/Users/eartigau/GitHubProjects/pca2d_preclean")
from pca2d import readers                                      # noqa: E402

C = 299792.458
OUT = "/Users/eartigau/tmp/diagnostic_harps/abso_check.pdf"
MODELS = "/Users/eartigau/tmp/pca2d_harps_v2/lbl_data/models"
WINDOWS = [(686.5, 689.5, "O2 bande B, saturee"),
           (626.5, 629.5, "O2 gamma"),
           (593.0, 595.0, "eau 594, propre a l'oeil"),
           (646.0, 652.0, "eau 649")]
NAME = "r.HARPS.2003-08-08T07:11:56.372_S2D_A.fits"


def main():
    inst = readers.instrument("HARPS", "ESO")
    inst.params.set("DATA_DIR", "/Users/eartigau/tmp/pca2d_harps_v2/lbl_data")
    spl_others, spl_water = readers._lbl_science(
        "tellu_clean").get_tapas_spl(inst, MODELS)

    raw = "data/HD10700_HARPS/" + NAME
    a = fits.getdata(raw, "SCIDATA")
    b = fits.getdata("data/HD10700_TC/" + NAME, "SCIDATA")
    w = fits.getdata(raw, "WAVEDATA_VAC_BARY") / 10.0
    head = fits.getheader(raw)
    airmass = float(head["HIERARCH ESO TEL AIRM START"])
    berv = float(head["HIERARCH ESO QC BERV"])
    doppler = np.sqrt((1 + berv / C) / (1 - berv / C))

    with PdfPages(OUT) as pdf:
        for lo, hi, label in WINDOWS:
            k = (w >= lo) & (w <= hi) & np.isfinite(a) & (a > 0)
            if k.sum() < 200:
                continue
            ww, aa, bb = w[k], a[k], b[k]
            o = np.argsort(ww)
            ww, aa, bb = ww[o], aa[o], bb[o]
            # LE REPERE DU SOL: l'atmosphere est la, pas au barycentre
            ground = ww / doppler
            good = ((np.clip(spl_others(ground), 0, 1) ** airmass)
                    * (np.clip(spl_water(ground), 0, 1) ** 0.4))
            # et le meme modele au MAUVAIS endroit, pour la comparaison
            wrong = ((np.clip(spl_others(ww), 0, 1) ** airmass)
                     * (np.clip(spl_water(ww), 0, 1) ** 0.4))

            fig, ax = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True)
            ax[0].plot(ww, aa / np.nanpercentile(aa, 95), "k-", lw=0.7,
                       label="brut / continu")
            ax[0].plot(ww, bb / np.nanpercentile(bb, 95), "C0-", lw=0.7,
                       label="nettoye par lbl_telluclean")
            ax[0].plot(ww, good, "C2-", lw=1.1, alpha=0.85,
                       label="TAPAS dans le repere du SOL (BERV %+.1f km/s)" % berv)
            ax[0].plot(ww, wrong, "C3--", lw=0.9, alpha=0.6,
                       label="TAPAS sur la grille barycentrique (ce que LBL fait)")
            ax[0].set_ylim(-0.05, 1.25)
            ax[0].set_ylabel("flux normalise")
            ax[0].legend(fontsize=7.5, loc="lower left")
            ax[0].set_title("%s : %.1f-%.1f nm, masse d'air %.2f"
                            % (label, lo, hi, airmass))

            ax[1].plot(ww, bb / aa, "C0-", lw=0.8,
                       label="division appliquee (nettoye / brut)")
            ax[1].plot(ww, 1.0 / np.clip(good, 1e-3, None), "C2-", lw=1.1,
                       alpha=0.85, label="1 / TAPAS au bon endroit")
            ax[1].set_yscale("log")
            ax[1].set_xlabel("longueur d'onde dans le vide (nm)")
            ax[1].set_ylabel("facteur")
            ax[1].legend(fontsize=7.5)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    print(OUT)


if __name__ == "__main__":
    main()
