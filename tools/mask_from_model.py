"""Un masque dont les POSITIONS viennent du modele et les VALEURS des donnees.

L'UTILISATEUR, le 2026-10-15: « on peut trouver le masque avec le template [du
modele], mais ne *pas* l'utiliser comme template, les raies n'ont pas
exactement les bonnes formes » et « le template de la vitesse elle-meme doit
etre celui fait avec les positions ».

CE QUE CHACUN APPORTE.
    le modele PHOENIX    OU sont les raies. Il est sans bruit, donc ses
                         extrema sont tous reels: 30835 raies contre 51113
                         extrema dans le gabarit mesure, soit 20000 de trop
                         dans celui-ci, qui sont du bruit pris pour des raies.
                         C'est la que naissait le line_snr median de 1.46.
    le gabarit mesure    A QUOI elles ressemblent. Les profils du modele ne
                         sont pas exactement les bons (abondances, profil
                         instrumental, granulation), et c'est lui qui sert de
                         reference de vitesse.

Donc: find_mask_lines sur le modele pour les bornes, puis depth, value et
line_snr relus dans le gabarit mesure a ces positions, par la formule meme de
LBL (general.py:3288-3306): depth = 1 - f / moyenne des deux voisins, et
line_snr = |depth * flux / rms|.
"""
import sys

import numpy as np
from astropy.io import fits
from astropy.table import Table

sys.path.insert(0, "/Users/eartigau/GitHubProjects/pca2d_preclean")


def build(model_template, measured_template, out_path, snr_min=20.0):
    from pca2d import readers

    inst = readers.instrument("HARPS", "ESO")
    general = readers._lbl_science("general")

    model = Table(fits.getdata(model_template, "TEMPLATE"))
    lines = general.find_mask_lines(inst, model)
    centre = np.asarray(lines["ll_mask_s"], dtype=float)

    measured = fits.getdata(measured_template, "TEMPLATE")
    mw = np.asarray(measured["wavelength"], dtype=float)
    mf = np.asarray(measured["flux_savgol_d0"], dtype=float)
    md2 = np.asarray(measured["flux_savgol_d2"], dtype=float)
    mrms = np.asarray(measured["rms"], dtype=float)
    keep = np.isfinite(mw) & np.isfinite(mf)
    mw, mf, md2, mrms = mw[keep], mf[keep], md2[keep], mrms[keep]
    order = np.argsort(mw)
    mw, mf, md2, mrms = mw[order], mf[order], md2[order], mrms[order]

    # l'echantillon du gabarit mesure le plus proche de chaque raie du modele,
    # et ses deux voisins, pour la profondeur au sens de LBL
    index = np.clip(np.searchsorted(mw, centre), 1, mw.size - 2)
    f_here = mf[index]
    neighbours = 0.5 * (mf[index - 1] + mf[index + 1])
    with np.errstate(invalid="ignore", divide="ignore"):
        depth = 1.0 - f_here / neighbours
        snr = f_here / mrms[index]
        line_snr = np.abs(depth * snr)

    table = Table()
    table["ll_mask_s"] = centre
    table["ll_mask_e"] = np.asarray(lines["ll_mask_e"], dtype=float)
    # le poids est la courbure MESUREE, pas celle du modele: c'est elle qui
    # dit ce que la raie fait a une vitesse dans les donnees
    table["w_mask"] = md2[index]
    table["value"] = f_here
    table["depth"] = depth
    table["line_snr"] = line_snr

    # MASK_SNR_MIN, comme LBL l'applique a son propre masque
    good = np.isfinite(line_snr) & (line_snr > snr_min) & np.isfinite(table["w_mask"])
    table = table[good]
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU(table, name="MASK")]).writeto(
        out_path, overwrite=True)
    print("%s" % out_path)
    print("  %d raies du modele -> %d au-dessus de line_snr %.0f"
          % (len(centre), len(table), snr_min))
    d = np.asarray(table["depth"], dtype=float)
    s = np.asarray(table["line_snr"], dtype=float)
    print("  depth: mediane %+.4f, positives %.1f%%"
          % (np.median(d), 100 * np.mean(d > 0)))
    print("  line_snr: mediane %.1f" % np.median(s))
    return table


if __name__ == "__main__":
    build("/Users/eartigau/tmp/harps_test/phoenix_template.fits",
          "/Users/eartigau/tmp/harps_test/lbl/templates/"
          "Template_HD10700_HARPS_HARPS_ESO.fits",
          "/Users/eartigau/tmp/harps_test/mask_from_model.fits")
