"""Un gabarit LBL bati sur le modele PHOENIX de Goettingen.

POURQUOI. Le masque de raies de LBL vient du GABARIT (general.find_mask_lines:
il garde tous les points ou la derivee change de signe, donc les maxima ET les
minima), et le modele PHOENIX ne sert qu'a fixer la vitesse systemique. Un
gabarit fait du modele donne donc un masque independant des donnees, ce qui
tranche: si les vitesses s'ameliorent, le gabarit mesure etait en cause; si
elles ne bougent pas, il ne l'etait pas.

CE QU'IL FAUT FOURNIR, d'apres find_mask_lines et load_template:
    wavelength, flux, rms            la base
    flux_savgol_d0/d1/d2/d3          USE_SAVGOL_TEMPLATE est vrai: les raies
                                     sont cherchees dans la derivee lissee
    les colonnes odd/even            le gabarit de LBL les porte; on y met la
                                     meme chose, le modele n'a pas de parite
"""
import sys

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.signal import savgol_filter

MODELS = "/Users/eartigau/tmp/harps_test/lbl/models"
C_KMS = 299792.458


def air_to_vacuum(wave_nm):
    """Les longueurs d'onde PHOENIX sont dans le VIDE, donc rien a faire.

    Gardee comme fonction pour que l'hypothese soit ecrite plutot que
    supposee: les grilles HiRes de Goettingen sont dans le vide, comme le
    WAVEDATA_VAC_BARY des S2D de l'ESO auquel ce gabarit sera compare.
    """
    return wave_nm


def build(out_path, teff_file, wave_min=398.0, wave_max=702.0, dv=0.5,
          resolution=115000.0):
    wave = fits.getdata("%s/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits" % MODELS)
    flux = fits.getdata("%s/%s" % (MODELS, teff_file))
    wave = air_to_vacuum(np.asarray(wave, dtype=float) / 10.0)   # A -> nm
    flux = np.asarray(flux, dtype=float)

    inside = (wave >= wave_min) & (wave <= wave_max) & np.isfinite(flux)
    wave, flux = wave[inside], flux[inside]
    order = np.argsort(wave)
    wave, flux = wave[order], flux[order]

    # LA RESOLUTION DE L'INSTRUMENT, et non celle du modele. Un masque fait de
    # raies infiniment fines ne tomberait pas sur les raies que HARPS mesure:
    # elles y sont elargies par le profil instrumental, et leurs maxima et
    # minima se deplacent quand deux raies voisines se melangent. On convolue
    # donc le modele a R = 115000 avant d'y chercher quoi que ce soit.
    step_kms = C_KMS * np.median(np.diff(np.log(wave)))
    fwhm_kms = C_KMS / resolution
    sigma_pix = (fwhm_kms / 2.3548) / step_kms
    half = int(np.ceil(4 * sigma_pix))
    x = np.arange(-half, half + 1)
    kernel = np.exp(-0.5 * (x / sigma_pix) ** 2)
    kernel /= kernel.sum()
    flux = np.convolve(flux, kernel, mode="same")

    # la grille magique de LBL, log-uniforme au pas demande
    n = int(np.log(wave_max / wave_min) / (dv / C_KMS))
    grid = wave_min * np.exp(np.arange(n) * dv / C_KMS)
    grid = grid[(grid >= wave[0]) & (grid <= wave[-1])]
    model = np.interp(grid, wave, flux)

    # NORMALISER PAR LE CONTINU, parce que le gabarit de LBL est un spectre
    # normalise (son flux median vaut 1.005 sur tau Ceti) et que find_mask_lines
    # divise par la colonne rms pour obtenir un line_snr. Le continu est pris
    # comme l'enveloppe superieure, au 95e centile sur 500 km/s.
    window = int(500.0 / dv) | 1
    from scipy.ndimage import percentile_filter
    continuum = percentile_filter(model, 95, size=window, mode="nearest")
    norm = model / continuum

    # un rms CONSTANT: le modele n'a pas de bruit, et c'est le point. Il est
    # pose a la valeur qui donne au gabarit mesure son line_snr typique, pour
    # que MASK_SNR_MIN garde le meme sens des deux cotes.
    rms = np.full_like(norm, 0.0041)

    width = int(round((C_KMS / resolution) / dv)) | 1
    width = max(width, 7)
    table = Table()
    table["wavelength"] = grid
    table["flux"] = norm
    table["eflux"] = rms
    table["rms"] = rms
    for parity in ("odd", "even"):
        table["flux_%s" % parity] = norm
        table["eflux_%s" % parity] = rms
        table["rms_%s" % parity] = rms
    for name, prefix in (("", "flux"), ("odd_", "flux_odd"), ("even_", "flux_even")):
        for deriv in range(4):
            column = savgol_filter(norm, width, 3, deriv=deriv)
            table["%s_savgol_d%d" % (prefix, deriv)] = column
    hdu = fits.BinTableHDU(table, name="TEMPLATE")
    primary = fits.PrimaryHDU()
    primary.header["ORIGIN"] = ("PHOENIX", "modele de Goettingen, pas les donnees")
    primary.header["MODEL"] = (teff_file[:60], "le fichier PHOENIX utilise")
    primary.header["RESOLUT"] = (resolution, "convolue a cette resolution")
    primary.header["NORMWIN"] = (500.0, "km/s, fenetre du continu au 95e centile")
    sci = fits.BinTableHDU(Table({"FILENAME": np.array(["PHOENIX"]),
                                  "BERV": np.array([0.0])}), name="SCI_TABLE")
    fits.HDUList([primary, hdu, sci]).writeto(out_path, overwrite=True)
    finite = np.isfinite(norm)
    print("%s" % out_path)
    print("  %d points, %.1f-%.1f nm au pas de %.2f km/s"
          % (len(grid), grid[0], grid[-1], dv))
    print("  flux normalise: mediane %.4f, 5e centile %.4f"
          % (np.median(norm[finite]), np.percentile(norm[finite], 5)))


if __name__ == "__main__":
    build(sys.argv[1], "lte05300-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits")
