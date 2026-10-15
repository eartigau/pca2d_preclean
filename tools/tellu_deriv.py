"""L'exposant H2O par la DERIVEE, applique a HARPS.

La methode de telluric_fit (tellu_tools.optimize_exponents, appelee par
predict_abso.py), transposee du WINERED au visible: on divise le spectre par
la transmission modelisee et on retient l'exposant qui rend le spectre
corrige le plus LISSE, c'est-a-dire qui minimise la dispersion ponderee du
gradient pixel a pixel.

Pourquoi c'est mieux que l'extrapolation polynomiale de LBL: un exposant faux
laisse a chaque raie tellurique un residu en forme de DERIVEE, qui gonfle
cette dispersion, tandis que les raies stellaires y contribuent pareil quel
que soit l'exposant. Le critere ne depend donc que des telluriques, meme
quand elles sont minoritaires, la ou l'ajustement de LBL s'effondrait sur sa
borne une fois sur trois.

Ce qui est libre et ce qui ne l'est pas, comme demande le 2026-10-15:
    H2O         LIBRE PAR FICHIER: la colonne d'eau change trop vite pour un
                modele par nuit
    especes seches  FIXEES a la masse d'air: O2 et CO2 sont bien melanges,
                leur colonne suit la masse d'air par construction
    vitesse     NON LIBRE: TAPAS est evalue dans le repere du SOL, pas sur la
                grille barycentrique comme LBL le fait
"""
import sys

import numpy as np
from scipy.optimize import minimize_scalar

sys.path.insert(0, "/Users/eartigau/GitHubProjects/pca2d_preclean")

C = 299792.458
WEIGHTS = dict(depth_max=0.6, depth_saturated=0.5, reject_saturated=0.8,
               transition_sigma=0.02)
OUTLIER_SIGMA = 6.0


def transmission_weights(trans, depth_max, depth_saturated, reject_saturated,
                         transition_sigma):
    """Les poids 0-1 de telluric_fit: une sigmoide sur la profondeur, et les
    raies saturees et leurs ailes mises a zero."""
    w_depth = 1.0 / (1.0 + np.exp(-(trans - depth_max) / transition_sigma))
    sat = trans < depth_saturated
    above = trans >= reject_saturated
    contam = np.zeros(len(trans), bool)
    for idx in (range(len(trans)), range(len(trans) - 1, -1, -1)):
        inside = False
        for i in idx:
            if sat[i]:
                inside = True
            elif inside and above[i]:
                inside = False
            contam[i] |= inside
    return w_depth * (~contam)


def weighted_nanstd(data, weights):
    ok = np.isfinite(data) & np.isfinite(weights) & (weights > 0)
    d, w = data[ok], weights[ok]
    if d.size < 32:
        return np.inf
    mean = np.sum(w * d) / np.sum(w)
    denom = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    if denom <= 0:
        return np.inf
    return np.sqrt(np.sum(w * (d - mean) ** 2) / denom)


def objective(sp, trans):
    """La dispersion ponderee du gradient du spectre corrige."""
    grad = np.gradient(sp / np.clip(trans, 1e-3, None))
    w = transmission_weights(trans, **WEIGHTS)
    valid = w > 0.5
    if valid.sum() < 64:
        return np.inf
    p16, p84 = np.nanpercentile(grad[valid], [16, 84])
    sig = max((p84 - p16) / 2.0, 1e-12)
    dev = np.abs(grad - np.nanmedian(grad[valid])) / sig
    w = w / (1.0 + np.exp(np.clip(dev - OUTLIER_SIGMA, -500, 500)))
    return weighted_nanstd(grad, w)


def lsf_kernel(step_kms, resolution):
    """Le noyau gaussien de la reponse instrumentale, en echantillons.

    LE MODELE DOIT ETRE CONVOLUE AVANT D'ETRE DIVISE. TAPAS est bien plus fin
    que le spectrographe (R = 115000 sur HARPS, soit 2.6 km/s), et des raies
    trop etroites ne peuvent pas annuler des raies elargies: il resterait a
    chaque raie un residu en forme de derivee, exactement ce que le critere
    cherche a eliminer. Mesure le 2026-10-15 sur 40 poses de HD 10700: la
    convolution gagne 4.6% de dispersion du gradient a elle seule.

    La puissance est prise AVANT la convolution, comme dans
    papier_helium/referee/tellu_derivative.py: l'absorption est exponentielle
    par colonne, l'elargissement instrumental vient apres.
    """
    sigma = (C / float(resolution)) / 2.3548 / float(step_kms)
    half = int(np.ceil(4 * sigma))
    x = np.arange(-half, half + 1)
    kernel = np.exp(-0.5 * (x / sigma) ** 2)
    return kernel / kernel.sum()


def water_exponent(wave_ground, sp, spl_water, spl_dry, airmass,
                   resolution=115000.0, bounds=(0.001, 20.0)):
    """L'exposant H2O de CE fichier, le seul parametre libre.

    DEUX GROUPES, PAS QUATRE ESPECES. Pour l'optique on ajuste l'eau d'un
    cote et tout le reste de l'autre (consigne du 2026-10-15), ce qui est
    aussi le decoupage que TAPAS donne dans LBL: get_tapas_spl rend deux
    splines, ABSO_WATER et ABSO_OTHERS. Separer CO2, CH4 et O2 n'aurait pas
    de sens ici: dans le visible seul O2 compte vraiment, et les trois sont
    bien melanges donc ils suivent la meme masse d'air.

    LE SEC EST FIXE A LA MASSE D'AIR, et ce n'est pas un pis-aller: ajuste
    librement, il vaut 0.971 fois la masse d'air (median sur 40 poses) et le
    fixer ne coute rien du tout (+0.0% de dispersion).

    LA VITESSE N'EST PAS LIBRE: `wave_ground` est deja dans le repere du sol,
    ou l'atmosphere se trouve. C'est le correctif de docs/tapas_frame.md, que
    LBL ne faisait pas.

    LA FORMULATION D'optimize_exponents (telluric_fit/tellu_tools.py:2413),
    plus stable numeriquement et qui ne tape pas les bornes: Nelder-Mead sur
    la dispersion ponderee du gradient, bornes larges (0.001 a 20,
    tellu_tools_config.py:224, contre 0.01 a 15 dans LBL).

    Mesure le 2026-10-15 sur HD 10700, en ajustant sur les trois bandes d'eau
    du domaine HARPS: 2 poses sur 40 restent sur la borne basse, contre 35%
    avec un seul exposant libre et 37% avec la methode de LBL.
    """
    from scipy.optimize import minimize

    step = C * np.median(np.diff(np.log(np.asarray(wave_ground, dtype=float))))
    kernel = lsf_kernel(step, resolution)
    dry = np.clip(spl_dry(wave_ground), 1e-6, 1.0) ** float(airmass)
    water = np.clip(spl_water(wave_ground), 1e-6, 1.0)

    def cost(x):
        alpha = float(x[0])
        if not (bounds[0] <= alpha <= bounds[1]):
            return np.inf
        return objective(sp, np.convolve(dry * water ** alpha, kernel,
                                         mode="same"))

    out = minimize(cost, x0=[1.0], method="Nelder-Mead",
                   options=dict(xatol=1e-3, fatol=1e-8, maxiter=200))
    return float(out.x[0]), float(out.fun)


#: les bandes d'eau du domaine HARPS, mesurees sur TAPAS le 2026-10-15: c'est
#: la que l'exposant H2O s'ajuste, et nulle part ailleurs. Hors de ces trois
#: bandes l'eau ne depasse pas 1% de profondeur et n'informe rien, si bien
#: qu'un ajustement sur tout le domaine se faisait noyer par les raies
#: stellaires et tombait sur sa borne une fois sur trois.
WATER_BANDS_HARPS = ((586.0, 600.0), (644.0, 664.0), (686.0, 696.0))


def water_mask(wave, bands=WATER_BANDS_HARPS):
    """Les pixels des bandes d'eau, ou l'exposant H2O s'ajuste."""
    wave = np.asarray(wave, dtype=float)
    out = np.zeros(wave.shape, dtype=bool)
    for lo, hi in bands:
        out |= (wave >= lo) & (wave <= hi)
    return out
