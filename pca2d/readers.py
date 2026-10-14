"""Lire un spectre de n'importe quel instrument, en passant par LBL.

« Tout ce qui passe dans le LBL peut être envoyé au PCA2D » (2026-10-14). LBL
sait déjà lire quinze combinaisons instrument/source: SPIRou et NIRPS en APERO,
CADC et ESO, HARPS et HARPS-N en ORIG, ESO et ESSP, ESPRESSO, CARMENES, SOPHIE,
CORALIE, MAROON-X, EXPRES, NEID, plus une classe générique. Réimplémenter cette
lecture ici reviendrait à la maintenir deux fois et à diverger; et puisque
pca2d a vocation à finir dans LBL, mieux vaut en emprunter dès maintenant le
chemin.

    read(path, instrument, data_source) -> le même dict que tfits.read_tfits

Ce que LBL donne nativement et ce qu'il faut en tirer:

    flux    inst.load_science_file(path), toujours 2D (ordres, pixels)
    wave    inst.get_wave_solution(path, data, header), en nm
    blaze   inst.load_blaze_from_science(...), normalisé par ordre
    meta    les mots-clés par leur nom LBL (KW_BERV, KW_MJDATE...), que chaque
            classe traduit pour son instrument

Ce que LBL ne donne pas, et ce qu'on en fait:

    recon   la transmission tellurique. APERO l'écrit dans le t.fits (Recon);
            HARPS et les autres ne l'ont pas, et LBL la calcule lui-même à
            partir de son modèle TAPAS (lbl/science/tellu_clean.py, le fichier
            tapas_lbl.fits qu'il télécharge). `tapas_transmission` ci-dessous
            l'évalue sur la grille du spectre, ce qui donne à pca2d le même
            objet que le Recon d'APERO.
    sky     l'émission du ciel. Aucun spectrographe visible ne la livre par
            pose; quality.max_sky_ratio est alors sans objet et le masque de
            ciel est laissé à l'instrument. None, et le reste du pipeline en
            tient compte.
    snr     le rapport signal sur bruit par ordre. Quand l'en-tête ne le porte
            pas, il est mesuré sur le flux lui-même (robuste, par ordre).

LES DEUX PARITÉS. pca2d ramène chaque pose à deux spectres, les ordres pairs et
les ordres impairs, parce que sur un échelle les ordres n et n+1 se recouvrent
alors que n et n+2 ne se touchent jamais: tout le spectre tient donc sans
redondance sur deux lignes. La règle est gardée pour tous les instruments
(décidé le 2026-10-14), mais elle est VÉRIFIÉE par `parity_is_safe`, qui mesure
le recouvrement réel et prévient quand l'hypothèse ne tient pas.
"""

from __future__ import annotations

import os
import warnings

import numpy as np

#: ce que LBL appelle les instruments, et la source de données associée, tels
#: que select.InstDict les épelle. Ce qui n'est pas ici n'est pas refusé: la
#: paire est passée telle quelle à LBL, qui dira ce qu'il en pense.
KNOWN = {
    "SPIROU": ("APERO", "CADC"),
    "NIRPS_HA": ("APERO", "CADC", "ESO"),
    "NIRPS_HE": ("APERO", "CADC", "ESO"),
    "HARPS": ("ORIG", "ESO", "ESSP"),
    "HARPSN": ("ORIG", "ESO", "ESSP"),
    "ESPRESSO": ("None",),
    "CARMENES": ("None",),
    "SOPHIE": ("None",),
    "CORALIE": ("None",),
    "MAROONX": ("RED", "BLUE"),
    "EXPRES": ("ORIG", "ESSP"),
    "NEID": ("ESSP",),
}
#: GENERIC existe comme classe dans LBL mais pas dans sa table instrument/source:
#: il se configure à la main, donc pca2d ne l'annonce pas ici. Il reste
#: utilisable en passant son nom directement à `instrument`.


def instrument(name, data_source, config_file=None):
    """La classe LBL qui lit cet instrument, instanciée.

    `parse=False` est ce qui permet de l'instancier depuis du code plutôt que
    depuis une ligne de commande; c'est déjà ainsi que lbl.check_profile
    vérifie un profil avant une campagne.
    """
    from lbl.instruments import select

    kwargs = dict(instrument=str(name), data_source=str(data_source))
    if config_file:
        kwargs["config_file"] = os.path.abspath(config_file)
    args = select.parse_args(
        ["INSTRUMENT", "DATA_SOURCE", "DATA_TYPE", "INPUT_FILE"],
        kwargs, __name__, parse=False)
    return select.load_instrument(args, plogger=None)


def tapas_transmission(inst, wave, model_dir, airmass=1.0):
    """La transmission tellurique sur `wave`, du modèle TAPAS de LBL.

    APERO écrit la transmission qu'il a divisée dans le t.fits (Recon), et
    pca2d s'en sert pour savoir où le spectre est aveugle. Les spectrographes
    visibles ne la livrent pas: LBL la modélise, et c'est ce modèle qu'on
    évalue ici, élevé à la masse d'air de la pose.

    Rend None quand le modèle n'est pas là, ce que l'appelant doit accepter:
    sans transmission, les coupures qui en dépendent sont simplement inactives.
    """
    try:
        from lbl.science import tellu_clean
    except Exception:                                          # noqa: BLE001
        return None
    try:
        splines = tellu_clean.get_tapas_spl(inst, model_dir)
    except Exception:                                          # noqa: BLE001
        return None
    if splines is None:
        return None
    # get_tapas_spl rend les splines des deux composantes de TAPAS, l'eau et
    # le reste; la transmission totale est leur produit, chacune à la masse
    # d'air de la pose
    out = np.ones_like(np.asarray(wave, dtype=float))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for spline in (splines if isinstance(splines, (list, tuple))
                       else [splines]):
            try:
                piece = np.asarray(spline(np.asarray(wave, dtype=float)),
                                   dtype=float)
            except Exception:                                  # noqa: BLE001
                return None
            out *= np.clip(piece, 0.0, 1.0) ** float(airmass)
    return out


def order_snr(flux, header_snr=None, n_orders=None):
    """Le rapport signal sur bruit par ordre.

    Celui de l'en-tête quand l'instrument l'écrit (APERO le fait, EXTSNxxx),
    mesuré sur le flux sinon: la médiane divisée par l'écart robuste des
    différences voisines, qui ne suppose rien sur le modèle de bruit.
    """
    if header_snr is not None and np.size(header_snr):
        return np.asarray(header_snr, dtype=float)
    flux = np.atleast_2d(np.asarray(flux, dtype=float))
    out = np.full(flux.shape[0], np.nan)
    for i, row in enumerate(flux):
        good = np.isfinite(row) & (row > 0)
        if good.sum() < 16:
            continue
        level = float(np.median(row[good]))
        step = np.diff(row[good])
        noise = 1.4826 * float(np.median(np.abs(step - np.median(step))))
        # la différence de deux voisins a sqrt(2) fois le bruit d'un seul
        out[i] = level / (noise / np.sqrt(2.0)) if noise > 0 else np.nan
    return out if n_orders is None else out[:n_orders]


def parity_is_safe(wave, tolerance=0.0):
    """(vrai si les ordres n et n+2 ne se recouvrent pas, le pire recouvrement).

    pca2d met les ordres pairs sur une ligne et les impairs sur l'autre, ce qui
    suppose que deux ordres de même parité ne mesurent jamais la même longueur
    d'onde. C'est vrai sur SPIRou et NIRPS; ailleurs il faut le vérifier plutôt
    que l'espérer, parce qu'un recouvrement ferait écrire deux mesures dans la
    même colonne de la grille.
    """
    wave = np.atleast_2d(np.asarray(wave, dtype=float))
    worst = 0.0
    for i in range(wave.shape[0] - 2):
        a, b = wave[i], wave[i + 2]
        if not (np.isfinite(a).any() and np.isfinite(b).any()):
            continue
        lo, hi = np.nanmax(a), np.nanmin(b)
        if np.isfinite(lo) and np.isfinite(hi):
            worst = max(worst, float(lo - hi))
    return worst <= float(tolerance), worst


def read(path, instrument_name, data_source, config_file=None,
         model_dir=None, inst=None):
    """Un spectre, dans la forme que pca2d attend, quel que soit l'instrument.

    Même dict que tfits.read_tfits: flux, wave, blaze, recon, sky, snr, meta.
    `inst` évite de réinstancier la classe LBL à chaque fichier d'une campagne.
    """
    from lbl.core import io as lbl_io                           # noqa: F401

    inst = inst or instrument(instrument_name, data_source, config_file)
    flux, header = inst.load_science_file(path)
    flux = np.atleast_2d(np.asarray(flux, dtype=float))
    wave = np.atleast_2d(np.asarray(
        inst.get_wave_solution(path, flux, header), dtype=float))
    try:
        blaze, _ = inst.load_blaze_from_science(path, flux, header,
                                                calib_directory=model_dir,
                                                normalize=True)
        blaze = np.atleast_2d(np.asarray(blaze, dtype=float))
    except Exception:                                           # noqa: BLE001
        blaze = np.ones_like(flux)

    def key(name, default=None):
        """Un mot-clé par son nom LBL, que chaque classe traduit pour son
        instrument: KW_BERV est BERV sur APERO et HIERARCH ESO QC BERV sur ESO."""
        try:
            return header[inst.params[name]]
        except Exception:                                       # noqa: BLE001
            return default

    airmass = key("KW_AIRMASS", 1.0)
    meta = {
        "filename": os.path.basename(path),
        "berv": float(key("KW_BERV", np.nan) or np.nan),
        "mjdate": float(key("KW_MJDATE", np.nan) or np.nan),
        "exptime": float(key("KW_EXPTIME", np.nan) or np.nan),
        "airmass": float(airmass if airmass is not None else np.nan),
        "object": str(key("KW_OBJNAME", "") or "").strip(),
        "instrument": str(instrument_name),
        "data_source": str(data_source),
    }
    meta["drsobjn"] = meta["object"]
    recon = (tapas_transmission(inst, wave, model_dir, meta["airmass"])
             if model_dir else None)
    return dict(flux=flux, wave=wave, blaze=blaze, recon=recon, sky=None,
                snr=order_snr(flux), meta=meta)


#: les instruments que pca2d lit avec son propre lecteur, parce que les t.fits
#: d'APERO portent des choses que LBL ne rend pas (le ciel dans OHLine, la
#: transmission dans Recon, le SNR par ordre dans l'en-tete) et que ce lecteur
#: les exploite. Tout le reste passe par LBL.
NATIVE = ("SPIROU", "NIRPS", "NIRPS_HA", "NIRPS_HE")


def reads_natively(config):
    """Faut-il lire ces spectres avec tfits, ou passer par LBL ?

    APERO est lu nativement (le t.fits porte le ciel et la transmission);
    HARPS, SOPHIE, ESPRESSO et les autres passent par LBL. `input.reader`
    tranche quand on veut forcer l'un ou l'autre.
    """
    asked = str((config.get("input") or {}).get("reader") or "auto").lower()
    if asked in ("tfits", "native"):
        return True
    if asked in ("lbl", "bridge"):
        return False
    instrument = str((config.get("input") or {}).get("instrument") or "").upper()
    source = str((config.get("lbl") or {}).get("data_source") or "").upper()
    # un t.fits d'APERO: notre lecteur en tire plus que LBL n'en donne
    return instrument in NATIVE and source in ("", "APERO", "CADC")


def read_any(path, config, inst=None):
    """Un spectre, par le lecteur que la configuration designe.

    C'est le seul point ou pca2d decide comment lire: cube, reconstruct et
    plotting passent par ici, donc ajouter un instrument ne touche qu'a
    `reads_natively`.
    """
    if reads_natively(config):
        from . import tfits as sptf
        return sptf.read_tfits(path)
    block = config.get("lbl") or {}
    name = block.get("instrument") or (config.get("input") or {}).get("instrument")
    source = block.get("data_source") or "ESO"
    models = (config.get("output") or {}).get("model_directory")
    return read(path, name, source, model_dir=models, inst=inst)
