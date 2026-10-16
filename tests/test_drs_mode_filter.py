"""Garder un seul mode du spectrographe, parce que deux modes sont deux
instruments.

Mesure le 2026-10-15 sur les 100 poses HARPS de HD 10700: 8 sont en mode EGGS,
la fibre de 1.4" a haute efficacite, et ce sont exactement les 8 dont
l'ajustement tellurique refuse tout exposant d'eau. Elles veulent le modele
TAPAS decale de +16.0 km/s la ou les poses HARPS le veulent a 0.0, et baisser
la resolution du modele de 115000 a 50000 n'y change rien: autre fibre, autre
solution en longueur d'onde.
"""
import numpy as np

from pca2d.config import DEFAULTS, cache_key, load_config
from pca2d.cube import _quality_reason


WAVE = np.linspace(500.0, 600.0, 64)
QUALITY = dict(DEFAULTS["quality"])


def meta(mode):
    return {"bjd": 2453718.5, "snr_band": 200.0, "median_flux": 1e4,
            "nan_fraction": 0.0, "berv": 10.0, "drs_mode": mode,
            "filename": "r.HARPS.2005-12-16T01:40:24_S2D_A.fits"}


def rejected(mode, quality):
    return _quality_reason(meta(mode), WAVE, quality)


def test_the_mode_asked_for_is_kept_and_the_other_is_not():
    quality = dict(QUALITY, drs_modes=["HARPS"])
    assert rejected("HARPS", quality) is None
    said = rejected("EGGS", quality)
    assert said and "EGGS" in said and "not one of the 1" in said


def test_several_modes_may_be_named():
    quality = dict(QUALITY, drs_modes=["HARPS", "EGGS"])
    assert rejected("HARPS", quality) is None
    assert rejected("EGGS", quality) is None
    assert rejected("HAM", quality) is not None


def test_whitespace_does_not_decide():
    """Les en-tetes FITS sont rembourres, et un mode rembourre est le meme."""
    quality = dict(QUALITY, drs_modes=["HARPS"])
    assert rejected("HARPS  ", quality) is None
    quality = dict(QUALITY, drs_modes=[" HARPS"])
    assert rejected("HARPS", quality) is None


def test_an_absent_mode_is_rejected_when_a_filter_is_set():
    """Le filtre est un choix explicite: si l'en-tete ne dit pas le mode, on
    ne peut pas affirmer que c'est le bon, et la raison le dit."""
    quality = dict(QUALITY, drs_modes=["HARPS"])
    said = rejected("", quality)
    assert said and "(unset)" in said


def test_no_filter_keeps_every_mode():
    """Le defaut, parce qu'un seul mode est le cas ordinaire et qu'un
    instrument qui n'ecrit pas cette cle ne doit pas perdre ses poses."""
    for value in (None, [], False):
        quality = dict(QUALITY, drs_modes=value)
        assert rejected("HARPS", quality) is None
        assert rejected("EGGS", quality) is None
        assert rejected("", quality) is None


def test_unset_does_not_re_key_the_cubes_already_built():
    """Le piege que la suite de tests a attrape: ajouter ce reglage avec un
    defaut None changeait la cle de TOUS les cubes existants."""
    base = load_config(None)
    plain = cache_key(base)
    base["quality"]["drs_modes"] = None
    assert cache_key(base) == plain, "unset ne change aucune cle"
    base["quality"]["drs_modes"] = ["HARPS"]
    keyed = cache_key(base)
    assert keyed != plain, "le filtre change quelles poses le cube contient"
    base["quality"]["drs_modes"] = ["HARPS", "EGGS"]
    assert cache_key(base) != keyed, "d'autres modes, un autre cube"


def test_the_harps_block_asks_for_harps_only():
    """La config LIVREE, pas seulement le mecanisme. load_config(None) ne lit
    que les defauts, donc il faut nommer le fichier du depot."""
    import os
    path = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                        "config.yaml")
    config = load_config(path, instrument="HARPS")
    assert config["quality"]["drs_modes"] == ["HARPS"]
    assert config["input"]["instrument"] == "HARPS"


def test_the_other_instruments_keep_every_mode():
    """Le filtre est un reglage de HARPS, pas une nouvelle regle generale:
    SPIRou et NIRPS n'ecrivent pas cette cle de la meme facon et perdraient
    toutes leurs poses."""
    import os
    path = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                        "config.yaml")
    for name in ("SPIROU", "NIRPS", "ESPRESSO"):
        config = load_config(path, instrument=name)
        assert config["quality"]["drs_modes"] is None, name
