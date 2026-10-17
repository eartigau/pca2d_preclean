"""Sur un cube empile par nuit, chaque pose recoit les coupures de SA nuit.

Avant le 2026-10-17, l'etape correct les cherchait sous le nom du fichier, et
la ligne d'une nuit porte celui de sa premiere pose: seule celle-ci les
recevait (GJ 1: 2658 echantillons coupes sur la premiere pose, 0 sur les
autres)."""

import os

import numpy as np
from astropy.table import Table

from pca2d.reconstruct import cut_key, night_row_names


def _cube(tmp_path):
    # deux nuits, deux parites chacune, la ligne nommee d'apres la 1re pose
    meta = Table({
        "filename": ["NIRPS.A1t.fits", "NIRPS.A1t.fits",
                     "NIRPS.B1t.fits", "NIRPS.B1t.fits"],
        "bjd": [2460000.70, 2460000.70, 2460001.80, 2460001.80],
    })
    meta.write(os.path.join(tmp_path, "meta.fits"))
    return str(tmp_path)


def test_each_night_is_named_after_its_row(tmp_path):
    names = night_row_names(_cube(tmp_path))
    assert names == {2460000: "NIRPS.A1t.fits", 2460001: "NIRPS.B1t.fits"}


def test_every_exposure_of_a_night_gets_the_nights_cuts(tmp_path):
    names = night_row_names(_cube(tmp_path))
    cuts = {"NIRPS.A1t.fits": np.array([True, False]),
            "NIRPS.B1t.fits": np.array([False, True])}
    # trois poses de la premiere nuit: la ligne de coefficients est celle de
    # la nuit, donc son BJD est celui de la nuit
    night_row = {"bjd": 2460000.70}
    for path in ("/d/NIRPS.A1t.fits", "/d/NIRPS.A2t.fits", "/d/NIRPS.A3t.fits"):
        assert cut_key(path, night_row, names) == "NIRPS.A1t.fits"
        assert cuts[cut_key(path, night_row, names)].tolist() == [True, False]
    assert cut_key("/d/NIRPS.B2t.fits", {"bjd": 2460001.80}, names) == "NIRPS.B1t.fits"


def test_unstacked_runs_keep_the_files_own_name():
    assert cut_key("/d/2811170t.fits", {"bjd": 2459890.9}, None) == "2811170t.fits"


def test_a_night_the_cube_does_not_know_falls_back_to_the_file(tmp_path):
    names = night_row_names(_cube(tmp_path))
    assert cut_key("/d/X.fits", {"bjd": 2460050.7}, names) == "X.fits"
