"""Les options de diagnostic de la nuit du 2026-10-17 (pca2d.diagnostics),
l'outil d'injection, et leur passage de la configuration aux etapes."""

import importlib.util
import os

import numpy as np
import pytest

from pca2d import diagnostics as dg

HERE = os.path.dirname(os.path.abspath(__file__))


def _inject_module():
    path = os.path.join(HERE, "..", "tools", "inject_planet.py")
    spec = importlib.util.spec_from_file_location("inject_planet", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ------------------------------------------------------------------ holdout
def test_holdout_alternates_nights_and_keeps_parities_together():
    # trois nuits, deux lignes (parites) chacune, dans le desordre
    bjd = np.array([60002.7, 60000.6, 60001.9, 60000.61, 60002.71, 60001.91])
    odd = dg.holdout_rows(bjd, "odd")
    even = dg.holdout_rows(bjd, "even")
    assert np.all(odd ^ even)                       # une partition
    nights = dg.night_of(bjd)
    for n in np.unique(nights):
        assert len(set(odd[nights == n])) == 1      # une nuit, un sort
    # nuits 60000, 60001, 60002 de rangs 0, 1, 2: la nuit du milieu est impaire
    assert odd[nights == 60001].all() and not odd[nights == 60000].any()


def test_holdout_refuses_other_words():
    with pytest.raises(ValueError):
        dg.holdout_rows([60000.7], "half")


# -------------------------------------------------------------- permutation
def test_derangement_never_serves_a_night_its_own():
    nights = np.repeat(np.arange(50), 3)            # 50 nuits de 3 poses
    source = dg.night_derangement(nights, seed=17)
    assert source.shape == nights.shape
    assert np.all(nights[source] != nights)
    assert source.min() >= 0 and source.max() < nights.size


def test_derangement_is_reproducible():
    nights = np.repeat(np.arange(20), 2)
    assert np.array_equal(dg.night_derangement(nights, 3),
                          dg.night_derangement(nights, 3))
    assert not np.array_equal(dg.night_derangement(nights, 3),
                              dg.night_derangement(nights, 4))


def test_derangement_needs_two_nights():
    with pytest.raises(ValueError):
        dg.night_derangement([5, 5, 5], 1)


# ------------------------------------------------------------------ support
def test_support_is_where_any_row_absorbs():
    trans = np.ones((200, 6))
    trans[150, 1] = 0.98                            # une pose, au-dela d'un bloc
    trans[:, 2] = 0.999
    trans[:, 3] = np.nan                            # rien de connu
    trans[7, 4] = np.nan
    trans[8, 4] = 0.5
    keep = dg.support_columns(trans, 0.995)
    assert keep.tolist() == [False, True, False, False, True, False]


def test_support_refuses_a_threshold_outside_the_unit_interval():
    with pytest.raises(ValueError):
        dg.support_columns(np.ones((2, 2)), 1.5)


# ---------------------------------------------------------------- injection
def test_injection_doppler_is_relativistic():
    inj = _inject_module()
    v = 30000.0
    beta = v / 299792458.0
    assert inj.doppler(v) == pytest.approx(np.sqrt((1 + beta) / (1 - beta)),
                                           rel=1e-15)
    # et pas le premier ordre: l'ecart est (v/c)^2 / 2
    assert abs(inj.doppler(v) - (1 + beta)) == pytest.approx(beta ** 2 / 2,
                                                             rel=1e-3)


def test_injection_velocity_is_the_sum_of_its_planets():
    inj = _inject_module()
    t = np.linspace(60000, 60400, 7)
    v = inj.planets_velocity(t, [(37.3, 3.0, 0.4), (365.25, 2.0, 1.1)])
    expect = (3.0 * np.sin(2 * np.pi * t / 37.3 + 0.4)
              + 2.0 * np.sin(2 * np.pi * t / 365.25 + 1.1))
    assert np.allclose(v, expect)


# ------------------------------------------------- configuration -> etapes
def test_options_reach_their_stage():
    from pca2d.cli import diagnostic_args
    from pca2d.reconstruct import parse_args as correct_args
    from pca2d.twoframe import CONFIGURABLE, parse_args as fit_args

    assert diagnostic_args({"correct": {}}) == []
    out = diagnostic_args({"correct": {"q_support": 0.995, "permute_b": 17}})
    assert out == ["--q-support", "0.995", "--permute-b", "17"]
    args = correct_args(["--fits", "x.fits", "--file", "y.fits", *out])
    assert args.q_support == 0.995 and args.permute_b == 17

    assert "holdout" in CONFIGURABLE
    assert fit_args(["--cube", "c", "--outdir", "o", "--holdout", "odd"]).holdout == "odd"


def test_new_knobs_do_not_rekey_the_cube():
    from pca2d.config import DEFAULTS, cache_key
    import copy
    cfg = copy.deepcopy(DEFAULTS)
    before = cache_key(cfg)
    cfg["twoframe"]["holdout"] = "odd"
    cfg["correct"]["q_support"] = 0.995
    cfg["correct"]["permute_b"] = 3
    assert cache_key(cfg) == before
