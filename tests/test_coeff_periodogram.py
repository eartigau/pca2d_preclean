"""The coefficient periodograms, including the fit that has no star coefficient.

n_star 0 is the nominal since 2026-09-12, and this figure asked the star
coefficients for the dates, so it failed on every run made with it.
"""
import numpy as np

from pca2d.figures.coeff_periodogram import main


def write_fit(path, n_star, n_earth, n=60):
    rng = np.random.default_rng(3)
    bjd = 2460000.0 + np.arange(n) * 3.0
    np.savez(path, bjd=bjd, rejected=np.zeros(n, bool),
             a=rng.normal(size=(n, n_star)), b=rng.normal(size=(n, n_earth)))


def test_a_fit_with_no_star_coefficient_still_draws_its_observer_ones(tmp_path):
    fit = str(tmp_path / "fit.npz")
    out = str(tmp_path / "p.pdf")
    write_fit(fit, n_star=0, n_earth=3)
    main(["--fit", fit, "--out", out, "--bootstrap", "0"])
    assert open(out, "rb").read(4) == b"%PDF"


def test_the_star_coefficients_are_drawn_when_there_are_some(tmp_path):
    fit = str(tmp_path / "fit.npz")
    out = str(tmp_path / "p.pdf")
    write_fit(fit, n_star=2, n_earth=3)
    main(["--fit", fit, "--out", out, "--bootstrap", "0", "--planets", "5.8"])
    assert open(out, "rb").read(4) == b"%PDF"


def test_a_fit_with_no_component_at_all_says_so_and_draws_nothing(tmp_path):
    fit = str(tmp_path / "fit.npz")
    out = str(tmp_path / "p.pdf")
    write_fit(fit, n_star=0, n_earth=0)
    main(["--fit", fit, "--out", out, "--bootstrap", "0"])
    import os
    assert not os.path.exists(out), "no panels is not an empty figure"
