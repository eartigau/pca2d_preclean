"""The star as one B-spline: carried by evaluation, updated exactly."""
import numpy as np
import pytest
from astropy.table import Table

from pca2d.splinestar import (SplineStar, coefficients, evaluate, solve_components,
                              values)
from pca2d.twoframe import LanczosShifter, main

N_PIX = 400


def _line(x0=200.0, sigma=5.0, depth=-0.3):
    x = np.arange(N_PIX, dtype=float)
    return depth * np.exp(-0.5 * ((x - x0) / sigma) ** 2)


def test_values_and_coefficients_are_inverse():
    v = _line()[None, :]
    assert np.allclose(values(coefficients(v))[:, 5:-5], v[:, 5:-5], atol=1e-12)


def test_no_shift_gives_the_spectrum_back():
    s = SplineStar(LanczosShifter(N_PIX, a=8, max_shift=20))
    v = _line()[None, :]
    out = s.carry(s.prepare(v), np.array([0.0]))[0]
    assert np.allclose(out[:, 5:-5], v[:, 5:-5], atol=1e-12)


def test_a_shift_moves_the_line_the_way_the_lanczos_kernel_does():
    lanczos = LanczosShifter(N_PIX, a=8, max_shift=20)
    s = SplineStar(lanczos)
    v = _line()[None, :]
    p = np.array([3.3])
    ours = s.carry(s.prepare(v), p)[0, 0]
    theirs = lanczos.carry(lanczos.prepare(v), p)[0, 0]
    assert abs(int(np.argmin(ours)) - 203) <= 1
    assert np.max(np.abs(ours - theirs)) < 2e-3 * 0.3


def test_the_data_side_is_still_the_lanczos_kernel():
    lanczos = LanczosShifter(N_PIX, a=8, max_shift=20)
    s = SplineStar(lanczos)
    rows = np.vstack([_line(), _line(150)])
    pix = np.array([2.4, -1.7])
    assert np.array_equal(s.rows(rows, pix), lanczos.rows(rows, pix))
    assert np.array_equal(s.adjoint(rows, pix), lanczos.adjoint(rows, pix))


def test_the_update_is_the_exact_least_squares_solution():
    r = np.random.default_rng(5)
    n_rows, pad = 10, 12
    truth = _line() + _line(100, 3, -0.2) + 0.01 * r.normal(size=N_PIX)
    truth_c = coefficients(truth[None, :])[0]
    delta = r.uniform(-6, 6, n_rows)
    a = r.uniform(0.5, 2.0, (n_rows, 1))
    w = r.uniform(0.5, 2.0, (n_rows, N_PIX))
    w[:, :15] = w[:, -15:] = 0.0
    padded = np.pad(truth_c[None, :], [(0, 0), (pad, pad)])
    data = a * evaluate(padded, delta, N_PIX, pad)[:, 0, :]
    data += 0.001 * r.normal(size=data.shape)
    # the same problem as dense matrices, solved by numpy
    normal, rhs = np.zeros((N_PIX, N_PIX)), np.zeros(N_PIX)
    for i in range(n_rows):
        E = np.zeros((N_PIX, N_PIX))
        for j in range(N_PIX):
            unit = np.zeros((1, N_PIX + 2 * pad))
            unit[0, pad + j] = 1.0
            E[:, j] = evaluate(unit, delta[i:i + 1], N_PIX, pad)[0, 0]
        normal += a[i, 0] ** 2 * E.T @ (w[i][:, None] * E)
        rhs += a[i, 0] * E.T @ (w[i] * data[i])
    diag = np.diag(normal)
    keep = diag > 1e-2 * np.median(diag[diag > 0])
    c = np.zeros(N_PIX)
    c[keep] = np.linalg.solve(normal[np.ix_(keep, keep)], rhs[keep])
    dense = values(c)[0]
    dense /= np.linalg.norm(dense)
    ours = solve_components(data.copy(), w, a, delta, pad)[0]
    assert np.allclose(ours, dense, atol=1e-7)
    inner = slice(30, -30)
    assert np.corrcoef(ours[inner], truth[inner])[0, 1] > 0.999


@pytest.fixture(scope="module")
def cube(tmp_path_factory):
    r = np.random.default_rng(7)
    n_rows, n_pix, dv = 48, 2000, 0.5
    x = np.arange(n_pix, dtype=float)
    grid = 1500.0 * np.exp(np.arange(n_pix) * dv / 299792.458)
    star = -0.3 * sum(np.exp(-0.5 * ((x - c) / 5.0) ** 2) for c in (300, 900, 1500))
    earth = -0.2 * sum(np.exp(-0.5 * ((x - c) / 4.0) ** 2) for c in (600, 1200))
    berv = np.linspace(-12.0, 12.0, n_rows // 2).repeat(2)
    data = np.zeros((n_rows, n_pix))
    for i in range(n_rows):
        data[i] = np.roll(star * (1.0 + 0.05 * r.normal()), int(round(-berv[i] / dv)))
        data[i] += earth * (1.0 + 0.3 * r.normal())
    data += 0.005 * r.normal(size=data.shape)
    path = tmp_path_factory.mktemp("cube")
    np.save(path / "grid.npy", grid)
    np.save(path / "data.npy", data)
    np.save(path / "sigma.npy", np.full_like(data, 0.005))
    meta = Table()
    meta["filename"] = ["%04dt.fits" % (i // 2) for i in range(n_rows)]
    meta["bjd"] = 2459000.0 + np.arange(n_rows) * 1.7
    meta["berv"] = berv
    meta["airmass"] = np.full(n_rows, 1.2)
    meta["snr_band"] = np.full(n_rows, 100.0)
    meta["exposure"] = np.arange(n_rows) // 2
    meta["parity"] = np.arange(n_rows) % 2
    meta.write(path / "meta.fits", overwrite=True)
    return str(path)


def _fit(cube, out, *extra):
    main(["--cube", cube, "--outdir", str(out), "--iters", "4", "-k", "1", "-j", "2",
          "--max-mad", "0", "--mean", "star", "--keep", "last", "--patience", "99",
          *extra])
    return np.load(out / "fit.npz")


def test_the_fit_runs_with_the_spline_star_and_fits_as_well(cube, tmp_path):
    grid_fit = _fit(cube, tmp_path / "grid", "--star-basis", "grid")
    spline_fit = _fit(cube, tmp_path / "spline", "--star-basis", "spline")
    assert str(spline_fit["star_basis"]) == "spline"
    assert np.all(np.isfinite(spline_fit["P"]))
    assert float(spline_fit["chi2_best"]) <= 1.02 * float(grid_fit["chi2_best"])
