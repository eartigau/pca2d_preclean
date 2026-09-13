"""A corrected file holds what panel 3 of the sequence figure shows.

Panel 3, "minus the OBSERVER block: what a corrected file holds", is the cube's
data minus the fit's per-parity mean minus the observer block, with every
sample the fit gave no weight to hidden. On 2026-09-10 the files were found to
keep the parity mean the panel removes (slanted stripes at 1267 nm on TOI2120),
and the panel was made the reference for the corrected spectra. This pins the
two together: what the correction divides out of an order is exactly what
panel 3 subtracts from that order's parity, and the samples it blanks are
exactly the ones the fit did not weight.
"""

import os

import numpy as np
import pytest
from astropy.table import Table

from pca2d.plotting import live_mask
from pca2d.reconstruct import (correction_on_grid, fit_weights_mask, load_model,
                               order_correction)
from pca2d.twoframe import fit_means, load_cube, main, mean_rows, row_parity

N_ROWS, N_PIX, DV = 12, 2000, 0.5


@pytest.fixture(scope="module")
def fitted(tmp_path_factory):
    r = np.random.default_rng(4)
    x = np.arange(N_PIX, dtype=float)
    grid = 1500.0 * np.exp(np.arange(N_PIX) * DV / 299792.458)
    star = -0.3 * sum(np.exp(-0.5 * ((x - c) / 5.0) ** 2) for c in (300, 900, 1500))
    earth = -0.2 * sum(np.exp(-0.5 * ((x - c) / 4.0) ** 2) for c in (600, 1200))
    berv = np.linspace(-12.0, 12.0, N_ROWS // 2).repeat(2)
    data = np.zeros((N_ROWS, N_PIX))
    for i in range(N_ROWS):
        data[i] = np.roll(star, int(round(-berv[i] / DV)))
        data[i] += earth * (1.0 + 0.3 * r.normal()) + 0.02 * (i % 2)
    data += 0.01 * r.normal(size=data.shape)
    data[3, 700:720] = np.nan              # a hole the fit cannot weight
    path = tmp_path_factory.mktemp("cube")
    np.save(path / "grid.npy", grid)
    np.save(path / "data.npy", data)
    np.save(path / "sigma.npy", np.full_like(data, 0.01))
    meta = Table()
    meta["filename"] = ["%04dt.fits" % (i // 2) for i in range(N_ROWS)]
    meta["bjd"] = 2459000.0 + np.arange(N_ROWS) * 1.7
    meta["berv"] = berv
    meta["airmass"] = np.full(N_ROWS, 1.2)
    meta["snr_band"] = np.full(N_ROWS, 100.0)
    meta["exposure"] = np.arange(N_ROWS) // 2
    meta["parity"] = np.arange(N_ROWS) % 2
    meta.write(path / "meta.fits", overwrite=True)
    out = tmp_path_factory.mktemp("fit")
    main(["--cube", str(path), "--outdir", str(out), "--iters", "3",
          "-k", "1", "-j", "2", "--max-mad", "0"])
    return str(path), out


def test_what_is_divided_out_is_what_panel_3_subtracts(fitted):
    cube, out = fitted
    fit = np.load(out / "fit.npz")
    model = load_model(str(out / "twoframe_components.fits"))
    _, _, _, meta = load_cube(cube, columns=np.arange(1))
    n, grid = len(meta), model["grid"]
    parity = row_parity(meta, n)
    means, group = fit_means(fit, meta, n, grid.size)
    # sequence.window_arrays, panel 3: data - offset - earth
    subtracted = np.asarray(mean_rows(means, group, n)) + fit["b"] @ fit["Q"]
    coeffs = {str(row["filename"]): row for row in model["coeffs"]}
    inner = slice(100, -100)
    for i in range(n):
        row = coeffs[str(meta["filename"][i])]
        correction, _, j = correction_on_grid(model, row, 0, model["n_earth"])
        values, live = order_correction(model, correction, j, int(parity[i]),
                                        grid[inner])
        assert live.any()
        assert np.allclose(values[live], subtracted[i][inner][live],
                           rtol=0, atol=1e-9), i


def test_the_samples_blanked_are_the_samples_the_fit_did_not_weight(fitted):
    cube, _ = fitted
    alive = fit_weights_mask(cube)
    _, _, w0, meta = load_cube(cube)
    parity = row_parity(meta, len(meta))
    for i in range(len(meta)):
        want = live_mask(w0[i:i + 1])[0]
        assert np.array_equal(alive[str(meta["filename"][i])][int(parity[i]) % 2],
                              want), i
    assert not alive["0001t.fits"][1][700:720].any(), "the hole must be blanked"


def test_panel_6_is_what_panel_3_took_out_and_it_shrinks(fitted):
    """Panel 6 shows the correction the files have divided out, and panel 3 is
    the data less exactly that, whether the correction was shrunk or not."""
    from pca2d.figures.sequence import load_context, window_arrays
    cube, out = fitted
    applied = {}
    for shrink in (False, True):
        ctx = load_context(cube, str(out / "fit.npz"), shrink=shrink)
        arr = window_arrays(ctx["cube"], ctx["fit"], ctx["means"], ctx["group"],
                            ctx["grid"], ctx["dv"], ctx["delta"],
                            float(ctx["grid"][1000]), 1.0,
                            templates=ctx["templates"], correct=ctx["correct"])
        home = arr["home"]
        ok = np.isfinite(home["given"]) & np.isfinite(home["applied"])
        assert np.allclose(home["corrected"][ok], (home["given"] - home["applied"])[ok],
                           atol=1e-5)
        applied[shrink] = np.abs(home["applied"][ok]).sum()
    assert applied[True] < applied[False], "shrinking did not take anything off"


def test_a_shrunk_file_divides_out_panel_6(fitted):
    """With --shrink, what the correct stage divides out of an order is panel 6
    as the figure computes it: from fit.npz's amplitudes and chi2 and the
    cube's weights, where the correct stage reads the components file's
    coefficient table and the weights in float32. Checked on a real TOI-2120
    file on 2026-09-11 (largest difference 8e-11 in ln f); this pins it."""
    from pca2d.reconstruct import shrink_model
    from pca2d.shrink import correction_basis
    cube, out = fitted
    fit = np.load(out / "fit.npz")
    fits_path = str(out / "twoframe_components.fits")
    model = load_model(fits_path)
    shrink_model(model, fits_path, cube)
    _, _, w0, meta = load_cube(cube)
    n, grid = len(meta), model["grid"]
    parity = row_parity(meta, n)
    chi2 = np.asarray(fit["chi2_red"], dtype=float)
    scale = float(np.median(chi2[np.isfinite(chi2) & (chi2 > 0)]))
    Q_figure, _ = correction_basis(fit["Q"], fit["b"], w0, scale, True)
    assert not np.allclose(Q_figure, fit["Q"]), "nothing was shrunk: the test is empty"
    means, group = fit_means(fit, meta, n, grid.size)
    panel6 = np.asarray(mean_rows(means, group, n)) + fit["b"] @ Q_figure
    coeffs = {str(row["filename"]): row for row in model["coeffs"]}
    inner = slice(100, -100)
    for i in range(n):
        row = coeffs[str(meta["filename"][i])]
        correction, _, j = correction_on_grid(model, row, 0, model["n_earth"])
        values, live = order_correction(model, correction, j, int(parity[i]),
                                        grid[inner])
        assert live.any()
        assert np.allclose(values[live], panel6[i][inner][live], rtol=0, atol=1e-7), i


def test_a_common_mask_hides_in_every_row_what_one_exposure_lost(fitted):
    """With correct --mask common a file blanks every sample any exposure left
    unweighted, so the panels hide those in every row, not only in the row that
    lost them."""
    from pca2d.figures.sequence import load_context, window_arrays
    cube, out = fitted
    hidden = {}
    for mode in ("exposure", "common"):
        ctx = load_context(cube, str(out / "fit.npz"), mask=mode)
        arr = window_arrays(ctx["cube"], ctx["fit"], ctx["means"], ctx["group"],
                            ctx["grid"], ctx["dv"], ctx["delta"],
                            float(ctx["grid"][710]), 1.0,
                            templates=ctx["templates"], correct=ctx.get("correct"))
        hidden[mode] = ~np.isfinite(arr["home"]["corrected"])
    assert hidden["common"].sum() > hidden["exposure"].sum(), \
        "the hole one exposure has is not hidden in the others"
    assert not (hidden["exposure"] & ~hidden["common"]).any(), \
        "a sample the per-exposure mask hides must stay hidden"


def test_the_common_mask_of_a_panel_is_taken_within_one_parity(fitted):
    """A wavelength is reached by ONE order parity: orders n and n+2 do not
    overlap. Intersected over every row, the mask asked a sample to be alive in
    rows that never cover it, and was empty everywhere the parities do not meet:
    the three H-band windows of the joint run vanished from the report.

    This calls window_arrays itself, because the first fix was written against a
    reimplementation of the rule and shipped with the grid shadowed by the loop
    variable, which no such test could see.
    """
    from pca2d.figures.sequence import load_context, window_arrays
    cube, out = fitted
    ctx = load_context(cube, str(out / "fit.npz"), mask="common")
    # a sample only one parity can see, as a wavelength only one order reaches
    w = np.load(os.path.join(cube, "sigma.npy"))
    w[0::2, 300:340] = 0.0
    np.save(os.path.join(cube, "sigma.npy"), w)
    ctx = load_context(cube, str(out / "fit.npz"), mask="common")
    arr = window_arrays(ctx["cube"], ctx["fit"], ctx["means"], ctx["group"],
                        ctx["grid"], ctx["dv"], ctx["delta"],
                        float(ctx["grid"][320]), 1.0,
                        templates=ctx["templates"], correct=ctx.get("correct"))
    assert arr is not None
    assert np.ndim(arr["grid"]) == 1 and arr["grid"].size > 10, \
        "the grid comes back as the grid, not as a group label"
    alive = np.isfinite(arr["home"]["corrected"])
    odd = np.asarray(ctx["group"]) % 2 == 1
    assert alive[odd].any(), "the parity that covers this window still has rows"


def test_the_figure_runs_as_the_bundle_runs_it(fitted, tmp_path):
    """The bundle runs sequence.py as a script, where a relative import fails:
    on 2026-09-11 one did, and the report came out without its sequence pages."""
    import os
    import subprocess
    import sys

    import pca2d.figures.sequence as seq
    cube, out = fitted
    grid = np.load(os.path.join(cube, "grid.npy"))
    pdf = tmp_path / "sequence.pdf"
    r = subprocess.run([sys.executable, seq.__file__, "--cube", cube,
                        "--fit", str(out / "fit.npz"),
                        "--windows", "%.3f:1" % grid[1000], "--out", str(pdf),
                        "--shrink", "--shrink-smooth", "--smooth-components", "1",
                        "--resolution", "70000"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-800:]
    assert pdf.exists()


def test_a_joint_cube_gets_one_page_per_star_and_one_of_them_all():
    """Three campaigns on one axis, ordered by BERV, interleave three stars'
    lines at three systemic velocities: the H-band pages of the first joint
    report were unreadable for it."""
    from pca2d.figures.sequence import pages_for

    plan = pages_for(np.array(["PROXIMA"] * 3 + ["GJ1"] * 2 + ["GJ3090"]))
    assert [label for _only, label in plan] == \
        ["PROXIMA", "GJ1", "GJ3090", "PROXIMA + GJ1 + GJ3090   (all)"]
    assert plan[0][0].tolist() == [True, True, True, False, False, False]
    assert plan[1][0].sum() == 2
    assert plan[-1][0] is None, "the last page keeps every row"

    solo = pages_for(np.array(["TOI2120"] * 4))
    assert solo == [(None, "TOI2120")], "one object, one page, named"
    assert pages_for(None) == [(None, None)], "a cube without the column"
    assert pages_for(np.array(["", ""])) == [(None, None)]
