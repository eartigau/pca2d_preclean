"""The star block as LBL's template and as its variability basis.

pca2d hands LBL a template of the star taken from the fit (the first star
component at its mean amplitude) and, past the first component, vectors in the
place of LBL's DTEMP temperature gradients, named STRPCA2..N. These pin their
content and their form: what LBL's writer and its RESPROJ reader are handed.
"""

import os

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from pca2d import lbltemplate as lt
from pca2d.grids import doppler


def _fit(k=3, m=500, rows=20, seed=1):
    r = np.random.default_rng(seed)
    P = np.linalg.qr(r.normal(size=(m, k)))[0].T
    a = r.normal(size=(rows, k)) + np.arange(k, 0, -1) * 10.0
    rejected = np.zeros(rows, dtype=bool)
    rejected[0] = True
    return {"a": a, "P": P, "rejected": rejected, "template": np.zeros(m)}


def test_the_star_is_each_component_at_its_mean_amplitude():
    fit = _fit()
    ln, abar, P = lt.mean_star(fit)
    assert np.allclose(abar, fit["a"][1:].mean(axis=0)), "rejected rows are left out"
    assert ln.shape == (2, 500)
    assert np.allclose(ln[0], abar @ P) and np.allclose(ln[1], abar @ P)


def test_a_per_parity_star_frame_mean_goes_into_its_own_parity():
    fit = _fit()
    fit["templates"] = np.vstack([np.full(500, 0.1), np.full(500, -0.2)])
    ln, abar, P = lt.mean_star(fit)
    assert np.allclose(ln[0] - abar @ P, 0.1) and np.allclose(ln[1] - abar @ P, -0.2)


def test_the_template_columns_are_what_lbl_writes():
    m = 200
    ln = np.vstack([np.full(m, -0.1), np.full(m, -0.1)])
    count = np.vstack([np.full(m, 10.0), np.full(m, 10.0)])
    count[1, :50] = 3                                  # odd orders barely saw these
    wsum = count * 4.0                                 # 1/sigma^2 = 4 per row
    cols = lt.template_columns(ln, count, wsum, np.array([10, 10]))
    assert set(lt.COLUMNS) <= set(cols)
    assert np.allclose(cols["flux_even"], np.exp(-0.1))
    assert np.isnan(cols["flux_odd"][:50]).all() and np.isfinite(cols["flux_odd"][50:]).all()
    assert np.allclose(cols["rms_even"], np.exp(-0.1) * 0.5), "flux x one exposure's sigma"
    assert np.allclose(cols["eflux"], cols["rms"])
    assert np.isfinite(cols["flux"]).all(), "the even orders saw everything"


class _Inst:
    """LBL's instrument, as far as write_template uses it, recording its props."""

    def calculate_savgol_template(self, dv_grid, flux_dict):
        return {"%s_savgol_d%d" % (k, d): v for k, v in flux_dict.items() for d in range(4)}

    def load_header(self, filename):
        return {"BERV": 12.0 if "a" in filename else 12.4}

    def get_berv(self, header):
        return header["BERV"] * 1000.0

    def populate_sci_table(self, filename, tdict, header, berv=0.0):
        tdict.setdefault("FILENAME", []).append(filename)
        tdict.setdefault("BERV", []).append(berv)
        return tdict

    def load_science_file(self, filename):
        class Header(dict):
            def __setitem__(self, key, value, comment=None):
                dict.__setitem__(self, key, (value, comment))
        return None, Header()

    def write_template(self, path, props, header, sci_table):
        self.written = (path, props, header, sci_table)


def test_write_template_hands_lbls_writer_everything_it_reads():
    inst = _Inst()
    m = 50
    cols = {key: np.ones(m) for key in lt.COLUMNS}
    lt.write_template(inst, "t.fits", np.linspace(1500, 1501, m), cols,
                      ["a.fits", "b.fits"], 500.0, {lt.PROVENANCE: ("1-3v", "run")})
    path, props, header, sci = inst.written
    # the keys LBL's Instrument.write_template reads out of props
    wanted = {"wavelength", "template_type", "template_coverage", "total_nobs_berv",
              "template_nobs", "savgol_fluxes"} | set(lt.COLUMNS)
    assert wanted <= set(props)
    assert "flux_odd_savgol_d3" in props["savgol_fluxes"]
    assert props["template_type"] == "LBL_SAVGOL"
    assert np.all(np.asarray(sci["VSYS"]) == 0) and len(sci["VSYS"]) == 2
    assert header[lt.PROVENANCE][0] == "1-3v"
    assert props["template_coverage"] == 1, "12.0 and 12.4 km/s share one km/s bin"


def test_the_strpca_tables_are_lbl_resproj_tables(tmp_path):
    fit = _fit(k=3)
    grid = np.exp(np.linspace(np.log(1500), np.log(1510), 500))
    flux = np.exp(-0.05 * np.abs(np.sin(np.arange(500))))
    flux[:20] = np.nan
    tables = lt.write_strpca(str(tmp_path), "TOI_PCA2D_2-3v", grid, fit["P"], flux, 21256.75)
    assert list(tables) == ["STRPCA2", "STRPCA3"], "one per component past the first"
    t = Table.read(tmp_path / tables["STRPCA2"])
    assert {"wavelength", "fractional_gradient"} <= set(t.colnames)
    # the rest frame, as LBL's doppler_shift(grid, v) = grid sqrt((1 - b) / (1 + b))
    beta = 21256.75 / 299792458.0
    assert np.allclose(t["wavelength"], grid * np.sqrt((1 - beta) / (1 + beta)), rtol=1e-12)
    assert np.allclose(t["wavelength"], grid / doppler(21.25675), rtol=1e-12)
    g = np.asarray(t["fractional_gradient"])
    assert np.allclose(g[20:], flux[20:] * fit["P"][1][20:])
    assert np.all(g[:20] == 0), "no NaN for LBL to interpolate across"


def test_a_template_lbl_wrote_is_not_ours(tmp_path):
    ours, theirs = tmp_path / "ours.fits", tmp_path / "theirs.fits"
    fits.PrimaryHDU(header=fits.Header([(lt.PROVENANCE, "1-3v")])).writeto(ours)
    fits.PrimaryHDU().writeto(theirs)
    assert lt.is_ours(str(ours)) and not lt.is_ours(str(theirs))
    assert not lt.is_ours(str(tmp_path / "absent.fits"))


def test_the_systemic_velocity_is_read_from_the_mask(tmp_path):
    path = tmp_path / "mask.fits"
    fits.PrimaryHDU(header=fits.Header([("SYSTVELO", 21256.75)])).writeto(path)
    assert lt.mask_systemic(str(path)) == pytest.approx(21256.75)


def _template_file(path, grid, flux, stamp="A"):
    header = fits.Header([(lt.PROVENANCE, "2-3v"), ("PCA2FIT", stamp)])
    table = fits.BinTableHDU(Table({"wavelength": grid, "flux": flux}), name="TEMPLATE")
    fits.HDUList([fits.PrimaryHDU(header=header), table]).writeto(path)


def test_the_template_goes_where_lbl_looks_and_never_over_lbls_own(tmp_path):
    grid = np.linspace(1500, 1501, 10)
    made = tmp_path / "star_template.fits"
    _template_file(made, grid, np.ones(10), stamp="A")
    slot = tmp_path / "templates" / "LBL_Template_X_spirou.fits"
    assert lt.place(str(made), str(slot)) == "copied" and lt.stamp(str(slot)) == "A"
    assert lt.place(str(made), str(slot)) == "same"
    made.unlink()
    _template_file(made, grid, np.ones(10), stamp="B")
    assert lt.place(str(made), str(slot)) == "replaced" and lt.stamp(str(slot)) == "B"
    theirs = tmp_path / "templates" / "LBL_Template_Y_spirou.fits"
    fits.PrimaryHDU().writeto(theirs)
    before = theirs.read_bytes()
    assert lt.place(str(made), str(theirs)) == "theirs"
    assert theirs.read_bytes() == before, "LBL's own template is left exactly as it was"
    assert not list((tmp_path / "templates").glob("*.part"))


def test_strpca_from_the_files_a_run_leaves(tmp_path):
    fit = _fit(k=2)
    np.savez(tmp_path / "fit.npz", P=fit["P"])
    grid = np.exp(np.linspace(np.log(1500), np.log(1510), 500))
    flux = np.full(500, 0.9)
    flux[:5] = np.nan
    _template_file(tmp_path / "t.fits", grid, flux)
    mask = tmp_path / "mask.fits"
    fits.PrimaryHDU(header=fits.Header([("SYSTVELO", -1234.5)])).writeto(mask)
    tables = lt.strpca_from(str(tmp_path / "fit.npz"), str(tmp_path / "t.fits"),
                            str(mask), str(tmp_path), "X_PCA2D_2-3v", run="2-3v")
    assert tables == {"STRPCA2": "STRPCA2_X_PCA2D_2-3v.fits"}
    t = Table.read(tmp_path / tables["STRPCA2"])
    assert np.allclose(t["wavelength"], grid / doppler(-1.2345), rtol=1e-12)
    assert np.allclose(np.asarray(t["fractional_gradient"])[5:], 0.9 * fit["P"][1][5:])


def test_strpca_needs_the_mask_and_a_fit_from_the_same_run(tmp_path):
    np.savez(tmp_path / "fit.npz", P=_fit(k=2)["P"])
    _template_file(tmp_path / "t.fits", np.linspace(1500, 1501, 499), np.ones(499))
    mask = tmp_path / "mask.fits"
    args = (str(tmp_path / "fit.npz"), str(tmp_path / "t.fits"), str(mask), str(tmp_path), "X")
    with pytest.raises(SystemExit, match="mask step"):
        lt.strpca_from(*args)
    fits.PrimaryHDU(header=fits.Header([("SYSTVELO", 0.0)])).writeto(mask)
    with pytest.raises(SystemExit, match="same run"):
        lt.strpca_from(*args)


def test_lbls_in_place_division_is_recognised_and_a_copy_is_not():
    aliased = ("        frac_diff_seg = diff_seg\n"
               "        frac_mean_rms = mean_rms\n"
               "        frac_diff_seg /= (b_ratio_seg * norm_seg)\n")
    copied = ("        frac_diff_seg = diff_seg / (b_ratio_seg * norm_seg)\n"
              "        frac_mean_rms = mean_rms / (b_ratio_seg * norm_seg)\n")
    assert lt._aliased(aliased)
    assert not lt._aliased(copied)
    assert not lt._aliased("        frac_diff_seg = diff_seg.copy()\n"
                           "        frac_diff_seg /= norm\n")


def test_the_stamp_changes_when_the_fit_does(tmp_path):
    path = tmp_path / "fit.npz"
    np.savez(path, a=np.zeros(3))
    first = lt.fit_stamp(str(path))
    np.savez(path, a=np.zeros(30))
    assert lt.fit_stamp(str(path)) != first
