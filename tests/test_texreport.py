"""The run report as a document, and the numbers it is made of.

On 2026-09-16 the report was called what it was: a compilation of matplotlib
outputs. It is a LaTeX document now, and these check what it says: which way
each number went, whose planets are whose, when velocities are older than the
fit, and that nothing about it can stop a run.
"""

import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest
import yaml
from astropy.table import Table

from pca2d import texreport as tr

HAVE_TEX = tr.find_pdflatex() is not None


@pytest.fixture(autouse=True)
def no_simbad(monkeypatch, tmp_path):
    """No network in a test: SIMBAD is 'not reached', and nothing is cached
    where a real report would find it."""
    from pca2d import simbad

    def offline(*_args, **_kw):
        raise OSError("no network in the tests")

    monkeypatch.setattr(simbad, "_get", offline)
    monkeypatch.setattr(simbad, "CACHE", str(tmp_path / "simbad_cache"))


# ----------------------------------------------------------------- text ---
def test_every_character_a_run_name_can_hold_prints():
    assert tr.tex("SMETHELLS_20") == r"SMETHELLS\_20"
    assert tr.tex("100% {x} & $y# ~") == \
        r"100\% \{x\} \& \$y\# \textasciitilde{}"
    assert tr.tex("(dF/dv)²") == r"(dF/dv)\textsuperscript{2}"
    # no em dash anywhere in what the report writes
    assert "—" not in tr.tex("a — b") and tr.tex("a — b") == "a - b"
    assert tr.tex("中") == "?", "unknown glyphs never stop the document"


def test_long_values_get_places_to_break():
    broken = tr.tex_break("_PCA2D_{tag}_GL406_3b7955")
    assert r"\_\allowbreak{}" in broken and r"\{tag\}" in broken


# ---------------------------------------------------------------- stats ---
def test_which_way_a_number_went():
    assert tr.verdict(10.0, 8.0) == "gain"
    assert tr.verdict(10.0, 12.0) == "loss"
    assert tr.verdict(10.0, 10.1) == "same"
    assert tr.verdict(0.40, 0.10, better="zero") == "gain"
    assert tr.verdict(-0.10, 0.40, better="zero") == "loss"
    assert tr.verdict(np.nan, 1.0) == "same"


def test_two_periods_are_one_peak_within_a_resolution_element():
    assert tr.same_period(5.38, 5.379, 900.0)
    assert not tr.same_period(5.38, 5.60, 900.0)
    assert not tr.same_period(None, 5.0, 900.0)


def test_outliers_do_not_scale_d2v():
    r = np.random.default_rng(3)
    x = r.normal(0, 1, 500)
    x[:5] = 1e6
    assert tr.inliers(x).sum() == 495
    assert tr.robust_sigma(x) == pytest.approx(1.0, rel=0.15)


def series(t, v, e, berv=None, d2v=None, label="x"):
    run = {"label": label, "t": t, "v": v, "e": e, "berv": berv,
           "d2v": d2v, "sd2v": None if d2v is None else np.full_like(d2v, 1e4),
           "mtime": 0.0}
    return run


def test_the_numbers_of_a_correction_that_removed_a_berv_bias():
    from pca2d import bervbias

    r = np.random.default_rng(5)
    t = np.sort(60000 + r.uniform(0, 600, 240))
    berv = 25 * np.sin(2 * np.pi * (t - 60000) / 365.25)
    planet = 5.0 * np.sin(2 * np.pi * t / 7.3)
    noise = r.normal(0, 2, t.size)
    e = np.full(t.size, 2.0)
    d2v = r.normal(0, 3e4, t.size)
    bias = bervbias.shape(berv, -8.0, bervbias.sigma_of(5.0))
    before = series(t, planet + noise + bias, e, berv, d2v, "delivered")
    after = series(t, planet + noise, e, berv, d2v, "corrected")
    numbers = tr.star_numbers(before, after, planets=[7.3])
    b, a = numbers["before"], numbers["after"]
    assert b["rms"] > a["rms"] and np.isfinite(numbers["removed"])
    assert b["bias_detected"] and not a["bias_detected"]
    assert b["bias_fwhm"] == pytest.approx(5.0, abs=1.5)
    assert "berv_r" not in b and "berv_slope" not in b, \
        "no straight line against BERV: it has no meaning here"
    assert a["planets"][0][0] == pytest.approx(5.0, abs=1.0), \
        "the planet is still there afterwards"
    better, worse, moved = tr.verdict_lines(numbers)
    assert "rms" in better and "the fitted BERV bias" in better
    assert not moved, "d2v and the planet were left alone"
    table = tr.summary_table("X", numbers)
    assert "BERV bias at its peak" in table and "slope" not in table
    assert r"$\Delta$BIC, no bias $-$ bias" in table
    assert b["bias_delta_bic"] > 10 > a["bias_delta_bic"]
    assert tr.bic_verdict(b["bias_delta_bic"], a["bias_delta_bic"]) == "gain"


# ------------------------------------------------------ whose planets ---
def test_a_joint_member_never_gets_the_first_members_planets(tmp_path):
    config = {"input": {"object": "PROXIMA+GJ1"},
              "target": {"planets": [11.18]},
              "provenance": {"command": "cli.py --objects PROXIMA,GJ1"}}
    assert tr.planets_of(str(tmp_path), config, "PROXIMA", joint=True) == [11.18]
    assert tr.planets_of(str(tmp_path), config, "GJ1", joint=True) == []
    # the member's own configuration beside the run is what counts
    (tmp_path / "cube_config_GJ1.yaml").write_text(
        yaml.safe_dump({"target": {"planets": [3.0]}}))
    assert tr.planets_of(str(tmp_path), config, "GJ1", joint=True) == [3.0]
    # and else its block in the configuration file the run was given
    layered = tmp_path / "config.yaml"
    layered.write_text(yaml.safe_dump(
        {"objects": {"TOI756": {"target": {"planets": [1.2392569, 149.4]}}}}))
    config["provenance"]["command"] = ("cli.py --objects PROXIMA,TOI756"
                                       " --config %s" % layered)
    assert tr.planets_of(str(tmp_path), config, "TOI756", joint=True) == \
        [1.2392569, 149.4]


# -------------------------------------------------- an old bound report ---
def test_an_old_bound_pdf_is_cut_along_its_bookmarks(tmp_path):
    from matplotlib.backends.backend_pdf import PdfPages
    from pypdf import PdfWriter

    pages = tmp_path / "pages.pdf"
    with PdfPages(pages) as pdf:
        for _ in range(6):
            fig = plt.figure()
            pdf.savefig(fig)
            plt.close(fig)
    writer = PdfWriter()
    writer.append(str(pages))
    for title, page in (("Front matter", 0), ("Variance", 2),
                        ("The sequence, step by step", 3),
                        ("What did not build", 5)):
        writer.add_outline_item(title, page)
    bound = tmp_path / "bound.pdf"
    with open(bound, "wb") as handle:
        writer.write(handle)
    kept = tr.split_bundle(str(bound), str(tmp_path / "report"))
    assert [(k["title"], k["pages"]) for k in kept] == \
        [("Variance", 1), ("The sequence, step by step", 2)]
    for k in kept:
        assert tr.pdf_pages(str(tmp_path / "report" / k["file"])) == k["pages"]


# ------------------------------------------------------- the document ---
def write_rdb(path, t, v, e, berv, d2v, dtemp=None):
    table = Table({"rjd": t, "vrad": v, "svrad": e, "d2v": d2v,
                   "sd2v": np.full_like(d2v, 1e4), "BERV": berv,
                   "EXTSN060": np.full_like(t, 80.0)})
    if dtemp is not None:
        table["DTEMP3500"] = dtemp
        table["sDTEMP3500"] = np.full_like(dtemp, 2.0)
    table.write(path, format="ascii.rdb", overwrite=True)


def a_run(tmp_path, star="TOI_756", with_rdb=True):
    """A run folder as the stages leave it, small enough to be a test."""
    root = tmp_path / "out"
    run = root / star / "0-3"
    (run / "report" / "figures").mkdir(parents=True)
    tree = tmp_path / "lbl"
    (tree / "lblrdb").mkdir(parents=True)
    config = {"input": {"object": star},
              "twoframe": {"n_star": 0, "n_earth": 3},
              "target": {"planets": [7.3]},
              "output": {"directory": str(root), "windows": ["1200:2"]},
              "lbl": {"directory": str(tree), "suffix": "_PCA2D_{tag}"},
              "provenance": {"commit": "0123456789abcdef", "dirty": False,
                             "command": "cli.py --object %s" % star}}
    (run / "resolved_config.yaml").write_text(yaml.safe_dump(config))
    fig = plt.figure(figsize=(4, 3))
    plt.plot([0, 1], [0, 1])
    fig.savefig(run / "report" / "figures" / "04-variance.pdf")
    plt.close(fig)
    tr.write_manifest(str(run / "report"), {
        "figures": [{"title": "Variance", "file": "figures/04-variance.pdf",
                     "pages": 1}],
        "failures": [{"what": "oh_residual.py", "why": "no OHLine",
                      "skipped": True}],
        "run": [["object", star], ["components", "0 star + 3 observer"]]})
    if with_rdb:
        r = np.random.default_rng(9)
        t = np.sort(60000 + r.uniform(0, 400, 120))
        berv = 20 * np.sin(2 * np.pi * (t - 60000) / 365.25)
        e = np.full(t.size, 2.0)
        noise = r.normal(0, 2, t.size)
        d2v = r.normal(0, 3e4, t.size)
        # a temperature that rotates at 7.3 d, and a telluric residual in the
        # delivered one that follows BERV
        spot = 8.0 * np.sin(2 * np.pi * t / 7.3) + r.normal(0, 2, t.size)
        write_rdb(tree / "lblrdb" / ("lbl_%s_%s.rdb" % (star, star)),
                  t, 1000 + noise + 0.5 * berv, e, berv, d2v,
                  spot + 0.4 * berv)
        name = star + "_PCA2D_0-3"
        write_rdb(tree / "lblrdb" / ("lbl_%s_%s.rdb" % (name, name)),
                  t, 1000 + noise, e, berv, d2v, spot)
    return run


@pytest.mark.skipif(not HAVE_TEX, reason="no pdflatex on this machine")
def test_a_run_with_velocities_becomes_a_document(tmp_path):
    run = a_run(tmp_path)
    out = tr.render(str(run))
    assert out == str(run / "TOI_756_0-3.pdf") and os.path.exists(out)
    assert tr.pdf_pages(out) >= 5
    tex = (run / "report" / "TOI_756_0-3.tex").read_text()
    for said in (r"\section{Summary}", r"\section{Velocities}",
                 r"\section{The run}", r"\section{The correction}",
                 r"\section{What did not build}", r"\tableofcontents",
                 r"TOI\_756", r"\gain{gain}", "K at 7.3000 d",
                 r"\listoffigures"):
        assert said in tex, said
    # every figure has a short name for the list, not its whole caption
    assert tex.count(r"\begin{figure}") == tex.count(r"\caption[")
    assert r"\caption[TOI\_756: the BERV bias, fitted]" in tex
    assert "<<" not in tex, "every placeholder of the template is filled"
    for kind in ("time", "berv", "corner", "d2v", "dtemp", "change",
                 "periods"):
        assert (run / "report" / "figures" / ("rv-toi-756-%s.pdf" % kind)).exists()
    assert "older than this run" not in tex
    # the temperature projection, both ways, and what the correction did to it
    assert r"DTEMP3500 robust sigma (K)" in tex
    assert "DTEMP3500's structure against $V_\\mathrm{tot}$" in tex
    # what was observed, from the headers, before the numbers
    assert r"TOI\_756: the observations, from the headers." in tex
    assert r"median SNR (EXTSN060) & 80.0" in tex
    assert r"calendar dates with $|V_\mathrm{tot}| < 4$ km/s, every year" in tex
    assert re.search(r"& \d+ [A-Z][a-z]{2} to \d+ [A-Z][a-z]{2}[^:&]*: \d+ of the"
                     r" 120 exposures \(\d+\.\d\\%\) fall on these dates", tex)
    assert r"median = " not in tex, "the median is on the figure's axis"
    # the date, and the time the PDF was written under it
    assert re.search(r"\\date\{\d{4}-\d{2}-\d{2}\\\\\n\{\\small written at"
                     r" \d{2}:\d{2}:\d{2} ", tex)
    assert "No temperature projection" not in tex


@pytest.mark.skipif(not HAVE_TEX, reason="no pdflatex on this machine")
def test_before_lbl_the_report_says_there_is_nothing_to_compare(tmp_path):
    run = a_run(tmp_path, with_rdb=False)
    assert tr.render(str(run))
    tex = (run / "report" / "TOI_756_0-3.tex").read_text()
    assert "no velocities to compare yet" in tex
    assert r"\section{Velocities}" not in tex


@pytest.mark.skipif(not HAVE_TEX, reason="no pdflatex on this machine")
def test_velocities_older_than_the_fit_are_said_to_be(tmp_path):
    run = a_run(tmp_path)
    np.savez(run / "fit.npz", a=np.zeros((2, 0)))
    future = os.path.getmtime(run / "fit.npz") + 3600
    os.utime(run / "fit.npz", (future, future))
    assert tr.render(str(run))
    tex = (run / "report" / "TOI_756_0-3.tex").read_text()
    assert "older than this run" in tex


def test_no_latex_leaves_the_bound_pdf_as_it_was(tmp_path, monkeypatch):
    run = a_run(tmp_path)
    (run / "TOI_756_0-3.pdf").write_bytes(b"%PDF bound")
    monkeypatch.setattr(tr, "find_pdflatex", lambda: None)
    assert tr.render(str(run)) is None
    assert (run / "TOI_756_0-3.pdf").read_bytes() == b"%PDF bound"


def test_a_report_that_fails_never_stops_the_run(monkeypatch):
    from pca2d import cli

    def broken(*_args, **_kw):
        raise RuntimeError("anything at all")

    monkeypatch.setattr(tr, "render", broken)
    said = []
    monkeypatch.setattr(cli, "log", lambda text, level="info": said.append(level))
    assert cli.render_report({"outdir": "/nowhere", "config": {}}) is None
    assert said == ["warn"]


def test_the_template_ships_with_the_package():
    import tomllib
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    with open(os.path.join(here, "pyproject.toml"), "rb") as handle:
        data = tomllib.load(handle)
    assert "templates/*.tex" in data["tool"]["setuptools"]["package-data"]["pca2d"]
    assert os.path.exists(tr.TEMPLATE)
    body = open(tr.TEMPLATE).read()
    for key in ("TITLE", "SUBTITLE", "DATE", "RUNHEAD", "ABSTRACT", "BODY"):
        assert "<<%s>>" % key in body
    assert "—" not in body and "---" not in body, "no em dash"


@pytest.mark.skipif(not HAVE_TEX, reason="no pdflatex on this machine")
def test_an_old_report_rewritten_is_kept_whole_first(tmp_path):
    """--run on a run from before the report replaces its bound PDF, the only
    copy of its figures: the old file is kept, so that can be undone."""
    from matplotlib.backends.backend_pdf import PdfPages
    from pypdf import PdfWriter

    run = a_run(tmp_path)
    os.remove(run / "report" / "manifest.json")
    pages = tmp_path / "pages.pdf"
    with PdfPages(pages) as pdf:
        for _ in range(3):
            fig = plt.figure()
            pdf.savefig(fig)
            plt.close(fig)
    writer = PdfWriter()
    writer.append(str(pages))
    writer.add_outline_item("Front matter", 0)
    writer.add_outline_item("Variance", 1)
    bound = run / "TOI_756_0-3.pdf"
    with open(bound, "wb") as handle:
        writer.write(handle)
    before = bound.read_bytes()
    assert tr.render(str(run))
    assert (run / "report" / "bound_before_the_report.pdf").read_bytes() == before
    assert bound.read_bytes() != before


def test_the_two_posteriors_are_drawn_on_the_same_axes(tmp_path):
    """Side by side, they are compared: one amp axis and one FWHM axis for
    both, whatever each posterior spans."""
    from pca2d import bervbias

    r = np.random.default_rng(8)
    berv = r.uniform(-25, 25, 150)
    e = np.full(150, 5.0)
    noise = r.normal(0, 5, 150)
    fits = {"before": bervbias.fit(berv, bervbias.shape(berv, -10, bervbias.sigma_of(5.0)) + noise,
                                   e, seed=1),
            "after": bervbias.fit(berv, noise, e, seed=2)}
    before = {"label": "delivered", "fits": fits}
    after = {"label": "corrected"}
    captured = []
    real = plt.Figure.savefig

    def keep(fig, *args, **kwargs):
        captured.append(fig)
        return real(fig, *args, **kwargs)

    import unittest.mock
    with unittest.mock.patch.object(plt.Figure, "savefig", keep):
        assert tr.figure_corner(before, after, str(tmp_path / "c.pdf"))
    joints = [ax for ax in captured[0].axes
              if ax.get_ylabel().startswith("amp")]
    assert len(joints) == 2
    assert joints[0].get_ylim() == joints[1].get_ylim()
    assert joints[0].get_xlim() == joints[1].get_xlim()


@pytest.mark.skipif(not HAVE_TEX, reason="no pdflatex on this machine")
def test_a_run_without_dtemp_says_so(tmp_path):
    """Every run before 2026-09-16: LBL was given no temperature table."""
    run = a_run(tmp_path)
    tree = tmp_path / "lbl" / "lblrdb"
    for path in tree.glob("*.rdb"):
        table = Table.read(path, format="ascii.rdb")
        table.remove_columns(["DTEMP3500", "sDTEMP3500"])
        table.write(path, format="ascii.rdb", overwrite=True)
    assert tr.render(str(run))
    tex = (run / "report" / "TOI_756_0-3.tex").read_text()
    assert "No temperature projection" in tex
    assert "DTEMP3500 robust sigma" not in tex
    assert not (run / "report" / "figures" / "rv-toi-756-dtemp.pdf").exists()


def test_the_dtemp_table_is_the_one_nearest_the_star():
    from pca2d.lbl import dtemp_table, runparams

    assert dtemp_table({}, 3304.0) == {"DTEMP3500": "temperature_gradient_3500.fits"}
    assert dtemp_table({"dtemp": "auto"}, 2700.0) == \
        {"DTEMP3000": "temperature_gradient_3000.fits"}
    assert dtemp_table({"dtemp": 5200}, 3000.0) == \
        {"DTEMP5000": "temperature_gradient_5000.fits"}
    assert dtemp_table({"dtemp": False}, 3300.0) == {}
    assert dtemp_table({}, None) == {}, "no Teff, nothing to choose by"
    with pytest.raises(SystemExit):
        dtemp_table({"dtemp": "hot"}, 3300.0)
    params = runparams({"lbl": {"steps": ["compute"]}}, "lbl", "NIRPS_HE",
                       "CADC", ["A", "A_PCA2D_0-3"], "c.yaml", teff=3717.0)
    assert params["RESPROJ_TABLES"] == \
        {"DTEMP3500": "temperature_gradient_3500.fits"}


def test_dtemp_stays_first_when_strpca_joins_it(tmp_path):
    """Only the first RESPROJ table is right in the LBL installed here, and
    DTEMP is the one both objects are compared on."""
    from pca2d.lbl import write_runner

    path = write_runner(str(tmp_path / "run_lbl.py"),
                        [("AFTER", "x", {"RESPROJ_TABLES": {"DTEMP3500": "t"}})],
                        {"fit": "f"}, "A", "2-3", str(tmp_path / "c.yaml"))
    body = open(path).read()
    assert 'dict(AFTER.get("RESPROJ_TABLES") or {}, **strpca_from(**STRPCA))' in body


def test_the_minus_signs_are_in_the_figures_text(tmp_path):
    """A viewer showed the posterior's amp axis as 20, 10, 0, 10, 20: the
    Unicode minus of a Type 3 font was drawn but not in the PDF's text. The
    report's figures embed a real font and write an ASCII minus, and the amp
    axis carries its sign."""
    from pypdf import PdfReader

    from pca2d import bervbias

    r = np.random.default_rng(8)
    berv = r.uniform(-25, 25, 150)
    e = np.full(150, 5.0)
    noise = r.normal(0, 5, 150)
    fits = {"before": bervbias.fit(berv, bervbias.shape(berv, -10, bervbias.sigma_of(5.0)) + noise,
                                   e, seed=1),
            "after": bervbias.fit(berv, noise, e, seed=2)}
    path = str(tmp_path / "corner.pdf")
    assert tr.figure_corner({"label": "delivered", "fits": fits},
                            {"label": "corrected"}, path)
    text = PdfReader(path).pages[0].extract_text()
    assert "-10" in text and "+10" in text
    fonts = PdfReader(path).pages[0]["/Resources"]["/Font"]
    assert all(f.get_object()["/Subtype"] != "/Type3" for f in fonts.values())


def test_the_bias_is_fitted_against_the_total_velocity():
    """V_tot = vrad/1000 - BERV: a star at +40 km/s has none of its
    exposures listed as close to the tellurics, since V_tot never comes near
    zero; one at +15 km/s whose velocities carry a bias in V_tot has it
    found there, where against BERV alone it would sit at +15."""
    from pca2d import bervbias

    r = np.random.default_rng(21)
    t = np.sort(60000 + r.uniform(0, 700, 200))
    berv = 25 * np.sin(2 * np.pi * (t - 60000) / 365.25)
    e = np.full(t.size, 2.0)
    v = 40000.0 + r.normal(0, 2, t.size)
    before = series(t, v, e, berv, None, "delivered")
    after = series(t, 40000.0 + r.normal(0, 2, t.size), e, berv, None,
                   "corrected")
    assert np.allclose(tr.vtot_of(before), v / 1000.0 - berv, atol=1e-6)
    numbers = tr.star_numbers(before, after)
    o = numbers["observations"]
    assert o["systemic"] == pytest.approx(40.0, abs=0.01)
    assert o["vtot"][0] > 4.0 and o["windows"] == [] and o["close_n"] == 0
    assert o["first"] <= o["last"] and o["nights"] == numbers["nights"]
    table = tr.observations_table("X", numbers)
    assert "never comes within 4 km/s of zero" in table
    assert "median SNR" not in table, "no SNR column, no SNR row"
    # V_tot never within 10 km/s of zero: nothing to fit, and said so
    assert "bias_peak" not in numbers["before"]
    assert numbers["before"]["bias_near"] < bervbias.MIN_NEAR

    vtot = 15.0 - berv
    v = 15000.0 + bervbias.shape(vtot, -8.0, bervbias.sigma_of(5.0)) \
        + r.normal(0, 2, t.size)
    before = series(t, v, e, berv, None, "delivered")
    after = series(t, 15000.0 + r.normal(0, 2, t.size), e, berv, None,
                   "corrected")
    numbers = tr.star_numbers(before, after)
    assert numbers["before"]["bias_detected"]
    assert not numbers["after"]["bias_detected"]
    assert numbers["before"]["bias_fwhm"] == pytest.approx(5.0, abs=1.5)


def test_the_close_dates_are_calendar_windows_every_year_alike():
    dates = np.array(["2023-06-26", "2023-06-28", "2024-06-17", "2025-07-02",
                      "2024-12-28", "2025-01-04", "2024-02-29"])
    windows = tr.calendar_windows(dates)
    assert [(tr.calendar_day(a), tr.calendar_day(b)) for a, b in windows] == \
        [("28 Dec", "4 Jan"), ("28 Feb", "28 Feb"), ("17 Jun", "2 Jul")], \
        "the year left out, a window across the new year kept whole"
    everything = np.array(["2026-06-20", "2022-01-02", "2024-03-15",
                           "2023-12-30", "2024-07-03"])
    assert list(tr.in_windows(everything, windows)) == \
        [True, True, False, True, False]
    assert tr.calendar_windows(np.array([], str)) == []


def test_the_periodograms_name_the_year_its_harmonics_and_the_month():
    periods = dict((name, p) for p, name in tr.REFERENCE_PERIODS)
    assert periods == pytest.approx({"1 yr": 365.25, "1/2 yr": 182.625,
                                     "1/3 yr": 121.75, "month": 29.53})


def test_the_periodogram_gives_the_power_of_each_false_alarm_level():
    r = np.random.default_rng(3)
    t = np.sort(60000 + r.uniform(0, 500, 150))
    y = r.normal(0, 1, t.size)
    found = tr.periodogram(t, y, np.ones_like(t))
    levels = found[5]
    assert len(levels) == len(tr.FAP_LEVELS) == 3
    assert np.all(np.diff(levels) > 0), "rarer false alarms need more power"
    assert 0 < levels[0] < levels[-1] < 1


def test_the_snr_quoted_is_the_measured_one_not_the_goal():
    from astropy.table import Table
    table = Table({"SNRGOAL": [150.0, 150.0], "EXTSN035": [80.0, 90.0],
                   "BERV": [1.0, 2.0]})
    name, values = tr.snr_column(table)
    assert name == "EXTSN035" and list(values) == [80.0, 90.0]
    assert tr.snr_column(Table({"SNRGOAL": [150.0]})) == (None, None)
    assert tr.snr_column(Table({"SNR": [42.0]}))[0] == "SNR"
