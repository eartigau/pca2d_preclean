"""The velocities a run measured belong in the run's own report.

Until 2026-09-15 a run ended with an rdb under lbl/lblrdb and a report that
never mentioned a velocity: the one number a correction is judged by was not in
the document the correction produced. These pin what the pages promise: every
star of the run gets one, joint runs included; the delivered spectra are on it
whenever lbl.before asked for them; and the two series are compared on the
exposures they share and not on whatever each of them happens to have.
"""

import os

import numpy as np
import pytest
from astropy.table import Table

from pca2d.rvpages import (load_series, numbers, on_common, rdb_path,
                           stars_of, velocity_pages)


def write_rdb(data_dir, name, rjd, vrad, svrad):
    path = rdb_path(data_dir, name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    t = Table()
    t["rjd"], t["vrad"], t["svrad"] = rjd, vrad, svrad
    t["filename"] = ["%.4f.fits" % x for x in rjd]
    t.write(path, format="ascii.rdb", overwrite=True)
    return path


def plan_for(tmp_path, objects, before=True, tag="0-7"):
    data_dir = str(tmp_path / "lbl")
    outdir = tmp_path / "out" / tag
    outdir.mkdir(parents=True)
    config = {"input": {"object": "+".join(objects)},
              "lbl": {"directory": data_dir, "before": before,
                      "suffix": "_PCA2D_{tag}"}}
    plan = {"config": config, "tag": tag, "outdir": str(outdir)}
    if len(objects) > 1:
        plan["members"] = [{"object": o} for o in objects]
    else:
        config["input"]["object"] = objects[0]
    return plan, data_dir


def campaign(seed, scatter):
    r = np.random.default_rng(seed)
    rjd = 59000.0 + np.repeat(np.arange(14.0), 3) + r.uniform(0, 0.2, 42)
    return rjd, r.normal(0, scatter, 42), np.full(42, 1.5)


def test_one_star_gets_a_page_with_both_series(tmp_path):
    plan, data_dir = plan_for(tmp_path, ["GJ1"])
    rjd, v, e = campaign(1, 8.0)
    write_rdb(data_dir, "GJ1", rjd, v, e)
    write_rdb(data_dir, "GJ1_PCA2D_0-7", rjd, v * 0.4, e)

    out = velocity_pages(plan)
    assert out and os.path.exists(out)
    from pypdf import PdfReader
    r = PdfReader(out)
    assert len(r.pages) == 2, "the numbers, then the one star's curves"
    text = r.pages[0].extract_text()
    assert "delivered" in text and "corrected" in text
    assert "gain in rms" in text


def test_every_star_of_a_joint_run_gets_its_own_page(tmp_path):
    """A joint fit corrects several stars off one observer basis; LBL measures
    each on its own, so one page each, not one page for the lot."""
    plan, data_dir = plan_for(tmp_path, ["PROXIMA", "GJ1", "GJ3090"])
    for k, star in enumerate(("PROXIMA", "GJ1", "GJ3090")):
        rjd, v, e = campaign(k, 5.0 + k)
        write_rdb(data_dir, star, rjd, v, e)
        write_rdb(data_dir, "%s_PCA2D_0-7" % star, rjd, v * 0.5, e)

    assert [s[0] for s in stars_of(plan)] == ["PROXIMA", "GJ1", "GJ3090"]
    out = velocity_pages(plan)
    from pypdf import PdfReader
    r = PdfReader(out)
    assert len(r.pages) == 4, "one numbers page, three stars"
    text = r.pages[0].extract_text()
    for star in ("PROXIMA", "GJ1", "GJ3090"):
        assert star in text


def test_without_lbl_before_the_page_is_the_corrected_alone(tmp_path):
    plan, data_dir = plan_for(tmp_path, ["GJ1"], before=False)
    rjd, v, e = campaign(3, 6.0)
    write_rdb(data_dir, "GJ1_PCA2D_0-7", rjd, v, e)

    out = velocity_pages(plan)
    from pypdf import PdfReader
    r = PdfReader(out)
    assert len(r.pages) == 2
    text = r.pages[0].extract_text()
    assert "corrected" in text and "delivered" not in text


def test_a_star_lbl_has_not_measured_is_left_out_not_fatal(tmp_path):
    plan, data_dir = plan_for(tmp_path, ["GJ1", "GJ3090"])
    rjd, v, e = campaign(4, 7.0)
    write_rdb(data_dir, "GJ1", rjd, v, e)
    write_rdb(data_dir, "GJ1_PCA2D_0-7", rjd, v * 0.6, e)

    out = velocity_pages(plan)
    from pypdf import PdfReader
    assert len(PdfReader(out).pages) == 2, "GJ1 only"

    empty, _ = plan_for(tmp_path / "none", ["NOBODY"])
    (tmp_path / "none").mkdir(exist_ok=True)
    assert velocity_pages(empty) is None, "nothing measured, nothing written"


def test_the_series_are_compared_on_shared_exposures_only(tmp_path):
    """A corrected set that lost its worst nights would read as an improvement
    made of nothing."""
    plan, data_dir = plan_for(tmp_path, ["GJ1"])
    rjd, v, e = campaign(5, 9.0)
    v[:6] = 60.0                       # six terrible exposures
    write_rdb(data_dir, "GJ1", rjd, v, e)
    write_rdb(data_dir, "GJ1_PCA2D_0-7", rjd[6:], v[6:], e[6:])

    before = load_series(data_dir, "GJ1", "delivered")
    after = load_series(data_dir, "GJ1_PCA2D_0-7", "corrected")
    assert before["stats"]["n"] == 42 and after["stats"]["n"] == 36
    cut, n_common = on_common([before, after])
    assert n_common == 36
    assert cut[0]["stats"]["n"] == cut[1]["stats"]["n"] == 36
    assert cut[0]["stats"]["rms"] == pytest.approx(cut[1]["stats"]["rms"]), \
        "same exposures, same velocities: no gain from dropping rows"
    assert before["stats"]["n"] == 42, "the uncut series is not modified"
    assert "gain in rms" in numbers("GJ1", cut, n_common)


def test_the_pages_are_appended_to_the_report_that_is_there(tmp_path):
    """The report is bound at the figures stage, long before LBL runs. These
    pages go on the end of it rather than into a second document nobody opens.
    """
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt

    plan, data_dir = plan_for(tmp_path, ["GJ1"])
    report = os.path.join(plan["outdir"], "GJ1_0-7.pdf")
    with PdfPages(report) as pdf:
        for _ in range(3):
            fig = plt.figure()
            pdf.savefig(fig)
            plt.close(fig)
    rjd, v, e = campaign(6, 5.0)
    write_rdb(data_dir, "GJ1", rjd, v, e)
    write_rdb(data_dir, "GJ1_PCA2D_0-7", rjd, v * 0.3, e)

    out = velocity_pages(plan)
    assert out == report, "the run's own report, not a new file"
    from pypdf import PdfReader
    r = PdfReader(report)
    assert len(r.pages) == 5, "the three it had, plus numbers and one star"
    assert not os.path.exists(report + ".rv.pdf")
    titles = [it.title for it in r.outline if not isinstance(it, list)]
    assert "The velocities LBL measured" in titles
