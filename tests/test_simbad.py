"""The star as SIMBAD and APERO describe it, without touching the network.

What is checked: which names are tried and in what order, where the header's
position comes from and at what epoch, that a star moving 10 arcsec a year is
compared at the right date, that SIMBAD being unreachable is an answer and
not an exception, and that APERO's PP_RV is read in m/s whatever its comment.
"""

import json
import os

import pytest

from pca2d import simbad
from pca2d import texreport as tr


def test_the_names_tried_start_with_aperos_own():
    header = {"PP_OBJNS": "WT 351", "PP_OBJN": "TOI756", "OBJECT": "TOI-756",
              "DRSOBJN": "TOI756"}
    assert simbad.candidates("TOI756", header) == ["WT 351", "TOI756",
                                                   "TOI-756"]
    # a Gaia identifier is the least ambiguous name of all
    spirou = {"GAIAID": 513299860904522752, "GAIADR": "Gaia DR2",
              "OBJECT": "TOI2120"}
    assert simbad.candidates("TOI2120", spirou)[0] == \
        "Gaia DR2 513299860904522752"
    # an instrument in the folder's name is not part of the star's
    assert "GL699" in simbad.candidates("GL699_SPIROU", {})


def test_positions_in_both_observatories_forms():
    assert simbad._sexagesimal("124825.216", True) == pytest.approx(192.10507, abs=1e-5)
    assert simbad._sexagesimal("-452814.148", False) == pytest.approx(-45.47060, abs=1e-5)
    assert simbad._sexagesimal("17:57:47.49", True) == pytest.approx(269.44788, abs=1e-5)
    assert simbad._sexagesimal("-0:30:00", False) == pytest.approx(-0.5)


def test_the_headers_position_and_its_epoch():
    apero = {"PP_RA": 281.71909, "PP_DEC": -62.177182, "PP_EPOCH": 2457206.0,
             "ESO TEL TARG ALPHA": 184652.552, "ESO TEL TARG DELTA": -621036.612}
    ra, dec, epoch, whose = simbad.header_position(apero)
    assert (ra, whose) == (281.71909, "APERO's")
    assert epoch == pytest.approx(2015.5, abs=0.01)
    eso = {"ESO TEL TARG ALPHA": 184652.552, "ESO TEL TARG DELTA": -621036.612,
           "ESO TEL TARG EPOCH": 2000.0}
    assert simbad.header_position(eso)[2:] == (2000.0, "the telescope's")
    cfht = {"OBJRA": "17:57:47.49", "OBJDEC": "4:44:49.6", "MJDATE": 58383.3}
    assert simbad.header_position(cfht)[2] == pytest.approx(2018.72, abs=0.01)
    assert simbad.header_position({}) is None


def test_a_fast_star_is_compared_where_it_was_that_night():
    """Barnard's star sat 194 arcsec from its J2000 position in a 2018 SPIRou
    header; carried along its proper motion, SIMBAD's position lands on it."""
    barnard = {"ra": 269.4520769586187, "dec": 4.693364966576667,
               "pmra": -801.551, "pmdec": 10362.394}
    cfht = {"OBJRA": "17:57:47.49", "OBJDEC": "4:44:49.6", "MJDATE": 58383.3}
    arcsec, epoch, whose = simbad.offset(barnard, cfht)
    assert arcsec < 2.0 and whose == "the telescope's"
    still = dict(barnard, pmra=0.0, pmdec=0.0)
    assert simbad.offset(still, cfht)[0] > 150


def test_no_network_is_an_answer(tmp_path, monkeypatch):
    def down(*_args, **_kw):
        raise OSError("no route to host")

    monkeypatch.setattr(simbad, "_get", down)
    found = simbad.star("TOI756", {"OBJECT": "TOI-756"}, cache=str(tmp_path))
    assert found["ok"] is False and "no route to host" in found["error"]
    assert not os.listdir(str(tmp_path)), "a failure is not cached"


def test_an_answer_is_kept_and_the_offset_worked_out_again(tmp_path,
                                                          monkeypatch):
    kept = {"ok": True, "main_id": "WT 351", "ra": 192.10506724,
            "dec": -45.47059666, "pmra": -216.502, "pmdec": 29.197,
            "offset": 999.0}
    with open(tmp_path / "TOI756.json", "w") as handle:
        json.dump(kept, handle)
    monkeypatch.setattr(simbad, "_get", lambda *a, **k: pytest.fail(
        "a fresh answer on disk is used as it is"))
    header = {"ESO TEL TARG ALPHA": 124825.216,
              "ESO TEL TARG DELTA": -452814.148, "ESO TEL TARG EPOCH": 2000.0}
    found = simbad.star("TOI756", header, cache=str(tmp_path))
    assert found["main_id"] == "WT 351"
    assert found["offset"] < 1.0, "measured against this header, not kept"


def test_aperos_velocity_is_in_metres_per_second_whatever_it_says(tmp_path):
    """NIRPS files label PP_RV [km/s] and hold GJ 1 at 25534: SIMBAD's
    25.534 km/s."""
    from astropy.io import fits

    folder = tmp_path / "data" / "GJ1"
    folder.mkdir(parents=True)
    science = fits.ImageHDU(name="FluxA")
    science.header["PP_OBJN"] = "GJ1"
    science.header["PP_RV"] = (25534.0, "The RV [km/s] used by the DRS")
    science.header["PP_RVS"] = "2018A&A...616A...7S"
    fits.HDUList([fits.PrimaryHDU(), science]).writeto(folder / "a_t.fits")
    header = tr.first_header({"input": {"directory": str(tmp_path / "data"),
                                        "object": "GJ1"}}, "GJ1")
    assert header["PP_RV_KMS"] == pytest.approx(25.534)
    table = tr.star_table("GJ1", {"ok": False, "error": "offline",
                                  "header": header}, [])
    assert "25.534" in table and r"\textbf{GJ1}" in table
    assert "offline" in table, "a star SIMBAD did not describe says why"
