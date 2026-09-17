"""The known planets, from the NASA Exoplanet Archive: asked by the right
identifiers, a TOI and its confirmed planet counted once, and a report that
still says something when the archive does not answer. No network here: the
archive's answers are canned."""

import numpy as np
import pytest

from pca2d import archive

FACTS = {"ok": True, "main_id": "TOI-782",
         "ids": ["TIC 429358906", "Gaia DR3 3518374197418907648",
                 "2MASS J12154108-1854365", "GJ 1154"]}

CONFIRMED = [{"pl_name": "TOI-782 b", "pl_letter": "b", "hostname": "TOI-782",
              "pl_orbper": "8.0240015", "pl_orbpererr1": "7.7e-06",
              "pl_tranmid": "2458577.04189", "pl_rvamp": "11.6",
              "pl_rvamperr1": "", "pl_orbeccen": "0.19"}]
TOIS = [{"toi": "782.01", "tid": "429358906", "pl_orbper": "8.0239811",
         "pl_orbpererr1": "", "pl_tranmid": "2458577.0421740",
         "tfopwg_disp": "CP"},
        {"toi": "782.02", "tid": "429358906", "pl_orbper": "2.5",
         "pl_orbpererr1": "", "pl_tranmid": "2458577.5",
         "tfopwg_disp": "PC"},
        {"toi": "782.03", "tid": "429358906", "pl_orbper": "1.1",
         "pl_orbpererr1": "", "pl_tranmid": "", "tfopwg_disp": "FP"}]


@pytest.fixture
def canned(monkeypatch, tmp_path):
    asked = []

    def answer(query):
        asked.append(query)
        return CONFIRMED if "pscomppars" in query else TOIS

    monkeypatch.setattr(archive, "_get", answer)
    monkeypatch.setattr(archive, "CACHE", str(tmp_path / "cache"))
    return asked


def test_the_star_is_asked_for_by_its_identifiers():
    ids = archive.identifiers(FACTS, folder="TOI782")
    assert ids["tic"] == "429358906" and ids["toi"] == "782"
    assert ids["gaia"] == "Gaia DR3 3518374197418907648"
    assert "Gl 1154" in ids["names"] and "GJ 1154" in ids["names"], \
        "the archive spells some GJ hosts Gl"
    assert archive.identifiers({}, folder="TOIM4508")["toi"] is None
    assert archive.identifiers({}, folder="TOI_756")["toi"] == "756"


def test_a_toi_and_its_confirmed_planet_are_one(canned):
    got = archive.planets(FACTS, folder="TOI782")
    assert got["ok"]
    names = [p["name"] for p in got["planets"]]
    assert names == ["TOI-782.02", "TOI-782 b"], \
        "shortest first, the confirmed name kept, the false positive out"
    b = got["planets"][1]
    assert b["k"] == pytest.approx(11.6) and b["source"] == "confirmed"
    assert b["also"] == ["TOI-782.01"]
    assert any("tic_id = 'TIC 429358906'" in q for q in canned)
    assert any("tid = 429358906" in q for q in canned)
    # asked once: the second call is the cache's
    archive.planets(FACTS, folder="TOI782")
    assert len(canned) == 2


def test_an_archive_that_does_not_answer_leaves_the_cache_or_nothing(
        canned, monkeypatch):
    archive.planets(FACTS, folder="TOI782")

    def down(_query):
        raise OSError("no network")

    monkeypatch.setattr(archive, "_get", down)
    kept = archive.planets(FACTS, folder="TOI782", refresh=True)
    assert kept["stale"] and len(kept["planets"]) == 2
    nothing = archive.planets({"ids": ["TIC 1"]}, folder="X")
    assert nothing["ok"] is False and nothing["planets"] == []
