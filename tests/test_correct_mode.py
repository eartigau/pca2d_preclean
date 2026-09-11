"""A fit of nights corrects every exposure, each with its own coefficients.

Proxima's 782 NIRPS exposures are fitted as 259 nights. The correct stage used
to pass --all, which corrects the files named in the fit's table, one per
night, with that night's coefficients: 254 corrected files for 782 exposures,
none of them with its own amplitudes. A cube of nights now switches it to
--by-night --refit.
"""

from astropy.table import Table

from pca2d.cli import correct_mode, nightly_stacked


def _cube(tmp_path, counts):
    path = tmp_path / "cube"
    path.mkdir(exist_ok=True)
    Table({"filename": ["f%d" % i for i in range(len(counts))],
           "n_exposures": counts}).write(path / "meta.fits", overwrite=True)
    return str(path)


def test_a_cube_of_nights_is_recognised(tmp_path):
    assert nightly_stacked(_cube(tmp_path, [3, 1, 2]))
    assert not nightly_stacked(_cube(tmp_path, [1, 1, 1]))
    assert not nightly_stacked(str(tmp_path / "absent"))


def test_every_exposure_of_a_night_is_refitted(tmp_path):
    plan = {"cube": _cube(tmp_path, [3, 2]), "written_config": "run/resolved_config.yaml"}
    assert correct_mode(plan) == ["--by-night", "--refit", "--config",
                                  "run/resolved_config.yaml"]
    plan["cube"] = _cube(tmp_path, [1, 1])
    assert correct_mode(plan) == ["--all"], "single exposures keep the fit's own rows"
