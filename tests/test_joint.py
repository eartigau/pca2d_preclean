"""Several stars, one observer basis: the joint cube and the labels it carries."""
import os

import numpy as np
import pytest
from astropy.table import Table

from pca2d import joint
from pca2d.reconstruct import group_of

N_PIX = 64


def _cube(path, rows, berv0=0.0, night0=0):
    os.makedirs(path, exist_ok=True)
    grid = 1500.0 * np.exp(np.arange(N_PIX) * 0.5 / 299792.458)
    np.save(os.path.join(path, "grid.npy"), grid)
    np.save(os.path.join(path, "data.npy"),
            np.arange(rows * N_PIX, dtype=np.float32).reshape(rows, N_PIX))
    np.save(os.path.join(path, "sigma.npy"), np.ones((rows, N_PIX), dtype=np.float32))
    meta = Table()
    meta["filename"] = ["%04dt.fits" % (night0 + i // 2) for i in range(rows)]
    meta["bjd"] = 2459000.0 + np.arange(rows) // 2
    meta["berv"] = berv0 + np.arange(rows) * 0.1
    meta["parity"] = np.arange(rows) % 2
    meta["exposure"] = np.arange(rows) // 2
    meta["object"] = ["?"] * rows
    meta.write(os.path.join(path, "meta.fits"), overwrite=True)
    return path


def test_the_rows_of_several_objects_sit_on_one_grid(tmp_path):
    a = _cube(str(tmp_path / "a"), 4)
    b = _cube(str(tmp_path / "b"), 6, berv0=10.0, night0=100)
    out = str(tmp_path / "joint")
    total = joint.build([a, b], ["ALPHA", "BETA"], out)
    assert total == 10
    data = np.load(os.path.join(out, "data.npy"))
    assert data.shape == (10, N_PIX)
    assert np.array_equal(data[:4], np.load(os.path.join(a, "data.npy")))
    assert np.array_equal(data[4:], np.load(os.path.join(b, "data.npy")))
    meta = Table.read(os.path.join(out, "meta.fits"))
    assert list(meta["object"]) == ["ALPHA"] * 4 + ["BETA"] * 6


def test_each_object_gets_its_own_star_spectra_and_keeps_the_parity(tmp_path):
    """The label is 2 * the object's index + the parity: the fit estimates one
    star spectrum per label, and `% 2` still gives the order parity, which
    means the same thing for every object."""
    a = _cube(str(tmp_path / "a"), 4)
    b = _cube(str(tmp_path / "b"), 4, berv0=10.0, night0=100)
    out = str(tmp_path / "joint")
    joint.build([a, b], ["ALPHA", "BETA"], out)
    meta = Table.read(os.path.join(out, "meta.fits"))
    parity = np.asarray(meta["parity"], dtype=int)
    assert list(parity) == [0, 1, 0, 1, 2, 3, 2, 3]
    assert list(parity % 2) == [0, 1] * 4
    exposure = np.asarray(meta["exposure"], dtype=int)
    assert exposure[0] == exposure[1], "the two rows of an exposure stay tied"
    assert not set(exposure[:4]) & set(exposure[4:]), "no id is shared"


def test_a_different_grid_is_refused(tmp_path):
    a = _cube(str(tmp_path / "a"), 2)
    b = _cube(str(tmp_path / "b"), 2)
    grid = np.load(os.path.join(b, "grid.npy"))
    np.save(os.path.join(b, "grid.npy"), grid * 1.001)
    with pytest.raises(SystemExit):
        joint.build([a, b], ["ALPHA", "BETA"], str(tmp_path / "joint"))


def test_an_object_is_told_which_star_spectra_are_its_own():
    objects = ["ALPHA", "BETA", "GAMMA"]
    assert [joint.star_group(objects, o) for o in objects] == [0, 2, 4]
    model = {"star_group": 4}
    assert group_of(model, 0, 6) == 4       # even orders of the third object
    assert group_of(model, 1, 6) == 5       # its odd ones
    assert group_of({}, 1, 2) == 1, "one object still picks its own parity"


def test_every_star_spectrum_gets_a_column_of_its_own(tmp_path):
    """Four groups, two objects: naming them all "all", as this did until
    2026-09-12, wrote one mean and dropped the rest, or refused the file for
    using a name twice. The correction reads them back by group."""
    from astropy.table import Table

    from pca2d.reconstruct import load_model
    from pca2d.twoframe import group_names, write_components_fits

    assert group_names(1) == ["all"]
    assert group_names(2) == ["even", "odd"]
    assert group_names(4) == ["g0_even", "g0_odd", "g1_even", "g1_odd"]
    assert len(set(group_names(6))) == 6

    m, groups = 16, 4
    grid = 1500.0 * np.exp(np.arange(m) * 0.5 / 299792.458)
    means = np.zeros((groups, m))
    templates = np.array([np.full(m, float(g + 1)) for g in range(groups)])
    table = Table({"filename": ["0001t.fits"], "berv": [0.0], "bjd": [2459000.0],
                   "rejected": [False], "chi2_red": [1.0], "b1": [1.0]})
    path = str(tmp_path / "twoframe_components.fits")
    write_components_fits(path, grid, np.zeros((0, m)), np.ones((1, m)),
                          np.zeros(m), means, np.zeros(0), np.ones(1),
                          1.0, 0.5, table, 0.5, templates=templates,
                          mean_mode="star")
    model = load_model(path)
    assert list(model["templates"]) == group_names(groups)
    for g in range(groups):
        model["star_group"] = g - g % 2
        from pca2d.reconstruct import template_for
        assert np.allclose(template_for(model, g % 2), templates[g]), g


def test_the_joint_cube_is_named_after_its_objects(tmp_path):
    assert joint.joint_name(["A", "B"]) == "A+B"
    one = joint.cube_path("cache", "tfits", "abc123", ["A", "B"])
    same = joint.cube_path("cache", "tfits", "abc123", ["A", "B"])
    other = joint.cube_path("cache", "tfits", "abc123", ["B", "A"])
    assert one == same and one != other, "the order of the objects is part of it"
    assert one.startswith(os.path.join("cache", "cube_tfits_abc123_j"))
