"""The periods a component must not vary at, and that they reach the figure.

The coefficient periodogram exists to answer one question: does any component of
the basis vary at the period of a known planet? A component that does subtracts
that planet out of the spectra, the velocities come back cleaner BECAUSE the
signal is gone, and nothing in the residuals looks wrong. Until 2026-09-12 the
figure was never given a single period, so the question was never put to it.
"""
import types

from pca2d.config import cache_key, load_config


def test_each_target_carries_its_published_periods(tmp_path):
    """Read against a data root that holds nothing: what a target's block says
    is a property of the configuration, and asking for it must not depend on a
    disk being mounted."""
    empty = str(tmp_path)
    for name, want in (("GJ3090", [2.85310198, 15.9407]),
                       ("TOI2120", [5.7998164]),
                       ("TOI4552", [0.30110032]),
                       ("TOI1452", [11.06201]),
                       ("PROXIMA", [11.18465, 5.12338])):
        config = load_config("config.yaml", object_name=name, data_dir=empty)
        assert config["target"]["planets"] == want, name
    quiet = load_config("config.yaml", object_name="GJ1", data_dir=empty)
    assert quiet["target"]["planets"] == [], "GJ 1 is the reference with none"


def test_knowing_a_period_does_not_rebuild_a_cube():
    """It describes the sky's owner, not the data: the cube is the same cube."""
    config = load_config("config.yaml", object_name="GJ3090")
    before = cache_key(config)
    config["target"]["planets"] = [1.0, 2.0, 3.0]
    config["target"]["prot"] = 90.0
    assert cache_key(config) == before


def test_a_joint_run_marks_every_members_planets(data_root):
    """One shared basis must not vary at ANY of the stars' periods, and the joint
    config is a copy of the first member's, which would have marked only its."""
    from pca2d.cli import joint_plan
    args = types.SimpleNamespace(
        objects=["PROXIMA", "GJ1", "GJ3090"], object="PROXIMA",
        config="config.yaml", data_dir=data_root("PROXIMA", "GJ1", "GJ3090"),
        out_dir=None, instrument="NIRPS",
        n_star=0, n_earth=None, windows=None, rebuild_cube=False, run_lbl=False,
        variant=None)
    plan = joint_plan(args, None)
    assert plan["config"]["target"]["planets"] == \
        sorted([2.85310198, 15.9407, 5.12338, 11.18465])


def test_the_bundle_hands_them_to_the_figure():
    """What is configured has to reach the command line, or the figure is a
    picture again."""
    from pca2d.figures.bundle import periodogram_args

    argv = periodogram_args("py", "p.py", "fit.npz", "out.pdf",
                            load_config("config.yaml", object_name="GJ3090"))
    assert argv[:6] == ["py", "p.py", "--fit", "fit.npz", "--out", "out.pdf"]
    assert argv[argv.index("--planets") + 1:] == ["2.85310198", "15.9407"], \
        "the archive's own digits, not a rounded period"

    quiet = periodogram_args("py", "p.py", "fit.npz", "out.pdf",
                             load_config("config.yaml", object_name="GJ1"))
    assert "--planets" not in quiet, "no planet, nothing to mark"
    assert "--prot" not in quiet

    spotted = load_config("config.yaml", object_name="GJ1")
    spotted["target"]["prot"] = 91.3
    assert "--prot" in periodogram_args("py", "p.py", "f", "o", spotted)
