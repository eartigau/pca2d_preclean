"""A variant is the nominal plus what it changes, and a run says which code made it."""
import os

import pytest
import yaml

from pca2d.cli import load_variant, name_variant, reused_fit
from pca2d.config import VARIANT_META, load_config
from pca2d.provenance import code_version, stamp

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _config(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text("general:\n"
                    "  highpass:\n    width_kms: 100\n"
                    "  twoframe:\n    iters: 16\n    star_basis: spline\n"
                    "  lbl:\n    suffix: _PCA2D_{tag}\n"
                    "  output:\n    directory: outputs\n")
    (tmp_path / "variants").mkdir()
    return path


def test_a_variant_wins_over_the_config_and_its_notes_stay_out(tmp_path):
    path = _config(tmp_path)
    (tmp_path / "variants" / "v1.yaml").write_text(
        "reuse_fit: base\nnote: a result\ntwoframe:\n  iters: 6\n  keep: last\n")
    variant = load_variant(str(path), "v1")
    cfg = load_config(str(path), variant=variant)
    assert cfg["twoframe"]["iters"] == 6 and cfg["twoframe"]["keep"] == "last"
    assert cfg["twoframe"]["star_basis"] == "spline", \
        "what the variant does not change is the config's"
    assert not any(k in cfg for k in VARIANT_META)


def test_a_variant_in_samples_keeps_its_window(tmp_path):
    """The series of 2026-09-11 ran with a 151-sample high pass; a variant that
    says so must get that cube, not the config's 100 km/s one."""
    path = _config(tmp_path)
    cfg = load_config(str(path), variant={"highpass": {"window": 151}})
    assert cfg["highpass"]["window"] == 151
    assert cfg["highpass"]["width_kms"] is None
    assert load_config(str(path))["highpass"]["window"] == 201


def test_a_variant_is_named_beside_the_nominal(tmp_path):
    path = _config(tmp_path)
    variant = {"twoframe": {"iters": 6}}
    cfg = load_config(str(path), variant=variant)
    assert name_variant(cfg, "v1", variant) == "outputs"
    assert cfg["output"]["directory"] == os.path.join("outputs", "_v1")
    assert cfg["lbl"]["suffix"] == "_PCA2D_{tag}_v1"
    assert cfg["variant"] == {"name": "v1"}


def test_a_variant_that_names_its_folder_and_object_keeps_them(tmp_path):
    path = _config(tmp_path)
    variant = {"output": {"directory": "outputs/_other"},
               "lbl": {"suffix": "_PCA2D_{tag}_x"}}
    cfg = load_config(str(path), variant=variant)
    name_variant(cfg, "v1", variant)
    assert cfg["output"]["directory"] == "outputs/_other"
    assert cfg["lbl"]["suffix"] == "_PCA2D_{tag}_x"


def test_no_variant_changes_nothing(tmp_path):
    cfg = load_config(str(_config(tmp_path)))
    before = yaml.safe_dump(cfg)
    name_variant(cfg, None, None)
    assert yaml.safe_dump(cfg) == before


def test_a_reused_fit_is_found_under_the_nominal_root():
    assert reused_fit("outputs", "star_spl", "TOI2120", "1-3") == \
        os.path.join("outputs", "_star_spl", "TOI2120", "1-3")
    assert reused_fit("outputs", "nominal", "TOI2120", "1-3") == \
        os.path.join("outputs", "TOI2120", "1-3")


def test_an_unknown_variant_stops_the_run(tmp_path):
    with pytest.raises(SystemExit):
        load_variant(str(_config(tmp_path)), "nope")


def test_the_code_version_is_a_commit_and_says_if_it_was_modified():
    v = code_version()
    s = stamp(v)
    if v["commit"]:
        assert len(v["commit"]) == 40
        assert s[:12] == v["commit"][:12]
        assert s.endswith("+") == bool(v["dirty"])
    else:
        assert s == "unknown"


def test_every_variant_in_the_repository_gives_what_it_sets():
    folder = os.path.join(REPO, "variants")
    names = sorted(f[:-5] for f in os.listdir(folder) if f.endswith(".yaml"))
    assert names
    config = os.path.join(REPO, "config.yaml")
    for name in names:
        variant = load_variant(config, name)
        cfg = load_config(config, variant=variant)
        for section, keys in variant.items():
            if section in VARIANT_META:
                continue
            for key, value in keys.items():
                assert cfg[section][key] == value, (name, section, key)


def test_a_variant_may_not_write_over_the_fit_it_reuses(tmp_path):
    """A correction-only variant borrows another run's fit. Writing its own
    corrected spectra into that run's folder would leave that run's velocities
    describing files that are no longer there: on 2026-09-25 `--out-dir
    outputs`, which overrides the variant's own `_<name>` folder, did exactly
    that to 19 of TOI2120 0-7's spectra."""
    from pca2d.cli import refuse_to_overwrite_the_reused_fit as refuse

    base = str(tmp_path / "outputs" / "TOI2120" / "0-7")
    # its own folder: nothing to say
    refuse(str(tmp_path / "outputs" / "_v1" / "TOI2120" / "0-7"), base, "v1",
           str(tmp_path / "outputs"))
    with pytest.raises(SystemExit) as stopped:
        refuse(base + "/", base, "0-7chi2", str(tmp_path / "outputs"))
    said = str(stopped.value)
    assert "0-7chi2" in said and base in said
    assert "--out-dir" in said and "_0-7chi2" in said, \
        "it has to say how to get the variant its own folder"
