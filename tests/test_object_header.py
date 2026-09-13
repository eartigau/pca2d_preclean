"""A folder may be named differently from what the headers say.

Barnard's star is `Gl699` in both pipelines' headers, and its two campaigns live
here as GL699_SPIROU and GL699_NIRPS, because a run is one instrument and the
two must not share a folder or an LBL object name. Matched on the folder's name,
every one of its spectra was skipped and the cube build ended with "every
exposure was rejected".
"""
from pca2d.config import cache_key, load_config
from pca2d.cube import object_key, same_object


def test_a_name_is_matched_on_what_identifies_it():
    assert object_key("Gl699") == object_key("GL699") == "GL699"
    assert object_key("TOI-2120") == object_key("TOI 2120") == "TOI2120"
    assert object_key("GL699_SPIROU") == "GL699SPIROU", \
        "which is NOT the same as GL699, hence the need for object_header"


def test_the_header_name_is_what_the_files_are_matched_on():
    meta = {"object": "Gl699", "drsobjn": "GL699"}
    assert same_object(meta, "Gl699")
    assert same_object(meta, "GL699"), "case and punctuation are not identity"
    assert not same_object(meta, "GL699_SPIROU"), \
        "the folder's name, which is why it cannot be the only thing tried"


def test_both_barnard_folders_name_the_header_they_match():
    for folder in ("GL699_SPIROU", "GL699_NIRPS"):
        config = load_config("config.yaml", object_name=folder,
                             instrument="SPIROU" if "SPIROU" in folder else "NIRPS")
        assert config["input"]["object"] == folder, "the folder is the object"
        assert config["input"]["object_header"] == "Gl699"
        assert same_object({"object": "Gl699"},
                           config["input"]["object_header"])


def test_an_unset_header_name_does_not_re_key_every_cube():
    """Adding the knob must not orphan cubes that took twenty minutes to build."""
    config = load_config("config.yaml", object_name="PROXIMA")
    assert config["input"]["object_header"] is None
    assert cache_key(config) == "13269e89fc73", \
        "the key PROXIMA's cube was built under, before the knob existed"
    config["input"]["object_header"] = "Proxima"
    assert cache_key(config) != "13269e89fc73", \
        "but a value describes which files are read, so it is hashed"


def test_the_header_name_is_part_of_the_object_and_not_of_the_build():
    """A joint run demands that its objects be built the same way, comparing
    their keys with the object taken out. Both NAMES of the object have to come
    out: left in, GL699_NIRPS had a different shared key and the four-object run
    would have refused the whole set."""
    from pca2d.cli import without_object

    keys = {}
    for name in ("PROXIMA", "GJ1", "GJ3090", "GL699_NIRPS"):
        config = load_config("config.yaml", object_name=name, instrument="NIRPS")
        keys[name] = cache_key(without_object(config))
    assert len(set(keys.values())) == 1, keys
    assert keys["GL699_NIRPS"] == keys["PROXIMA"], \
        "a folder named otherwise than its headers is built like any other"
