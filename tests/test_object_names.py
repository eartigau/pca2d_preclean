"""A file belongs to the object asked for whatever case its OBJECT was typed in.

OBJECT is what somebody typed into an observing form. The NIRPS files of
Proxima say 'Proxima', APERO calls them PROXIMA (DRSOBJN), and the folder and
the --object were PROXIMA: an exact comparison skipped all 89 of them, one
line each, with nothing to say what it had been expecting.
"""

from pca2d.cube import object_key, same_object


def test_case_and_punctuation_do_not_make_another_object():
    assert object_key("Proxima") == object_key("PROXIMA") == "PROXIMA"
    assert object_key("TOI-2120") == object_key("TOI 2120") == object_key("toi2120")
    assert object_key(None) == ""


def test_the_proxima_files_are_proxima():
    assert same_object({"object": "Proxima"}, "PROXIMA")
    assert same_object({"object": "TOI2120"}, "TOI-2120")


def test_apero_s_own_name_is_accepted_when_object_was_typed_otherwise():
    assert same_object({"object": "Proxima Cen", "drsobjn": "PROXIMA"}, "PROXIMA")
    assert not same_object({"object": "Proxima Cen"}, "PROXIMA")


def test_another_target_in_the_folder_is_still_skipped():
    assert not same_object({"object": "GL699", "drsobjn": "GL699"}, "PROXIMA")
    assert not same_object({"object": ""}, "PROXIMA")
