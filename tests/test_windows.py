"""A window written centre:width is a window, not a base-60 number.

PyYAML follows YAML 1.1, where an unquoted `1267:2` is 1267 * 60 + 2 = 76022.
Three of the eight windows in config.yaml were read that way, landed tens of
thousands of nanometres outside the grid, and were skipped by the figure
scripts in a subprocess whose output nobody sees on success: the PDF simply
had five windows instead of eight. The ones that survived were the ones whose
centre happened to carry a decimal point.
"""

import yaml

from pca2d.config import check_windows, read_yaml


def test_a_window_stays_the_text_it_was_written_as():
    doc = read_yaml("windows: [1200.3:2, 1220:4, 1267:2, 2450:5, 2450:0.5]")
    assert doc["windows"] == ["1200.3:2", "1220:4", "1267:2", "2450:5",
                              "2450:0.5"]


def test_the_trap_is_real_with_the_ordinary_loader():
    """So nobody later swaps read_yaml back for yaml.safe_load as a tidy-up."""
    doc = yaml.safe_load("windows: [1267:2, 2450:0.5]")
    assert doc["windows"] == [76022, 147000.5]


def test_nothing_else_in_the_real_config_changes_type():
    """Only the base-60 windows may differ between the two loaders.

    config.yaml quotes its windows since 2026-09-11, so nothing differs at all
    any more; the claim that stays is that the package's loader changes the type
    of nothing else. That it reads an unquoted 1220:4 as text is tested above.
    """
    def walk(a, b, path=""):
        if isinstance(a, dict) and isinstance(b, dict):
            for key in set(a) | set(b):
                yield from walk(a.get(key), b.get(key), "%s.%s" % (path, key))
        elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
            for i, (x, y) in enumerate(zip(a, b)):
                yield from walk(x, y, "%s[%d]" % (path, i))
        elif a != b or type(a) is not type(b):
            yield path, a, b

    with open("config.yaml") as handle:
        old = yaml.safe_load(handle)
    with open("config.yaml") as handle:
        new = read_yaml(handle)
    differences = list(walk(old, new))
    for path, before, after in differences:
        assert path.startswith(".general.output.windows"), (
            "%s changed from %r to %r" % (path, before, after))
        assert isinstance(before, int) and isinstance(after, str)


def test_ordinary_numbers_are_still_numbers():
    doc = read_yaml("a: 2\nb: 0.5\nc: -0.5\nd: 16.0\ne: .inf\nf: true\ng: null"
                    "\nh: 0x1f\ni: 1_000")
    assert doc == {"a": 2, "b": 0.5, "c": -0.5, "d": 16.0, "e": float("inf"),
                   "f": True, "g": None, "h": 31, "i": 1000}


def test_a_window_no_figure_can_show_is_said_out_loud(capsys):
    domain = {"wave_min": 955.0, "wave_max": 2500.0}
    check_windows(["1267:2", 76022, "3100:2", "not a window"], domain)
    out = capsys.readouterr().out
    assert "output.windows: 1267:2" not in out, "a good window was complained about"
    assert "output.windows: 76022" in out
    assert "YAML 1.1" in out, "the number did not get the base-60 explanation"
    assert "output.windows: 3100:2" in out and "outside the domain" in out
    assert "not centre:width" in out


def test_the_sequence_figure_has_the_panel_without_the_star():
    """Under panel 3, the star block taken out instead of the observer one."""
    import inspect

    from pca2d.figures import sequence

    source = inspect.getsource(sequence)       # the panels live in window_arrays
    three = source.index("3. minus the OBSERVER block")
    four = source.index("4. minus the STAR block")
    five = source.index("5. minus everything")
    assert three < four < five
