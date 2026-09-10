"""The timestamp is on every printed line, including the ones a wrap makes.

A message longer than the terminal used to run off the edge or wrap wherever
the terminal chose, leaving continuation text with no timestamp in front of
it. Now it is wrapped to the window, never inside a word, and each line carries
the stamp, so the convention `YYMMDD HH:MM:SS.SS | message` holds line by line.
"""

import re

from pca2d import logger

STAMP = re.compile(r"^\d{6} \d{2}:\d{2}:\d{2}\.\d{2} \| ")
LONG = ("the star is never exactly where the BERV alone would put it: it has a"
        " planet, it has activity, and the carry operator leaves a residual of"
        " its own")


def printed(capsys):
    return capsys.readouterr().out.rstrip("\n").split("\n")


def test_a_long_message_wraps_with_the_stamp_on_every_line(monkeypatch, capsys):
    monkeypatch.setattr(logger, "_terminal_width", lambda: 70)
    logger.log(LONG)
    lines = printed(capsys)
    assert len(lines) > 1, "a message wider than the window was not wrapped"
    assert all(STAMP.match(line) for line in lines), "a line lost its stamp"
    assert len({STAMP.match(line).group(0) for line in lines}) == 1, (
        "one message, so one timestamp, repeated")
    assert all(len(line) <= 70 for line in lines)
    assert " ".join(STAMP.sub("", l) for l in lines) == LONG, "text was lost"


def test_nothing_is_wrapped_without_a_window(monkeypatch, capsys):
    """A log file has no width, and a grep looks for whole lines."""
    monkeypatch.setattr(logger, "_terminal_width", lambda: None)
    logger.log(LONG * 3)
    assert len(printed(capsys)) == 1


def test_a_path_or_a_name_is_never_cut(monkeypatch, capsys):
    monkeypatch.setattr(logger, "_terminal_width", lambda: 60)
    path = "outputs/TOI-2120/2-7v/corrected/2811170t_0-7_with_a_long_name.fits"
    logger.log("writing the corrected spectrum of TOI-2120 (2-7v) to %s now" % path)
    out = capsys.readouterr().out
    assert path in out, "the path was split and can no longer be copied"
    assert "TOI-2120" in out and "(2-7v)" in out


def test_an_indented_message_stays_indented_when_it_wraps(monkeypatch, capsys):
    monkeypatch.setattr(logger, "_terminal_width", lambda: 60)
    logger.log("  iter 3  R2=0.548569  left/raw=0.3730  clipped=0.756%"
               "  cond(A) med/max 2.9/9.8  shift 160.4 m/s rms")
    lines = printed(capsys)
    assert len(lines) > 1
    assert all(STAMP.sub("", line).startswith("  ") for line in lines)


def test_a_window_too_narrow_to_help_is_not_used(monkeypatch, capsys):
    """Wrapping into a 10-character column is worse than not wrapping."""
    monkeypatch.setattr(logger, "_terminal_width", lambda: 30)
    logger.log(LONG)
    assert len(printed(capsys)) == 1


# ------------------------------------------------ and the bars around them ---
def test_a_message_asks_the_bars_to_step_aside(monkeypatch, capsys):
    """Printed straight, a message landed on the end of the bar being drawn."""
    calls = []

    class Recorder:
        @staticmethod
        def write(text, file=None):
            calls.append(text)

    monkeypatch.setattr(logger, "_tqdm", Recorder)
    logger.log("iter 3  R2=0.548569")
    assert len(calls) == 1 and "iter 3  R2=0.548569" in calls[0], (
        "log did not go through tqdm.write, so a live bar would be drawn over")


def test_a_bar_inside_a_sweep_says_which_sweep():
    from pca2d import progress

    progress.set_label("sweep 3")
    try:
        assert progress._labelled("coefficients") == "sweep 3: coefficients"
    finally:
        progress.set_label(None)
    assert progress._labelled("coefficients") == "coefficients"


def test_there_is_no_bar_over_the_sweeps():
    """One bar per phase of a sweep, and none over all sixteen of them."""
    import inspect

    from pca2d import twoframe

    source = inspect.getsource(twoframe.main)
    assert 'desc="fitting"' not in source
    assert "sweeps." not in source
