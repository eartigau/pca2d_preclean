"""The window in four languages, and the buttons that choose between them.

English and French were there first; Spanish (with the Colombian flag) and
Portuguese (with Portugal's) came on 2026-09-16. No window is opened: a second
one would write its settings over those of the window already running.
"""

import re

import pytest

from pca2d.gui import (ALL_OPTIONS, EN, FR, LANGUAGES, STAGES, TEXTS, App,
                       text)
from pca2d.gui_es import ES
from pca2d.gui_pt import PT

PLACEHOLDER = re.compile(r"%[-+ 0#]*\d*(?:\.\d+)?[sdfgr%]")
FLAGS = {"en": "\U0001F1EC\U0001F1E7", "fr": "\U0001F1EB\U0001F1F7",
         "es": "\U0001F1E8\U0001F1F4", "pt": "\U0001F1F5\U0001F1F9"}


@pytest.mark.parametrize("name,table", [("es", ES), ("pt", PT)])
def test_every_english_text_is_translated_with_its_placeholders(name, table):
    assert set(EN) <= set(table), sorted(set(EN) - set(table))
    for key, english in EN.items():
        assert PLACEHOLDER.findall(table[key]) == PLACEHOLDER.findall(english), \
            (name, key)
        assert table[key].count("\n") == english.count("\n"), (name, key)
        assert (re.findall(r"<[^>]+>", table[key]).__len__()
                == re.findall(r"<[^>]+>", english).__len__()), (name, key)
        if key.startswith("help_"):
            assert len(table[key]) > 40, ("too short to explain", name, key)
    for key, _path, _kind in ALL_OPTIONS:
        assert table.get("opt_" + key), (name, key)
    for stage in STAGES:
        assert len(table["help_stage_" + stage]) > 40
    assert not any("—" in value for value in table.values()), "no em dash"


def test_four_languages_each_named_with_its_own_flag():
    assert LANGUAGES == ("en", "fr", "es", "pt")
    assert TEXTS == {"en": EN, "fr": FR, "es": ES, "pt": PT}
    for code in LANGUAGES:
        assert TEXTS[code]["lang"].startswith(FLAGS[code] + " "), code
    assert ES["lang"] == "\U0001F1E8\U0001F1F4 Español", "Colombia's flag"
    assert PT["lang"] == "\U0001F1F5\U0001F1F9 Português", "Portugal's flag"
    assert text("es", "run") == "Ejecutar" and text("pt", "run") == "Executar"
    assert text("de", "run") == "Run", "an unknown language falls back"


def test_portuguese_is_portugal_s():
    """ficheiro, not arquivo: the flag on the button is Portugal's."""
    words = set(re.findall(r"\w+", " ".join(PT.values()).lower()))
    assert "ficheiro" in words and "arquivo" not in words
    assert "registo" in words and "tela" not in words


def window_in(lang):
    window = App.__new__(App)
    said = []

    class Widget:
        def configure(self, **kwargs):
            said.append(kwargs.get("text"))

    window.labels = [(Widget(), "run", "text")]
    window.tabs = []
    window._draw_headings = lambda: None
    window._state = lambda: None
    window._sync = lambda: None
    window.lang = lang
    return window, said


def test_the_window_goes_round_all_four():
    window, said = window_in("en")
    for expected in ("fr", "es", "pt", "en"):
        window.switch_language()
        assert window.lang == expected
    assert said == ["Lancer", "Ejecutar", "Executar", "Run"]
    window.set_language("xx")
    assert window.lang == "en"


def test_the_buttons_are_the_other_languages():
    """At the top, one button per language the window is not in: click the
    one you want."""
    window, _said = window_in("es")
    made = []

    class Bar:
        def winfo_children(self):
            return list(made)

    class Button:
        def __init__(self, parent, text, command):
            self.text, self.command = text, command
            made.append(self)

        def pack(self, **_kwargs):
            pass

        def destroy(self):
            made.remove(self)

    class Ttk:
        pass

    Ttk.Button = Button
    window.ttk = Ttk
    window.lang_bar = Bar()
    window._tip = lambda widget, key: None
    window._draw_languages()
    assert [b.text for b in made] == [EN["lang"], FR["lang"], PT["lang"]]
    french = next(b for b in made if b.text == FR["lang"])
    french.command()
    assert window.lang == "fr"
    assert [b.text for b in made] == [EN["lang"], ES["lang"], PT["lang"]], \
        "drawn again: the others are now English, Spanish and Portuguese"


def test_each_language_writes_its_own_log_lines():
    window, _said = window_in("pt")
    line = App._line(window, "log_ended", 3)
    assert line.endswith("| a execução terminou, código de saída 3\n")
    window.lang = "es"
    assert App._line(window, "log_ended", 3).endswith(
        "| la ejecución terminó, código de salida 3\n")
