"""Refuser un arbre LBL dont un dossier est un lien MORT.

Le 2026-10-15, outputs/lbl/ portait dix liens vers un disque externe qui
n'etait pas monte. LBL s'est arrete au bout de 1.9 seconde sur une
FileExistsError en creant son dossier de log, et comme les quatre etapes
precedentes avaient reussi, la chaine avait l'air d'etre allee au bout.

Le piege: Path.exists() SUIT le lien, donc un lien casse "n'existe pas" et
io.make_dir de LBL part le creer; mais le nom est bien pris, et mkdir echoue.
"""
import os

import pytest

from pca2d.cli import _check_lbl_tree


def test_a_tree_with_no_links_passes(tmp_path):
    for name in ("templates", "masks", "lblrv"):
        (tmp_path / name).mkdir()
    _check_lbl_tree(str(tmp_path))          # ne leve pas


def test_a_live_link_passes(tmp_path):
    real = tmp_path / "ailleurs"
    real.mkdir()
    os.symlink(str(real), str(tmp_path / "templates"))
    _check_lbl_tree(str(tmp_path))


def test_a_dead_link_stops_the_run(tmp_path):
    os.symlink("/Volumes/pas_monte/lbl/templates",
               str(tmp_path / "templates"))
    with pytest.raises(SystemExit) as caught:
        _check_lbl_tree(str(tmp_path))
    assert caught.value.code == 2


def test_a_missing_directory_is_not_an_error(tmp_path):
    """L'arbre peut ne pas encore exister: ce n'est pas ce qu'on teste ici,
    et le faire echouer empecherait un premier run."""
    _check_lbl_tree(str(tmp_path / "pas_la"))


def test_every_dead_link_is_named(tmp_path, capsys):
    """Nommer le lien ET sa cible, puisque c'est ce que l'usager doit
    reparer; LBL, lui, ne nommait meme pas la cause."""
    for name in ("templates", "masks"):
        os.symlink("/Volumes/pas_monte/lbl/" + name, str(tmp_path / name))
    (tmp_path / "science").mkdir()
    with pytest.raises(SystemExit):
        _check_lbl_tree(str(tmp_path))
    said = capsys.readouterr().out
    assert "templates" in said and "masks" in said
    assert "/Volumes/pas_monte" in said
