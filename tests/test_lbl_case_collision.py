"""Nommer les produits LBL qu'une collision de CASSE fait relire.

Le 2026-10-15, un run sur HD10700_tc avait laisse un template, un masque et
un .rdb. Le run sur HD10700_TC les a pris pour les siens, le systeme de
fichiers de macOS ne distinguant pas la casse. Le pire etant le .rdb: LBL
venait de calculer 99 vitesses fraiches, et il a compile un tableau de 74
lignes a partir du fichier de l'autre run, sans rien signaler.
"""
import os

from pca2d.lbl import _warn_case_collisions


def tree(tmp_path, **files):
    """Un arbre LBL avec les fichiers demandes, par dossier."""
    for sub, names in files.items():
        folder = tmp_path / sub
        folder.mkdir(parents=True, exist_ok=True)
        for name in names:
            (folder / name).write_text("")
    return str(tmp_path)


def test_a_product_of_this_run_says_nothing(tmp_path, capsys):
    root = tree(tmp_path, templates=["Template_HD10700_TC_HARPS_ESO.fits"])
    _warn_case_collisions(root, ["HD10700_TC"])
    assert "ONLY BY CASE" not in capsys.readouterr().out


def test_a_product_differing_only_by_case_is_named(tmp_path, capsys):
    root = tree(tmp_path, templates=["Template_HD10700_tc_HARPS_ESO.fits"])
    _warn_case_collisions(root, ["HD10700_TC"])
    said = capsys.readouterr().out
    assert "ONLY BY CASE" in said
    assert "Template_HD10700_tc_HARPS_ESO.fits" in said
    assert "HD10700_TC" in said and "HD10700_tc" in said


def test_a_different_object_says_nothing(tmp_path, capsys):
    """Les autres campagnes de l'arbre ne sont pas des collisions."""
    root = tree(tmp_path, templates=["Template_HD10700_N2_HARPS_ESO.fits",
                                     "Template_PROXIMA_HARPS_ESO.fits"])
    _warn_case_collisions(root, ["HD10700_TC"])
    assert "ONLY BY CASE" not in capsys.readouterr().out


def test_every_named_folder_is_looked_at(tmp_path, capsys):
    """Le .rdb est le plus dangereux: il fait compiler un tableau faux."""
    root = tree(tmp_path,
                templates=["Template_HD10700_tc_HARPS_ESO.fits"],
                masks=["LBL_Mask_HD10700_tc_full_harps.fits"],
                lblrdb=["lbl_HD10700_tc_HD10700_tc.rdb"],
                lblreftable=["ref_table_LBL_Mask_HD10700_tc_full_harps.csv"])
    _warn_case_collisions(root, ["HD10700_TC"])
    said = capsys.readouterr().out
    for sub in ("templates", "masks", "lblrdb", "lblreftable"):
        assert sub in said, sub


def test_both_object_names_are_checked(tmp_path, capsys):
    """La serie brute ET la serie corrigee."""
    root = tree(tmp_path, lblrdb=["lbl_hd10700_tc_pca2d_0-2_x.rdb"])
    _warn_case_collisions(root, ["HD10700_TC", "HD10700_TC_PCA2D_0-2"])
    assert "ONLY BY CASE" in capsys.readouterr().out


def test_a_missing_tree_is_not_an_error(tmp_path):
    _warn_case_collisions(str(tmp_path / "pas_la"), ["HD10700_TC"])
