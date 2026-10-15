"""Qui gagne quand deux ordres de MEME parité mesurent la même longueur d'onde.

pca2d met les ordres pairs sur une ligne et les impairs sur l'autre, ce qui
suppose que deux ordres de même parité ne se touchent jamais. Vrai sur SPIRou
et HARPS; FAUX sur ESPRESSO, dont les 170 ordres sont 85 ordres physiques
coupés en deux tranches: à l'intérieur d'une tranche, deux ordres consécutifs
se recouvrent de 3,4 nm (mesuré le 2026-10-14 sur un S2D de tau Ceti), dont
1,8 nm survivent à l'érosion des bords.

Sans arbitrage, les deux écriraient dans la même colonne de la grille. La
règle: celui qui porte le plus de signal la garde, l'autre passe à NaN.
"""

import numpy as np

from pca2d.tfits import _overlap_blaze


def zone():
    """Une grille, un ordre dessus, et le masque qui le dit tout bon."""
    grid = np.linspace(500.0, 510.0, 1000)
    wave = np.linspace(500.0, 506.0, 200)
    return grid, wave, np.ones(200, dtype=bool)


def test_le_flux_tranche_quand_le_blaze_est_plat():
    """Un S2D de l'ESO a un blaze qui vaut UN partout, parce que la classe de
    LBL ne déblaze pas: il ne départage rien. C'est alors le flux médian qui
    dit lequel des deux ordres est au milieu et lequel est sur son bord.
    Mesuré sur ESPRESSO, zone commune aux ordres 0 et 2: 16036 contre 11749."""
    grid, wave, ok = zone()
    plat = np.ones(200)
    fort = _overlap_blaze(wave, plat, ok, grid, 0, 500, np.full(200, 1000.0))
    faible = _overlap_blaze(wave, plat, ok, grid, 0, 500, np.full(200, 100.0))
    assert fort > faible
    assert fort == 1000.0 and faible == 100.0


def test_le_blaze_repond_seul_quand_il_dit_quelque_chose():
    """Sur un t.fits d'APERO le blaze décrit la réponse de l'ordre et suffit:
    le flux ne doit pas s'en mêler, sans quoi une raie profonde ferait perdre
    un ordre qui est pourtant au milieu de sa course."""
    grid, wave, ok = zone()
    courbe = np.linspace(0.2, 1.0, 200)
    avec_peu = _overlap_blaze(wave, courbe, ok, grid, 0, 500, np.full(200, 1.0))
    avec_beaucoup = _overlap_blaze(wave, courbe, ok, grid, 0, 500,
                                   np.full(200, 1e6))
    assert avec_peu == avec_beaucoup
    assert 0.2 < avec_peu < 1.0


def test_une_zone_vide_ne_designe_personne():
    """Hors de l'ordre il n'y a rien à comparer, et rendre un nombre ferait
    gagner ou perdre un ordre sur une zone qu'il ne mesure pas."""
    grid, wave, ok = zone()
    plat = np.ones(200)
    # la zone [900, 950] de la grille est au-delà de 506 nm
    assert _overlap_blaze(wave, plat, ok, grid, 900, 950, np.ones(200)) is None
    # et une zone à l'envers non plus
    assert _overlap_blaze(wave, plat, ok, grid, 500, 400, np.ones(200)) is None


def test_les_deux_ordres_sont_mesures_sur_la_MEME_zone():
    """Le piège qu'il a fallu corriger: comparer un ordre moyenné sur toute sa
    largeur à un autre pris sur son bord ferait gagner le premier à tous les
    coups. resample_exposure garde donc de quoi redemander au précédent ce
    qu'il porte sur la zone du recouvrement, et non ce qu'il portait ailleurs."""
    import inspect

    from pca2d import tfits

    source = inspect.getsource(tfits.resample_exposure)
    # le precedent est reinterroge avec SES longueurs d'onde et SON masque
    assert "last_order[parity] = (w_shift, ok, order)" in source
    assert "p_shift, p_ok, p_order = previous" in source
    assert "_overlap_blaze(p_shift, blaze[p_order], p_ok, grid," in source
