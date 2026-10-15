"""Une solution en longueur d'onde qui recule.

Le rééchantillonnage est un spline, qui exige des abscisses croissantes. La
plupart des fichiers les ont; pas tous. Mesuré le 2026-10-14: 5 des 60 poses
ESPRESSO de tau Ceti ont, au bord bleu de leur ordre le plus bleu, 177 pixels
dont la longueur d'onde DÉCROÎT avant de repartir (recul de 0,00026 nm contre
un pas médian de 0,00062, donc un bout de solution qui ne tient pas, pas une
oscillation numérique). scipy s'arrêtait sur « x must be increasing » et la
campagne entière s'arrêtait avec lui.
"""

import numpy as np

from pca2d.preprocess import register


def test_un_ordre_qui_recule_se_reechantillonne_quand_meme():
    """Trier remet chaque pixel à la longueur d'onde que sa solution lui
    donne, ce qui est la seule lecture honnête de ce qu'elle dit."""
    wave = np.linspace(500.0, 501.0, 200)
    # les dix premiers reculent, comme sur ESPRESSO
    wave[:10] = wave[10] - np.arange(10, 0, -1) * 1e-4
    wave[:10] = wave[:10][::-1]
    assert not np.all(np.diff(wave) > 0), "le cas à tester doit bien reculer"

    values = np.ones_like(wave)
    good = np.ones_like(wave, dtype=bool)
    grid = np.linspace(500.2, 500.8, 50)
    out, out_good = register(wave, values, good, 0.0, 0.0, grid)
    assert out_good.any()
    # un signal constant doit ressortir constant
    assert np.allclose(out[out_good], 1.0, atol=1e-6)


def test_deux_pixels_a_la_meme_longueur_d_onde_ne_comptent_qu_une_fois():
    """Un spline refuse deux abscisses égales, et le second pixel ne mesure
    rien que le premier n'ait déjà mesuré."""
    wave = np.linspace(500.0, 501.0, 100)
    wave[50] = wave[49]
    values = np.ones_like(wave)
    good = np.ones_like(wave, dtype=bool)
    grid = np.linspace(500.2, 500.8, 30)
    out, out_good = register(wave, values, good, 0.0, 0.0, grid)
    assert out_good.any()
    assert np.allclose(out[out_good], 1.0, atol=1e-6)


def test_un_ordre_deja_croissant_n_est_pas_touche():
    """Le tri ne doit rien changer là où il n'y a rien à trier: SPIRou, NIRPS
    et HARPS passent tous par ici à chaque pixel de chaque pose."""
    wave = np.linspace(1600.0, 1601.0, 300)
    rng = np.random.default_rng(7)
    values = 1.0 + 0.01 * rng.normal(size=wave.size)
    good = np.ones_like(wave, dtype=bool)
    grid = np.linspace(1600.2, 1600.8, 80)

    out, out_good = register(wave, values, good, 12.0, 0.0, grid)
    # le meme appel, avec des entrees explicitement triees: bit pour bit
    order = np.argsort(wave)
    same, same_good = register(wave[order], values[order], good[order],
                               12.0, 0.0, grid)
    assert np.array_equal(out, same)
    assert np.array_equal(out_good, same_good)


def test_les_pixels_qui_reculent_sont_retires_du_masque():
    """order_quality les rejette comme il rejette un flux négatif: le tri
    permet au spline d'exister, le masque dit qu'ils ne mesurent rien."""
    import inspect

    from pca2d import tfits

    source = inspect.getsource(tfits.order_quality)
    assert "rising" in source
    assert "np.diff(wave, axis=1) > 0" in source
