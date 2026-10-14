"""Lire n'importe quel instrument en passant par LBL.

« Tout ce qui passe dans le LBL peut être envoyé au PCA2D » (2026-10-14). Rien
ici ne lit un spectre: ce sont les règles du pont qui sont testées, pour que la
suite tourne sans données et sans réseau. La lecture elle-même est vérifiée à
la main sur de vrais fichiers, et ce qu'elle a donné est consigné ci-dessous.
"""
import numpy as np
import pytest

from pca2d import readers


def test_les_instruments_de_lbl_sont_ceux_que_pca2d_annonce():
    """La table KNOWN doit rester celle de LBL, sinon elle ment."""
    select = pytest.importorskip("lbl.instruments.select",
                                 reason="LBL est dans environment.yml")
    for name, sources in readers.KNOWN.items():
        assert name in select.InstDict, name
        for source in sources:
            assert source in select.InstDict[name], "%s/%s" % (name, source)
    # les quatre que l'utilisateur a nommés le 2026-10-14
    for name in ("HARPS", "HARPSN", "SOPHIE", "ESPRESSO"):
        assert name in readers.KNOWN, name


def test_les_parites_sont_verifiees_et_non_supposees():
    """pca2d met les ordres pairs sur une ligne et les impairs sur l'autre, ce
    qui suppose que deux ordres de MEME parité ne se recouvrent jamais. Mesuré
    le 2026-10-14: vrai sur HARPS (0.000 nm), FAUX sur SPIRou dans le bleu, où
    14 paires se recouvrent jusqu'à 5.4 nm (ordres 0 à 13, 980-1150 nm), alors
    que tfits.py affirmait le contraire. D'où cette vérification."""
    # trois ordres qui ne se recouvrent pas deux à deux
    wave = np.array([[500.0, 510.0], [505.0, 515.0], [511.0, 521.0]])
    ok, worst = readers.parity_is_safe(wave)
    assert ok and worst <= 0

    # l'ordre 0 déborde sur l'ordre 2: c'est le cas SPIRou
    wave = np.array([[980.0, 986.2], [993.0, 998.8], [980.7, 991.0]])
    ok, worst = readers.parity_is_safe(wave)
    assert not ok
    assert worst == pytest.approx(5.5, abs=0.1)


def test_le_snr_est_mesure_quand_l_en_tete_ne_le_porte_pas():
    """APERO écrit EXTSNxxx, HARPS non: il faut alors le mesurer sur le flux."""
    rng = np.random.default_rng(3)
    flux = 1000.0 + rng.normal(scale=10.0, size=(2, 2000))
    snr = readers.order_snr(flux)
    assert snr.shape == (2,)
    # 1000 / 10 = 100, à la précision de l'estimateur robuste
    assert np.all(np.abs(snr - 100) < 15), snr
    # celui de l'en-tête a toujours priorité
    given = np.array([42.0, 43.0])
    assert np.allclose(readers.order_snr(flux, header_snr=given), given)
    # un ordre vide ne fabrique pas de valeur
    assert np.isnan(readers.order_snr(np.full((1, 4), np.nan))[0])


def test_le_contrat_de_lecture_est_celui_des_t_fits():
    """read() doit rendre exactement ce que tfits.read_tfits rend, sinon le
    reste du pipeline ne peut pas s'en servir sans le savoir."""
    import inspect
    source = inspect.getsource(readers.read)
    for key in ("flux", "wave", "blaze", "recon", "sky", "snr", "meta"):
        assert "%s=" % key in source, key
    # recon et sky peuvent manquer: aucun spectrographe visible ne les livre,
    # et l'appelant doit l'accepter plutôt que de planter
    assert "sky=None" in source
