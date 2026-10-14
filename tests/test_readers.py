"""Lire n'importe quel instrument en passant par LBL.

« Tout ce qui passe dans le LBL peut être envoyé au PCA2D » (2026-10-14). Rien
ici ne lit un spectre: ce sont les règles du pont qui sont testées, pour que la
suite tourne sans données et sans réseau. La lecture elle-même est vérifiée à
la main sur de vrais fichiers, et ce qu'elle a donné est consigné ci-dessous.
"""
import numpy as np
import pytest

from pca2d import readers


def tmp_fits(hdus):
    """Un fichier FITS jetable, dans le dossier temporaire de pytest."""
    import tempfile
    from astropy.io import fits

    path = tempfile.mktemp(suffix=".fits")
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return path


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


def test_le_temps_et_la_pose_ne_viennent_pas_de_la_table_de_lbl():
    """Mesuré le 2026-10-14 sur un S2D HARPS: la classe HARPS de LBL fait
    pointer KW_MJDATE *et* KW_EXPTIME sur le même mot-clé, HIERARCH ESO QC BJD,
    parce que LBL n'a besoin que d'un axe temporel. pca2d s'en sert autrement:
    une pose de 60 s lue comme 2452859.8 fausserait toute pondération, et un
    BJD pris pour un MJD décalerait la campagne de 2.4 millions de jours."""
    import inspect
    source = inspect.getsource(readers.read)
    assert 'plain("MJD-OBS")' in source
    assert 'plain("EXPTIME")' in source
    # un BJD qui aurait échappé au mot-clé standard est ramené en MJD
    assert "2400000.5" in source


def test_la_transmission_n_est_jamais_none():
    """order_quality et _order_sigma indexent recon comme un tableau: un None
    les ferait planter. Sans modèle, la réponse honnête est une transmission
    de 1 (« aucune absorption connue ici »), et meta dit laquelle des deux."""
    import inspect
    source = inspect.getsource(readers.read)
    assert "np.ones_like(flux)" in source
    assert "telluric_model" in source


def test_le_dossier_de_donnees_du_depot_ne_doit_pas_masquer_lbl():
    """Ce dépôt contient un dossier `lbl/` de DONNÉES (calib, masks, models,
    science...). Dès que sa racine est sur sys.path, `import lbl` en fait un
    paquet-espace-de-noms: lbl.science est masqué alors que lbl.instruments ne
    l'est pas, donc l'instrument se charge et la transmission non, sans que
    rien ne le dise. Mesuré le 2026-10-14."""
    pytest.importorskip("lbl.instruments.select",
                        reason="LBL est dans environment.yml")
    assert readers._lbl_root() is not None, "le vrai paquet LBL reste introuvable"
    module = readers._lbl_science("tellu_clean")
    assert module is not None, "lbl.science reste masqué par le dossier de données"
    assert "lbl/science/tellu_clean.py" in str(getattr(module, "__file__", ""))


def test_l_extension_du_flux_est_retrouvee_et_non_devinee():
    """Le fichier corrigé est une copie de l'original dont on remplace le flux,
    donc il faut savoir laquelle des extensions le porte. Le nom change d'un
    instrument à l'autre (FluxAB sur SPIRou, SCIDATA sur un S2D de l'ESO) et
    LBL ne le nomme pas: sa classe par défaut prend la première extension qui
    porte des données. On la retrouve donc par la forme, comme astropy l'a
    choisie, et None quand rien ne correspond: l'appelant refuse alors
    d'écrire plutôt que d'écrire au mauvais endroit."""
    from astropy.io import fits

    flux = np.arange(12, dtype=float).reshape(3, 4)
    path = tmp_fits([fits.PrimaryHDU(),
                     fits.ImageHDU(flux, name="SCIDATA"),
                     fits.ImageHDU(np.zeros((3, 4)), name="ERRDATA")])
    assert readers.flux_extension(path, flux) == "SCIDATA"
    # une forme qu'aucune extension ne porte: rien, et surtout pas un nom
    assert readers.flux_extension(path, np.zeros((7, 7))) is None


def test_les_sentinelles_deviennent_des_nan():
    """Un en-tête qui n'a pas mesuré quelque chose l'écrit par une valeur
    convenue: -1 pour un seeing que l'ESO n'a pas (vu le 2026-10-14 sur un S2D
    HARPS), -9999.9 pour un mât météo en panne. La prendre pour un nombre
    mettrait un seeing négatif dans un diagnostic."""
    assert -1.0 in readers.SENTINELS
    for value in (-9999.0, -9999.9, -999.0):
        assert value in readers.SENTINELS


def test_le_contrat_de_meta_est_celui_de_read_tfits():
    """cube.py lit meta["bjd"], meta["berv"], meta["object"]... sans savoir
    d'où vient le spectre. Une clé manquante n'est pas une lecture dégradée,
    c'est un KeyError au milieu d'une campagne: mesuré le 2026-10-14, le pont
    s'est arrêté sur meta["bjd"] après avoir lu 50 fichiers."""
    import inspect
    from pca2d import tfits

    source = inspect.getsource(readers.read)
    for name in list(tfits.META_KEYS) + ["object", "drsobjn", "airmass"]:
        assert '"%s"' % name in source, name
