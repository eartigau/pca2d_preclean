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


def test_le_berv_vient_de_lbl_et_non_de_l_en_tete():
    """La question n'est pas « quel était le mouvement barycentrique », que
    l'en-tête donne toujours, mais « reste-t-il à l'appliquer ». Un t.fits de
    SPIRou est dans le repère de l'observateur; un S2D de l'ESO porte
    WAVEDATA_VAC_BARY, DÉJÀ corrigée. LBL a tranché classe par classe, donc on
    lui pose la question.

    Mesuré le 2026-10-14 sur HD 10700: en prenant le BERV de l'en-tête, les
    poses au BERV le plus négatif sortaient à 1.3 et 4 km/s du groupe et le
    rms de LBL était de 1311 m/s."""
    select = pytest.importorskip("lbl.instruments.select",
                                 reason="LBL est dans environment.yml")
    del select

    class Zero:
        """Une classe qui dit « déjà corrigé », comme HARPS_ESO."""
        def get_berv(self, header):
            return 0.0

    class FromHeader:
        """Une classe qui lit l'en-tête, comme SPIROU_CADC. LBL rend des m/s."""
        def get_berv(self, header):
            return 23742.8

    assert readers.berv_of(Zero(), {}, 23.7428) == 0.0
    assert readers.berv_of(FromHeader(), {}, None) == pytest.approx(23.7428)

    class Broken:
        def get_berv(self, header):
            raise RuntimeError("pas implémenté")

    # une classe qui ne répond pas laisse l'en-tête parler, et NaN sinon
    assert readers.berv_of(Broken(), {}, 12.5) == pytest.approx(12.5)
    assert np.isnan(readers.berv_of(Broken(), {}, None))


def test_defaire_la_correction_barycentrique_est_relativiste_et_reversible():
    """`input.to_observer` rend le spectre au repère de l'observateur, pour
    que l'ajustement à deux repères retrouve le contraste dont il vit. Sur
    HARPS, qui livre une grille déjà barycentrique, le BERV vu par pca2d
    valait zéro pour toutes les poses et la correction faisait passer le rms
    de LBL de 11.2 à 51.6 m/s; avec cette option, 14.2 m/s (2026-10-14).

    Le Doppler est relativiste: le premier ordre se tromperait de 1.5 m/s à
    30 km/s, ce qui compte ici."""
    c = 299792.458
    for v in (23.7428, -27.5, 0.0):
        beta = v / c
        assert readers._doppler(v) == pytest.approx(
            np.sqrt((1 + beta) / (1 - beta)), rel=1e-12)
    # et le premier ordre N'EST PAS la réponse
    assert readers._doppler(30.0) != pytest.approx(1 + 30.0 / c, rel=1e-12)

    import inspect
    source = inspect.getsource(readers.read)
    # la grille est divisée par le Doppler: barycentrique -> observateur
    assert "wave / _doppler(header_berv)" in source
    # et le BERV rendu redevient celui qu'il reste à appliquer
    assert "berv = float(header_berv)" in source


def test_defaire_la_correction_barycentrique_est_relativiste():
    """`input.to_observer` rend le spectre au repère de l'observateur, pour
    que l'ajustement à deux repères retrouve le contraste dont il vit. Sur
    HARPS, qui livre une grille déjà barycentrique, le BERV vu par pca2d
    valait zéro pour toutes les poses et la correction faisait passer le rms
    de LBL de 11.2 à 51.6 m/s; avec cette option, 14.2 (2026-10-14).

    Le Doppler est celui de grids, relativiste: le premier ordre se
    tromperait de 1.5 m/s à 30 km/s, ce qui compte ici."""
    from pca2d.grids import doppler

    for v in (23.7428, -27.5, 0.0):
        assert readers._doppler(v) == pytest.approx(float(doppler(v)), rel=1e-12)
    # et le premier ordre N'EST PAS la réponse
    assert readers._doppler(30.0) != pytest.approx(1 + 30.0 / 299792.458, rel=1e-12)

    import inspect
    source = inspect.getsource(readers.read)
    # la grille est divisée par le Doppler: barycentrique -> observateur
    assert "wave / _doppler(header_berv)" in source
    # et le BERV rendu redevient celui qu'il reste à appliquer
    assert "berv = float(header_berv)" in source


def test_to_observer_refuse_d_inventer_un_berv():
    """Défaire une correction barycentrique sans savoir de combien serait pire
    que de ne rien faire: le spectre sortirait décalé d'une valeur inconnue.
    On s'arrête plutôt que de deviner."""
    import inspect
    source = inspect.getsource(readers.read)
    assert "input.to_observer demande de defaire" in source
    assert "raise ValueError" in source


def test_to_observer_et_le_masque_a_la_majorite_vont_ensemble():
    """`input.to_observer` ramène la lecture au repère du sol, pour que
    l'ajustement travaille dans les deux repères comme sur SPIRou et NIRPS;
    seule la grille de départ diffère, et c'est pour cela qu'on la ramène à
    leur convention.

    IL A UN COÛT, ET CE COÛT A UN REMÈDE. Dans le repère du sol chaque pose
    couvre un intervalle décalé de son BERV, et `correct.mask: common` exige
    l'accord de TOUTES: l'intersection de 85 intervalles décalés ne garde que
    39,7 % des échantillons contre 61,5 % sur la grille livrée. Mesuré le
    2026-10-15 sur HD 10700, 64 poses communes:

        livrés                              1,56 m/s
        0-7 BERV nul                        1,89
        0-7 to_observer, unanimité          5,19   (43,0 % de NaN)
        0-7 to_observer, 75 % des poses     2,44   (17,0 % de NaN)
        0-7 to_observer, 50 % des poses     2,76   (14,2 % de NaN)

    Donc: un instrument qui active `to_observer` doit abaisser `mask_share`,
    sans quoi l'intersection mange ce que le levier rapporte.
    """
    import yaml

    with open("config.yaml", encoding="utf-8") as stream:
        loaded = yaml.safe_load(stream)
    table = loaded.get("instruments") or {}
    nominal = ((loaded.get("general") or {}).get("correct") or {})
    for name, block in table.items():
        if not (block.get("input") or {}).get("to_observer"):
            continue
        share = ((block.get("correct") or {}).get("mask_share")
                 if block.get("correct") else None)
        if share is None:
            share = nominal.get("mask_share", 1.0)
        mode = ((block.get("correct") or {}).get("mask")
                if block.get("correct") else None) or nominal.get("mask")
        assert mode != "common" or float(share) < 1.0, (
            "%s active to_observer avec un masque commun a l'unanimite:"
            " l'intersection de poses decalees coute plus que le levier ne"
            " rapporte (5,19 m/s contre 2,44 a 75%%)" % name)
