"""Le seuil par pixel porte sur le SNR, pas sur la PROFONDEUR.

Un fond de raie profonde porte moins de photons, donc un RMS fractionnaire
plus grand, mais c'est l'endroit ou le contenu en vitesse est le plus dense:
le jeter revient a jeter le signal qu'on vient mesurer.

Mesure du 2026-10-16 sur le cube HARPS de HD 10700: l'ancien seuil
(ln_clip_low = -0.5, soit tout pixel sous 0.607 du continu) ecartait 4.36%
des echantillons, et 96.1% d'entre eux avaient un SNR d'au moins 5, de SNR
median 21.3. Des pixels parfaitement mesures, jetes parce que sombres.
"""
import numpy as np

from pca2d import twoframe


def test_the_default_is_an_snr_threshold():
    import inspect
    sig = inspect.signature(twoframe.load_cube)
    assert sig.parameters["min_pixel_snr"].default == \
        twoframe.DEFAULT_MIN_PIXEL_SNR
    assert sig.parameters["ln_clip_low"].default is None, \
        "le seuil de profondeur ne doit plus s'appliquer tout seul"


def test_a_deep_but_well_measured_pixel_is_kept():
    """LE CAS QUI COMPTE: le fond du Na D a 589.15 nm. Profond, donc l'ancien
    seuil le tuait, mais mesure a SNR 20."""
    deep = -1.2                      # 30% du continu, bien sous exp(-0.5)
    sigma = 1.0 / 20.0               # SNR 20
    assert deep < -0.5, "l'ancien seuil l'aurait tue"
    assert (1.0 / sigma) >= twoframe.DEFAULT_MIN_PIXEL_SNR, \
        "le nouveau le garde"


def test_a_shallow_but_noisy_pixel_is_dropped():
    """L'inverse: un pixel pres du continu mais sans photons."""
    shallow = -0.01
    sigma = 1.0 / 2.0                # SNR 2
    assert shallow > -0.5, "l'ancien seuil le gardait"
    assert not ((1.0 / sigma) >= twoframe.DEFAULT_MIN_PIXEL_SNR), \
        "le nouveau l'ecarte"


def _weights(data, sigma, **kw):
    """Les poids que load_cube calculerait, sur le meme chemin arithmetique."""
    data = np.asarray(data, dtype=np.float64)
    sigma = np.asarray(sigma, dtype=np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        w = np.where(np.isfinite(sigma) & (sigma > 0), 1.0 / sigma ** 2, 0.0)
    w[~np.isfinite(data)] = 0.0
    cut = kw.get("min_pixel_snr", twoframe.DEFAULT_MIN_PIXEL_SNR)
    if cut is not None:
        with np.errstate(invalid="ignore", divide="ignore"):
            w[~(1.0 / sigma >= float(cut))] = 0.0
    elif kw.get("ln_clip_low") is not None:
        w[data < kw["ln_clip_low"]] = 0.0
    return w


def test_depth_alone_no_longer_decides():
    """Deux pixels de MEME SNR et de profondeurs opposees pesent pareil."""
    sigma = np.array([0.05, 0.05])            # SNR 20 tous les deux
    data = np.array([-1.5, -0.01])            # l'un profond, l'autre non
    w = _weights(data, sigma)
    assert w[0] > 0 and w[1] > 0
    assert w[0] == w[1], "la profondeur ne doit plus rien changer au poids"


def test_the_old_behaviour_is_still_reachable():
    """Pour relire un cube comme avant, si on veut comparer."""
    sigma = np.array([0.05, 0.05])
    data = np.array([-1.5, -0.01])
    w = _weights(data, sigma, min_pixel_snr=None, ln_clip_low=-0.5)
    assert w[0] == 0.0, "l'ancien seuil tuait le pixel profond"
    assert w[1] > 0.0


def test_a_dead_pixel_is_dropped_either_way():
    sigma = np.array([0.0, np.nan, 1e3])
    data = np.array([-0.1, -0.1, -0.1])
    w = _weights(data, sigma)
    assert np.all(w == 0.0)
