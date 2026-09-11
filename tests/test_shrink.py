"""Observer components are kept where the data detect them and nowhere else."""
import numpy as np

from pca2d.shrink import shrink_factors


def _problem(seed=2, rows=300, cols=4000):
    """Two components with a real feature each at a few columns, pure noise
    elsewhere, every row with its own noise level; Q estimated per column by
    weighted least squares, as the fit's observer update does."""
    r = np.random.default_rng(seed)
    b = r.normal(size=(rows, 2)) * [3.0, 1.0] + [5.0, 0.0]
    sigma = r.uniform(0.5, 3.0, rows)[:, None] * np.ones((1, cols))
    w = 1.0 / sigma ** 2
    truth = np.zeros((2, cols))
    truth[0, 100:110] = 0.2           # a clearly detected feature in component 1
    truth[1, 2000:2010] = 0.3         # and one in component 2, elsewhere
    data = b @ truth + sigma * r.normal(size=(rows, cols))
    fisher = np.einsum("rc,rj,rk->cjk", w, b, b)
    rhs = np.einsum("rc,rj,rc->cj", w, b, data)
    Q = np.linalg.solve(fisher, rhs[..., None])[..., 0].T
    return Q, b, w


def test_a_detected_feature_is_kept_and_noise_is_dropped():
    Q, b, w = _problem()
    s = shrink_factors(Q, b, w)
    assert s[0, 100:110].min() > 0.9, "a significant feature was shrunk"
    assert s[1, 2000:2010].min() > 0.9
    noise = np.ones(Q.shape[1], dtype=bool)
    noise[90:120] = noise[1990:2020] = False
    for j in (0, 1):
        assert np.median(s[j, noise]) == 0.0
        assert s[j, noise].mean() < 0.25


def test_each_component_is_judged_on_its_own():
    """At the columns of component 1's feature, component 2 has nothing: it is
    dropped there while component 1 is kept."""
    Q, b, w = _problem()
    s = shrink_factors(Q, b, w)
    assert s[0, 100:110].mean() > 0.9
    assert s[1, 100:110].mean() < 0.5


def test_a_row_with_almost_no_weight_changes_almost_nothing():
    Q, b, w = _problem()
    loud = w.copy()
    loud[0] *= 1e-8                    # the first row made enormously noisy
    s_all = shrink_factors(Q, b, w)
    s_loud = shrink_factors(Q, b, loud)
    s_without = shrink_factors(Q, np.delete(b, 0, axis=0), np.delete(w, 0, axis=0))
    assert np.allclose(s_loud, s_without, atol=1e-6)
    assert not np.allclose(s_all, s_without, atol=1e-6)


def test_a_larger_reduced_chi2_shrinks_more():
    Q, b, w = _problem()
    assert shrink_factors(Q, b, w, chi2_scale=2.0).sum() < shrink_factors(Q, b, w).sum()
