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


def _weak_line_problem(seed=3, rows=300, cols=3000, depth=0.03):
    """One component, a weak line 9 samples wide, pure noise elsewhere."""
    r = np.random.default_rng(seed)
    b = r.normal(5.0, 1.0, (rows, 1))
    sigma = r.uniform(0.5, 3.0, rows)[:, None] * np.ones((1, cols))
    w = 1.0 / sigma ** 2
    x = np.arange(cols, dtype=float)
    truth = depth * np.exp(-0.5 * ((x - 1500) / (9 / 2.3548)) ** 2)[None, :]
    data = b @ truth + sigma * r.normal(size=(rows, cols))
    fisher = (b[:, 0] ** 2) @ w
    Q = ((b[:, 0][:, None] * w * data).sum(axis=0) / fisher)[None, :]
    return Q, b, w, truth


def test_smoothing_the_significance_judges_a_line_over_its_width():
    """A line too weak per pixel is punched with holes pixel by pixel; judged
    over a resolution element it is kept along its whole width, and the noise,
    whose averaged z^2 gathers near 1, stays mostly dropped."""
    Q, b, w, _ = _weak_line_problem()
    line = slice(1496, 1505)
    noise = np.r_[0:1400, 1600:3000]
    raw = shrink_factors(Q, b, w)[0]
    smoothed = shrink_factors(Q, b, w, smooth_fwhm=9)[0]
    assert raw[line].min() == 0.0, "the test line is meant to be weak per pixel"
    assert smoothed[line].min() > 0.3, "a hole was left inside the line"
    assert smoothed[line].mean() > raw[line].mean()
    assert np.mean(smoothed[noise] == 0.0) > 0.4
    assert smoothed[noise].max() < raw[noise].max()


def test_a_smoothed_component_gains_significance_and_the_others_are_untouched():
    from pca2d.shrink import component_variance, smoothed_components
    Q, b, w, truth = _weak_line_problem()
    two = np.vstack([Q, Q[:, ::-1]])
    # an independent second amplitude: an affine copy of the first would make
    # the Fisher matrix nearly singular and every variance huge
    b2 = np.hstack([b, np.random.default_rng(11).normal(0.0, 3.0, b.shape)])
    var = component_variance(b2, w)
    out, factors = smoothed_components(two, var, [0], 9)
    assert np.array_equal(out[1], two[1]), "an unlisted component was changed"
    assert np.all(factors[1] == 1.0)
    line = slice(1496, 1505)
    assert factors[0][line].min() > 0.8, "the line was not kept along its width"
    assert factors[0][1498:1503].min() > 0.9
    near = slice(1400, 1600)
    assert (np.corrcoef(out[0][near], truth[0][near])[0, 1]
            > np.corrcoef(two[0][near], truth[0][near])[0, 1]), (
        "smoothing and shrinking did not bring the component closer to the truth")


def test_the_correct_stage_is_told_what_the_config_asks():
    from pca2d.cli import shrink_args
    cfg = {"correct": {"shrink": True, "shrink_smooth": True, "smooth_components": [2, 3]},
           "twoframe": {"resolution": 70000}}
    assert shrink_args(cfg) == ["--shrink", "--shrink-smooth", "--smooth-components",
                                "2,3", "--resolution", "70000.0"]
    assert shrink_args({"correct": {}, "twoframe": {"resolution": 70000}}) == []
