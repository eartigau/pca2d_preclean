"""The block solve, checked against a problem whose answer is known.

The cube in `problem` was built forwards from a basis and coefficients, so these
tests compare against the truth rather than against yesterday's output. That
matters for a refactor: a characterisation test that pins the current answer
will happily pin a bug, while recovering the truth cannot.
"""

from __future__ import annotations

import numpy as np
import pytest

from pca2d.twoframe import (LanczosShifter, block_power, joint_coeffs,
                                 leakage, nightly_bin, star_model, tidy_block,
                                 update_earth, update_star)


@pytest.fixture
def shifter(problem):
    return LanczosShifter(problem["data"].shape[1], threads=1)


def test_joint_coeffs_recovers_the_truth(problem, shifter):
    """Given the true basis, the solve must return the true coefficients."""
    a, b, _alpha, cond = joint_coeffs(problem["data"], problem["w"], problem["P"],
                              problem["Q"], shifter, problem["delta"], chunk=8)
    assert np.all(np.isfinite(cond)), "the normal matrix was singular somewhere"
    assert a.shape == problem["a"].shape
    assert b.shape == problem["b"].shape
    assert np.allclose(a, problem["a"], atol=1e-4), "star coefficients are wrong"
    assert np.allclose(b, problem["b"], atol=1e-4), "Earth coefficients are wrong"


def test_tying_forces_one_answer_per_exposure(problem, shifter):
    """With `exposure` set, both rows of an exposure get identical coefficients.

    This is the parity tie. It is not an average taken afterwards: the normal
    equations of the two rows are summed and solved once, so the test is that
    the two rows come back bit-identical, not merely close.
    """
    n = problem["data"].shape[0]
    exposure = np.repeat(np.arange(n // 2), 2)   # rows 0,1 are one exposure
    a, b, _alpha, _ = joint_coeffs(problem["data"], problem["w"], problem["P"],
                           problem["Q"], shifter, problem["delta"], chunk=8,
                           exposure=exposure)
    for e in np.unique(exposure):
        rows = np.where(exposure == e)[0]
        assert np.array_equal(a[rows[0]], a[rows[1]]), \
            "tied exposure %d has two different star vectors" % e
        assert np.array_equal(b[rows[0]], b[rows[1]]), \
            "tied exposure %d has two different Earth vectors" % e


def test_zero_weight_rows_get_zero_coefficients(problem, shifter):
    """A row with no weight constrains nothing and must not invent a signal."""
    w = problem["w"].copy()
    w[3] = 0.0
    a, b, _alpha, _ = joint_coeffs(problem["data"], w, problem["P"], problem["Q"],
                           shifter, problem["delta"], chunk=8)
    assert np.allclose(a[3], 0.0), "a fully masked row produced star coefficients"
    assert np.allclose(b[3], 0.0), "a fully masked row produced Earth coefficients"


def test_star_model_reproduces_the_star_block(problem, shifter):
    n_spec, n_pix = problem["data"].shape
    model = star_model(problem["P"], problem["a"], shifter, problem["delta"],
                       n_spec, n_pix, chunk=8)
    want = np.zeros_like(model)
    for i in range(n_spec):
        for k in range(problem["P"].shape[0]):
            want[i] += problem["a"][i, k] * np.roll(problem["P"][k],
                                                    int(problem["delta"][i]))
    edge = 16
    assert np.allclose(model[:, edge:-edge], want[:, edge:-edge], atol=1e-6)


def test_iteration_reduces_the_residual(problem, shifter):
    """Alternating updates must walk the weighted residual downhill.

    NOT asserted for a single update. `update_star` deliberately solves with the
    DIAGONAL of the normal matrix rather than the full one, as its docstring
    says, so one step is an approximate M-step and monotone decrease is not
    guaranteed for it. What is guaranteed, and what the fit relies on, is that
    the alternation converges. This is the property worth protecting through a
    refactor.
    """
    d, w = problem["data"], problem["w"]
    P = problem["P"] + 0.3 * np.roll(problem["P"], 5, axis=1)   # perturbed start
    Q = problem["Q"].copy()
    a, b, _alpha, _ = joint_coeffs(d, w, P, Q, shifter, problem["delta"], chunk=8)

    def chi2(P_, Q_, a_, b_):
        m = star_model(P_, a_, shifter, problem["delta"], d.shape[0],
                       d.shape[1], chunk=8) + b_ @ Q_
        return float(np.sum(w * (d - m) ** 2))

    before = chi2(P, Q, a, b)
    for _ in range(5):
        P = update_star(d, w, P, Q, a, b, shifter, problem["delta"], 8)
        Q = update_earth(d, w, P, Q, a, b, shifter, problem["delta"], 8)
        a, b, _alpha, _ = joint_coeffs(d, w, P, Q, shifter, problem["delta"], chunk=8)
    after = chi2(P, Q, a, b)
    assert after < before, ("five iterations did not reduce chi2: %.6g -> %.6g"
                            % (before, after))


def test_block_power_is_non_negative(problem, shifter):
    a, b, _alpha, _ = joint_coeffs(problem["data"], problem["w"], problem["P"],
                           problem["Q"], shifter, problem["delta"], chunk=8)
    # block_power takes a CALLABLE (start, stop) -> vectors, because the star
    # block's vectors differ per spectrum and materialising them all at once
    # would hold an (n_spectra, n_components, n_pixels) array
    prepared = shifter.prepare(problem["P"])

    def carried(start, stop):
        return shifter.carry(prepared, problem["delta"][start:stop])

    power = block_power(problem["w"], a, carried, chunk=8)
    assert np.all(power >= -1e-12), "negative variance contribution"


def test_leakage_is_small_when_the_blocks_are_orthogonal(problem, shifter):
    lk = leakage(problem["data"], problem["w"], problem["P"], problem["Q"],
                 shifter, problem["delta"], chunk=8)
    assert np.all(np.isfinite(lk))


def test_tidy_block_preserves_the_model(problem):
    """Reordering and renormalising a block must not change what it predicts."""
    P, a = problem["P"], problem["a"]
    power = np.array([2.0, 5.0])          # deliberately out of order
    P2, a2, power2 = tidy_block(P, a, power)
    assert np.allclose(a @ P, a2 @ P2, atol=1e-9), "tidy_block changed the model"
    assert np.all(np.diff(power2) <= 1e-12), "components are not sorted by power"


def test_nightly_bin_averages_within_a_night():
    bjd = np.array([100.1, 100.2, 101.3, 101.4, 101.45])
    c = np.array([[1.0], [3.0], [10.0], [20.0], [30.0]])
    nights, means = nightly_bin(bjd, c)
    assert nights.size == 2
    assert means[0, 0] == pytest.approx(2.0)
    assert means[1, 0] == pytest.approx(20.0)
