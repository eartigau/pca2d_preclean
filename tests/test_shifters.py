"""The shift operators, which everything else in the fit rests on.

The star basis is carried into each spectrum's frame by these, so an error here
is invisible in the residuals and fatal in the coefficients. Two properties are
worth more than any pinned number and are checked first:

* shifting by zero changes nothing;
* the adjoint really is the adjoint, <S x, y> == <x, S* y>, which is what makes
  the normal equations correct rather than merely plausible.
"""

from __future__ import annotations

import numpy as np
import pytest

from pca2d.twoframe import FourierShifter, LanczosShifter


@pytest.fixture(params=["lanczos", "fft"])
def shifter(request, problem):
    n = problem["data"].shape[1]
    return (LanczosShifter(n, threads=1) if request.param == "lanczos"
            else FourierShifter(n))


def test_zero_shift_is_identity(shifter, problem):
    x = problem["P"]
    out = shifter.rows(x, np.zeros(x.shape[0]))
    assert np.allclose(out, x, atol=1e-9), "a shift of zero moved the data"


def test_adjoint_is_the_adjoint(shifter, problem, rng):
    """<S x, y> == <x, S* y>, to machine precision.

    This is the property the block solve depends on: the star update solves
    normal equations built from S and S*, and if they are not a genuine
    adjoint pair the solution is the answer to a different problem.
    """
    n = problem["data"].shape[1]
    pix = np.array([3.0, -2.0, 0.0, 7.0])
    x = rng.normal(size=(pix.size, n))
    y = rng.normal(size=(pix.size, n))
    lhs = float(np.sum(shifter.rows(x, pix) * y))
    rhs = float(np.sum(x * shifter.adjoint(y, pix)))
    assert lhs == pytest.approx(rhs, rel=1e-8, abs=1e-8)


def test_integer_shift_matches_roll(problem):
    """On a band-limited row an integer shift is exactly numpy.roll.

    Only asserted for Lanczos: the FFT shifter wraps a phase ramp around the
    whole array and is checked by the adjoint and identity tests instead.
    """
    n = problem["data"].shape[1]
    s = LanczosShifter(n, threads=1)
    x = problem["P"]
    for d in (-3.0, 1.0, 5.0):
        got = s.rows(x, np.full(x.shape[0], d))
        want = np.array([np.roll(row, int(d)) for row in x])
        edge = 16
        assert np.allclose(got[:, edge:-edge], want[:, edge:-edge], atol=1e-6), \
            "integer shift by %+d does not match numpy.roll" % d


def test_carry_matches_rows(problem):
    """carry() is the optimised path; it must agree with rows() exactly."""
    n = problem["data"].shape[1]
    s = LanczosShifter(n, threads=1)
    pix = problem["delta"][:6]
    prepared = s.prepare(problem["P"])
    carried = s.carry(prepared, pix)
    for i, d in enumerate(pix):
        want = s.rows(problem["P"], np.full(problem["P"].shape[0], d))
        assert np.allclose(carried[i], want, atol=1e-9), \
            "carry and rows disagree at row %d" % i


def test_threads_do_not_change_the_answer(problem):
    """The thread count is a speed knob and must not be a physics knob."""
    n = problem["data"].shape[1]
    prepared_1 = LanczosShifter(n, threads=1)
    prepared_4 = LanczosShifter(n, threads=4)
    pix = problem["delta"][:8]
    one = prepared_1.carry(prepared_1.prepare(problem["P"]), pix)
    four = prepared_4.carry(prepared_4.prepare(problem["P"]), pix)
    assert np.array_equal(one, four), "the answer depends on the thread count"


def test_diag_normal_is_positive(shifter, problem):
    """The diagonal of S* W S cannot be negative for non-negative weights."""
    d = shifter.diag_normal(problem["w"], problem["delta"])
    assert np.all(d >= -1e-12), "negative diagonal in the normal matrix"
