"""One velocity per exposure, fitted so the observer block cannot have it.

The star is never exactly where the BERV alone would put it: it has a planet,
it has activity, and the carry operator leaves a residual of its own. That
residual has one shape, the derivative of the star, and an observer block with
several free vectors and an amplitude per exposure will describe it. Dividing
that block back out then writes a velocity into the corrected spectrum, which
is measurable: on TOI2120 it explained 41% of the variance of what the
correction did to LBL's velocities.

These tests build a problem where the answer is known: a star spectrum given a
shift of a few hundred m/s per exposure, on top of an observer-frame pattern
that has nothing to do with it.
"""

import numpy as np
import pytest

from pca2d.twoframe import (LanczosShifter, joint_coeffs, star_model,
                            velocity_column)

DV = 0.5                      # km/s per sample, as the magic grid is built
N_PIX, N_SPEC = 2048, 40


@pytest.fixture(scope="module")
def shifted():
    """A star with a per-exposure velocity, plus an unrelated observer signal."""
    r = np.random.default_rng(11)
    x = np.arange(N_PIX, dtype=float)

    # a star: a few absorption lines wide enough for a derivative to mean
    # something at this sampling
    star = np.zeros(N_PIX)
    for centre in (300, 700, 1150, 1600, 1900):
        star -= np.exp(-0.5 * ((x - centre) / 6.0) ** 2)
    star -= star.mean()

    # an observer-frame pattern, fixed in the observer frame by construction
    earth = np.zeros(N_PIX)
    for centre in (450, 980, 1420):
        earth -= np.exp(-0.5 * ((x - centre) / 4.0) ** 2)
    earth -= earth.mean()
    earth /= np.linalg.norm(earth)

    berv = np.linspace(-15.0, 15.0, N_SPEC)          # km/s
    truth_v = 300.0 * np.sin(np.arange(N_SPEC) * 0.7)   # m/s, the star's own
    depth = 1.0 + 0.3 * r.normal(size=N_SPEC)          # observer amplitude

    shifter = LanczosShifter(N_PIX, a=8, threads=1)
    P = (star / np.linalg.norm(star))[None, :]
    delta = -berv / DV                                # pixels, as the fit uses
    a = np.full((N_SPEC, 1), float(np.linalg.norm(star)))

    # the star at its BERV, then displaced by its own velocity: first order,
    # which is exactly what the fitted column claims to describe
    base = star_model(P, a, shifter, delta, N_SPEC, N_PIX, chunk=8)
    alpha_true = truth_v / 1000.0 / DV                # m/s -> pixels
    data = base + alpha_true[:, None] * np.gradient(base, axis=1)
    data += depth[:, None] * earth[None, :]
    w = np.ones_like(data)
    return dict(data=data, w=w, P=P, Q=earth[None, :], shifter=shifter,
                delta=delta, a=a, truth_v=truth_v, depth=depth, base=base)


def test_the_fitted_shift_is_the_velocity_that_was_put_in(shifted):
    a, b, alpha, _ = joint_coeffs(
        shifted["data"], shifted["w"], shifted["P"], shifted["Q"],
        shifted["shifter"], shifted["delta"], chunk=8,
        velocity_from=shifted["a"])
    got = alpha * DV * 1000.0                          # pixels -> m/s
    err = got - shifted["truth_v"]
    assert np.std(err) < 15.0, (
        "recovered velocity is off by %.1f m/s rms" % np.std(err))
    assert abs(np.corrcoef(got, shifted["truth_v"])[0, 1] - 1) < 1e-3


def test_without_the_column_the_observer_block_takes_the_velocity(shifted):
    """The failure this exists to prevent, made to happen on purpose."""
    plain = joint_coeffs(shifted["data"], shifted["w"], shifted["P"],
                         shifted["Q"], shifted["shifter"], shifted["delta"],
                         chunk=8)
    with_v = joint_coeffs(shifted["data"], shifted["w"], shifted["P"],
                          shifted["Q"], shifted["shifter"], shifted["delta"],
                          chunk=8, velocity_from=shifted["a"])
    truth = shifted["depth"]
    b_plain = plain[1][:, 0]
    b_vel = with_v[1][:, 0]
    err_plain = np.std(b_plain - truth)
    err_vel = np.std(b_vel - truth)
    assert err_vel < err_plain, (
        "the observer amplitude is no better with the velocity term: %.4f"
        " against %.4f" % (err_vel, err_plain))


def test_the_column_is_the_derivative_of_the_model_not_of_the_data(shifted):
    """A derivative of one noisy spectrum is a derivative of its noise."""
    r = np.random.default_rng(3)
    clean = shifted["base"][0]
    noisy = clean + 0.05 * r.normal(size=clean.size)

    from_model = velocity_column(clean)
    from_data = velocity_column(noisy)

    # measured where the star has no line, so the true derivative is zero and
    # whatever is left is what the column would contribute as pure noise
    flat = np.abs(clean - np.median(clean)) < 0.01 * np.ptp(clean)
    assert flat.sum() > 200, "the test spectrum has no line-free region"
    quiet_model = np.std(from_model[flat])
    quiet_data = np.std(from_data[flat])
    signal = np.std(from_model[~flat])

    assert quiet_data > 10 * quiet_model, (
        "a derivative of one noisy spectrum should be far noisier than the"
        " model's: %.3g against %.3g" % (quiet_data, quiet_model))
    assert quiet_data > 0.2 * signal, (
        "and comparable to the signal it is supposed to carry: %.3g against"
        " %.3g" % (quiet_data, signal))


def test_a_shift_of_zero_is_found_to_be_zero(shifted):
    """No velocity in, no velocity out: the term must not invent one."""
    data = shifted["base"] + shifted["depth"][:, None] * shifted["Q"]
    _, _, alpha, _ = joint_coeffs(data, shifted["w"], shifted["P"],
                                  shifted["Q"], shifted["shifter"],
                                  shifted["delta"], chunk=8,
                                  velocity_from=shifted["a"])
    assert np.max(np.abs(alpha * DV * 1000.0)) < 5.0


# ------------------------------------------------- the scaling of the solve --
def test_equilibrating_changes_the_conditioning_and_not_the_answer():
    """Both halves of the claim, since only one of them is obvious.

    The velocity column arrives with a norm an order of magnitude away from the
    orthonormal ones, which costs two decades of condition number without any
    two columns being collinear. Scaling it is exact arithmetic, so the answer
    must not move.
    """
    from pca2d.twoframe import equilibrated_solve

    r = np.random.default_rng(5)
    n = 6
    M = r.normal(size=(400, n))
    M[:, -1] *= 300.0                      # one column on a different scale
    amat = M.T @ M
    truth = r.normal(size=n)
    bvec = amat @ truth

    got, cond_scaled = equilibrated_solve(amat, bvec)
    plain = np.linalg.solve(amat + 1e-12 * np.trace(amat) * np.eye(n), bvec)

    assert np.allclose(got, truth, rtol=1e-8, atol=1e-8), "the answer moved"
    assert np.allclose(got, plain, rtol=1e-6, atol=1e-6), (
        "equilibrating is a change of scale, not a change of problem")
    assert cond_scaled < np.linalg.cond(amat) / 100, (
        "conditioning barely improved: %.3g against %.3g"
        % (cond_scaled, np.linalg.cond(amat)))


def test_an_unconstrained_column_still_solves(shifted):
    """A component whose support is fully masked has a zero diagonal."""
    from pca2d.twoframe import equilibrated_solve

    amat = np.diag([4.0, 0.0, 9.0])
    bvec = np.array([8.0, 0.0, 27.0])
    got, cond = equilibrated_solve(amat, bvec)
    assert np.isfinite(got).all()
    assert np.allclose(got[[0, 2]], [2.0, 3.0])
    assert got[1] == 0.0, "an unconstrained column must come back as zero"


def test_the_error_matrix_keeps_the_shape_its_callers_slice(shifted):
    """The velocity sigma must not ride along as an extra column.

    It did once. Every caller slices [:, n_star:] for the observer block, so
    the extra column became an eighth error bar on seven components and the
    plotting died on a broadcast, sixteen minutes of fit later.
    """
    from pca2d.twoframe import coefficient_errors

    a, b, alpha, _ = joint_coeffs(
        shifted["data"], shifted["w"], shifted["P"], shifted["Q"],
        shifted["shifter"], shifted["delta"], chunk=8,
        velocity_from=shifted["a"])
    formal, scaled, chi2_red, sig_vel = coefficient_errors(
        shifted["data"], shifted["w"], shifted["P"], shifted["Q"],
        shifted["shifter"], shifted["delta"], a, b, chunk=8, alpha=alpha)

    n_star, n_earth = shifted["P"].shape[0], shifted["Q"].shape[0]
    assert formal.shape[1] == scaled.shape[1] == n_star + n_earth
    assert scaled[:, n_star:].shape[1] == n_earth, "an eighth error bar is back"
    assert sig_vel is not None and sig_vel.shape == (shifted["data"].shape[0],)
    assert np.all(sig_vel[np.isfinite(sig_vel)] >= 0)


def test_without_the_term_nothing_extra_comes_back(shifted):
    from pca2d.twoframe import coefficient_errors

    a, b, _, _ = joint_coeffs(shifted["data"], shifted["w"], shifted["P"],
                              shifted["Q"], shifted["shifter"],
                              shifted["delta"], chunk=8)
    formal, scaled, _, sig_vel = coefficient_errors(
        shifted["data"], shifted["w"], shifted["P"], shifted["Q"],
        shifted["shifter"], shifted["delta"], a, b, chunk=8)
    assert sig_vel is None
    assert formal.shape[1] == shifted["P"].shape[0] + shifted["Q"].shape[0]


# ------------------------------------------- where the shift may be measured --
def test_the_column_is_silent_outside_the_mask(shifted):
    from pca2d.twoframe import velocity_column

    row = shifted["base"][0]
    mask = np.zeros(row.size)
    mask[500:1500] = 1.0
    column = velocity_column(row, mask)
    assert np.all(column[:500] == 0) and np.all(column[1500:] == 0)
    assert np.any(column[500:1500] != 0)


def test_a_shift_is_still_recovered_from_a_masked_column(shifted):
    """Restricting where it is measured must not bias what it measures."""
    from pca2d.twoframe import joint_coeffs

    mask = np.zeros(N_PIX)
    mask[200:1800] = 1.0          # drops two of the five lines' wings
    _, _, alpha, _ = joint_coeffs(
        shifted["data"], shifted["w"], shifted["P"], shifted["Q"],
        shifted["shifter"], shifted["delta"], chunk=8,
        velocity_from=shifted["a"], velocity_mask=mask)
    got = alpha * DV * 1000.0
    assert np.std(got - shifted["truth_v"]) < 25.0, (
        "masked recovery is off by %.1f m/s rms" % np.std(got - shifted["truth_v"]))


def test_clean_columns_reads_the_cube_and_asks_for_the_bad_nights(tmp_path):
    """A column clean on average and absorbed on bad nights is not clean."""
    from pca2d.twoframe import clean_columns

    n_rows, n_cols = 20, 6
    trans = np.ones((n_rows, n_cols))
    trans[:, 1] = 0.2                       # always absorbed
    trans[:, 2] = 0.99                      # always clean
    trans[:, 3] = 0.99                      # clean except on the worst nights
    trans[:3, 3] = 0.4
    trans[:, 4] = 0.90                      # never quite clean enough
    np.save(tmp_path / "trans.npy", trans)

    mask = clean_columns(str(tmp_path), threshold=0.95)
    assert mask is not None
    assert mask[0] and mask[2], "a clean column was rejected"
    assert not mask[1] and not mask[4], "an absorbed column was kept"
    assert not mask[3], (
        "a column absorbed on 15% of the nights passed; the cut is meant to ask"
        " about the bad nights, not the median one")


def test_no_transmission_in_the_cube_is_not_an_error(tmp_path):
    from pca2d.twoframe import clean_columns

    assert clean_columns(str(tmp_path)) is None
    assert clean_columns(str(tmp_path / "does_not_exist")) is None
