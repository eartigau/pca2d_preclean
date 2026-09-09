"""Shared fixtures: one small, fixed two-frame problem the tests can lean on.

Everything here is deterministic. No file is read, no network is touched, and
the whole suite runs in seconds, because a regression net that takes a minute
gets skipped and a net that gets skipped is not a net.
"""

from __future__ import annotations

import numpy as np
import pytest

# Small enough to be instant, large enough that a shift of several pixels and a
# five-component basis are not degenerate.
N_PIX = 512
N_SPEC = 24
N_STAR = 2
N_EARTH = 2


@pytest.fixture(scope="session")
def rng():
    return np.random.default_rng(20260901)


@pytest.fixture(scope="session")
def problem():
    """A synthetic star-plus-Earth cube with a KNOWN answer.

    Built forwards from a basis and coefficients, so a solver can be checked
    against the truth rather than against its own previous output. The star
    block is shifted by a per-spectrum offset and the Earth block is not, which
    is the whole geometry of the method in miniature.
    """
    r = np.random.default_rng(7)
    grid = np.arange(N_PIX, dtype=float)

    def smooth(n):
        """Band-limited rows, so the shift operator has something it can move."""
        f = r.normal(size=(n, N_PIX))
        k = np.exp(-0.5 * (np.arange(N_PIX) - N_PIX / 2) ** 2 / 6.0 ** 2)
        return np.array([np.convolve(row, k, mode="same") for row in f])

    P = smooth(N_STAR)
    Q = smooth(N_EARTH)
    a = r.normal(size=(N_SPEC, N_STAR))
    b = r.normal(size=(N_SPEC, N_EARTH))
    # integer offsets keep the truth exact under any correct shifter
    delta = r.integers(-4, 5, size=N_SPEC).astype(float)

    data = np.zeros((N_SPEC, N_PIX))
    for i in range(N_SPEC):
        for k in range(N_STAR):
            data[i] += a[i, k] * np.roll(P[k], int(delta[i]))
        data[i] += b[i] @ Q

    w = np.ones_like(data)
    # The truth above uses a circular np.roll; a local interpolating operator
    # cannot reproduce that within its own half-width of an edge, and neither
    # can any real spectrograph. Zero the margin so the comparison is made
    # where both descriptions agree, which is also where the data lives.
    edge = 24
    w[:, :edge] = 0.0
    w[:, -edge:] = 0.0
    # a realistic weight map has holes; put some in, in different places per row
    for i in range(0, N_SPEC, 3):
        start = edge + (i * 7) % (N_PIX - 2 * edge - 20)
        w[i, start:start + 20] = 0.0
    return dict(grid=grid, data=data, w=w, P=P, Q=Q, a=a, b=b, delta=delta)
