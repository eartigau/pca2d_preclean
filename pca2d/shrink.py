"""Each observer component kept at a column only as far as the data detect it there.

The correction divides b_n . Q out of every file, and Q was estimated from the
same noisy spectra it corrects: its error at a column is common to every
exposure, a fixed pattern in the observer's frame that lands somewhere else in
the star's at every BERV. Where a column holds nothing the data can detect,
the correction adds that pattern and removes nothing. So each component is
shrunk at each column by the positive-part James-Stein factor of its own
significance there:

    F_i = sum_n w_ni b_n b_n^T        the components' Fisher matrix at column i
    sigma(Q_ji)^2 = [F_i^-1]_jj        every row weighted by its own noise
    s_ji = max(0, 1 - 1/z_ji^2),       z_ji = Q_ji / sigma(Q_ji)

zero where |z| <= 1, close to one where |z| >> 1, and continuous in between.
The factor belongs to the column and is the same for every exposure, so a
pixel's time series stays one series; each exposure still gets its own
amount of correction through its b_n. Per component rather than one factor
per column, since at one column the water can be significant and the oxygen
not. The pixels themselves keep their weight: where the correction goes to
zero, the file keeps its flux as delivered.
"""

import numpy as np


def shrink_factors(Q, b, w, chi2_scale=1.0):
    """(J, M) factors in [0, 1] for the observer components Q (J, M).

    `b` holds the rows' amplitudes (rows, J) and `w` their weights (rows, M),
    one row per cube row. `chi2_scale`, the fit's reduced chi2, scales the
    weights so the errors are the ones the residual shows rather than the
    formal ones. Columns no row weights get zero.
    """
    Q = np.asarray(Q, dtype=float)
    b = np.asarray(b, dtype=float)
    n_comp, n_cols = Q.shape
    fisher = np.zeros((n_cols, n_comp, n_comp))
    for j in range(n_comp):
        for k in range(j, n_comp):
            f = (b[:, j] * b[:, k]) @ w / float(chi2_scale)
            fisher[:, j, k] = f
            fisher[:, k, j] = f
    diag = fisher[:, np.arange(n_comp), np.arange(n_comp)]
    ok = np.all(diag > 0, axis=1)
    out = np.zeros((n_comp, n_cols))
    if not ok.any():
        return out
    sub = fisher[ok]
    # a ridge far below any real curvature, so a column where two components
    # are nearly degenerate still inverts
    sub += 1e-12 * np.trace(sub, axis1=1, axis2=2)[:, None, None] * np.eye(n_comp)
    var = np.linalg.inv(sub)[:, np.arange(n_comp), np.arange(n_comp)]
    z2 = Q[:, ok].T ** 2 / var
    out[:, ok] = np.clip(1.0 - 1.0 / np.maximum(z2, 1e-300), 0.0, 1.0).T
    return out
