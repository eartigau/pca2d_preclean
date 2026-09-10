"""Every Doppler shift is relativistic: sqrt((1 + v/c) / (1 - v/c)).

The first-order 1 + v/c is off by (v/c)**2 / 2, 1.5 m/s at a 30 km/s BERV,
the size of the signals the corrected spectra are measured for.
"""

import pathlib
import re

import numpy as np

from pca2d.grids import (C_KMS, doppler, log_shift, pixel_shift,
                         shift_velocity, velocity)


def test_there_and_back_is_the_identity():
    v = np.array([-31.0, -2.5, 0.0, 0.001, 18.8, 31.0])
    assert np.allclose(doppler(v) * doppler(-v), 1.0, rtol=0, atol=1e-15)
    assert np.allclose(velocity(log_shift(v)), v, rtol=1e-12, atol=1e-12)
    assert np.allclose(shift_velocity(pixel_shift(v, 0.5), 0.5), v,
                       rtol=1e-12, atol=1e-12)


def test_it_is_the_log_of_the_ratio_and_composes_relativistically():
    a, b = 25.0, -12.0
    assert np.isclose(np.log(doppler(a)), log_shift(a), rtol=1e-12)
    both = (a + b) / (1 + a * b / C_KMS ** 2)
    assert np.isclose(doppler(a) * doppler(b), doppler(both), rtol=0, atol=1e-15)


def test_the_first_order_form_is_off_by_one_and_a_half_metres_per_second():
    v = 30.0
    error_ms = (doppler(v) - (1 + v / C_KMS)) * C_KMS * 1000.0
    assert np.isclose(error_ms, 0.5 * (v / C_KMS) ** 2 * C_KMS * 1000.0, rtol=1e-3)
    assert 1.4 < error_ms < 1.6


def test_the_pixel_shift_is_atanh_over_the_log_step():
    dv = 0.5
    assert np.isclose(pixel_shift(30.0, dv),
                      np.arctanh(30.0 / C_KMS) / (dv / C_KMS), rtol=1e-12)
    # on a log grid, the plain v/dv was already this to a millionth of a pixel;
    # it is the wavelength factor of the registration that was not
    assert abs(pixel_shift(30.0, dv) - 30.0 / dv) < 1e-6


def test_no_module_uses_the_first_order_form():
    root = pathlib.Path(__file__).resolve().parents[1] / "pca2d"
    first_order = re.compile(r"1(?:\.0)?\s*[+-]\s*[\w\[\]\"'.()]+\s*/\s*C_KMS")
    found = ["%s:%d" % (p.relative_to(root), n) for p in root.rglob("*.py")
             for n, line in enumerate(p.read_text().splitlines(), 1)
             if first_order.search(line)]
    assert not found, found
