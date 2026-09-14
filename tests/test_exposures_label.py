"""A nightly-stacked fit must not call its nights exposures.

GL699 on SPIRou has 1976 spectra over 458 nights. The cube was stacked, so the
fit has 916 rows, one per night and parity, and the correlation matrix was
titled "458 exposures": a reader counts 1976 files on disk and concludes the
figure dropped three quarters of the campaign. The count was right, the noun
was wrong, and the noun is the whole message.
"""

import numpy as np
from astropy.table import Table

from pca2d.twoframe import _per_row, exposures_label


def _names(n_rows):
    """One file name per exposure, twice over: even and odd order parities."""
    return np.array(["f%03d.fits" % (i // 2) for i in range(n_rows)])


def test_single_exposures_are_called_exposures():
    names = _names(8)
    assert exposures_label(names, np.ones(8, bool)) == "4 exposures"
    # an all-ones n_exposures is a cube that was NOT stacked, and must read the
    # same as no column at all
    assert exposures_label(names, np.ones(8, bool), np.ones(8)) == "4 exposures"


def test_stacked_nights_say_nights_and_how_many_spectra():
    # four nights of 5, 3, 6 and 2 exposures: 16 spectra, 8 rows
    per_night = [5, 3, 6, 2]
    names = _names(8)
    per_row = np.repeat(per_night, 2).astype(float)
    label = exposures_label(names, np.ones(8, bool), per_row)
    assert label == "4 nights of 16 spectra", label


def test_a_night_is_counted_once_not_once_per_parity():
    """Both rows of a night carry that night's count; summing them doubles it."""
    names = _names(4)
    per_row = np.array([7.0, 7.0, 3.0, 3.0])
    assert exposures_label(names, np.ones(4, bool), per_row) == "2 nights of 10 spectra"


def test_the_keep_mask_selects_which_nights_are_counted():
    names = _names(6)
    per_row = np.array([4.0, 4.0, 9.0, 9.0, 2.0, 2.0])
    keep = np.array([True, True, False, False, True, True])
    assert exposures_label(names, keep, per_row) == "2 nights of 6 spectra"


def test_one_parity_kept_still_counts_the_whole_night():
    """A window reached by one parity only keeps one row of each night."""
    names = _names(6)
    per_row = np.array([4.0, 4.0, 9.0, 9.0, 2.0, 2.0])
    keep = np.array([True, False, True, False, True, False])
    assert exposures_label(names, keep, per_row) == "3 nights of 15 spectra"


def test_without_names_only_rows_can_be_counted():
    assert exposures_label(None, np.ones(8, bool)) == "8 rows"


def test_per_row_reads_the_column_and_ignores_an_unstacked_one():
    stacked = Table({"filename": ["a", "a"], "n_exposures": [3, 3]})
    assert _per_row(stacked) is not None
    flat = Table({"filename": ["a", "a"], "n_exposures": [1, 1]})
    assert _per_row(flat) is None, "all ones means the rows ARE exposures"
    old = Table({"filename": ["a", "a"]})
    assert _per_row(old) is None, "a cube from before stacking has no column"


def test_per_row_accepts_an_npz_archive():
    """fit.npz exposes `files`, not `colnames`; both have to work."""

    class Archive(dict):
        @property
        def files(self):
            return list(self)

    assert _per_row(Archive(n_exposures=np.array([4.0, 4.0]))) is not None
    assert _per_row(Archive(bjd=np.zeros(2))) is None
