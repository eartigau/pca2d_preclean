"""A template made from a smoothed star is smoothed the same way, and says so."""
import numpy as np

import pca2d.lbltemplate as lt


def test_the_template_reads_how_the_fit_smoothed_its_star():
    assert lt.star_fwhm({"star_fwhm": np.array(5)}) == 5
    assert lt.star_fwhm({"star_fwhm": np.array(0)}) == 0
    assert lt.star_fwhm({}) == 0, "a fit from before the key is unsmoothed"


def test_a_smoothed_fit_gets_a_stamp_of_its_own(tmp_path):
    plain, smoothed = tmp_path / "plain.npz", tmp_path / "smoothed.npz"
    np.savez(plain, star_fwhm=0)
    np.savez(smoothed, star_fwhm=5)
    assert not lt.template_stamp(str(plain)).endswith(" s5")
    assert lt.template_stamp(str(smoothed)).endswith(" s5")
