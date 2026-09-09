"""Two-frame weighted PCA of NIRPS spectra.

M components anchored in the stellar rest frame and N in the observer frame,
fitted together, so the fit decides which frame a spectral feature belongs to
rather than the analyst deciding by choosing a registration.

    build     t.fits -> a cached, registered cube          (pca2refs-cube)
    twoframe  the M+N fit                                  (pca2refs-fit)
    reconstruct  rebuild or correct one exposure           (pca2refs-apply)
    wizard    validate a configuration against the data, and fill it in

Input is APERO `t.fits`, order by order. The single-frame pipeline this grew
out of is in `_obsolete/`. Documentation: https://eartigau.github.io/pca2refs/
"""
__version__ = "1.0.0"
