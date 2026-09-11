# Variants

A variant is the nominal configuration, `config.yaml`, plus what it changes,
written in one file here: `variants/<name>.yaml`. It runs with

    pca2d-preclean --object TOI2120 --variant <name>

Its products go to `outputs/_<name>/<object>/<M-N>/` and its LBL object is
`<object>_PCA2D_<M-N>_<name>`, unless the file sets `output.directory` or
`lbl.suffix` itself. A file with `reuse_fit: <other>` changes the correction
only: it takes the fit of `<other>` (`reuse_fit: nominal` for the nominal's),
and the cube, fit and figures stages are not run. `note` records the result.
Neither of those two keys configures anything.

Every run writes `resolved_config.yaml` beside its outputs, with the commit of
the code that ran it under `provenance`. The corrected spectra carry that
commit as `PCA2GIT`, the LBL template too, and `fit.npz` as `pca2d_code`.

## TOI-2120 series of 2026-09-11

Every run of the series used cube a92a1b985cdd, with the high pass of that
day, 151 samples (75 km/s), and 6 sweeps keeping the last. Each file repeats
those settings, so it reproduces its run on its own. LBL on 316 exposures,
in m/s:

| variant | what it changes | rms | robust | nightly |
|---|---|---:|---:|---:|
| (delivered) | no correction | 47.71 | 34.99 | |
| s5 | per-parity offset (mean: offset), grid star | 31.74 | 23.73 | |
| star | per-parity star spectrum, grid star | 20.36 | 16.57 | 18.92 |
| star_shr | star, correction shrunk (correction only) | 20.35 | 16.16 | 18.90 |
| star_sg | star smoothed to 1 resolution element, grid | 24.14 | 19.38 | 22.61 |
| star_sg05 | star smoothed to 1/2 element, grid | 22.34 | 17.76 | |
| star_spl | spline star | 20.99 | 17.28 | 19.55 |
| spl_shr | spline, correction shrunk (correction only) | 21.05 | 16.16 | 19.60 |
| spl_shrsig | spl_shr, significance smoothed | 21.17 | 16.39 | 19.73 |
| spl_shrQ23 | spl_shr, observer components 2 and 3 smoothed | 21.00 | 16.61 | 19.54 |
| spl_ss025 | spline star smoothed to 1/4 element, shrunk | 21.48 | 15.84 | 19.98 |
| spl_ss025s | spl_ss025, LBL template smoothed like the star | pending | | |
| spl_ss05 | spline star smoothed to 1/2 element, shrunk | 21.91 | 16.45 | 20.35 |

- **Templates built before the smoothing reached them:** star_sg, star_sg05 and spl_ss025 were measured against LBL templates built before the star smoothing propagated to the template.
- **Correction-only variants:** star_shr and the spl_*shr* variants were run by a scratch script against a copy of their base's LBL template and masks. With `--variant`, the lbl stage builds the same template from the same fit.
- **spl_ss025s:** it wrote into spl_ss025's folder, whose `resolved_config.yaml` now carries its LBL suffix.
- **spl_shrQ23:** it smooths observer components, which are not to be smoothed; it is kept as a record only.
