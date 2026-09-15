# pca2d-preclean

Two-frame weighted PCA of echelle spectra, and the front end of an LBL run.

One basis travels with the star and one stands still on the detector. Both are
fitted to every exposure at once, the one that stands still is taken out, and
what is left goes to [LBL](https://github.com/njcuk9999/lbl) together with a
template of the star made by the same fit.

```
pca2d-preclean --object TOI-2120
```

runs all of it, from the t.fits to LBL's velocities.

## What it buys

TOI-2120, 316 SPIRou exposures over 452 days, one star and three observer
components, and LBL on the same exposures before and after:

| spectra | rms | robust sigma | nightly rms | median error |
| --- | ---: | ---: | ---: | ---: |
| as delivered | 47.7 m/s | 35.0 m/s | 46.7 m/s | 7.29 m/s |
| corrected, the star fixed (0-3) | 17.4 m/s | 14.1 m/s | 15.4 m/s | 7.93 m/s |

The robust sigma is 1.4826 times the MAD, and the nightly rms is that of the
weighted nightly means. `docs/make_figures.py` draws the velocities from the
rdb files and prints these numbers.

It is not free everywhere, and what decides it is measured
(`paper/pca2d.pdf`). On four campaigns, against the delivered spectra on the
exposures both series share:

| target | instrument | delivered | corrected |
| --- | --- | ---: | ---: |
| TOI-2120 | SPIRou | 47.7 m/s | 17.4 m/s |
| GJ 1 | NIRPS | 2.70 m/s | 2.61 m/s |
| Proxima | NIRPS | 2.86 m/s | 3.04 m/s |
| TOI-4552 | NIRPS | 12.6 m/s | 14.2 m/s |

The gain is large where telluric residuals dominate and absent where the
delivered velocities are already good. What separates them is the barycentric
span of the exposures a fit is built on: two 14-night slices of the TOI-4552
campaign, at the same signal-to-noise, spanning 42.8 and 0.7 km/s of BERV, go
from 8.5 to 4.5 m/s of robust scatter and from 9.1 to 21.8 respectively. Where
the star does not move against the sky the two frames are one frame, the
observer block takes up stellar structure, and dividing it out removes part of
the star. Fitting several stars against one observer basis is the way out
(`docs/joint_fit.md`).

## The model

Each exposure becomes the log of its flux less a Savitzky-Golay of it, 100
km/s wide (`highpass.width_kms`, 201 samples of 0.5 km/s; runs saved before
2026-09-11 used 151 samples, 75 km/s), on one log-uniform grid in the
observer's frame, as two rows: its even orders and its odd ones. Then

```
y_n  =  S_n P a_n  +  v_n d(S_n P a_n)/dv  +  Q b_n  +  m_p
```

* `S_n` carries the star basis `P` into exposure `n` by its barycentric
  velocity. On a log grid a Doppler shift is a translation, of
  atanh(v/c) / (dv/c) samples (relativistic, always), and `S_n` is an exact
  Lanczos operator for it rather than an interpolation that would smear what it
  touched.
* `P a_n` is the star. `Q b_n` is what stands still on the detector: telluric
  absorption, OH emission, the instrument.
* `v_n` is one velocity per exposure, times the derivative of the star as
  reconstructed. Without it the observer block learns the star's own motion,
  and the correction divides that out of the flux. With it the motion is
  fitted, and stays in.
* `m_p` is one offset per order parity, in the observer's frame. Consecutive
  orders overlap, and at a given wavelength one parity samples near an order's
  centre and the other near its edge, where the resolution is not the same.
  It is taken out once before the fit and put back into the correction.

The coefficients of both blocks are solved jointly for each exposure and the
bases one at a time, by block coordinate descent, until chi2 turns over. The
two rows of an exposure share its coefficients: they are one measurement.

The weights come from the noise that can be measured, not from a photon model.
Past about 2200 nm on SPIRou the thermal background dominates. Those photons
were counted, so they arrive with their own Poisson noise, and the pipeline
subtracts their mean and not their variance: the noise stays after the signal
it belonged to is gone. A photon sigma computed from the flux that is left
describes a spectrum nobody recorded, and calls the noisiest part of the array
the quietest, by a factor of fifty. The sample-to-sample scatter is measured
from what is actually there, and the weights take whichever of the two is
larger.

`--mean iterate` gives each frame its own mean per parity instead,
re-estimated at every sweep. It converges on synthetic data and on a 70 nm
slice of TOI-2120, and did not on the whole cube, so `offset` is the default.

## What a corrected file holds

```
f_corrected  =  f * exp(-(Q b_n + m_p))
```

on the file's own pixels, and NaN wherever the fit gave the sample no weight or
the sky outshines the star. That is panel 3 of the report's sequence figure,
sample for sample: the two differ only by the continuum, which the file keeps
and the panel's high pass takes off. `tests/test_panel3_is_the_correction.py`
keeps them together.

The star block and the velocity term stay in the flux; they are what LBL is
about to measure. Everything else in the t.fits, the wavelength solution and
the blaze included, is copied through, so LBL reads a corrected file with the
same class as the original. A corrected file is named for what was divided out
of it: `2811170t_0-3.fits` is no star component and three observer ones.

The primary header says what was done. `PCA2xxxx` belongs to the run as a
whole; `PCASTR` and `PCAOBS` are the two frames, with `_N` the number of
amplitude cards that follow and `_D` how many of those components were divided
out:

```
PCA2SKYR=                  4.0 / sky/flux above this was set to NaN
PCA2SKYN=                    3 / samples removed as sky-dominated
PCA2REF =                    T / two-frame PCA correction applied
PCA2BERV=   2.6686512682488392 / km/s used to carry the star basis
PCA2NPIX=               180801 / samples corrected
PCA2MEAN=                    T / observer-frame parity mean divided out
PCA2WNAN=                 6557 / samples the fit gave no weight, set to NaN
PCA2REJ =                    F / exposure was MAD-rejected
PCASTR_V=    77.03979913296997 / m/s, star shift the fit absorbed
PCASTR_N=                    1 / star-frame comps listed, PCASTR01..01
PCASTR_D=                    0 / star-frame comps divided out of the flux
PCASTR01=   44.223559850801585 / star-frame comp 1 amplitude, left in the flux
PCAOBS_N=                    3 / observer-frame comps listed, PCAOBS01..03
PCAOBS_D=                    3 / observer-frame comps divided out of the flux
PCAOBS01=  -18.060037548051366 / observer-frame comp 1 amplitude, divided out
```

Every amplitude the fit has is written, not only the ones divided out: an rdb
knows only what LBL measured, so the header is the one place a component can be
lined up with its exposure. Two digits, since `PCASTR001` would be nine
characters and a HIERARCH card.

## LBL, built in

The `lbl` stage puts two objects side by side in one LBL tree, from the same
instrument profile:

```
lbl/science/TOI-2120/               symlinks to the spectra as delivered
lbl/science/TOI-2120_PCA2D_2-3v/    symlinks to what this run corrected
```

The name carries the run's tag because LBL takes a whole science folder, and two
corrections landing in one would be measured as one series without a word.

**The corrected object's template is the fit's star.** The first star
component at its mean amplitude is a high-passed template of the star, fitted
to every exposure at once in the barycentric frame with the atmosphere already
described by the other block. Each order parity gets its own: the star block
plus that parity's weighted mean of what the fit left, since LBL measures an
order against the template of its parity and the two parities see a line at
different resolutions. The stage writes it in LBL's template format, by
LBL's own writer (`pca2d/lbltemplate.py`), as `star_template.fits` beside the
run's outputs, and copies it to where LBL looks for that object's template, so
LBL's template step finds it there and has nothing to do. A template LBL built
itself is never replaced. Through the same 75 km/s high pass, it and LBL's own
template for the same spectra agree with a correlation of 0.93 and a slope of
0.94 (0.99 and 0.99 in H).

**Star components past the first are variability indicators.** LBL projects
every line's residual on RESPROJ tables; its DTEMP tables are temperature
gradients of model spectra. With two or more star components the stage writes
the others as the same kind of table, `STRPCA2` to `STRPCAn`, measured on the
star itself: `fractional_gradient` is the template's flux times the component,
so LBL's `STRPCA2` for an exposure is that exposure's `a2` less its mean, measured
line by line. LBL evaluates these tables in the star's rest frame, which its
mask step measures, so the corrected object runs its mask first, the tables are
written, and then its velocities are measured with them; the rdb gets a
`STRPCA2` and an `sSTRPCA2` column. The LBL installed here divides the residual
in place for each table, so a second table would be projected on a residual
divided twice: `STRPCA2` is right, and the stage warns when a fit has three star
components or more.

On TOI-2120 with two star components (2-3v), LBL's `STRPCA2` follows the fit's
own `a2` exposure by exposure with a correlation of 0.969 and a slope of 0.81
(`docs/figures/strpca_2-3v.svg`): the tables carry what they should. That fit is
also the warning that goes with them. Its `a2` follows the barycentric velocity
(-0.55), the seeing (-0.58) and the sun's elevation (-0.49), which is a star
component modelling the sky, and its velocity term absorbed shifts of 232 m/s
rms that follow the barycentric velocity. The correction leaves that in the
flux, and LBL's velocities on the 2-3v spectra scatter by 67 m/s against 27 m/s
for 1-3v. A second star component is a variability indicator once the
correlation table says it belongs to the star.

Beside the run's outputs:

| file | what it is |
| --- | --- |
| `lbl_config.yaml` | LBL's configuration in LBL's own keys; it refuses any other, which is a reason to have it written rather than typed |
| `run_lbl.py` | an ordinary LBL wrap script, one runparams dict per object (`BEFORE`, `AFTER`, and `STRPCA` when there is one), to be read and run by hand |
| `star_template.fits` | the fit's star template, made again only when the fit changes |

`lbl.run: true` in `config.yaml`, or `--run-lbl`, has the stage run LBL rather
than only prepare it. Which LBL instrument a spectrograph is comes from its
block in `config.yaml`: LBL calls NIRPS `NIRPS_HA` or `NIRPS_HE` by the mode it
was observed in, and the wrong one raises no error, it returns velocities from
another instrument's profile. The effective temperature LBL needs for its mask
is read from the spectra (`OBJTEMP`), unless `lbl.teff` says otherwise.

## Installing it

```
git clone https://github.com/eartigau/pca2d_preclean.git
cd pca2d_preclean
conda env create -f environment.yml
conda activate pca2d-preclean
```

One alias, written once, so that a new terminal is one word from being ready:

```
echo "alias pca2d='source $PWD/pca2d.sh'" >> ~/.zshrc
```

`pca2d` then leaves whatever conda environments are active, however many deep,
enters this one, and prints what there is to run: the window, a solo run, a
joint run, a dry run, and the tests. Leaving rather than activating on top is
the point, since the two environments that matter here both hold a version of
LBL. Run the file instead of sourcing it and it says so, with the alias to keep.

If `conda env create` ends in `CondaEnvException: Pip failed`, with
`git version did not run successfully` and `xcrun: error: unable to load
libxcrun` above it, the environment is x86_64 under Rosetta on an Apple Silicon
Mac and pip cannot run Apple's git stub. The page has the three lines that
confirm it and the three ways out:
<https://eartigau.github.io/pca2d_preclean/#trouble>.

That environment holds **both codes**, this one and LBL. Having to deactivate
one to run the other is how a t.fits gets measured by the wrong version of
something.

It is why the versions are nailed down rather than floored. LBL pins its
dependencies exactly, `numpy==2.3.3`, `astropy==7.2.0`, `scipy==1.17.0` and the
rest, and asks for python >=3.12,<3.13; `environment.yml` repeats those pins so
that conda installs them and pip finds them already satisfied. Only `lbl`
itself, which is not on PyPI, this package in place, and one PyPI-only
dependency come through pip. `astropy-base` and `matplotlib-base` rather than
the metapackages: nothing here imports pyarrow or bqplot, and every entry point
calls `matplotlib.use("Agg")` before it draws. `pyproject.toml` keeps this
package's own looser floor, python 3.10, for anyone installing it on its own.

## A window, if the command line is one thing too many

```
pca2d-gui
```

lists the objects of the data root with their median SNR and exposure time, a
box to tick each one, runs one of them or several together, refuses two
instruments at once, explains what every choice implies when the pointer rests
on it, in English or in French, writes the command it is about to run so it can
be copied into a terminal, and shows the run's log as it comes. What it read of
a data root is remembered under `~/.pca2d/`, so a folder opens at once the
second time and only new files are read. `docs/gui.md`.

## Several stars, one observer basis

```
pca2d-preclean --objects PROXIMA,GJ1,GJ3090 --n-star 0
```

The atmosphere and the instrument belong to the night, not to the target, so
several campaigns can be fitted against a single observer basis while each star
keeps its own spectrum per order parity. No basis can follow one star when
three of them, with different barycentric coverage and systemic velocities,
constrain it.

Nights in common are NOT the criterion, and few of them are an advantage: the
observer components are a parasite that is always present, what the fit measures
is their PATTERN, and nights no other star of the group saw widen the range of
conditions that pattern is determined over. Nor is comparable brightness a
requirement, since the weights are 1/sigma^2 and absolute. The first run,
PROXIMA + GJ 1 + GJ 3090, came back worse for GJ 1 than its solo fit (3.04
against 2.61 robust), and the cause found on 2026-09-13 was the common mask,
intersected over every object so that a sample lost by one night of the faintest
star was blanked in every exposure of the others. That is fixed, and the joint
result has to be measured again before anything is concluded from it. The window
reports the nights a selection shares as information. `docs/joint_fit.md`.

## Which choices are nominal, and what each one was worth

`docs/options.md` is the inventory: every parameter that is a choice, the
measurement that settled it, what is still open, and what is a candidate for
removal once understood. It is the document to read before changing a default.

## Running it

Two roots, and a run reads from one and writes to the other. Spectra go under
the input root, one folder per target, and the only thing ever written there is
one `pca2d_index.csv` per campaign folder, the log of what has been read of that
campaign, which travels with it (`docs/gui.md`). The
instrument is read from the `INSTRUME` keyword of the files, never chosen on a
command line, because reading the wrong extension raises no error: it returns
different photons.

```
data/                      <- input.directory
  TOI-2120/
    2811170t.fits ...
outputs/                   <- output.directory
  TOI-2120/
    1-3v/
      resolved_config.yaml, fit.npz, twoframe_components.fits
      TOI-2120_1-3v.pdf
      corrected/
      star_template.fits, lbl_config.yaml, run_lbl.py
```

The tag is `<star>-<observer>` components, and `v` when the velocity term was
fitted. The five stages, each announcing how long it took:

| stage | what it leaves |
| --- | --- |
| `cube` | every spectrum on one log-uniform grid, in `cache/`; about twenty minutes the first time, reused after |
| `fit` | the two bases and their coefficients |
| `figures` | **one** multipage PDF: the resolved parameters, then every plot |
| `correct` | the corrected t.fits, in `corrected/` |
| `lbl` | the star template, both objects staged, LBL's config and run script, and LBL run if asked |

On 316 SPIRou exposures, fit, figures and correction take twenty to thirty
minutes, and LBL about half an hour per object.

`output.fits_directory` keeps a run's products elsewhere, an external disk for
instance: the run folder becomes one link to the same path under it, made
before anything is written, and so do LBL's folders except `lbl/science`, which
is made of links an exFAT disk cannot hold. A run stops rather than write
locally when that disk is not mounted. `cache/` stays put.

Useful flags: `--dry-run` resolves everything and touches nothing, `--stages
cube,fit` runs part of it, `--n-star`/`--n-earth` override the component
counts, `--rebuild-cube` ignores the cache, `--run-lbl` runs LBL.

## Configuring it

One file, `config.yaml`, in three layers merged in order:

```yaml
general:        # everything that does not depend on the instrument
instruments:    # what does: extensions, wavelength range, bands, LBL profile
  SPIROU: ...
  NIRPS: ...
objects:        # anything one target needs and the others do not
```

**Adding a spectrograph is adding a block under `instruments`**, not editing
any code. The block names which FITS extension holds the flux, the wavelength
solution, the blaze, the telluric reconstruction and the sky model, plus the
wavelength range, the photometric bands and LBL's name for it.

Two conventions worth knowing before changing anything:

* The wavelength solution comes from the **extension**, never from header
  polynomials. The blaze comes from the blaze extension of the same file.
* A reported MAD is `median(|x - median x|)` with no 1.4826 factor. Where a
  sigma-equivalent is wanted the factor is applied at the point of use and
  named a sigma there.

## The grid

The destination grid is log-uniform, `lambda_i = wave0 * exp(i * dv / c)`, so a
Doppler shift is a pure translation of an integer-plus-fraction number of
samples. That is what lets the star basis be carried into each exposure's frame
by an exact Lanczos operator.

The step is `domain.dv`, and `domain.smart_dv: true` measures it from the data
instead: the finest pixel step in the first spectrum, sampled at 70% of its own
width. SPIRou pixels are about 2.27 km/s, so that is a grid of about 1.58 km/s
against the 0.5 the config asks for otherwise: three times fewer samples and a
fit faster by about as much. **With `smart_dv` on, `dv` is not read at all**,
and the run says so where it announces the domain. The high pass is a width
in km/s (`highpass.width_kms`), so it covers the same velocity on any grid.
What smart_dv does not do is rescale the one window still counted in samples,
`weights.empirical_noise_box`; the run says what it now covers in km/s.

## Testing

```
./check.sh
```

489 tests in about thirteen seconds, with no spectrum and no network
access, so a clone with no `data/` at all runs them. They pin
the things that have gone wrong here: the adjoint identity the block solve
depends on, the parity tie producing identical coefficients, the corrected file
being panel 3, a file's OBJECT matching its object whatever case it was typed
in, the rejection threshold being in robust sigmas rather than MADs,
LBL's own reader accepting the config, the runparams and the template written
for it, the STRPCA tables being written between the mask and the velocities,
and a header card that used to be tested against a dict while the real rows
were FITS records.

## Access

The repository is private for now. Contact **Étienne Artigau** for access.
