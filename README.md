# pca2d-preclean

Two-frame weighted PCA of echelle spectra. One basis anchored in the **star's**
rest frame, one anchored in the **observer's**, fitted jointly against the same
data, and then only the observer's half is divided back out. What is fixed to
the sky comes off; what belongs to the star stays.

The point is what a velocity code sees afterwards.

```
pca2d-preclean --object TOI-2120
```

That is the whole interface.

## What it does

A spectrum contains two things that move differently. The star's lines shift
with the barycentric velocity; the atmosphere's do not. Over a year of
observing, the barycentric velocity sweeps by tens of km/s, which is enough to
separate them if the decomposition is told that one basis travels and the other
does not.

So the model is

    y_n  =  S_n P a_n  +  Q b_n

where `S_n` is an exact translation by the exposure's own barycentric velocity,
`P` is the basis that follows the star and `Q` the one that stands still on the
detector. Both are solved for together, by block coordinate descent, under
weights that come from the noise you can measure rather than the noise a photon
model claims.

That last part is not fastidiousness. Past about 2200 nm on SPIRou the thermal
background dominates the detector. Those are photons and they were counted, so
they arrive with their own Poisson noise, and the pipeline subtracts their MEAN
and not their variance: the noise stays in the data after the signal it belonged
to has gone. A photon sigma computed from the flux that is left describes a
spectrum nobody recorded, and it declares the noisiest part of the array the
quietest, by a factor of fifty. The sample-to-sample scatter is measured from
what is actually there and cannot make that mistake, so the weights use
whichever of the two is larger.

Nothing is subtracted before the fit except one instrumental offset. There is
no median template and no mean spectrum taken out in front, because whatever
comes out before the fit is outside the model, and the correction removes
components: a term subtracted early would be fitted by nothing and removed from
nothing.

## Installing it

```
git clone https://github.com/eartigau/pca2d_preclean.git
cd pca2d_preclean
conda env create -f environment.yml
conda activate pca2d-preclean
```

That environment holds **both codes**: `pca2d-preclean` and
[LBL](https://github.com/njcuk9999/lbl), which is what the corrected spectra
exist to be fed to. Having to deactivate one to run the other is how a t.fits
gets measured by the wrong version of something.

It is why the versions are nailed down rather than floored. LBL pins its
dependencies exactly, `numpy==2.3.3`, `astropy==7.2.0`, `scipy==1.17.0` and the
rest, and asks for python >=3.12,<3.13; `environment.yml` repeats those pins so
that **conda** installs them and pip finds them already satisfied. Only `lbl`
itself, which is not on PyPI, this package in place, and one PyPI-only
dependency come through pip, so conda and pip never fight over numpy and
`conda env update` later cannot clobber what pip put there. `astropy-base` and
`matplotlib-base` rather than the metapackages: nothing here imports pyarrow or
bqplot, and every entry point calls `matplotlib.use("Agg")` before it draws, so
the GUI toolkits would be installed to be never loaded.

`pyproject.toml` keeps this package's own honest floor of 3.10 and its looser
bounds, for anyone installing it on its own.

After that, `pca2d-preclean`, `lbl_find`, `lbl_setup`, `lbl_demo` and
`lbl_reset` are all on the path, and `./check.sh` is 83 tests in under a
second.

## Running it

Two roots, and a run reads from one and writes to the other. Spectra go under
the input root, one folder per target, and nothing is ever written there, so it
can be a shared or read-only archive. The instrument is read from the
`INSTRUME` keyword of the files themselves, never chosen on a command line,
because reading the wrong extension raises no error: it returns different
photons.

```
data/                      <- general.input.directory, the input root
  TOI-2120/
    2811170t.fits
    2811171t.fits
    ...
outputs/                   <- general.output.directory, the output root
  TOI-2120/
    2-7/
      resolved_config.yaml
      fit.npz, twoframe_components.fits
      TOI-2120_2-7.pdf
      corrected/
```

Then

```
pca2d-preclean --object TOI-2120
```

runs four stages, each announcing how long it took:

| stage | what it leaves |
| --- | --- |
| `cube` | every spectrum on one log-uniform grid, in `cache/` |
| `fit` | the two bases and their coefficients, in `outputs/<object>/<M>-<N>/` |
| `figures` | **one** multipage PDF: the resolved parameters, then every plot |
| `correct` | the observer block divided out, as t.fits, in `outputs/<object>/<M>-<N>/corrected/` |
| `lbl` | both sets of spectra set up for LBL, delivered and corrected, and the two files that run it |

Both roots live in `config.yaml`, since a run is a config and an object and
nothing else. `--data-dir` and `--out-dir` override them for one run, which is
what a scratch disk or a second machine needs. `cache/` is neither: it holds
rebuildable intermediates, and it stays put so that changing where a run's
products go does not orphan a cube that took twenty minutes.

Useful flags: `--dry-run` resolves everything and touches nothing, `--stages
cube,fit` runs part of it, `--n-star`/`--n-earth` override the component
counts, `--rebuild-cube` ignores the cache.

## Measuring it

Correcting a spectrum is worth what the velocity is worth, and the only honest
way to know whether it helped is to measure both. So the `lbl` stage never sets
LBL up on the corrected files alone: it puts **two objects** side by side in one
LBL tree, from the same instrument profile, each building its own template.

```
lbl/science/TOI-2120/               -> symlinks to the spectra as delivered
lbl/science/TOI-2120_PCA2D_2-7/     -> symlinks to what this run corrected
```

Symlinks, not copies: the spectra already exist twice and a third copy buys
nothing. The name carries the component counts because LBL globs a science
folder, so a 2-7 and a 3-5 correction landing in one folder would be measured
as a single series with nothing said about it.

Beside the run's other outputs it leaves the two files LBL needs, and they are
both meant to be read and edited:

| file | what it is |
| --- | --- |
| `lbl_config.yaml` | LBL's configuration in LBL's own keys, which `lbl_compute --config` reads. Every key is one LBL knows; it refuses any other, which is a good reason to have it written rather than typed |
| `run_lbl.py` | an ordinary LBL wrap script, with both objects in the `rparams` dict LBL users already know |

Running LBL is hours, so the stage prepares and stops there, and says what to
run. `lbl.run: true` in the config, or `--run-lbl`, has it run instead.

Which LBL instrument a spectrograph is comes from its block in `config.yaml`,
not from the header: LBL calls NIRPS `NIRPS_HA` or `NIRPS_HE` by the mode it
was observed in, and picking the wrong one raises no error, it returns
velocities from another instrument's profile. `lbl.teff` is worth filling in,
since LBL chooses the stellar model its mask comes from by effective
temperature.

## Configuring it

One file, `config.yaml`, in three layers merged in order:

```yaml
general:        # everything that does not depend on the instrument
instruments:    # what does: extensions, wavelength range, bands
  SPIROU: ...
  NIRPS: ...
objects:        # anything one target needs and the others do not
```

**Adding a spectrograph is adding a block under `instruments`**, not editing
any code. The block names which FITS extension holds the flux, the wavelength
solution, the blaze, the telluric reconstruction and the sky model, plus the
wavelength range and the photometric bands.

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
by an exact Lanczos operator instead of an interpolation that would smear
whatever it touched.

The step is `domain.dv`, and `domain.smart_dv: true` measures it from the data
instead of taking it on faith: the finest pixel step in the first spectrum,
sampled at 70% of its own width, which is all a resampling has to guarantee.
SPIRou pixels are about 2.27 km/s, so that is a grid of about 1.58 km/s against
the 0.5 the config asks for otherwise: three times fewer samples, three times
less memory, and a fit that is faster by about as much.

**With `smart_dv` on, `dv` is not read at all.** It is overwritten and not
combined with anything, so `dv: 0.5`, `dv: 4.0` and `dv: null` beside a true
`smart_dv` are the same run down to the cube's cache key. The run says so on
the line where it announces the domain, and `dv: null` is the honest way to
write it. What smart_dv does not do is rescale the windows that are counted in
samples, `highpass.window` and `weights.empirical_noise_box`; the run says what
each of them now covers in km/s and leaves that decision where it belongs.

Echelle orders are split by parity into two rows per exposure. Consecutive
orders overlap, and at a given wavelength one parity samples near an order
centre and the other near an edge, where the resolution is not the same; a
single mean over both leaves that difference in the residual, coherent and
large. It is the one thing taken out before the fit.

## Testing

```
./check.sh
```

83 tests, under a second, no file and no network access. They pin the things
that have gone wrong here: the adjoint identity the block solve depends on, the
parity tie producing bit-identical coefficients, the rejection threshold being
in robust sigmas rather than MADs, the isolated-sample rule, the memory
decision behind nightly coadding, the object name being joined to the input
root exactly once, that the corrected object carries the component counts
that made it so two runs cannot be measured as one, and that neither
command-line tool returns a value to `sys.exit`.

## Access

The repository is private for now. Contact **Étienne Artigau** for access.
