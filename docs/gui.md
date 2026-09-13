# The window

    pca2d-gui

For anyone who should not have to read a command line before their first run.
It does nothing the command line cannot: every setting is one key of
`config.yaml` or one flag of `pca2d-preclean`, and the command it is about to
run is written out so that it can be copied into a terminal.

## What is on it

- **A data root**, shown in full: `data` and `config.yaml` mean different
  folders from different working directories, and here `data` is a folder of
  links onto a shared disk. Nothing is ever written in the data root, which can
  be a read-only archive.
- **An output root**, proposed as `corrected` beside the data root when the
  field is empty: the corrected spectra are a copy of the campaign, tens of
  gigabytes, and belong on the disk the campaign is already on.
- **The targets**, each with how many spectra it has, **its median SNR**, **its
  median exposure time**, **its magnitude with the band it is in** and its
  instrument. The rows are grouped and tinted by instrument, one box per
  instrument shows or hides its targets, and clicking a column heading sorts by
  it, again to reverse, the object heading a third time for the default order.
  What is not known sorts last, never first.
- **A box per target.** Tick one for a solo run. **Tick several and they are
  fitted together against one observer basis**, each keeping its own star
  spectrum per order parity (`docs/joint_fit.md`). Click the box, double-click
  the row, or press the space bar; `all` and `none` do the whole list. They must
  come from ONE instrument: two ticked instruments refuse to run, because a run
  is one domain, one grid and one set of extensions, all read from the
  instrument.
- **The settings that change a result**, and only those: how many components in
  each frame, the static part, the star basis, the velocity term, the number of
  sweeps, the shrinkage, which samples a corrected file blanks, the high pass,
  the grid step, nightly coadding, and whether LBL is run. Hovering one writes
  what it means at the bottom of the window, and after a moment in a box.
- **The stages.** cube, fit, figures, correct, lbl: any subset, in that order.
- **A variant**, from `variants/`: the nominal plus the lines that variant
  changes. What was measured with it is in `variants/README.md`.
- **The command**, in full, before anything runs.
- **The output**, streamed as the run prints it, in the colours it uses in a
  terminal: green for progress, blue for a number, orange for something
  skipped, red for what stops a run. The window's own lines are in the same
  form, `YYMMDD HH:MM:SS.SS | message`, and in the window's language.
- **The barycentric coverage** of whatever is ticked, as a histogram in bins of
  3 km/s, stacked by star, with three numbers under it: what the campaign
  covers, its span, and what its ecliptic latitude ALLOWS. That last one is
  fixed by the sky, `|BERV| <= 29.78 cos(beta)`, and it decides whether the
  method can work on a target at all: TOI-1452 at +80.5 degrees can never span
  more than 9.8 km/s. A coverage short of the span means gaps, which more nights
  fill; a span short of the possible means a young campaign; a small possible
  means the wrong target. While the scan is still reading, the panel says so
  rather than looking final.
- **English or French**, one button, everything included: the labels, the
  explanations and the window's own log lines.

## The index of a data root

The first time a data root is opened, one header of every spectrum is read: the
per-order extraction SNR APERO wrote (`EXTSNxxx`), `EXPTIME`, `BERV`, the
coordinates, the magnitude under whichever keyword that instrument writes, and
`INSTRUME`, which is read from the FILE and never from a configuration. The answers go into
ONE file per root under `~/.pca2d/scans/`, named after the root with a digest of
its full path, so two roots whose last folder is called `science` cannot share
one index.

Every later visit reads that index first and draws the list at once, then checks
the root: a file is taken as it was if its size and modification time are
unchanged, and only what was added or replaced is read. The names and the file
counts go up before a single header is read, and each campaign's numbers arrive
every ten spectra, marked with a tilde until all of them are in. Those ten are
spread across the campaign rather than taken in order, since file names sort by
date and the first ten of a campaign are one night's weather: on GJ 1 they gave
a signal-to-noise of 134 against the campaign's 163.

A field added to the index re-reads the files that lack it, and only those,
while the rows keep showing everything else they already knew. **Rescan** is that same
check, for after copying new spectra in. A file that cannot be read is recorded
as unreadable rather than re-read at every visit. Deleting an index costs the
few seconds of one scan, nothing else.

## The LBL window

**LBL settings...** opens the `lbl:` block on its own: what is measured (the
delivered spectra, the corrected ones, or both), what the corrected object is
called next to the delivered one, which template it is measured against, where
LBL's tree is, which of its steps run, and whether the spectra get there as
symlinks or copies. That is the equivalent of the wrapper, and it is the step
that produces velocities.

## The buttons that write something

- **Export YAML** saves the settings that differ from the configuration as a
  variant file, which is how one run is made reproducible. A variant that
  repeats the nominal says nothing, and the window refuses to write one.
- **Save as defaults** writes those same settings into `config.yaml` itself, as
  what every later run starts from. It edits the lines it has to and nothing
  else, so the comments that say which measurement chose each value survive; a
  `yaml.safe_dump` of the parsed document would have deleted all of them. It
  asks first, and shows exactly which keys it is about to change.
- **Save log** writes what the window has shown.

## What it does not do

It does not resolve the configuration itself, does not touch the data, and holds
no state of its own beyond the last settings and what was ticked, which it
remembers in `~/.pca2d_gui.json`, and the per-root index under `~/.pca2d/`.
Everything else belongs to `pca2d-preclean`, which it runs as a subprocess and
can stop.
