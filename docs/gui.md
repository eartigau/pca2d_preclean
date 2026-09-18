# The window

    pca2d-gui

For anyone who should not have to read a command line before their first run.
It does nothing the command line cannot: every setting is one key of
`config.yaml` or one flag of `pca2d-preclean`, and the command it is about to
run is written out so that it can be copied into a terminal.

## Three pages

Everything at once was a wall. The window is six pages, and each has room to be
read: **targets**, what is being reduced; **settings**, how; **LBL**, what the
step that produces velocities is asked; **analysis**, the command, the
buttons that start and end it, and everything it says; **runs**, what has been
run before; and **cleanup**, what all of it has left on the disks. The banner
and the language button stay above them. What follows is where each thing is.

The corners are round, which ttk does not do: no theme here draws one, so the
tabs, the buttons, the fields and the tick boxes are pictures of themselves,
nine-patch images that ttk stretches along the middle and leaves alone at the
corners (`App._round`). They are drawn eight times over and shrunk with a
Lanczos filter, since ImageDraw has no antialiasing and a corner drawn at its
final size is a staircase. The selected tab is the same ground as the page under
it, so it reads as that page's own edge rather than a card floating over it.
Without PIL none of it happens and the window is the square one it was.

## What is on it

- **A data root**, shown in full: `data` and `config.yaml` mean different
  folders from different working directories, and here `data` is a folder of
  links onto a shared disk. The only thing ever written in a data root is one
  `pca2d_index.csv` per campaign folder, the log below; the spectra are never
  touched, and a campaign whose folder refuses the file, an archive mounted
  read-only, stays exactly as it was.
- **An output root**, proposed as `corrected` beside the data root when the
  field is empty: the corrected spectra are a copy of the campaign, tens of
  gigabytes, and belong on the disk the campaign is already on.
- **A products disk**, optional (`output.fits_directory`): each run folder
  becomes one link to it, so the corrected spectra, the fit and the figures land
  there rather than on the internal disk. It is a path on ONE machine, so
  `config.yaml` no longer carries one and the window **empties the field
  whenever what it names is not there**, whether that came from the
  configuration or from what the window itself remembered. Empty is a decision
  and the command says it, `--no-fits-dir`: everything stays under the output
  root, which works.
- **On a first run** the two roots are EMPTY and the configuration is the
  `config.yaml` that came with the installation, found beside the package
  rather than resolved against whatever folder the window was started from.
  Where the spectra are and where their copies go are choices about somebody's
  disks, and a window that opens with a plausible path in those fields invites
  a run against a folder nobody picked. So the status line says `pick a data
  root` for as long as there is none, the log says what such a folder holds,
  and the command box says `tick a target` instead of showing a command with no
  `--object` in it.
- **The targets**, each with how many spectra it has, **its median SNR**, **its
  median exposure time**, **its magnitude with the band it is in** and its
  instrument. The rows are grouped and tinted by instrument, one box per
  instrument shows or hides its targets, and clicking a column heading sorts by
  it, again to reverse, the object heading a third time for the default order.
  What is not known sorts last, never first. A campaign still being read
  carries a disc beside its tick box, filling clockwise from twelve o'clock in
  twelve steps as its files come in, and gone once they all are, so a finished
  list is not a column of symbols. Drawn rather than written, because a
  character is the size of the name beside it.
- **A box per target.** Tick one for a solo run. **Tick several and they are
  fitted together against one observer basis**, each keeping its own star
  spectrum per order parity (`docs/joint_fit.md`). Click the box, double-click
  the row, or press the space bar; `all` and `none` do the whole list. They must
  come from ONE instrument: two ticked instruments refuse to run, because a run
  is one domain, one grid and one set of extensions, all read from the
  instrument.
- **The settings that change a result**, and only those: how many components in
  each frame, the velocity term, the number of sweeps, the shrinkage, the **correction fit
  metric** as two buttons, `F` and `(dF/dv)²`, the second being the nominal and
  weighting every sample by the star's own derivative there, the high pass, the
  grid step and nightly coadding. **Each travels to the run under its own
  flag**: until 2026-09-15 only the two component counts did, so a value typed
  here was shown, changed, and then ignored by the run, which used the
  configuration's own. What is settled is not among them,
  and is fixed in the code (`config.DEFAULTS`): the star is one cubic B-spline
  (`star_basis`), the static part is one star spectrum per order parity
  (`mean: star`), and a sample lost by one exposure is blanked in all of them,
  so the campaign carries one set of lines (`correct.mask: common`). The older
  paths are still read from a config or a variant file, for redoing the runs
  that were made on them. Whether LBL is run is on the LBL page, since the
  stages already say whether the lbl step happens at all. Hovering one writes
  what it means at the bottom of the window, and after a moment in a box.
- **The stages, as the chain they are:** `cube → fit → correct → lbl`, with the
  arrows drawn, because that is the order they happen in and the order they
  depend on each other in. Unticking one unticks everything after it, since
  nothing downstream has its input any more, and ticking `lbl` ticks `correct`,
  which is what it measures. Ticking does not pull the chain the other way: the
  fit is always redone when it is asked for, so correcting again with the fit
  that is already there has to stay sayable. **figures** sits apart on the same
  line, with no arrow to it: it draws what the fit left, and nothing waits for a
  drawing.
- **No variant picker.** The parameters converged, and a second set of settings
  offered beside the settings is a window that contradicts itself.
  `pca2d-preclean --object X --variant NAME` still runs one, and Export YAML
  still writes one, which is what `variants/README.md` reproduces from.
- **The command**, in full, before anything runs.
- **The output**, streamed as the run prints it, in the colours it uses in a
  terminal: green for progress, blue for a number, orange for something
  skipped, red for what stops a run. The window's own lines are in the same
  form, `YYMMDD HH:MM:SS.SS | message`, and in the window's language.
- **The barycentric coverage** of whatever is ticked, as a histogram in bins of
  3 km/s, stacked by star, counted on a numbered y axis, filling whatever height
  the page has beside the list and redrawn as the window is resized, with three numbers under it: what the campaign
  covers, its span, and what its ecliptic latitude ALLOWS. That last one is
  fixed by the sky, `|BERV| <= 29.78 cos(beta)`, and it decides whether the
  method can work on a target at all: TOI-1452 at +80.5 degrees can never span
  more than 9.8 km/s. A coverage short of the span means gaps, which more nights
  fill; a span short of the possible means a young campaign; a small possible
  means the wrong target. While the scan is still reading, the panel says so
  rather than looking final.
- **The signal-to-noise against the barycentric velocity**, one dot per
  spectrum in its star's colour, under the histogram. The histogram says which
  velocities a campaign covers; this says what it covers them WITH. The fit
  weighs a spectrum by 1/sigma^2, so a range whose far end is covered by a
  campaign's worst nights is not the range the fit really sees, and the two
  frames come apart less well than the coverage promises.
- **English, French, Spanish or Portuguese**, two letters each at the top, the
  other languages than the one the window is in: everything is included, the
  labels, the explanations and the window's own log lines. No flag, and no
  emoji anywhere else either: a colour emoji needs a font that carries it, and
  a machine without that font draws a box, or nothing (WSL, 2026-09-17). The
  ticks of the target list are `[x]` and `[ ]` for the same reason.

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

**The log beside the data.** The same answers are written to
`pca2d_index.csv` **inside each campaign's own folder**, one line per spectrum
with the keywords read from it, appended every ten files and rewritten when that
campaign is finished. Inside the folder rather than at the top of the root
because campaigns are copied and rsynced one at a time: a log that stays behind
is a campaign read again at the other end. That is what a second window reads:
open the same root anywhere else, on this machine or another one mounting the
same disk, and it starts from what has already been read. The header row is the
schema, so a log written by a version with other columns is ignored rather than
half believed, and every line is still checked against its file's size and
modification time before it is used.

A field added to the index re-reads the files that lack it, and only those,
while the rows keep showing everything else they already knew. **Rescan** is that same
check, for after copying new spectra in, and the window does it for itself:
every ten seconds it lists the data root, one folder listing per campaign with
nothing opened, and reads what appeared. A campaign copied in while the window
sits open turns up on its own, rather than the next time somebody thinks to
press Rescan. A file count is compared with the PREVIOUS look and never with the
list, since the list counts what the index holds and a spectrum the scan could
not read is missing from it for good. A file that cannot be read is recorded
as unreadable rather than re-read at every visit. Deleting an index costs the
few seconds of one scan, nothing else.

## The LBL page

The `lbl:` block, between the settings and the analysis: what is measured (the
delivered spectra, the corrected ones, or both), what the corrected object is
called next to the delivered one, which template it is measured against, where
LBL's tree is, which of its steps run, whether LBL is actually run, and whether
the spectra get there as symlinks or copies. That is the equivalent of the
wrapper, and it is the step that produces velocities. It had a window of its own
until the pages existed, which was one window too many: these are settings like
the others, they are simply LBL's rather than the fit's.

## The runs page

Every run the output root holds, read from the runs themselves rather than from
anything the window remembers: each writes `resolved_config.yaml` into its own
folder before its first stage, with every setting as it resolved it, the command
it was given word for word, the code it ran and the time it started
(`pca2d/runs.py`). A run launched from a terminal, or by somebody else on the
same disk, is therefore in the list too, and one that is going appears as soon
as it has started.

The list gives the run's hash first (the six characters its folder and its PDF
carry, the hash of its command: two runs of the same targets and the same
counts are told apart by it), then the time it started, its targets, its
counts, its components and whether its compilation PDF is there; picking one fills the table under it with
every option that run was given, and the two buttons open its PDF or its folder.
A double-click opens the PDF.

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
- **Quit**, in the banner beside the language button, where every page can
  reach it: it ends the window, which is not a property of whichever page
  happens to be open. It closes the window. The
  settings are written at every change, so nothing here is lost by leaving; a
  run is a subprocess of the window and is stopped by leaving, so it asks first.
  The window's own close box asks the same question.

## What it does not do

It does not resolve the configuration itself, does not touch the data, and holds
no state of its own beyond the last settings and what was ticked, which it
remembers in `~/.pca2d_gui.json`, and the per-root index under `~/.pca2d/`.
Everything else belongs to `pca2d-preclean`, which it runs as a subprocess and
can stop.
