# The window

    pca2d-gui

For anyone who should not have to read a command line before their first run.
It does nothing the command line cannot: every setting is one key of
`config.yaml` or one flag of `pca2d-preclean`, and the command it is about to
run is written out so that it can be copied into a terminal.

## What is on it

- **A data root.** The objects it holds appear in the list, with how many
  spectra each has and, read from the first file that opens, which instrument
  took them.
- **The objects.** Select one for a solo run. **Select several and they are
  fitted together against one observer basis**, each keeping its own star
  spectrum per order parity (`docs/joint_fit.md`). They must come from one
  instrument; the window says so before the run refuses them.
- **The settings that change a result**, and only those: how many components in
  each frame, the static part, the star basis, the velocity term, the number of
  sweeps, the shrinkage, which samples a corrected file blanks, the high pass,
  the grid step, nightly coadding, and whether LBL is run. Hovering one writes
  what it means at the bottom of the window.
- **The stages.** cube, fit, figures, correct, lbl: any subset, in that order.
- **A variant**, from `variants/`: the nominal plus the lines that variant
  changes. What was measured with it is in `variants/README.md`.
- **The command**, in full, before anything runs.
- **The output**, streamed as the run prints it, in the colours it uses in a
  terminal: green for progress, blue for a number, orange for something
  skipped, red for what stops a run. The stage it is in shows at the right.

## The two buttons that write something

- **Export YAML** saves the settings that differ from the configuration as a
  variant file. A variant that repeats the nominal says nothing, and the window
  refuses to write one.
- **Save log** writes what the window has shown.

## What it does not do

It does not resolve the configuration itself, does not touch the data, and
holds no state of its own beyond the last settings, which it remembers in
`~/.pca2d_gui.json`. Everything else belongs to `pca2d-preclean`, which it runs
as a subprocess and can stop.
