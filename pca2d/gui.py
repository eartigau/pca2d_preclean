"""A window to run the pipeline, for people who should not have to read the
command line first.

    pca2d-gui

It shows the objects that are in the data root, lets one of them be run on its
own or several of them together against one observer basis, exposes the
handful of settings that change a result, writes the command it is about to
run so that it can be copied into a terminal, and streams the run's own log
into the window with the colours it would have in a terminal.

Every item explains itself on hover, in English or in French, and says what
the choice implies rather than only what it is called. Nothing here decides
anything: every option maps to one key of config.yaml or one flag of
`pca2d-preclean`, the command is shown before it runs, and the settings can be
exported as a variant file, which is how a run is made reproducible
(variants/README.md).
"""
from __future__ import annotations

import glob
import json
import os
import queue
import re
import shutil
import subprocess
import sys
import threading
import time

from . import scan
from .logger import stamp

HOME_STATE = os.path.expanduser("~/.pca2d_gui.json")
ANSI = re.compile(r"\033\[(\d+)m")
#: One palette, taken from the family the rest of these tools belong to
#: (themes_outils/DESIGN_SYSTEM.md): the same sky-blue accent and the same dark
#: ground, so that going from one of them to another feels like one place. The
#: surfaces here are light, because this window is mostly a table of numbers,
#: and the log alone is dark, because it is a terminal stream and reads better
#: that way.
BG = "#eef2f8"          # the window's ground
SURFACE = "#ffffff"     # cards, the table, the entries
INK = "#16202e"         # text
MUTED = "#5d6b81"       # secondary text
ACCENT = "#0b6ba8"      # the family's #62c2ff, darkened to sit on white
ACCENT_SOFT = "#d8ecfa"
LINE = "#c9d4e4"
LOG_BG = "#0d1826"      # the family's --bg
LOG_INK = "#dfe8f5"
#: the logger's colours, and how this window paints them, on the dark log
LEVELS = {"32": ("info", "#5fd78a"), "34": ("value", "#62c2ff"),
          "33": ("warn", "#ffb454"), "31": ("error", "#ff6b6b")}
STAGES = ("cube", "fit", "figures", "correct", "lbl")
#: what the window can change, and the config key each one is
OPTIONS = [
    ("n_star", "twoframe.n_star", "int"),
    ("n_earth", "twoframe.n_earth", "int"),
    ("mean", "twoframe.mean", ("star", "offset", "full", "iterate")),
    # no star_basis here: the star is one cubic B-spline, decided in the code
    # (config.DEFAULTS). The grid path is still read from a variant file, for
    # redoing the runs that were made on it.
    ("velocity_term", "twoframe.velocity_term", "bool"),
    ("iters", "twoframe.iters", "int"),
    ("shrink", "correct.shrink", "bool"),
    ("mask", "correct.mask", ("common", "exposure", "none")),
    ("width_kms", "highpass.width_kms", "float"),
    ("dv", "domain.dv", "float"),
    ("nightly_stack", "input.nightly_stack", ("auto", "true", "false")),
    ("run", "lbl.run", "bool"),
]
#: the `lbl:` block, in its own window: what the wrapper would have been asked
OPTIONS_LBL = [
    ("lbl_prepare", "lbl.prepare", "bool"),
    ("lbl_before", "lbl.before", "bool"),
    ("lbl_after", "lbl.after", "bool"),
    ("lbl_star_template", "lbl.star_template", "bool"),
    ("lbl_strpca", "lbl.strpca", "bool"),
    ("lbl_directory", "lbl.directory", "text"),
    ("lbl_suffix", "lbl.suffix", "text"),
    ("lbl_teff", "lbl.teff", "text"),
    ("lbl_template", "lbl.template", "text"),
    ("lbl_steps", "lbl.steps", "list"),
    ("lbl_link", "lbl.link", ("symlink", "copy")),
]
#: every option the window can change, wherever it is shown
ALL_OPTIONS = OPTIONS + OPTIONS_LBL
#: a target is run or not run, and the list says which with a box
CHECKED, UNCHECKED = "☑", "☐"
#: shown where an instrument is not known YET, as against not known at all
UNREAD = "…"
#: one tint per instrument, on the row's background. A joint fit is one
#: instrument, so which instrument a target belongs to is the first thing the
#: list has to make obvious; the colours are pale enough that the ticked box and
#: the selection still read over them.
INSTRUMENT_TINT = {"NIRPS": "#eaf3ff", "SPIROU": "#fff1e6"}
OTHER_TINTS = ("#eefaf0", "#f6eeff", "#fdf6e3", "#f0f0f0")
#: one per STAR, for the bars of the coverage histogram. Not per instrument: a
#: pool of four NIRPS campaigns stacked in one blue says nothing about which
#: star fills which bin, which is the whole reason for stacking them. A
#: colour-blind safe categorical palette, and it survives being printed.
STAR_COLOURS = ("#0072b2", "#d55e00", "#009e73", "#cc79a7", "#56b4e9",
                "#e69f00", "#332288", "#882255")

EN = {
    "subtitle": "two-frame precleaning, then LBL. Pick a data root, pick"
                " objects, look at the command, run it.",
    "data_dir": "data root", "config": "config",
    "out_dir": "output root (optional)",
    "browse": "Browse", "rescan": "Rescan",
    "objects": "objects", "settings": "settings", "stages": "stages",
    "variant": "variant", "command": "the command this runs", "output": "output",
    "col_object": "object", "col_files": "files", "col_instrument": "instrument",
    "col_snr": "SNR", "col_exptime": "exp (s)", "col_mag": "mag",
    'run_name': 'run name',
    'timeline': 'when the ticked campaigns were observed',
    'timeline_none': 'tick a target to see when it was observed',
    'timeline_note': 'keeping %s to %s: %d exposures of %d',
    'whole': 'all of it',
    'min_rjd': 'from',
    'max_rjd': 'to',
    'exists': '⚠ this run already exists',
    'help_run_name': 'A name for this run. Its products go to <output root>/_NAME/ and its LBL object is <object>_PCA2D_<M-N>_NAME, so two runs of the same targets at different settings never write into one folder nor under one LBL name, where LBL would measure the mixture without a word. Empty is the nominal path. It is proposed from the dates when they are set, and can be anything.',
    'help_dates': 'Keep only the exposures between these two dates, in reduced Julian date (BJD - 2400000). What is excluded is neither fitted nor corrected. It is how a campaign is cut to a season, which the barycentric coverage sometimes asks for: two 14-night slices of TOI-4552 at the same signal-to-noise differ by a factor of sixty in coverage and the correction changes sign between them.',
    "berv": "barycentric coverage of the ticked targets",
    "berv_none": "tick a target to see its barycentric coverage",
    "berv_wait": "no barycentric velocity read yet",
    "berv_building": "UNDER CONSTRUCTION: still reading",
    "berv_note": "effective coverage %.0f km/s, span %.0f, at most %s possible"
                 " for this target's ecliptic latitude; bins of %.0f km/s,"
                 " %d exposures",
    "help_berv":
        "How much of the barycentric range the ticked campaigns actually cover,"
        " in bins of 3 km/s, one colour per instrument and stacked. This is the"
        " quantity that decides whether the correction helps: the observer block"
        " is identifiable only because the star moves through it. Two 14-night"
        " slices of TOI-4552 at the same signal-to-noise, one spanning"
        " 42.8 km/s and the other 0.7, went from a gain of a factor two to a"
        " loss of a factor 2.4. The EFFECTIVE coverage counts only the bins that"
        " hold an exposure, so a campaign observed at two extremes and nowhere"
        " between is not credited with the range between them. The third number"
        " is what the sky ALLOWS: a target's |BERV| never exceeds"
        " 29.78 cos(ecliptic latitude) km/s, so TOI-1452 at +80.5 deg can only"
        " ever span 9.8 km/s and no amount of observing will separate the two"
        " frames for it. A coverage short of the span means gaps, which more"
        " nights can fill; a span short of the possible means the campaign is"
        " young; a small possible means the target is the wrong one.",
    "run": "Run", "stop": "Stop", "dry": "Dry run", "export": "Export YAML...",
    "savelog": "Save log...", "openout": "Open outputs",
    "savedefaults": "Save as defaults...", "lblwin": "LBL settings...",
    "all": "all", "none": "none",
    "idle": "idle", "running": "running", "lang": "Français",
    "lbl_title": "LBL settings", "close": "Close",
    "nothing_export": "every setting is the configuration's own, so a variant"
                      " file would say nothing.",
    "export_title": "save these settings as a variant",
    "mixed_title": "two instruments",
    "mixed": "%s are not from one instrument (%s).\n\nOne run is one domain, one"
             " grid and one set of extensions, all read from the instrument, so"
             " objects from two spectrographs cannot be fitted together. Tick"
             " the ones from a single instrument.",
    "defaults_title": "save as defaults",
    "defaults_ask": "Write these values into %s as the defaults for every"
                    " run?\n\n%s\n\nThe comments in the file are kept.",
    "log_index": "index of this data root: %s",
    "log_out_proposed": "output root proposed, beside the data: %s",
    "log_scan_start": "reading the data root %s",
    "log_scan_done": "%s: %d objects, %d spectra read, %d already known, %d gone",
    "log_scan_none": "no folder with spectra under %s",
    "log_busy": "still reading the data root: wait for it to finish",
    "log_no_root": "not a folder: %s",
    "log_no_object": "no object ticked: tick at least one in the list",
    "log_mixed": "two instruments ticked (%s): a run is one instrument",
    "log_mixed_refused": "refusing to run: %s are from two instruments (%s)",
    "log_nights": "%s: %d nights in common, out of %s",
    "log_nights_thin":
        "%s: %d nights in common out of %s. Few is an ADVANTAGE here: the"
        " observer components are a parasite that is always present, what is"
        " being measured is its pattern, and nights no other star saw widen the"
        " range of conditions that pattern is measured over.",
    "log_command": "running: %s",
    "log_ended": "the run ended, exit code %d",
    "log_stopping": "asking the run to stop",
    "log_export": "variant written: %s",
    "log_defaults": "defaults written into %s: %s",
    "log_nothing": "nothing to write: every setting is the configuration's own",
    "log_saved_log": "log written: %s",
    "log_no_outputs": "nothing written there yet: %s",
    "log_failed": "could not start it: %s",
    "log_no_command": "nothing was installed, so there is nothing to run",
    "log_installing": "installing the command into this environment, from %s",
    "no_command_title": "the command is not installed",
    "no_command":
        "pca2d-preclean cannot be found for\n\n  %s\n\nInstall it now from"
        " %s?\n\nThat runs `pip install -e . --no-deps` there, which puts the"
        " command in this interpreter's own bin and touches nothing else.",
    "no_command_still": "still not found after installing",
    "opt_n_star": "star components", "opt_n_earth": "observer components",
    "opt_mean": "static part",
    "opt_velocity_term": "fit a velocity per exposure",
    "opt_iters": "sweeps at most",
    "opt_shrink": "divide only what is significant",
    "opt_mask": "samples a corrected file blanks",
    "opt_width_kms": "high pass (km/s)", "opt_dv": "grid step (km/s)",
    "opt_nightly_stack": "coadd each night", "opt_run": "run LBL (hours)",
    "opt_lbl_prepare": "write LBL's tree",
    "opt_lbl_before": "measure the delivered spectra",
    "opt_lbl_after": "measure the corrected spectra",
    "opt_lbl_star_template": "use our star as LBL's template",
    "opt_lbl_strpca": "extra star components as RESPROJ",
    "opt_lbl_directory": "LBL data directory", "opt_lbl_suffix": "corrected name",
    "opt_lbl_teff": "effective temperature", "opt_lbl_template": "template file",
    "opt_lbl_steps": "steps", "opt_lbl_link": "spectra into LBL as",
    "help_lbl_prepare":
        "Write lbl_config.yaml and run_lbl.py beside the run's outputs and put"
        " both sets of spectra into LBL's science folders. Off, the correction"
        " is still written and nothing is prepared for a velocity measurement.",
    "help_lbl_before":
        "Measure the DELIVERED spectra as well, as their own object in the same"
        " LBL tree. Off, the correction is measured against nothing and the"
        " result cannot be read as better or worse than doing nothing.",
    "help_lbl_after":
        "Measure the CORRECTED spectra. Off with `before` on, LBL measures only"
        " the delivered ones, which is how a reference series is built once and"
        " then reused by every later run.",
    "help_lbl_star_template":
        "Hand LBL the star spectrum this fit built, instead of letting LBL build"
        " its own from the corrected spectra. OFF is the nominal: measured on"
        " Proxima, our template made LBL fit lines 20% wider and doubled the"
        " error per exposure, 0.97 to 1.67 m/s, for 3.34 m/s of rms against"
        " 3.04. Whatever our template lacks, LBL pays for per line.",
    "help_lbl_strpca":
        "With two or more star components, the ones past the first are handed to"
        " LBL as RESPROJ tables, the way its own DTEMP gradients are, so the rdb"
        " carries each one's projection per exposure and a correlation can be"
        " looked for rather than assumed absent.",
    "help_lbl_directory":
        "LBL's DATA_DIR: its own tree, shared by every object and every run, so"
        " it sits beside the outputs rather than inside one run's folder. Point"
        " it at an existing LBL tree and that tree is used as it is.",
    "help_lbl_suffix":
        "What the corrected object is called next to the delivered one:"
        " TOI-2120 and TOI-2120_PCA2D_2-7. `{tag}` is the run's component"
        " counts, and leaving it out makes two runs write their corrected"
        " spectra into ONE folder, where LBL measures the mixture silently.",
    "help_lbl_teff":
        "The effective temperature LBL is told, which is how it picks its line"
        " list. `auto` reads it from the header, and a number overrides it.",
    "help_lbl_template":
        "A template FILE for LBL to use, by path, instead of the one it would"
        " build. Empty is the nominal: LBL builds a template from the spectra it"
        " is measuring, which is the series each set is measured against.",
    "help_lbl_steps":
        "Which of LBL's steps to run, in order: template, mask, compute,"
        " compile. Fewer is for picking up a tree that already has the earlier"
        " ones, never for skipping work a later step needs.",
    "help_lbl_link":
        "How the spectra get into LBL's science folders. `symlink` is the"
        " nominal: a campaign is tens of gigabytes and LBL only reads them."
        " `copy` is for a filesystem that cannot hold a link, an exFAT disk"
        " above all.",
    "help_col_snr":
        "The median over this target's exposures of the per-order extraction SNR"
        " the pipeline wrote in each file, so a bright target and a faint one can"
        " be told apart before anything is run. Read from the headers alone and"
        " remembered, so a folder opens at once the second time.",
    "help_col_exptime":
        "The median exposure time of this target's files, in seconds. With the"
        " SNR and the file count, it says what kind of campaign this is: many"
        " short exposures of a bright star, or few long ones of a faint one.",
    "help_col_mag":
        "Click a column heading to sort by it, again to reverse, and the object"
        " heading twice to go back to instrument-then-name. Sorting by THIS one"
        " compares bands that are not the same: NIRPS writes J and SPIRou H, so"
        " a J of 5.3 and an H of 10.5 are ordered as numbers and not as"
        " brightnesses.\n\n"
        "The target's brightness as ITS OWN pipeline recorded it, with the band"
        " it is in: NIRPS writes the J magnitude, SPIRou writes H, and the two"
        " differ by about a magnitude on an M dwarf, so the band is shown rather"
        " than assumed. Read from the headers, never from a catalogue.",
    "help_check":
        "Tick a target to run it. Several ticked are fitted TOGETHER against one"
        " observer basis. Click the box, double-click the row, or press the space"
        " bar. They must all come from the same instrument.",
    "help_all_button": "Tick every target in the list, or untick every one.",
    "help_instrument_filter":
        "Show or hide the targets of one instrument. A run is ONE instrument, so"
        " hiding the others is the quickest way to a selection that can actually"
        " be run; the rows are grouped and tinted by instrument for the same"
        " reason. Hiding an instrument unticks its targets, and nothing is"
        " deleted: tick the box again and they come back as they were.",
    "help_savedefaults_button":
        "Write the settings that DIFFER from the configuration into config.yaml"
        " itself, as the defaults every later run starts from. The file's"
        " comments are kept, since they are the measurements that chose each"
        " value. A variant file leaves the nominal alone; this changes it.",
    "help_lblwin_button":
        "The LBL block in its own window: what is measured, how the corrected"
        " object is named, which template it is measured against, which of LBL's"
        " steps run. That is the step that produces velocities.",
    "help_data_dir":
        "The input ROOT, not one object's folder: the objects below are its"
        " subfolders. Nothing is ever written in it.",
    "help_config":
        "config.yaml: everything the window does not show. Three layers are"
        " merged in it: what is general, what belongs to the spectrograph, and"
        " what belongs to one target.",
    "help_out_dir":
        "Where a run writes. Proposed as `corrected` BESIDE the data root, since"
        " the corrected spectra are a copy of the campaign, tens of gigabytes,"
        " and belong on the disk the campaign is already on. Empty uses the"
        " configuration's own output root. A run puts its"
        " resolved configuration, its report, its corrected spectra and its LBL"
        " folders under <root>/<object>/<M>-<N>/.",
    "help_rescan":
        "Read the data root again. What was read before is remembered in an index"
        " under your home folder, never in the data root, so only files that were"
        " added or replaced are read: use this after copying new spectra in.",
    "help_objects":
        "One object ticked: a solo run. SEVERAL: they are fitted together against ONE"
        " observer basis, each keeping its own star spectrum per order parity."
        " The atmosphere and the instrument are shared, the stars are not, so a"
        " basis fitted on several stars cannot follow any one of them. They must"
        " come from the same instrument.",
    "help_n_star":
        "Components fitted in the STAR's rest frame, on top of the star spectrum"
        " itself, which is always taken out with a coefficient of exactly one."
        " 0 is the nominal: a free amplitude in front of a term in a logarithm"
        " is an exponent on the flux, and an exponent on a star's mean spectrum"
        " describes no star. Dropping it took TOI-4552 from 14.6 to 11.2 m/s of"
        " robust scatter, and TOI-2120 from 14.4 to 14.1.",
    "help_n_earth":
        "Components fitted in the OBSERVER's frame: the atmosphere and the"
        " instrument. Three is the nominal. More freedom describes the sky"
        " better and takes more of the star with it wherever the two frames are"
        " degenerate, which is what a narrow BERV coverage does.",
    "help_mean":
        "The part of the model that has no amplitude of its own: the mean"
        " spectrum. There is one per ORDER PARITY, since even and odd orders"
        " see a wavelength at different resolutions, and the question is which"
        " frame it lives in. `star`, the nominal: one spectrum per parity in the"
        " STAR's frame, a BERV-binned median taken out once before any"
        " component with a coefficient of exactly 1, and no observer-frame mean,"
        " so the correction divides out the observer block alone."
        " `offset`: no star spectrum, and only the part of the observer-frame"
        " mean that DIFFERS between parities, the shared part left in for the"
        " observer block to describe; it goes back into the correction."
        " `full`: the whole observer-frame mean per parity, shared part"
        " included. `iterate`: one mean per parity in EACH frame, re-estimated"
        " at every sweep. Measured: TOI-2120 20.4 m/s with `star` against 31.7"
        " with `offset`; on Proxima the observer-frame mean alone injected"
        " 46 m/s and the observer block alone 48, their sum 19. `iterate`"
        " converges on synthetic data and did not on a whole campaign.",
    "help_velocity_term":
        "Fits one velocity per exposure beside the components, to keep the"
        " star's own motion out of the observer block. It is fitted and written"
        " to the corrected files, never divided out. Measured on TOI-4552 it"
        " changed nothing: 16.4 against 15.9 m/s.",
    "help_iters":
        "Sweeps at most. A sweep solves every exposure's amplitudes, then"
        " updates each basis against the other's residual. The fit keeps the"
        " best iterate and stops when chi2 has clearly turned over, so this is a"
        " ceiling, not a duration.",
    "help_shrink":
        "Divides each observer component out only where the data detect it, at"
        " each column, with every exposure counted together: a pattern at half a"
        " sigma in each of N exposures is detected at about 0.5 sqrt(N) and is"
        " kept, while a column where the component is noise is left alone."
        " Without it TOI-4552 lost 3 m/s more: 18.1 against 15.1.",
    "help_mask":
        "A sample the fit gave no weight cannot be corrected, so it is blanked."
        " `common`, the nominal, blanks in every exposure every sample that any"
        " exposure lost, so the campaign carries ONE set of lines. `exposure`"
        " blanks each exposure's own, which cost 1.8 m/s on TOI-4552 with"
        " nothing divided out at all: a line set that moves from epoch to epoch"
        " is scatter. `none` keeps the delivered flux there, uncorrected.",
    "help_width_kms":
        "The width of the Savitzky-Golay filter that removes the continuum, in"
        " km/s, so that it treats a line the same way whatever the grid step. It"
        " is part of the cube: changing it builds a new one.",
    "help_dv":
        "The grid step, in km/s. The grid is uniform in log wavelength, which"
        " makes a Doppler shift an exact translation. Finer is heavier: the cube"
        " and the fit both scale with it. Changing it builds a new cube.",
    "help_nightly_stack":
        "Coadd the exposures of one night before fitting. A memory decision and"
        " never a modelling one: every exposure is registered on its own BERV"
        " first, so nothing is smeared, and the correction still solves each"
        " exposure's own amplitudes against the fixed basis. `auto` weighs the"
        " fit's footprint against the memory the config allows.",
    "help_run":
        "Run LBL once the spectra are corrected. That is the step that produces"
        " velocities, and it is hours. Off, everything LBL needs is still"
        " written and the run says how to start it by hand.",
    "help_stage_cube":
        "Reads every spectrum once onto the common grid, with its weights. A few"
        " minutes, and cached: a second run with the same settings reuses it.",
    "help_stage_fit":
        "The two-frame decomposition. The long one: it prints an R2 after every"
        " sweep.",
    "help_stage_figures":
        "One PDF: the campaign against time, the parameters, then one page per"
        " wavelength window showing every step of the model and what the"
        " correction divides out.",
    "help_stage_correct":
        "Writes the corrected t.fits, one per exposure, with the observer block"
        " divided out and the unweighted samples blanked.",
    "help_stage_lbl":
        "Hands both sets of spectra to LBL, the delivered ones and the corrected"
        " ones, as two objects in one tree, so the velocities can be compared"
        " rather than believed.",
    "help_variant":
        "A file of variants/: the nominal configuration plus the few lines that"
        " variant changes. What each one was measured to give is written in"
        " variants/README.md, and a run made with one is reproducible by name.",
    "help_command":
        "Exactly what the Run button will execute. Copy it into a terminal and"
        " it does the same thing: the window is a wrapper, not a second way of"
        " doing things.",
    "help_run_button":
        "Starts the command above. The output appears below as it is printed,"
        " and the run can be stopped.",
    "help_stop_button":
        "Asks the run to stop. What it has already written stays: a cube, a fit"
        " and corrected files are all reusable, and a stage that did not finish"
        " simply runs again next time.",
    "help_dry_button":
        "Resolves everything, prints the plan, touches nothing. The honest way"
        " to see which cube a set of options points at before spending an hour"
        " on it.",
    "help_export_button":
        "Writes the settings that DIFFER from the configuration as a variant"
        " file, which is how a run is made reproducible. A file that repeats the"
        " nominal says nothing, and is refused.",
    "help_savelog_button": "Writes what the window has shown to a file.",
    "help_openout_button": "Opens the output root in the file browser.",
    "help_lang": "Switch the window between English and French.",
    "help_apero":
        "APERO, the pipeline that reduced every spectrum this window reads"
        " (Cook et al. 2022, PASP 134, 114509). It is what wrote the extensions,"
        " the wavelength solution, the per-order signal-to-noise and the"
        " barycentric velocities used here; this package starts from its t.fits"
        " and never re-derives any of it.",
    "help_log":
        "The run's own log, in the colours a terminal would give it: green for"
        " progress, blue for a number, orange for something skipped, red for"
        " what stops a run.",
}

FR = {
    "subtitle": "prénettoyage à deux référentiels, puis LBL. Choisir un dossier"
                " de données, des objets, regarder la commande, la lancer.",
    "data_dir": "dossier de données", "config": "configuration",
    "out_dir": "dossier de sortie (optionnel)",
    "browse": "Parcourir", "rescan": "Relire",
    "objects": "objets", "settings": "réglages", "stages": "étapes",
    "variant": "variante", "command": "la commande qui sera lancée",
    "output": "sortie",
    "col_object": "objet", "col_files": "fichiers", "col_instrument": "instrument",
    "col_snr": "SNR", "col_exptime": "pose (s)", "col_mag": "mag",
    'run_name': 'nom du passage',
    'timeline': 'quand les campagnes cochées ont été observées',
    'timeline_none': 'cocher une cible pour voir quand elle a été observée',
    'timeline_note': 'on garde du %s au %s : %d poses sur %d',
    'whole': 'tout',
    'min_rjd': 'du',
    'max_rjd': 'au',
    'exists': '⚠ ce passage existe déjà',
    'help_run_name': "Un nom pour ce passage. Ses produits vont dans <racine de sortie>/_NOM/ et son objet LBL est <objet>_PCA2D_<M-N>_NOM, pour que deux passages des mêmes cibles à des réglages différents n'écrivent jamais dans un même dossier ni sous un même nom LBL, où le LBL mesurerait le mélange sans un mot. Vide, c'est le chemin nominal. Il est proposé à partir des dates quand elles sont fixées, et peut être n'importe quoi.",
    'help_dates': "Ne garder que les poses entre ces deux dates, en jour julien réduit (BJD - 2400000). Ce qui est exclu n'est ni ajusté ni corrigé. C'est ainsi qu'on coupe une campagne en saisons, ce que la couverture barycentrique réclame parfois : deux tranches de 14 nuits de TOI-4552 au même SNR diffèrent d'un facteur soixante en couverture, et la correction y change de signe.",
    "berv": "couverture en BERV des cibles cochées",
    "berv_none": "cocher une cible pour voir sa couverture en BERV",
    "berv_wait": "aucune vitesse barycentrique encore lue",
    "berv_building": "EN COURS DE CONSTRUCTION : lecture en cours",
    "berv_note": "couverture effective %.0f km/s, étendue %.0f, au plus %s"
                 " possible vu la latitude écliptique ; bins de %.0f km/s,"
                 " %d poses",
    "help_berv":
        "Quelle part de la gamme barycentrique les campagnes cochées couvrent"
        " réellement, en bins de 3 km/s, une couleur par instrument et empilées."
        " C'est la quantité qui décide si la correction aide : le bloc"
        " observateur n'est identifiable que parce que l'étoile s'y déplace."
        " Deux tranches de 14 nuits de TOI-4552 au même SNR, l'une couvrant"
        " 42,8 km/s et l'autre 0,7, passent d'un gain d'un facteur deux à une"
        " perte d'un facteur 2,4. La couverture EFFECTIVE ne compte que les bins"
        " qui contiennent une pose : une campagne observée à deux extrêmes et"
        " nulle part entre les deux n'est pas créditée de l'intervalle. Le"
        " troisième nombre est ce que le ciel AUTORISE : le |BERV| d'une cible"
        " ne dépasse jamais 29,78 cos(latitude écliptique) km/s, donc TOI-1452,"
        " à +80,5°, ne pourra jamais couvrir plus de 9,8 km/s et aucune quantité"
        " d'observations n'y séparera les deux référentiels. Une couverture"
        " inférieure à l'étendue signale des trous, que des nuits combleront ;"
        " une étendue inférieure au possible, une campagne jeune ; un possible"
        " faible, une cible mal choisie.",
    "run": "Lancer", "stop": "Arrêter", "dry": "Essai à blanc",
    "export": "Exporter le YAML...", "savelog": "Enregistrer le journal...",
    "openout": "Ouvrir les sorties",
    "savedefaults": "Enregistrer comme défauts...", "lblwin": "Réglages LBL...",
    "all": "tout", "none": "rien",
    "idle": "au repos", "running": "en cours", "lang": "English",
    "lbl_title": "réglages LBL", "close": "Fermer",
    "nothing_export": "tous les réglages sont ceux de la configuration : un"
                      " fichier de variante ne dirait rien.",
    "export_title": "enregistrer ces réglages comme variante",
    "mixed_title": "deux instruments",
    "mixed": "%s ne viennent pas du même instrument (%s).\n\nUn passage, c'est un"
             " seul domaine, une seule grille et un seul jeu d'extensions, tous"
             " lus dans l'instrument : des objets de deux spectrographes ne"
             " peuvent pas être ajustés ensemble. Cochez ceux d'un seul"
             " instrument.",
    "defaults_title": "enregistrer comme défauts",
    "defaults_ask": "Écrire ces valeurs dans %s comme défauts de tous les"
                    " passages ?\n\n%s\n\nLes commentaires du fichier sont"
                    " conservés.",
    "log_index": "index de ce dossier de données : %s",
    "log_out_proposed": "racine de sortie proposée, à côté des données : %s",
    "log_scan_start": "lecture du dossier de données %s",
    "log_scan_done": "%s : %d objets, %d spectres lus, %d déjà connus, %d disparus",
    "log_scan_none": "aucun dossier contenant des spectres sous %s",
    "log_busy": "lecture du dossier de données en cours : attendre la fin",
    "log_no_root": "ce n'est pas un dossier : %s",
    "log_no_object": "aucun objet coché : en cocher au moins un dans la liste",
    "log_mixed": "deux instruments cochés (%s) : un passage, c'est un instrument",
    "log_mixed_refused": "passage refusé : %s viennent de deux instruments (%s)",
    "log_nights": "%s : %d nuits en commun, sur %s",
    "log_nights_thin":
        "%s : %d nuits en commun sur %s. En avoir peu est un AVANTAGE ici : les"
        " composantes observateur sont un parasite toujours présent, ce qu'on"
        " mesure est son motif, et des nuits qu'aucune autre étoile n'a vues"
        " élargissent la gamme de conditions sur laquelle ce motif est mesuré.",
    "log_command": "lancement : %s",
    "log_ended": "passage terminé, code de sortie %d",
    "log_stopping": "demande d'arrêt du passage",
    "log_export": "variante écrite : %s",
    "log_defaults": "défauts écrits dans %s : %s",
    "log_nothing": "rien à écrire : tous les réglages sont ceux de la configuration",
    "log_saved_log": "journal écrit : %s",
    "log_no_outputs": "rien n'y est encore écrit : %s",
    "log_failed": "impossible de le lancer : %s",
    "log_no_command": "rien n'a été installé, il n'y a donc rien à lancer",
    "log_installing": "installation de la commande dans cet environnement, depuis %s",
    "no_command_title": "la commande n'est pas installée",
    "no_command":
        "pca2d-preclean est introuvable pour\n\n  %s\n\nL'installer"
        " maintenant depuis %s ?\n\nCela lance `pip install -e . --no-deps`"
        " là-bas, ce qui place la commande dans le bin de cet interpréteur et"
        " ne touche à rien d'autre.",
    "no_command_still": "toujours introuvable après l'installation",
    "opt_n_star": "composantes stellaires",
    "opt_n_earth": "composantes observateur",
    "opt_mean": "partie statique",
    "opt_velocity_term": "ajuster une vitesse par pose",
    "opt_iters": "itérations au plus",
    "opt_shrink": "ne diviser que le significatif",
    "opt_mask": "échantillons blanchis",
    "opt_width_kms": "passe-haut (km/s)", "opt_dv": "pas de grille (km/s)",
    "opt_nightly_stack": "empiler chaque nuit", "opt_run": "lancer LBL (heures)",
    "opt_lbl_prepare": "écrire l'arbre du LBL",
    "opt_lbl_before": "mesurer les spectres livrés",
    "opt_lbl_after": "mesurer les spectres corrigés",
    "opt_lbl_star_template": "notre étoile comme gabarit du LBL",
    "opt_lbl_strpca": "composantes stellaires en RESPROJ",
    "opt_lbl_directory": "dossier de données du LBL",
    "opt_lbl_suffix": "nom du corrigé",
    "opt_lbl_teff": "température effective", "opt_lbl_template": "fichier gabarit",
    "opt_lbl_steps": "étapes", "opt_lbl_link": "spectres vers LBL en",
    "help_lbl_prepare":
        "Écrire lbl_config.yaml et run_lbl.py à côté des sorties du passage et"
        " déposer les deux jeux de spectres dans les dossiers science du LBL."
        " Désactivé, la correction est quand même écrite et rien n'est préparé"
        " pour une mesure de vitesse.",
    "help_lbl_before":
        "Mesurer aussi les spectres LIVRÉS, comme objet distinct dans le même"
        " arbre LBL. Désactivé, la correction est mesurée contre rien et le"
        " résultat ne peut pas se lire comme meilleur ou pire que ne rien faire.",
    "help_lbl_after":
        "Mesurer les spectres CORRIGÉS. Désactivé avec `before` activé, le LBL ne"
        " mesure que les livrés : c'est ainsi qu'une série de référence est"
        " construite une fois puis réutilisée par tous les passages suivants.",
    "help_lbl_star_template":
        "Confier au LBL le spectre stellaire que cet ajustement a construit, au"
        " lieu de le laisser bâtir le sien à partir des spectres corrigés."
        " DÉSACTIVÉ est le nominal : mesuré sur Proxima, notre gabarit faisait"
        " ajuster au LBL des raies 20 % plus larges et doublait l'erreur par pose,"
        " de 0,97 à 1,67 m/s, pour 3,34 m/s de dispersion contre 3,04. Tout ce qui"
        " manque à notre gabarit, le LBL le paie raie par raie.",
    "help_lbl_strpca":
        "À partir de deux composantes stellaires, celles après la première sont"
        " confiées au LBL comme tables RESPROJ, comme le sont ses propres"
        " gradients DTEMP, pour que le rdb porte la projection de chacune par pose"
        " et qu'une corrélation puisse être cherchée au lieu d'être supposée"
        " absente.",
    "help_lbl_directory":
        "Le DATA_DIR du LBL : son arbre à lui, partagé par tous les objets et"
        " tous les passages, donc placé à côté des sorties plutôt que dans le"
        " dossier d'un passage. Pointé sur un arbre LBL existant, cet arbre est"
        " utilisé tel quel.",
    "help_lbl_suffix":
        "Comment s'appelle l'objet corrigé à côté du livré : TOI-2120 et"
        " TOI-2120_PCA2D_2-7. `{tag}` est le nombre de composantes du passage, et"
        " l'omettre fait écrire à deux passages leurs spectres corrigés dans UN"
        " seul dossier, où le LBL mesure le mélange sans un mot.",
    "help_lbl_teff":
        "La température effective annoncée au LBL, qui lui sert à choisir sa"
        " liste de raies. `auto` la lit dans l'en-tête, un nombre l'impose.",
    "help_lbl_template":
        "Un FICHIER gabarit à imposer au LBL, par son chemin, au lieu de celui"
        " qu'il construirait. Vide est le nominal : le LBL bâtit un gabarit à"
        " partir des spectres qu'il mesure, et c'est ce gabarit que chaque jeu"
        " mesure.",
    "help_lbl_steps":
        "Lesquelles des étapes du LBL lancer, dans l'ordre : template, mask,"
        " compute, compile. En retirer sert à reprendre un arbre qui possède déjà"
        " les précédentes, jamais à sauter un travail dont une étape suivante a"
        " besoin.",
    "help_lbl_link":
        "Comment les spectres arrivent dans les dossiers science du LBL."
        " `symlink` est le nominal : une campagne pèse des dizaines de"
        " gigaoctets et le LBL ne fait que les lire. `copy` est pour un système"
        " de fichiers incapable de porter un lien, un disque exFAT avant tout.",
    "help_col_snr":
        "La médiane, sur les poses de cette cible, du SNR d'extraction par ordre"
        " que le pipeline a écrit dans chaque fichier : une cible brillante et une"
        " faible se distinguent avant de rien lancer. Lu dans les seuls en-têtes"
        " et mémorisé, donc un dossier s'ouvre aussitôt la deuxième fois.",
    "help_col_exptime":
        "Le temps de pose médian des fichiers de cette cible, en secondes. Avec le"
        " SNR et le nombre de fichiers, il dit quel genre de campagne c'est :"
        " beaucoup de poses courtes d'une étoile brillante, ou peu de longues"
        " d'une faible.",
    "help_col_mag":
        "Cliquer un en-tête de colonne pour trier dessus, encore pour inverser,"
        " et deux fois sur l'en-tête des objets pour revenir à instrument puis"
        " nom. Trier sur CELLE-CI compare des bandes différentes : NIRPS écrit J"
        " et SPIRou H, donc un J de 5,3 et un H de 10,5 sont ordonnés comme des"
        " nombres, pas comme des éclats.\n\n"
        "L'éclat de la cible tel que SON pipeline l'a noté, avec la bande où il"
        " est mesuré : NIRPS écrit la magnitude J, SPIRou écrit H, et les deux"
        " diffèrent d'environ une magnitude sur une naine M ; la bande est donc"
        " affichée plutôt que supposée. Lu dans les en-têtes, jamais dans un"
        " catalogue.",
    "help_check":
        "Cocher une cible pour la traiter. Plusieurs cochées sont ajustées"
        " ENSEMBLE contre une seule base observateur. Cliquer la case,"
        " double-cliquer la ligne, ou appuyer sur la barre d'espace. Elles doivent"
        " toutes venir du même instrument.",
    "help_all_button":
        "Cocher toutes les cibles de la liste, ou les décocher toutes.",
    "help_instrument_filter":
        "Afficher ou masquer les cibles d'un instrument. Un passage, c'est UN"
        " instrument : masquer les autres est le chemin le plus court vers une"
        " sélection qui peut vraiment être lancée, et les lignes sont groupées"
        " et teintées par instrument pour la même raison. Masquer un instrument"
        " décoche ses cibles, et rien n'est effacé : recochez la case et elles"
        " reviennent telles quelles.",
    "help_savedefaults_button":
        "Écrire les réglages qui DIFFÈRENT de la configuration dans config.yaml"
        " lui-même, comme défauts dont partira tout passage ultérieur. Les"
        " commentaires du fichier sont conservés, puisqu'ils sont les mesures qui"
        " ont choisi chaque valeur. Un fichier de variante laisse le nominal"
        " intact ; ceci le change.",
    "help_lblwin_button":
        "Le bloc LBL dans sa propre fenêtre : ce qui est mesuré, comment l'objet"
        " corrigé est nommé, contre quel gabarit il est mesuré, lesquelles des"
        " étapes du LBL tournent. C'est l'étape qui produit les vitesses.",
    "help_data_dir":
        "La RACINE des données, pas le dossier d'un objet : les objets listés"
        " en dessous en sont les sous-dossiers. Rien n'y est jamais écrit.",
    "help_config":
        "config.yaml : tout ce que la fenêtre ne montre pas. Trois couches y"
        " sont fusionnées, ce qui est général, ce qui appartient au"
        " spectrographe, et ce qui appartient à une cible.",
    "help_out_dir":
        "Où le passage écrit. Proposé comme `corrected` À CÔTÉ du dossier de"
        " données, puisque les spectres corrigés sont une copie de la campagne,"
        " des dizaines de gigaoctets, et ont leur place sur le disque où la"
        " campagne est déjà. Vide, c'est la racine de sortie de la configuration. Un passage y dépose sa configuration résolue, son"
        " rapport, ses spectres corrigés et ses dossiers LBL, sous"
        " <racine>/<objet>/<M>-<N>/.",
    "help_rescan":
        "Relire le dossier de données. Ce qui a déjà été lu est mémorisé dans un"
        " index sous votre dossier personnel, jamais dans le dossier de données,"
        " donc seuls les fichiers ajoutés ou remplacés sont relus : à utiliser"
        " après y avoir copié des spectres.",
    "help_objects":
        "Un objet coché : passage solo. PLUSIEURS : ils sont ajustés ensemble contre"
        " UNE seule base observateur, chacun gardant son spectre stellaire par"
        " parité d'ordre. L'atmosphère et l'instrument sont communs, les étoiles"
        " non, donc une base ajustée sur plusieurs étoiles ne peut suivre aucune"
        " d'elles. Ils doivent venir du même instrument.",
    "help_n_star":
        "Composantes ajustées dans le référentiel de l'ÉTOILE, en plus du"
        " spectre stellaire lui-même, toujours retiré avec un coefficient de 1"
        " exactement. 0 est le nominal : une amplitude libre devant un terme en"
        " logarithme est un exposant sur le flux, et un exposant sur le spectre"
        " moyen d'une étoile ne décrit aucune étoile. La retirer a fait passer"
        " TOI-4552 de 14,6 à 11,2 m/s en robuste, et TOI-2120 de 14,4 à 14,1.",
    "help_n_earth":
        "Composantes ajustées dans le référentiel de l'OBSERVATEUR :"
        " l'atmosphère et l'instrument. Trois est le nominal. Plus de liberté"
        " décrit mieux le ciel et emporte davantage d'étoile partout où les deux"
        " référentiels sont dégénérés, ce que produit une couverture en BERV"
        " étroite.",
    "help_mean":
        "La partie du modèle qui n'a pas d'amplitude propre : le spectre moyen."
        " Il y en a un par PARITÉ D'ORDRE, puisque les ordres pairs et impairs"
        " voient une longueur d'onde à des résolutions différentes, et la"
        " question est dans quel référentiel il vit. `star`, le nominal : un"
        " spectre par parité dans le référentiel de l'ÉTOILE, une médiane"
        " groupée en BERV retirée une fois avant toute composante avec un"
        " coefficient exactement égal à 1, et aucune moyenne en référentiel"
        " observateur ; la correction ne divise donc que le bloc observateur."
        " `offset` : pas de spectre stellaire, et seulement la part de la"
        " moyenne observateur qui DIFFÈRE entre parités, la part commune restant"
        " dans les données pour le bloc observateur ; elle revient dans la"
        " correction. `full` : toute la moyenne observateur par parité, part"
        " commune comprise. `iterate` : une moyenne par parité dans CHAQUE"
        " référentiel, réestimée à chaque itération. Mesuré : TOI-2120 20,4 m/s"
        " avec `star` contre 31,7 avec `offset` ; sur Proxima la moyenne"
        " observateur seule injectait 46 m/s et le bloc observateur seul 48,"
        " leur somme 19. `iterate` converge sur des données synthétiques et n'a"
        " pas convergé sur une campagne entière.",
    "help_velocity_term":
        "Ajuste une vitesse par pose à côté des composantes, pour garder le"
        " mouvement propre de l'étoile hors du bloc observateur. Elle est"
        " ajustée et écrite dans les fichiers corrigés, jamais divisée. Mesurée"
        " sur TOI-4552, elle n'a rien changé : 16,4 contre 15,9 m/s.",
    "help_iters":
        "Nombre maximal d'itérations. Une itération résout les amplitudes de"
        " chaque pose, puis met à jour chaque base contre le résidu de l'autre."
        " L'ajustement garde la meilleure itération et s'arrête quand le chi2"
        " s'est clairement retourné : c'est donc un plafond, pas une durée.",
    "help_shrink":
        "Ne divise chaque composante observateur que là où les données la"
        " détectent, colonne par colonne, toutes les poses comptées ensemble :"
        " un motif à un demi-sigma dans chacune des N poses est détecté à"
        " environ 0,5 racine(N) et il est gardé, tandis qu'une colonne où la"
        " composante n'est que du bruit est laissée intacte. Sans elle,"
        " TOI-4552 perdait 3 m/s de plus : 18,1 contre 15,1.",
    "help_mask":
        "Un échantillon auquel l'ajustement n'a donné aucun poids ne peut pas"
        " être corrigé, il est donc blanchi. `common`, le nominal, blanchit dans"
        " toutes les poses tout échantillon qu'une pose a perdu : la campagne"
        " porte alors UN seul jeu de raies. `exposure` blanchit celui de chaque"
        " pose, ce qui a coûté 1,8 m/s sur TOI-4552 sans rien diviser du tout :"
        " un jeu de raies qui bouge d'une époque à l'autre est de la dispersion."
        " `none` y laisse le flux livré, non corrigé.",
    "help_width_kms":
        "La largeur du filtre de Savitzky-Golay qui retire le continu, en km/s,"
        " pour qu'il traite une raie de la même façon quel que soit le pas de"
        " grille. Elle fait partie du cube : la changer en construit un autre.",
    "help_dv":
        "Le pas de la grille, en km/s. La grille est uniforme en logarithme de"
        " longueur d'onde, ce qui fait d'un décalage Doppler une translation"
        " exacte. Plus fin est plus lourd : le cube et l'ajustement grandissent"
        " avec. La changer construit un autre cube.",
    "help_nightly_stack":
        "Coaddition des poses d'une nuit avant l'ajustement. C'est une décision"
        " de mémoire, jamais de modélisation : chaque pose est d'abord recalée"
        " sur sa propre BERV, donc rien n'est étalé, et la correction résout"
        " toujours les amplitudes propres de chaque pose contre la base fixée."
        " `auto` compare l'empreinte de l'ajustement à la mémoire permise.",
    "help_run":
        "Lancer LBL une fois les spectres corrigés. C'est l'étape qui produit"
        " les vitesses, et elle dure des heures. Désactivée, tout ce dont LBL a"
        " besoin est quand même écrit, et le passage dit comment le lancer.",
    "help_stage_cube":
        "Lit chaque spectre une fois sur la grille commune, avec ses poids."
        " Quelques minutes, et mis en cache : un second passage aux mêmes"
        " réglages le réutilise.",
    "help_stage_fit":
        "La décomposition à deux référentiels. La longue : elle affiche un R2"
        " après chaque itération.",
    "help_stage_figures":
        "Un PDF : la campagne en fonction du temps, les paramètres, puis une"
        " page par fenêtre de longueur d'onde montrant chaque étape du modèle et"
        " ce que la correction divise.",
    "help_stage_correct":
        "Écrit les t.fits corrigés, un par pose, le bloc observateur divisé et"
        " les échantillons sans poids blanchis.",
    "help_stage_lbl":
        "Confie les deux jeux de spectres à LBL, les livrés et les corrigés,"
        " comme deux objets d'un même arbre, pour que les vitesses se comparent"
        " au lieu de se croire.",
    "help_variant":
        "Un fichier de variants/ : la configuration nominale plus les quelques"
        " lignes que cette variante change. Ce que chacune a donné est écrit"
        " dans variants/README.md, et un passage fait avec l'une d'elles est"
        " reproductible par son nom.",
    "help_command":
        "Exactement ce que le bouton Lancer exécutera. Copiée dans un terminal,"
        " elle fait la même chose : la fenêtre est une enveloppe, pas une"
        " seconde façon de faire.",
    "help_run_button":
        "Lance la commande ci-dessus. La sortie apparaît en dessous au fil de"
        " l'écriture, et le passage peut être arrêté.",
    "help_stop_button":
        "Demande l'arrêt. Ce qui est déjà écrit reste : un cube, un ajustement"
        " et des fichiers corrigés se réutilisent, et une étape inachevée"
        " recommencera simplement la prochaine fois.",
    "help_dry_button":
        "Résout tout, affiche le plan, ne touche à rien. La façon honnête de"
        " voir quel cube un jeu d'options désigne avant d'y passer une heure.",
    "help_export_button":
        "Écrit les réglages qui DIFFÈRENT de la configuration dans un fichier de"
        " variante, ce qui rend un passage reproductible. Un fichier qui répète"
        " le nominal ne dit rien, et il est refusé.",
    "help_savelog_button": "Écrit dans un fichier ce que la fenêtre a montré.",
    "help_openout_button":
        "Ouvre la racine de sortie dans le navigateur de fichiers.",
    "help_lang": "Bascule la fenêtre entre l'anglais et le français.",
    "help_apero":
        "APERO, le pipeline qui a réduit tous les spectres que cette fenêtre"
        " lit (Cook et al. 2022, PASP 134, 114509). C'est lui qui a écrit les"
        " extensions, la solution en longueur d'onde, le SNR par ordre et les"
        " vitesses barycentriques employées ici ; ce paquet part de ses t.fits"
        " et n'en redérive rien.",
    "help_log":
        "Le journal du passage, dans les couleurs qu'un terminal lui donnerait :"
        " vert pour la progression, bleu pour une valeur, orange pour ce qui est"
        " sauté, rouge pour ce qui arrête un passage.",
}

TEXTS = {"en": EN, "fr": FR}


def text(lang, key, default=None):
    """The label or the explanation, in the window's language."""
    return TEXTS.get(lang, EN).get(key, default if default is not None else key)


def objects_in(root, pattern="*t.fits"):
    """[(name, number of files)] for every object folder under the data root.

    The folder listing alone, with nothing opened: what the list can show
    before the index has been read. The numbers beside each name (instrument,
    SNR, exposure time) come from pca2d.scan, which remembers them.
    """
    return [(name, len(files))
            for name, files in scan.objects_of(root, pattern).items()]


def absolute(path):
    """A path as it will be read, in full: `~` expanded, relative made absolute.

    The window showed `data` and `config.yaml` as written, which say nothing
    about WHERE a run will read from: the same two words mean a different folder
    from a different working directory, and the data here live on a shared disk
    reached through a link. What is shown is now what is used.
    """
    if not path:
        return ""
    return os.path.abspath(os.path.expanduser(str(path)))


def corrected_dir(data_root):
    """Where to propose putting the run's products, given where it reads.

    Beside the data and named `corrected`: the corrected spectra are a copy of
    the campaign, tens of gigabytes, and they belong on the disk the campaign is
    already on rather than on whatever disk the window happened to start from.
    A root already called `corrected` is left alone rather than nested inside
    itself.
    """
    full = absolute(data_root)
    if not full:
        return ""
    parent, name = os.path.split(full.rstrip(os.sep))
    if name.lower() == "corrected":
        return full
    return os.path.join(parent or os.sep, "corrected")


def preclean_argv():
    """How to run `pca2d-preclean` from HERE, or None if it cannot be found.

    Beside the interpreter first: the window is started by the entry point
    installed in an environment's bin, and its sibling is the command, whatever
    PATH the window inherited. A window started from another shell had none
    ("No such file or directory: 'pca2d-preclean'", 2026-09-13) although the
    command was installed all along. Then PATH, then the module, which works
    wherever the package imports.
    """
    beside = os.path.join(os.path.dirname(os.path.abspath(sys.executable)),
                          "pca2d-preclean")
    if os.path.exists(beside) and os.access(beside, os.X_OK):
        return [beside]
    found = shutil.which("pca2d-preclean")
    if found:
        return [found]
    try:
        import pca2d.cli                                        # noqa: F401
    except Exception:                                           # noqa: BLE001
        return None
    return [sys.executable, "-m", "pca2d.cli"]


def instruments_of(rows, names):
    """The instruments the named objects were read from, the unknown aside.

    More than one of them stops a run: the domain, the grid and the extensions
    to read all come from the instrument, so objects from two spectrographs
    have no common grid to be fitted on. An object whose instrument could not
    be read says '?' and is not counted as a second one, since it is not
    evidence of anything.
    """
    return {(rows.get(name) or {}).get("instrument") for name in names} - \
        {"?", "", None}


def build_command(state):
    """The `pca2d-preclean` command a set of choices means, as a list.

    One object or several, the stages asked for, the variant if one is chosen,
    and the two counts, which are flags rather than config keys so that a run
    says in its own name what it fitted.
    """
    argv = ["pca2d-preclean"]
    objects = list(state.get("objects") or [])
    if len(objects) > 1:
        argv += ["--objects", ",".join(objects)]
    elif objects:
        argv += ["--object", objects[0]]
    if state.get("config"):
        argv += ["--config", state["config"]]
    if state.get("data_dir"):
        argv += ["--data-dir", state["data_dir"]]
    if state.get("out_dir"):
        argv += ["--out-dir", state["out_dir"]]
    if str(state.get("run_name") or "").strip():
        argv += ["--name", str(state["run_name"]).strip()]
    for flag, key in (("--min-rjd", "min_rjd"), ("--max-rjd", "max_rjd")):
        value = str(state.get(key) or "").strip()
        if value:
            argv += [flag, value]
    if state.get("variant") and state["variant"] != "(none)":
        argv += ["--variant", state["variant"]]
    if state.get("n_star") not in (None, ""):
        argv += ["--n-star", str(state["n_star"])]
    if state.get("n_earth") not in (None, ""):
        argv += ["--n-earth", str(state["n_earth"])]
    stages = [s for s in STAGES if state.get("stage_" + s)]
    if stages and len(stages) != len(STAGES):
        argv += ["--stages", ",".join(stages)]
    if state.get("dry_run"):
        argv.append("--dry-run")
    return argv


def variant_yaml(state, defaults=None):
    """The settings that differ from the configuration, as a variant file.

    Only what was changed: a variant is the nominal plus its own lines, and a
    file that repeats the nominal says nothing (variants/README.md).
    """
    defaults = defaults or {}
    out = {}
    for key, path, kind in ALL_OPTIONS:
        if key not in state or state[key] in (None, ""):
            continue
        value = state[key]
        if kind == "int":
            value = int(value)
        elif kind == "float":
            value = float(value)
        elif kind == "bool":
            value = bool(value)
        elif kind == "list":
            value = [part.strip() for part in str(value).replace(",", " ").split()]
        elif kind == "text":
            value = str(value)
        elif value in ("true", "false"):
            value = value == "true"
        section, name = path.split(".")
        if defaults.get(section, {}).get(name, object()) == value:
            continue
        out.setdefault(section, {})[name] = value
    return out


def _read_state():
    try:
        with open(HOME_STATE) as handle:
            return json.load(handle)
    except Exception:                                         # noqa: BLE001
        return {}


def _write_state(state):
    try:
        with open(HOME_STATE, "w") as handle:
            json.dump(state, handle, indent=1)
    except OSError:
        pass


class Tip:
    """What an item means, in a small window, after a moment on it.

    The window, and nothing else. This used to write the same text into the
    status line as well, truncated, which made that line flicker with every
    mouse movement and, worse, wiped the scan's progress from it: leaving an
    item reset the line to "idle" while a scan was still reading files. The
    status line reports what the window is DOING; the floating window is what
    explains an item.
    """

    def __init__(self, app, widget, key):
        #: a string, or a callable returning one: the list's explanation depends
        #: on which column the pointer is over
        self.app, self.widget, self.key = app, widget, key
        self.window = None
        self.after = None
        widget.bind("<Enter>", self.enter, add="+")
        widget.bind("<Leave>", self.leave, add="+")
        widget.bind("<ButtonPress>", self.leave, add="+")

    def says(self):
        return self.app.t(self.key() if callable(self.key) else self.key)

    def enter(self, _event=None):
        self.after = self.widget.after(500, self.show)

    def show(self):
        import tkinter as tk
        if self.window is not None:
            return
        x = self.widget.winfo_rootx() + 20
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 8
        self.window = tk.Toplevel(self.widget)
        self.window.wm_overrideredirect(True)
        self.window.wm_geometry("+%d+%d" % (x, y))
        tk.Label(self.window, text=self.says(), justify="left",
                 background=LOG_BG, foreground=LOG_INK, relief="flat",
                 borderwidth=0, wraplength=460,
                 font=self.app.fonts.get("small", ("Helvetica", 11)),
                 padx=10, pady=8).pack()

    def leave(self, _event=None):
        if self.after is not None:
            self.widget.after_cancel(self.after)
            self.after = None
        if self.window is not None:
            self.window.destroy()
            self.window = None


def chip(canvas, x, y, text, ink, fill, font, pad=5, outline=None):
    """A label on its own ground, centred at (x, y), drawn over whatever is there.

    A warning written straight onto a histogram is read against the bars: at
    the density of a stacked campaign it is barely legible, and it is worst
    exactly where there is most data to warn about. The text goes down first,
    its box is measured from it rather than guessed from a character count, and
    the box is lowered underneath it.
    """
    item = canvas.create_text(x, y, text=text, fill=ink, font=font)
    x0, y0, x1, y1 = canvas.bbox(item)
    box = canvas.create_rectangle(x0 - pad, y0 - pad + 2, x1 + pad, y1 + pad - 2,
                                  fill=fill, outline=outline or fill)
    canvas.tag_lower(box, item)
    return item, box


class App:
    """The window itself."""

    def __init__(self, root):
        import tkinter as tk
        from tkinter import ttk

        self.tk, self.ttk = tk, ttk
        self.root = root
        self.proc = None
        self.lines = queue.Queue()
        self.saved = _read_state()
        self.lang = self.saved.get("lang", "en")
        self.labels = []          # (widget, key, how) to relabel on a switch
        self._cfg = None
        self.names = {}           # tree item -> object name
        self.rows = {}            # object name -> what the index knows of it
        self.checked = set(self.saved.get("checked") or [])
        self.index = {}
        self.scanning = False
        # None = the default order, by instrument then name
        self.sort_column = self.saved.get("sort_column") or None
        self.sort_reverse = bool(self.saved.get("sort_reverse"))
        self.lbl_window = None
        root.title("pca2d-preclean")
        # a 13-inch laptop has about 800 points of usable height, and the log
        # at the bottom is the part that was falling off the screen
        root.geometry("1280x800")
        root.minsize(1000, 620)
        style = ttk.Style()
        # clam rather than the native aqua: aqua ignores most colour options, so
        # a window styled under it stays the grey it was born with
        style.theme_use("clam" if "clam" in style.theme_names() else "default")
        self._style(style)
        root.configure(background=BG)

        self.vars = {}
        self.status = ttk.Label(root, text=self.t("idle"), style="Hint.TLabel")
        self._build_top(root)
        panes = ttk.Panedwindow(root, orient="horizontal")
        panes.pack(fill="x", expand=False, padx=14)
        left, right = ttk.Frame(panes), ttk.Frame(panes)
        panes.add(left, weight=3)
        panes.add(right, weight=4)
        self._build_objects(left)
        self._build_options(right)
        # one panel on each side: stacked on one, they made that column half as
        # tall again as the other and pushed the log off a laptop screen
        self._build_berv(left)
        self._build_time(right)
        self._build_command(root)
        self._build_log(root)
        self._propose_out()      # on opening, not only when the data root moves
        self.refresh_objects()
        self._drain_id = self.root.after(80, self._drain)
        self.root.protocol("WM_DELETE_WINDOW", self._close)

    # ---- the look -----------------------------------------------------
    def _style(self, style):
        """Every widget class, once, so nothing is styled at its call site."""
        body = ("Helvetica Neue" if sys.platform == "darwin" else "Helvetica")
        mono = ("Menlo" if sys.platform == "darwin" else "DejaVu Sans Mono")
        self.fonts = {"body": (body, 12), "small": (body, 11),
                      "head": (body, 15, "bold"), "mono": (mono, 11)}
        style.configure(".", background=BG, foreground=INK,
                        font=self.fonts["body"], borderwidth=0)
        style.configure("TFrame", background=BG)
        style.configure("TLabel", background=BG, foreground=INK)
        style.configure("Head.TLabel", font=self.fonts["head"], foreground=ACCENT)
        style.configure("Hint.TLabel", foreground=MUTED, font=self.fonts["small"])
        style.configure("TLabelframe", background=BG, bordercolor=LINE,
                        relief="solid", borderwidth=1)
        style.configure("TLabelframe.Label", background=BG, foreground=ACCENT,
                        font=(body, 11, "bold"))
        style.configure("TCheckbutton", background=BG, foreground=INK,
                        focuscolor=BG)
        style.map("TCheckbutton", background=[("active", BG)])
        style.configure("TEntry", fieldbackground=SURFACE, background=SURFACE,
                        bordercolor=LINE, lightcolor=LINE, darkcolor=LINE,
                        insertcolor=INK, padding=4)
        style.configure("TCombobox", fieldbackground=SURFACE, background=SURFACE,
                        bordercolor=LINE, arrowcolor=ACCENT, padding=3)
        style.map("TCombobox", fieldbackground=[("readonly", SURFACE)])
        style.configure("TButton", background=SURFACE, foreground=INK,
                        bordercolor=LINE, focuscolor=BG, padding=(10, 5),
                        relief="solid", borderwidth=1)
        style.map("TButton",
                  background=[("pressed", ACCENT_SOFT), ("active", ACCENT_SOFT)],
                  foreground=[("disabled", "#9aa6b8")])
        # the one action this window exists for
        style.configure("Run.TButton", background=ACCENT, foreground="#ffffff",
                        bordercolor=ACCENT, font=(body, 12, "bold"),
                        padding=(16, 6))
        style.map("Run.TButton",
                  background=[("pressed", "#08557f"), ("active", "#0d7cc2"),
                              ("disabled", "#9dbdd4")],
                  foreground=[("disabled", "#eef4f8")])
        style.configure("Treeview", background=SURFACE, fieldbackground=SURFACE,
                        foreground=INK, rowheight=23, bordercolor=LINE,
                        font=self.fonts["body"])
        style.configure("Treeview.Heading", background=BG, foreground=MUTED,
                        font=(body, 11, "bold"), relief="flat", padding=(4, 5))
        style.map("Treeview.Heading", background=[("active", ACCENT_SOFT)],
                  foreground=[("active", ACCENT)])
        style.map("Treeview", background=[("selected", ACCENT_SOFT)],
                  foreground=[("selected", INK)])
        style.configure("TPanedwindow", background=BG)
        style.configure("Vertical.TScrollbar", background=BG, troughcolor=BG,
                        bordercolor=BG, arrowcolor=MUTED)

    # ---- language -----------------------------------------------------
    def t(self, key):
        return text(self.lang, key)

    def _state(self):
        """The status line says what the window is doing, and only that.

        Left alone while a scan is reading files: its progress is the most
        useful thing the line can hold, and a language switch or a finished
        subprocess must not wipe it.
        """
        if not self.scanning:
            self.status.configure(text=self.t("running" if self.proc else "idle"))

    def _register(self, widget, key, how="text"):
        self.labels.append((widget, key, how))
        return widget

    def _tip(self, widget, key):
        Tip(self, widget, key)
        return widget

    def switch_language(self):
        self.lang = "fr" if self.lang == "en" else "en"
        for widget, key, how in self.labels:
            try:
                if how == "text":
                    widget.configure(text=self.t(key))
                elif how == "heading":
                    pass                  # redrawn together, below
            except Exception:                                 # noqa: BLE001
                pass
        self._draw_headings()
        self._state()
        if self.lbl_window is not None:
            self.lbl_window.title(self.t("lbl_title"))
        self._sync()

    # ---- widgets ------------------------------------------------------
    def _build_top(self, parent):
        ttk, tk = self.ttk, self.tk
        frame = ttk.Frame(parent)
        frame.pack(fill="x", padx=14, pady=(12, 8))
        # APERO's own logo, then the name: every spectrum this window reads was
        # reduced by APERO, so the pipeline signs the top left corner. The file
        # is in the package, never fetched at run time.
        head = ttk.Frame(frame)
        head.grid(row=0, column=0, sticky="w")
        logo = self._logo("apero_logo.png", height=22)
        if logo is not None:
            self.logo_label = ttk.Label(head, image=logo, background=BG)
            self.logo_label.image = logo      # or the garbage collector eats it
            self.logo_label.pack(side="left", padx=(0, 8))
            self._tip(self.logo_label, "help_apero")
        ttk.Label(head, text="pca2d", style="Head.TLabel").pack(side="left")
        self._register(ttk.Label(frame, style="Hint.TLabel", wraplength=640,
                                 justify="left",
                                 text=self.t("subtitle")), "subtitle").grid(
            row=0, column=1, columnspan=2, sticky="w", padx=10)
        button = ttk.Button(frame, text=self.t("lang"),
                            command=self.switch_language, width=10)
        button.grid(row=0, column=3, sticky="e")
        self._register(button, "lang")
        self._tip(button, "help_lang")
        # shown in full, since a relative path means a different folder from a
        # different working directory and these data are reached through a link
        for i, (key, default) in enumerate((
                ("data_dir", absolute(self.saved.get("data_dir", "data"))),
                ("config", absolute(self.saved.get("config", "config.yaml"))),
                ("out_dir", self.saved.get("out_dir", "")))):
            label = ttk.Label(frame, text=self.t(key))
            label.grid(row=i + 1, column=0, sticky="w", pady=2)
            self._register(label, key)
            var = tk.StringVar(value=default)
            self.vars[key] = var
            entry = ttk.Entry(frame, textvariable=var, width=70)
            entry.grid(row=i + 1, column=1, sticky="we", padx=6)
            self._tip(entry, "help_" + key)
            self._tip(label, "help_" + key)
            var.trace_add("write", lambda *_: self._sync())
            if key == "data_dir":
                # proposed when the data root is settled, not at every keystroke
                entry.bind("<FocusOut>", lambda _e: self._propose_out())
                entry.bind("<Return>", lambda _e: self._propose_out())
            browse = ttk.Button(frame, text=self.t("browse"),
                                command=lambda k=key: self._browse(k))
            browse.grid(row=i + 1, column=2)
            self._register(browse, "browse")
        # the run's name, proposed and editable, with what it would overwrite
        row = ttk.Frame(frame)
        row.grid(row=4, column=0, columnspan=4, sticky="we", pady=(6, 0))
        label = ttk.Label(row, text=self.t("run_name"))
        label.pack(side="left")
        self._register(label, "run_name")
        self._tip(label, "help_run_name")
        self.vars["run_name"] = tk.StringVar(value=self.saved.get("run_name", ""))
        entry = ttk.Entry(row, textvariable=self.vars["run_name"], width=26)
        entry.pack(side="left", padx=(6, 14))
        self._tip(entry, "help_run_name")
        self.vars["run_name"].trace_add("write", lambda *_: self._sync())
        for key in ("min_rjd", "max_rjd"):
            self.vars[key] = tk.StringVar(value=self.saved.get(key, ""))
            self.vars[key].trace_add("write", lambda *_: self._sync())
        self.exists = ttk.Label(row, style="Hint.TLabel", text="")
        self.exists.pack(side="left", padx=(6, 0))
        rescan = ttk.Button(frame, text=self.t("rescan"),
                            command=self.refresh_objects)
        rescan.grid(row=1, column=3, padx=4)
        self._register(rescan, "rescan")
        self._tip(rescan, "help_rescan")
        frame.columnconfigure(1, weight=1)

    def _build_objects(self, parent):
        ttk = self.ttk
        box = ttk.Labelframe(parent, text=self.t("objects"))
        box.pack(fill="both", expand=True, pady=4)
        self._register(box, "objects")
        self.tree = ttk.Treeview(box, columns=("files", "snr", "exptime", "mag",
                                               "instrument"),
                                 show="tree headings", selectmode="extended",
                                 height=7)
        self.headings = (("#0", "col_object"), ("files", "col_files"),
                         ("snr", "col_snr"), ("exptime", "col_exptime"),
                         ("mag", "col_mag"), ("instrument", "col_instrument"))
        for column, _key in self.headings:
            self.tree.heading(column, command=lambda c=column: self._sort_by(c))
        self._draw_headings()
        self.tree.column("#0", width=180)
        self.tree.column("files", width=52, anchor="e")
        self.tree.column("snr", width=56, anchor="e")
        self.tree.column("exptime", width=62, anchor="e")
        self.tree.column("mag", width=66, anchor="e")
        self.tree.column("instrument", width=80)
        self.tree.pack(fill="both", expand=True, padx=6, pady=(6, 2))
        # the box is in the first column, so a click on it toggles and a click
        # on the name still selects the row the usual way
        self.tree.bind("<Button-1>", self._clicked)
        self.tree.bind("<Double-1>", self._double)
        self.tree.bind("<space>", lambda _e: self._toggle(self.tree.selection()))
        self.tree.bind("<Motion>", self._hover_column, add="+")
        # the list holds five different things, so what it explains follows the
        # column the pointer is over rather than being one text for all of them
        self._tip(self.tree, self._tree_key)
        bar = ttk.Frame(box)
        bar.pack(fill="x", padx=6, pady=(0, 6))
        for key, value in (("all", True), ("none", False)):
            button = ttk.Button(bar, text=self.t(key), width=6,
                                command=lambda v=value: self._check_all(v))
            button.pack(side="left", padx=(0, 4))
            self._register(button, key)
            self._tip(button, "help_all_button")
        # one box per instrument found, filled in when the root is read: a run
        # is one instrument, so hiding the others is the first thing anybody
        # does before picking targets
        self.filters = ttk.Frame(bar)
        self.filters.pack(side="left", padx=(10, 0))
        self.instrument_vars = {}
        self.count = ttk.Label(bar, style="Hint.TLabel", text="")
        self.count.pack(side="right")

    def _build_berv(self, parent):
        """The barycentric coverage of what is ticked, as a histogram.

        The quantity that decides whether the correction helps at all, drawn
        where the choice is made. Plain canvas rather than a plotting library:
        it is redrawn on every tick and must cost nothing.
        """
        ttk = self.ttk
        box = ttk.Labelframe(parent, text=self.t("berv"))
        box.pack(fill="x", padx=6, pady=(0, 6))
        self._register(box, "berv")
        self.berv_canvas = self.tk.Canvas(box, height=86, highlightthickness=0,
                                          background=SURFACE)
        self.berv_canvas.pack(fill="x", padx=6, pady=(4, 2))
        self.berv_note = ttk.Label(box, style="Hint.TLabel", text="")
        self.berv_note.pack(anchor="w", padx=8, pady=(0, 4))
        self._tip(self.berv_canvas, "help_berv")
        self._tip(self.berv_note, "help_berv")
        self.berv_canvas.bind("<Configure>", lambda _e: self._draw_berv())

    def _build_time(self, parent):
        """When the ticked campaigns were observed, and what to keep of them.

        A reduced Julian date says nothing to anybody, so the dates are drawn
        and read in the calendar, and chosen with two sliders rather than typed.
        """
        ttk = self.ttk
        box = ttk.Labelframe(parent, text=self.t("timeline"))
        box.pack(fill="x", padx=6, pady=(0, 6))
        self._register(box, "timeline")
        self.time_canvas = self.tk.Canvas(box, height=80, highlightthickness=0,
                                          background=SURFACE)
        self.time_canvas.pack(fill="x", padx=6, pady=(4, 2))
        self.time_canvas.bind("<Configure>", lambda _e: self._draw_time())
        self._tip(self.time_canvas, "help_dates")
        sliders = ttk.Frame(box)
        sliders.pack(fill="x", padx=6, pady=(0, 2))
        self.scales = {}
        for key in ("min_rjd", "max_rjd"):
            label = ttk.Label(sliders, text=self.t(key), width=4)
            label.pack(side="left")
            self._register(label, key)
            # each slider starts at ITS end of the campaign, so the pair opens
            # on the whole of it. Both starting at 0 put the upper bound on the
            # first night, which reads as a window holding nothing.
            scale = ttk.Scale(sliders, from_=0.0, to=1.0,
                              value=0.0 if key == "min_rjd" else 1.0,
                              command=lambda v, k=key: self._slide(k, v))
            scale.pack(side="left", fill="x", expand=True, padx=(4, 12))
            self._tip(scale, "help_dates")
            self.scales[key] = scale
        self.time_note = ttk.Label(box, style="Hint.TLabel", text="")
        self.time_note.pack(anchor="w", padx=8, pady=(0, 4))
        self._tip(self.time_note, "help_dates")
        reset = ttk.Button(sliders, text=self.t("whole"), width=8,
                           command=self._reset_dates)
        reset.pack(side="left")
        self._register(reset, "whole")
        self._tip(reset, "help_dates")

    def _reset_dates(self):
        """Back to the whole campaign: no window, the nominal path."""
        self._sliding = True
        for key, scale in getattr(self, "scales", {}).items():
            scale.set(0.0 if key == "min_rjd" else 1.0)
            self.vars[key].set("")
        self._sliding = False
        self._draw_time()
        self._sync()

    def _slide(self, key, value):
        """A slider moved: 0..1 of the span the ticked campaigns cover."""
        if getattr(self, "_sliding", False):
            return
        span = getattr(self, "_time_span", None)
        if not span:
            return
        lo, hi = span
        rjd = lo + float(value) * (hi - lo)
        # the ends mean "no bound", so the nominal path stays reachable
        at_end = (key == "min_rjd" and float(value) <= 0.001) or \
                 (key == "max_rjd" and float(value) >= 0.999)
        self.vars[key].set("" if at_end else "%.2f" % rjd)
        self._draw_time()

    @staticmethod
    def handle_fraction(text, lo, hi, key):
        """Where a slider handle belongs for the bound `text`, as 0..1.

        An empty bound means "no bound", which is that slider's OWN end: 0 for
        the lower one and 1 for the upper. Both at 0 is an upper bound on the
        first night, a window holding nothing.
        """
        default = 0.0 if key == "min_rjd" else 1.0
        try:
            rjd = float(str(text).strip())
        except (TypeError, ValueError):
            return default
        if hi <= lo:
            return default
        return min(1.0, max(0.0, (rjd - lo) / (hi - lo)))

    def _place_handles(self, lo, hi):
        """Put the handles where the kept bounds are.

        The bounds are remembered between sessions and the handles were not, so
        a saved window came back as a date in the field and a handle at the end
        of the campaign, each contradicting the other.
        """
        if getattr(self, "_sliding", False):
            return
        self._sliding = True
        try:
            for key, scale in getattr(self, "scales", {}).items():
                scale.set(self.handle_fraction(self.vars[key].get(), lo, hi, key))
        finally:
            self._sliding = False

    def _draw_time(self):
        """One dot per exposure, in its star's colour, and the window kept."""
        canvas = getattr(self, "time_canvas", None)
        if canvas is None:
            return
        canvas.delete("all")
        names = self.picked()
        width = max(int(canvas.winfo_width()), 60)
        height = max(int(canvas.winfo_height()), 40)
        times = scan.exposure_times(self.index, names) if names else {}
        if not times:
            canvas.create_text(width // 2, height // 2,
                               text=self.t("timeline_none"), fill="#888",
                               font=("Helvetica", 10))
            self.time_note.configure(text="")
            self._time_span = None
            return
        every = [t for values in times.values() for t in values]
        lo, hi = min(every), max(every)
        if hi - lo < 1.0:
            hi = lo + 1.0
        self._time_span = (lo, hi)
        self._place_handles(lo, hi)
        pad, top, foot = 10, 8, 22
        rows = max(1, len(times))
        band = (height - top - foot) / rows

        def x_of(rjd):
            return pad + (rjd - lo) / (hi - lo) * (width - 2 * pad)

        # what is kept, as a lit band behind the dots
        keep_lo = self._bound("min_rjd", lo)
        keep_hi = self._bound("max_rjd", hi)
        canvas.create_rectangle(x_of(keep_lo), top - 4, x_of(keep_hi),
                                height - foot + 4, fill=ACCENT_SOFT, outline="")
        for i, (name, values) in enumerate(sorted(times.items())):
            y = top + band * (i + 0.5)
            colour = self._star_colour(name)
            step = max(1, len(values) // 700)     # a canvas is not a plot library
            for t in values[::step]:
                x = x_of(t)
                inside = keep_lo <= t <= keep_hi
                canvas.create_line(x, y - 4, x, y + 4,
                                   fill=colour if inside else "#d6dbe4")
            canvas.create_text(pad, y - 9, text=name, anchor="w", fill=colour,
                               font=("Helvetica", 8))
        # the calendar, which is what anybody reads
        canvas.create_line(pad, height - foot + 6, width - pad, height - foot + 6,
                           fill="#bbb")
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            when = scan.rjd_to_date(lo + frac * (hi - lo))
            x = pad + frac * (width - 2 * pad)
            canvas.create_line(x, height - foot + 6, x, height - foot + 10,
                               fill="#888")
            canvas.create_text(min(max(x, 22), width - 22), height - foot + 17,
                               text=when.strftime("%Y-%m"), fill="#555",
                               font=("Helvetica", 8))
        kept = sum(1 for t in every if keep_lo <= t <= keep_hi)
        self.time_note.configure(
            text=self.t("timeline_note")
            % (scan.rjd_to_date(keep_lo).strftime("%d %b %Y"),
               scan.rjd_to_date(keep_hi).strftime("%d %b %Y"), kept, len(every)))

    def _bound(self, key, default):
        try:
            return float(self.vars[key].get())
        except (TypeError, ValueError):
            return default

    #: what each column of the list is, for the explanation that follows the
    #: pointer: the box, the counts, then the two numbers read from the headers
    COLUMN_HELP = {"#0": "help_check", "#1": "help_objects", "#2": "help_col_snr",
                   "#3": "help_col_exptime", "#4": "help_col_mag",
                   "#5": "help_objects"}

    #: which row value each column sorts on, and whether it is a number
    SORT_KEY = {"#0": ("object", False), "files": ("files", True),
                "snr": ("snr", True), "exptime": ("exptime", True),
                "mag": ("mag", True), "instrument": ("instrument", False)}

    def _draw_headings(self):
        """The column names, with an arrow on the one the list is sorted by."""
        for column, key in self.headings:
            arrow = ""
            if column == getattr(self, "sort_column", None):
                arrow = "  ▼" if self.sort_reverse else "  ▲"
            self.tree.heading(column, text=self.t(key) + arrow)

    def _sort_by(self, column):
        """Sort by this column, and reverse it if it is already the one.

        Clicking the object column a second time goes back to the default,
        which is by instrument and then by name: a run is one instrument, so
        that grouping is what the list is for most of the time.
        """
        if column == getattr(self, "sort_column", None):
            if self.sort_reverse:
                self.sort_column, self.sort_reverse = None, False
            else:
                self.sort_reverse = True
        else:
            self.sort_column, self.sort_reverse = column, False
        self._draw_headings()
        self._fill(list(self.rows.values()))

    def _sorted(self, rows):
        """The rows in the order the list should show them.

        By instrument then name unless a column was clicked. A value that is
        missing sorts LAST either way: a blank is not a small number, and a
        target whose scan has not reached it should not head the list.
        """
        column = getattr(self, "sort_column", None)
        if column is None:
            # an instrument not read yet goes last here too, rather than first
            # because "?" precedes "nirps" in the alphabet
            return sorted(rows, key=lambda r: (
                (0, (r["instrument"] or "").lower())
                if r.get("instrument") and r["instrument"] != "?" else (1, ""),
                r["object"].lower()))
        field, numeric = self.SORT_KEY.get(column, ("object", False))
        reverse = bool(getattr(self, "sort_reverse", False))

        def key(row):
            value = row.get(field)
            if value is None or value == "":
                # last whichever way round, so the sort is reversed on the
                # values and never on what is not known
                return (1, 0.0 if numeric else "")
            if numeric:
                return (0, -float(value) if reverse else float(value))
            text = str(value).lower()
            return (0, text)

        out = sorted(rows, key=key)
        if not numeric and reverse:
            known = [r for r in out if r.get(field) not in (None, "")]
            out = known[::-1] + [r for r in out if r.get(field) in (None, "")]
        return out

    def _tree_key(self):
        return self.COLUMN_HELP.get(getattr(self, "_column", "#0"), "help_check")

    def _hover_column(self, event):
        """Which column the pointer is over, so the floating window explains
        THAT one. Nothing is written anywhere until it opens."""
        self._column = self.tree.identify_column(event.x)

    # ---- the ticks -----------------------------------------------------
    def _clicked(self, event):
        """A click on the box toggles; anywhere else selects, as usual."""
        item = self.tree.identify_row(event.y)
        if item and self.tree.identify_column(event.x) == "#0" and event.x <= 28:
            self._toggle([item])
            return "break"
        return None

    def _double(self, event):
        item = self.tree.identify_row(event.y)
        if item:
            self._toggle([item])
            return "break"
        return None

    def _toggle(self, items):
        for item in items:
            name = self.names.get(item)
            if name is None:
                continue
            if name in self.checked:
                self.checked.discard(name)
            else:
                self.checked.add(name)
            self._draw_check(item, name)
        self._after_ticks()

    def _check_all(self, on):
        for item, name in self.names.items():
            if on:
                self.checked.add(name)
            else:
                self.checked.discard(name)
            self._draw_check(item, name)
        self._after_ticks()

    def _draw_check(self, item, name):
        glyph = CHECKED if name in self.checked else UNCHECKED
        self.tree.item(item, text="%s  %s" % (glyph, name))

    def _after_ticks(self):
        """What the ticks mean, said as soon as they change rather than at Run."""
        names = self.picked()
        self.count.configure(text="%d / %d" % (len(names), len(self.names)))
        instruments = self.instruments(names)
        # said when it becomes true, not at every tick that keeps it true
        if len(instruments) > 1 and instruments != getattr(self, "_warned", None):
            self._say("log_mixed", ", ".join(sorted(instruments)), level="warn")
        self._warned = instruments if len(instruments) > 1 else None
        self._nights(names)
        self._draw_berv()
        self._draw_time()
        self._sync()

    def _nights(self, names):
        """How many nights the ticked campaigns share, said when it changes.

        The criterion that decided the first joint run and the one a list of
        names cannot show: see docs/joint_fit.md.
        """
        if len(names) < 2:
            self._shared = None
            return
        common, counts = scan.shared_nights(self.index, names)
        if not counts or len(counts) < 2:
            return
        key = (tuple(sorted(names)), len(common))
        if key == getattr(self, "_shared", None):
            return
        self._shared = key
        joined = ", ".join(str(c) for c in counts)
        # a fact worth knowing, and when the overlap is small it is good news,
        # not a fault: the observer components are a parasite that is always
        # there, so a pattern measured over a WIDER range of conditions is
        # better determined (the user, 2026-09-13)
        if len(common) < 0.2 * min(counts):
            self._say("log_nights_thin", " + ".join(names), len(common), joined,
                      level="value")
        else:
            self._say("log_nights", " + ".join(names), len(common), joined,
                      level="value")

    def picked(self):
        """The ticked objects, in the order the list shows them."""
        return [self.names[item] for item in self.tree.get_children()
                if self.names.get(item) in self.checked]

    def instruments(self, names):
        return instruments_of(self.rows, names)

    def _build_options(self, parent):
        ttk, tk = self.ttk, self.tk
        box = ttk.Labelframe(parent, text=self.t("settings"))
        box.pack(fill="both", expand=True, pady=4)
        self._register(box, "settings")
        grid = ttk.Frame(box)
        grid.pack(fill="both", expand=True, padx=8, pady=6)
        for i, (key, path, kind) in enumerate(OPTIONS):
            row, col = i % 4, (i // 4) * 3
            label = ttk.Label(grid, text=self.t("opt_" + key))
            label.grid(row=row, column=col, sticky="w", pady=2)
            self._register(label, "opt_" + key)
            default = self.saved.get(key, self._config_default(path, kind))
            if kind == "bool":
                var = tk.BooleanVar(value=bool(default))
                widget = ttk.Checkbutton(grid, variable=var)
            elif isinstance(kind, tuple):
                var = tk.StringVar(value=str(default))
                widget = ttk.Combobox(grid, textvariable=var, values=list(kind),
                                      width=10, state="readonly")
            else:
                var = tk.StringVar(value="" if default is None else str(default))
                widget = ttk.Entry(grid, textvariable=var, width=12)
            widget.grid(row=row, column=col + 1, sticky="w", padx=(6, 18))
            var.trace_add("write", lambda *_: self._sync())
            self.vars[key] = var
            self._tip(widget, "help_" + key)
            self._tip(label, "help_" + key)
        run = ttk.Frame(box)
        run.pack(fill="x", padx=8, pady=(0, 8))
        stages = ttk.Label(run, text=self.t("stages"))
        stages.grid(row=0, column=0, sticky="w")
        self._register(stages, "stages")
        for i, stage in enumerate(STAGES):
            var = tk.BooleanVar(value=self.saved.get("stage_" + stage, True))
            self.vars["stage_" + stage] = var
            box_ = ttk.Checkbutton(run, text=stage, variable=var,
                                   command=self._sync)
            box_.grid(row=0, column=i + 1, padx=3)
            self._tip(box_, "help_stage_" + stage)
        variant = ttk.Label(run, text=self.t("variant"))
        variant.grid(row=1, column=0, sticky="w", pady=4)
        self._register(variant, "variant")
        self.vars["variant"] = tk.StringVar(value=self.saved.get("variant",
                                                                 "(none)"))
        self.variant_box = ttk.Combobox(run, textvariable=self.vars["variant"],
                                        values=self._variants(), width=18,
                                        state="readonly")
        self.variant_box.grid(row=1, column=1, columnspan=3, sticky="w")
        self.vars["variant"].trace_add("write", lambda *_: self._sync())
        self._tip(self.variant_box, "help_variant")
        self._tip(variant, "help_variant")

    def _build_command(self, parent):
        ttk, tk = self.ttk, self.tk
        box = ttk.Labelframe(parent, text=self.t("command"))
        box.pack(fill="x", padx=14, pady=8)
        self._register(box, "command")
        self.command = tk.Text(box, height=2, wrap="word",
                               font=self.fonts["mono"], background=SURFACE,
                               foreground=INK, relief="flat", padx=8, pady=6,
                               highlightthickness=1,
                               highlightbackground=LINE, highlightcolor=LINE)
        self.command.pack(fill="x", padx=6, pady=6)
        self._tip(self.command, "help_command")
        bar = ttk.Frame(parent)
        bar.pack(fill="x", padx=14, pady=(2, 6))
        for key, command, attr, tip in (
                ("run", self.start, "run_button", "help_run_button"),
                ("stop", self.stop, "stop_button", "help_stop_button"),
                ("dry", lambda: self.start(dry=True), None, "help_dry_button"),
                ("lblwin", self.open_lbl, None, "help_lblwin_button"),
                ("export", self.export, None, "help_export_button"),
                ("savedefaults", self.save_defaults, None,
                 "help_savedefaults_button"),
                ("savelog", self.save_log, None, "help_savelog_button"),
                ("openout", self.open_outputs, None, "help_openout_button")):
            button = ttk.Button(bar, text=self.t(key), command=command,
                                style="Run.TButton" if key == "run"
                                else "TButton")
            button.pack(side="left", padx=(0 if key == "run" else 5, 0))
            self._register(button, key)
            self._tip(button, tip)
            if attr:
                setattr(self, attr, button)
        self.stop_button.configure(state="disabled")
        self.status.pack(side="right")

    def _build_log(self, parent):
        ttk, tk = self.ttk, self.tk
        box = ttk.Labelframe(parent, text=self.t("output"))
        box.pack(fill="both", expand=True, padx=14, pady=(4, 12))
        self._register(box, "output")
        self.log = tk.Text(box, wrap="word", height=3, font=self.fonts["mono"],
                           background=LOG_BG, foreground=LOG_INK,
                           insertbackground=LOG_INK, relief="flat",
                           padx=8, pady=6, highlightthickness=0)
        bar = ttk.Scrollbar(box, command=self.log.yview)
        self.log.configure(yscrollcommand=bar.set)
        bar.pack(side="right", fill="y")
        self.log.pack(fill="both", expand=True, padx=6, pady=6)
        # A long line wraps instead of being cut off, and its continuation is
        # indented to where the message starts, under `| `, so the stamps stay a
        # column of their own and a wrapped sentence still reads as one event.
        # lmargin2 is exactly that: the indent of every line of a paragraph but
        # the first.
        import tkinter.font as tkfont
        indent = tkfont.Font(font=self.fonts["mono"]).measure(
            "%s | " % stamp())
        for _code, (name, colour) in LEVELS.items():
            self.log.tag_configure(name, foreground=colour, lmargin2=indent)
        self.log.tag_configure("plain", foreground=LOG_INK, lmargin2=indent)
        self._tip(self.log, "help_log")

    # ---- state --------------------------------------------------------
    def _config_default(self, path, kind):
        section, name = path.split(".")
        value = (self._config().get(section) or {}).get(name)
        if value is None and kind == "bool":
            return False
        return value

    def _config(self):
        if self._cfg is None:
            try:
                from .config import load_config
                self._cfg = load_config(self.vars["config"].get()
                                        if "config" in self.vars else "config.yaml")
            except Exception:                                 # noqa: BLE001
                self._cfg = {}
        return self._cfg

    def _variants(self):
        folder = os.path.join(os.path.dirname(os.path.abspath(
            self.vars["config"].get() if "config" in self.vars
            else "config.yaml")), "variants")
        names = sorted(os.path.splitext(os.path.basename(p))[0]
                       for p in glob.glob(os.path.join(folder, "*.yaml")))
        return ["(none)"] + names

    def state(self):
        out = {"objects": self.picked(), "lang": self.lang}
        for key, var in self.vars.items():
            out[key] = var.get()
        return out

    def _sync(self, *_args):
        state = self.state()
        self.command.delete("1.0", "end")
        self.command.insert("1.0", " ".join(build_command(state)))
        self._warn_exists(state)
        keep = {k: v for k, v in state.items() if k != "objects"}
        keep["checked"] = sorted(self.checked)
        keep["sort_column"] = getattr(self, "sort_column", None)
        keep["sort_reverse"] = bool(getattr(self, "sort_reverse", False))
        _write_state(keep)

    def _warn_exists(self, state):
        """Say, beside the name, when this run has already been made.

        Not a refusal: re-running is how a fit is redone. But finding somebody
        else's fit under your own name, silently, is how two experiments become
        one set of numbers.
        """
        label = getattr(self, "exists", None)
        if label is None:
            return
        names = state.get("objects") or []
        if not names:
            label.configure(text="")
            return
        root = (state.get("out_dir")
                or (self._config().get("output") or {}).get("directory")
                or "outputs")
        named = str(state.get("run_name") or "").strip()
        tag = "%s-%s" % (state.get("n_star") or 0, state.get("n_earth") or 3)
        where = os.path.join(absolute(root), "_" + named if named else "")
        if len(names) > 1:
            where = os.path.join(where, "joint", "+".join(names), tag)
        else:
            where = os.path.join(where, names[0], tag)
        if os.path.exists(os.path.join(where, "fit.npz")):
            label.configure(text=self.t("exists"), foreground="#b26a00")
        else:
            label.configure(text="", foreground=MUTED)

    def refresh_objects(self):
        """Show what is remembered of this data root, then go and check it.

        The index is read first and drawn at once, so a folder on a disk that
        has to spin up is not waited for; the scan that follows reads only the
        files that were added or replaced, and redraws when it is done.
        """
        root = self.vars["data_dir"].get()
        self._cfg = None
        self.variant_box.configure(values=self._variants())
        if self.scanning:
            # one scan at a time: two threads walking the same index would
            # overwrite each other's answers
            self._say("log_busy", level="warn")
            return
        self.index = scan.load(root)
        self._fill(scan.summaries(self.index))
        if not os.path.isdir(root):
            self._say("log_no_root", root, level="warn")
            return
        self.scanning = True
        self._say("log_scan_start", root)
        self._say("log_index", scan.index_path(root), level="value")
        threading.Thread(target=self._scan, args=(root,), daemon=True).start()

    def _fill(self, rows):
        """Draw one row per object, keeping whatever was ticked.

        Grouped by instrument and alphabetical within a group, because a run is
        one instrument: the targets that can go together are then adjacent, and
        each row is tinted with its instrument's colour so the grouping survives
        a glance. Instruments whose box is unticked are not drawn at all.
        """
        self.rows = {row["object"]: row for row in rows}
        self._instrument_boxes(sorted({r.get("instrument") for r in rows
                                       if r.get("instrument")
                                       and r.get("instrument") != "?"}))
        shown = self._sorted([r for r in rows
                              if self._instrument_on(r.get("instrument"))])
        self.tree.delete(*self.tree.get_children())
        self.names = {}
        for row in shown:
            name = row["object"]
            tint = self._tint(row.get("instrument") or "?")
            item = self.tree.insert("", "end", values=self._values(row),
                                    tags=(tint,))
            self.names[item] = name
            self._draw_check(item, name)
        # an object that is gone, or whose instrument is hidden, is not run
        self.checked &= set(self.names.values())
        self.count.configure(text="%d / %d" % (len(self.picked()),
                                               len(self.names)))
        self._sync()

    def _draw_berv(self):
        """Redraw the coverage histogram for whatever is ticked."""
        canvas = getattr(self, "berv_canvas", None)
        if canvas is None:
            return
        canvas.delete("all")
        names = self.picked()
        width = max(int(canvas.winfo_width()), 50)
        height = max(int(canvas.winfo_height()), 40)
        if not names:
            canvas.create_text(width // 2, height // 2, text=self.t("berv_none"),
                               fill="#888", font=("Helvetica", 10))
            self.berv_note.configure(text="")
            return
        edges, counts, summary = scan.berv_coverage(self.index, names)
        if not counts:
            canvas.create_text(width // 2, height // 2, text=self.t("berv_wait"),
                               fill="#b26a00", font=("Helvetica", 10))
            self.berv_note.configure(text="")
            return
        pad, foot = 6, 16
        top = height - foot
        tallest = max(1, max(int(c.max()) for c in counts.values() if c.size))
        n_bins = len(edges) - 1
        step = (width - 2 * pad) / max(n_bins, 1)

        def x_of(v):
            return pad + (v - edges[0]) / max(edges[-1] - edges[0], 1e-9) * (
                width - 2 * pad)

        # what the SKY allows this selection, |BERV| <= 29.78 cos(beta): the
        # axis is the whole solar system, and a target at a high ecliptic
        # latitude can never fill it however long it is observed. Drawn under
        # the bars so an empty stretch inside the limits reads as a gap to be
        # filled, and one outside them as nothing anybody can do.
        reach = (summary or {}).get("possible")
        if reach:
            canvas.create_rectangle(x_of(-reach / 2.0), pad, x_of(reach / 2.0),
                                    top, fill="#eef3f7", outline="")
            for edge in (-reach / 2.0, reach / 2.0):
                canvas.create_line(x_of(edge), pad, x_of(edge), top,
                                   fill="#9fb8c9", dash=(2, 2))
        for i in range(n_bins):
            x0 = pad + i * step
            bottom = top
            for name in names:
                c = counts.get(name)
                if c is None or not c[i]:
                    continue
                # stacked, each star in the colour its row has in the list
                h = (top - pad) * c[i] / tallest
                canvas.create_rectangle(x0, bottom - h, x0 + max(step - 1, 1),
                                        bottom,
                                        fill=self._star_colour(name),
                                        outline="")
                bottom -= h
        canvas.create_line(pad, top, width - pad, top, fill="#bbb")
        # a legend, or the colours say nothing: one square and one name per star,
        # in the order they are stacked
        x = pad + 2
        for name in names:
            if name not in counts:
                continue
            canvas.create_rectangle(x, pad, x + 8, pad + 8,
                                    fill=self._star_colour(name), outline="")
            label = canvas.create_text(x + 11, pad + 4, text=name, anchor="w",
                                       fill="#444", font=("Helvetica", 8))
            x = canvas.bbox(label)[2] + 10
            if x > width - 60:
                break
        # every 10 km/s: on a fixed axis the labels are the ruler the bars are
        # read against, and two numbers at the ends are not a ruler
        for value in (-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0):
            if not edges[0] <= value <= edges[-1]:
                continue
            x = x_of(value)
            canvas.create_line(x, top, x, top + 3, fill="#888")
            canvas.create_text(x, top + 9, text="%+.0f" % value, fill="#555",
                               font=("Helvetica", 8))
        # said while the scan is still reading, since the histogram is then
        # drawn from part of the campaign and would otherwise look final
        if self._berv_partial(names):
            # below the legend, never over it, and on its own ground: the bars
            # it is warning about are what made it hard to read
            chip(canvas, width // 2, pad + 24, self.t("berv_building"),
                 ink="#b26a00", fill="#fdf1dd", outline="#e0a94a",
                 font=("Helvetica", 10, "bold"))
        possible = summary.get("possible")
        self.berv_note.configure(
            text=self.t("berv_note")
            % (summary["effective"], summary["span"],
               "%.0f" % possible if possible else "?",
               scan.BERV_BIN, summary["n"]))

    def _berv_partial(self, names):
        """Whether some spectrum of a ticked object has not been read yet."""
        if self.scanning:
            return True
        for name in names:
            row = self.rows.get(name) or {}
            known = len(((self.index.get("objects") or {}).get(name) or {})
                        .get("files") or {})
            if known < int(row.get("files") or 0):
                return True
        return False

    def _star_colour(self, name):
        """One colour per STAR, stable for as long as the window is open.

        Not per instrument: a pool of four NIRPS campaigns stacked in one blue
        shows nothing about which star fills which bin, and which star fills
        which bin is the whole reason for stacking them.
        """
        self._star_colours = getattr(self, "_star_colours", {})
        if name not in self._star_colours:
            self._star_colours[name] = STAR_COLOURS[len(self._star_colours)
                                                    % len(STAR_COLOURS)]
        return self._star_colours[name]

    def _logo(self, name, height=26):
        """A PNG from the package's assets, scaled to `height`, or None.

        None rather than a raised exception: a missing or unreadable image is a
        decoration that did not appear, never a window that did not open.
        """
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "assets", name)
        if not os.path.exists(path):
            return None
        try:
            from PIL import Image, ImageTk
            image = Image.open(path)
            scale = height / float(image.height)
            return ImageTk.PhotoImage(
                image.resize((max(1, int(image.width * scale)), height),
                             Image.LANCZOS))
        except Exception:                                       # noqa: BLE001
            try:
                photo = self.tk.PhotoImage(file=path)
                step = max(1, int(round(photo.height() / float(height))))
                return photo.subsample(step, step)
            except Exception:                                   # noqa: BLE001
                return None

    def _tint(self, instrument):
        """The row colour of an instrument, made once and kept.

        An instrument not read yet gets NO colour: a tint says "this
        instrument", and a row whose instrument is still unknown must not be
        painted as though it were a third one. It was, in the pale yellow of the
        fallback palette, which in this project is the colour of a missing
        sample (2026-09-13).
        """
        self._tints = getattr(self, "_tints", {})
        if not instrument or instrument in ("?", UNREAD):
            return ""
        tag = "inst_%s" % instrument
        if tag not in self._tints:
            colour = INSTRUMENT_TINT.get(instrument.upper())
            if colour is None:
                colour = OTHER_TINTS[len(self._tints) % len(OTHER_TINTS)]
            self._tints[tag] = colour
            self.tree.tag_configure(tag, background=colour)
        return tag

    def _instrument_on(self, instrument):
        var = self.instrument_vars.get(instrument or "?")
        return True if var is None else bool(var.get())

    def _instrument_boxes(self, instruments):
        """One box per instrument the data root holds, all ticked to begin with.

        Rebuilt only when the set of instruments changes, so that boxes somebody
        just unticked are not silently ticked again by a rescan.
        """
        if list(self.instrument_vars) == list(instruments):
            return
        for child in self.filters.winfo_children():
            child.destroy()
        keep = dict(self.instrument_vars)
        self.instrument_vars = {}
        for name in instruments:
            var = self.tk.BooleanVar(
                value=bool(keep[name].get()) if name in keep else True)
            box = self.ttk.Checkbutton(
                self.filters, text=name, variable=var,
                command=lambda: self._fill(list(self.rows.values())))
            box.pack(side="left", padx=(0, 6))
            self._tip(box, "help_instrument_filter")
            self.instrument_vars[name] = var

    def _values(self, row, approximate=False):
        """One row's columns.

        A number not read yet is blank, never a zero. The tilde marks what is
        still an ESTIMATE from the first few spectra and will sharpen as the rest
        are read: the median signal-to-noise and the median exposure time, which
        are medians over exposures. The magnitude is not one of them. It is a
        property of the star, the same in the first file as in the last, so a
        tilde on it would promise a precision that has nothing to gain.
        """
        tilde = "~" if approximate else ""
        mag = ("" if row.get("mag") is None else
               "%s=%.1f" % (row.get("mag_band") or "?", row["mag"]))
        return (row["files"],
                "" if row.get("snr") is None else "%s%.0f" % (tilde, row["snr"]),
                "" if row.get("exptime") is None
                else "%s%.0f" % (tilde, row["exptime"]),
                mag,
                # "?" claims the file was read and said nothing; UNREAD says
                # the scan has not got there yet, which is the usual case
                (row.get("instrument") if row.get("instrument") not in (None, "?")
                 else (UNREAD if approximate or not row.get("snr") else "?")))

    def _listed(self, counts):
        """The names and the counts, before a single header has been read."""
        rows = []
        for name in sorted(counts):
            row = scan.summary(self.index, name)
            row["files"] = counts[name]    # the folder is the truth for the count
            rows.append(row)
        self._fill(rows)

    def _one(self, name, known, total):
        """This object's numbers so far, in place, marked if they are partial.

        `known` spectra of `total` have been read. Fewer than all of them is an
        estimate and says so with a tilde; the count shown stays the folder's.
        """
        row = scan.summary(self.index, name)
        row["files"] = total
        self.rows[name] = row
        for item, shown in self.names.items():
            if shown == name:
                self.tree.item(item, values=self._values(row, known < total))
                break
        # The coverage panel and the timeline are drawn from what has been
        # READ, so a campaign that has just been finished changes both of them,
        # and the banner saying they are partial has to go when it stops being
        # true. Only when the campaign is done and it is ticked: every ten
        # spectra would be a redraw of both panels a few hundred times.
        if known >= total and name in self.picked():
            self._draw_berv()
            self._draw_time()

    def _scan(self, root):
        """Read what changed, off the main thread, touching no widget.

        Everything this thread has to say goes into the queue the main thread
        drains: Tk may only be called from the thread running its loop, and a
        widget poked from here raises or, worse, does not.
        """
        def on_file(name, done, total):
            self.lines.put(("status", "%s  %d/%d" % (name, done, total)))

        def on_listed(counts):
            # the names and the counts cost one folder listing, so they are
            # shown before any header is read: a first scan of a campaign on a
            # shared disk is minutes, and an empty list says nothing meanwhile
            self.lines.put(("listed", counts))

        saved = [0.0]

        def on_object(name, known, total):
            # saved when an object is done, and otherwise at most every 20 s: a
            # scan interrupted halfway keeps what it read, without writing the
            # whole index 220 times on the way
            now = time.time()
            if known >= total or now - saved[0] > 20.0:
                scan.save(self.index, root)
                saved[0] = now
            self.lines.put(("object", name, known, total))
        try:
            index, tally = scan.update(root, index=self.index, on_file=on_file,
                                       on_listed=on_listed, on_object=on_object)
            path = scan.save(index, root)
        except OSError as exc:
            self.lines.put(("line", self._line("log_failed", exc), "error"))
            self.scanning = False
            return
        self.lines.put(("scanned", root, index, tally, path))

    def _scanned(self, root, index, tally, path):
        self.scanning = False
        self.index = index
        rows = scan.summaries(index)
        self._fill(rows)
        # the scan is over, so both panels are final. Nothing else redraws them
        # until a tick changes or the window is resized, and "UNDER
        # CONSTRUCTION: still reading" sat on a finished histogram until then
        self._draw_berv()
        self._draw_time()
        self._state()
        if not rows:
            self._say("log_scan_none", root, level="warn")
            return
        self._say("log_scan_done", os.path.basename(path or root), len(rows),
                  tally["read"], tally["kept"], tally["gone"], level="value")

    # ---- running ------------------------------------------------------
    def start(self, dry=False):
        from tkinter import messagebox
        if self.proc is not None:
            return
        state = self.state()
        if not state["objects"]:
            self._say("log_no_object", level="error")
            return
        state["dry_run"] = dry
        # an error, not a warning: one run is one instrument, and a joint fit of
        # two spectrographs has no grid to be fitted on
        instruments = self.instruments(state["objects"])
        if len(instruments) > 1:
            self._say("log_mixed_refused", ", ".join(state["objects"]),
                      ", ".join(sorted(instruments)), level="error")
            messagebox.showerror(self.t("mixed_title"),
                                 self.t("mixed") % (", ".join(state["objects"]),
                                                    ", ".join(sorted(instruments))))
            return
        argv = build_command(state)
        # shown as `pca2d-preclean ...`, since that is what to paste into a
        # terminal, but RUN through whatever path actually holds it here
        found = preclean_argv()
        if found is None and not self._offer_install():
            return
        found = found or preclean_argv()
        if found is None:
            self._say("log_failed", self.t("no_command_still"), level="error")
            return
        self._say("log_command", " ".join(argv), level="value")
        argv = found + argv[1:]
        env = dict(os.environ, PCA2D_COLOUR="1", PYTHONUNBUFFERED="1")
        try:
            self.proc = subprocess.Popen(
                argv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1, env=env,
                cwd=os.path.dirname(os.path.abspath(self.vars["config"].get()))
                or None)
        except OSError as exc:
            self._say("log_failed", exc, level="error")
            self.proc = None
            return
        self.run_button.configure(state="disabled")
        self.stop_button.configure(state="normal")
        self.status.configure(text=self.t("running"))
        threading.Thread(target=self._reader, daemon=True).start()

    def _offer_install(self):
        """Offer to install the command when it cannot be found at all.

        `pip install -e .` in the folder holding the configuration, which is the
        repository: it puts the entry point in this interpreter's own bin, where
        preclean_argv looks first. Asked, never done behind anyone's back.
        """
        from tkinter import messagebox
        here = os.path.dirname(os.path.abspath(self.vars["config"].get()))
        if not messagebox.askyesno(self.t("no_command_title"),
                                   self.t("no_command") % (sys.executable, here)):
            self._say("log_no_command", level="error")
            return False
        self._say("log_installing", here, level="value")
        try:
            done = subprocess.run([sys.executable, "-m", "pip", "install", "-e",
                                   ".", "--no-deps"], cwd=here,
                                  capture_output=True, text=True, timeout=600)
        except Exception as exc:                                # noqa: BLE001
            self._say("log_failed", exc, level="error")
            return False
        for line in (done.stdout or "").splitlines()[-3:]:
            self._write("   %s\n" % line, "plain")
        if done.returncode:
            self._say("log_failed", (done.stderr or "").strip()[-200:],
                      level="error")
            return False
        return True

    def _reader(self):
        for line in self.proc.stdout:
            self.lines.put(line)
        code = self.proc.wait()
        self.proc = None
        # through the queue, not straight to the window: the last lines the run
        # printed are still in it, and "it ended" belongs after them
        self.lines.put(("line", self._line("log_ended", code),
                        "info" if not code else "error"))
        self.lines.put(("finished",))

    def _finished(self):
        self.run_button.configure(state="normal")
        self.stop_button.configure(state="disabled")
        self._state()

    def stop(self):
        if self.proc is not None:
            self._say("log_stopping", level="warn")
            self.proc.terminate()

    def _drain(self):
        """Everything the threads produced, applied here on the main thread.

        A thread that has something to show puts it in this queue and touches
        nothing: one place where widgets change, and the run's output and the
        window's own lines stay in the order they happened.
        """
        while True:
            try:
                item = self.lines.get_nowait()
            except queue.Empty:
                break
            if isinstance(item, str):
                self._write(item)
            elif item[0] == "line":
                self._write(item[1], item[2])
            elif item[0] == "status":
                self.status.configure(text=item[1])
            elif item[0] == "listed":
                self._listed(item[1])
            elif item[0] == "object":
                self._one(*item[1:])
            elif item[0] == "scanned":
                self._scanned(*item[1:])
            elif item[0] == "finished":
                self._finished()
        self._drain_id = self.root.after(80, self._drain)

    # ---- what the window itself says ----------------------------------
    def _line(self, key, *args):
        """One of the window's own lines, in the project's convention:
        `YYMMDD HH:MM:SS.SS | message`, and in the window's language."""
        body = self.t(key)
        try:
            body = body % args if args else body
        except (TypeError, ValueError):
            # a translation whose placeholders drifted must not stop the window
            body = "%s %s" % (body, " ".join(str(a) for a in args))
        return "%s | %s\n" % (stamp(), body)

    def _say(self, key, *args, **kwargs):
        self._write(self._line(key, *args), kwargs.get("level", "info"))

    def _write(self, line, forced=None):
        tag = forced or "plain"
        if forced is None:
            for code in ANSI.findall(line):
                if code in LEVELS:
                    tag = LEVELS[code][0]
                    break
        body = ANSI.sub("", line)
        if body.startswith("\r"):
            self.log.delete("end-2l", "end-1l")
            body = body.lstrip("\r")
        self.log.insert("end", body, tag)
        self.log.see("end")
        if "| stage " in body:
            self.status.configure(text=body.split("| ", 1)[-1].strip()[:60])

    # ---- the buttons that write things --------------------------------
    def export(self):
        from tkinter import filedialog, messagebox

        import yaml
        body = variant_yaml(self.state(), self._config())
        if not body:
            messagebox.showinfo("pca2d", self.t("nothing_export"))
            return
        path = filedialog.asksaveasfilename(
            title=self.t("export_title"),
            initialdir=os.path.join(os.path.dirname(os.path.abspath(
                self.vars["config"].get())), "variants"),
            defaultextension=".yaml", initialfile="mine.yaml")
        if not path:
            return
        head = ("# written by pca2d-gui: the nominal configuration plus these\n"
                "# lines. Run it with: pca2d-preclean --object NAME --variant %s\n"
                % os.path.splitext(os.path.basename(path))[0])
        with open(path, "w") as handle:
            handle.write(head)
            yaml.safe_dump(body, handle, sort_keys=False, default_flow_style=False)
        self._say("log_export", path, level="value")
        self.variant_box.configure(values=self._variants())

    def save_defaults(self):
        """Write what was changed into config.yaml itself, comments and all.

        A variant leaves the nominal alone and is the reproducible way to run
        something once; this is for when a choice has been settled and should
        be what every later run starts from. The file's comments are the
        measurements that chose each value, so they are kept (config.update_file).
        """
        from tkinter import messagebox

        from .config import update_file
        body = variant_yaml(self.state(), self._config())
        values = {"%s.%s" % (section, key): value
                  for section, block in body.items()
                  for key, value in block.items()}
        if not values:
            self._say("log_nothing", level="warn")
            messagebox.showinfo("pca2d", self.t("nothing_export"))
            return
        path = os.path.abspath(self.vars["config"].get())
        shown = "\n".join("  %s: %s" % (name, values[name])
                          for name in sorted(values))
        if not messagebox.askyesno(self.t("defaults_title"),
                                   self.t("defaults_ask") % (path, shown)):
            return
        try:
            written = update_file(path, values)
        except (SystemExit, OSError) as exc:
            self._say("log_failed", exc, level="error")
            messagebox.showerror("pca2d", str(exc))
            return
        self._cfg = None                     # the configuration has moved
        self._say("log_defaults", path, ", ".join(sorted(written)), level="value")

    def save_log(self):
        from tkinter import filedialog
        path = filedialog.asksaveasfilename(title=self.t("savelog"),
                                            defaultextension=".log")
        if path:
            with open(path, "w") as handle:
                handle.write(self.log.get("1.0", "end"))
            self._say("log_saved_log", path, level="value")

    def open_outputs(self):
        """Show the output root in the file browser: a run leaves one PDF and a
        folder of corrected spectra, and they are easier to find by looking."""
        root = self.vars["out_dir"].get() or (
            (self._config().get("output") or {}).get("directory") or "outputs")
        path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(
            self.vars["config"].get())), root))
        if not os.path.isdir(path):
            self._say("log_no_outputs", path, level="warn")
            return
        opener = ("open" if sys.platform == "darwin"
                  else "explorer" if os.name == "nt" else "xdg-open")
        subprocess.Popen([opener, path])

    def open_lbl(self):
        """The `lbl:` block in its own window: the wrapper's own questions.

        Built once and hidden rather than destroyed, so that the language switch
        keeps finding its labels and the window comes back as it was left.
        """
        ttk, tk = self.ttk, self.tk
        if self.lbl_window is not None:
            self.lbl_window.deiconify()
            self.lbl_window.lift()
            return
        window = self.lbl_window = tk.Toplevel(self.root)
        window.title(self.t("lbl_title"))
        window.geometry("620x380")
        head = ttk.Label(window, text=self.t("lbl_title"), style="Head.TLabel")
        head.pack(anchor="w", padx=12, pady=(12, 0))
        self._register(head, "lbl_title")
        note = ttk.Label(window, style="Hint.TLabel", wraplength=580,
                         justify="left", text=self.t("help_lblwin_button"))
        note.pack(anchor="w", padx=12, pady=(2, 8))
        self._register(note, "help_lblwin_button")
        grid = ttk.Frame(window)
        grid.pack(fill="both", expand=True, padx=12)
        for i, (key, path, kind) in enumerate(OPTIONS_LBL):
            row, col = i % 6, (i // 6) * 2
            label = ttk.Label(grid, text=self.t("opt_" + key))
            label.grid(row=row, column=col, sticky="w", pady=3)
            self._register(label, "opt_" + key)
            default = self.saved.get(key, self._config_default(path, kind))
            if kind == "list":
                default = ", ".join(default or []) if isinstance(default, list) \
                    else default
            if kind == "bool":
                var = tk.BooleanVar(value=bool(default))
                widget = ttk.Checkbutton(grid, variable=var)
            elif isinstance(kind, tuple):
                var = tk.StringVar(value=str(default))
                widget = ttk.Combobox(grid, textvariable=var, values=list(kind),
                                      width=12, state="readonly")
            else:
                var = tk.StringVar(value="" if default is None else str(default))
                widget = ttk.Entry(grid, textvariable=var, width=24)
            widget.grid(row=row, column=col + 1, sticky="w", padx=(8, 20))
            var.trace_add("write", lambda *_: self._sync())
            self.vars[key] = var
            self._tip(widget, "help_" + key)
            self._tip(label, "help_" + key)
        close = ttk.Button(window, text=self.t("close"), command=window.withdraw)
        close.pack(anchor="e", padx=12, pady=10)
        self._register(close, "close")
        window.protocol("WM_DELETE_WINDOW", window.withdraw)

    def _propose_out(self):
        """Offer a place for the run's products, beside the data it reads.

        Only when the field is empty: a proposal, not a decision. Said in the
        log, because a folder that is about to receive tens of gigabytes should
        not appear in a box without a word.
        """
        if self.vars["out_dir"].get().strip():
            return
        proposed = corrected_dir(self.vars["data_dir"].get())
        if proposed:
            self.vars["out_dir"].set(proposed)
            self._say("log_out_proposed", proposed, level="value")

    def _browse(self, key):
        from tkinter import filedialog
        if key == "config":
            path = filedialog.askopenfilename(title=self.t("config"),
                                              filetypes=[("YAML", "*.yaml")])
        else:
            path = filedialog.askdirectory(title=self.t(key))
        if path:
            self.vars[key].set(absolute(path))
            if key == "data_dir":
                self._propose_out()
                self.refresh_objects()

    def _close(self):
        if self.proc is not None:
            self.proc.terminate()
        # cancel the pending drain, or it fires into a window that is gone
        after = getattr(self, "_drain_id", None)
        if after is not None:
            try:
                self.root.after_cancel(after)
            except Exception:                                   # noqa: BLE001
                pass
        self.root.destroy()


def main(argv=None):
    import tkinter as tk

    root = tk.Tk()
    App(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    sys.exit(main())
