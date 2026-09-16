"""A window to run the pipeline, for people who should not have to read the
command line first.

    pca2d-gui

It shows the objects that are in the data root, lets one of them be run on its
own or several of them together against one observer basis, exposes the
handful of settings that change a result, writes the command it is about to
run so that it can be copied into a terminal, and streams the run's own log
into the window with the colours it would have in a terminal.

Every item explains itself on hover, in English, French, Spanish or Portuguese,
and says what
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

import numpy as np

from . import scan
from .logger import stamp

HOME_STATE = os.path.expanduser("~/.pca2d_gui.json")
#: the family that carries the colour emoji here. Tk falls back to the default
#: font when a family is unknown, so a machine without it shows the character
#: in black and white rather than failing to draw anything.
EMOJI_FONT = ("Apple Color Emoji" if sys.platform == "darwin"
              else "Segoe UI Emoji" if os.name == "nt" else "Noto Color Emoji")
#: EVERY escape sequence, not only the ones this window paints with. LBL
#: colours its own output (`\033[92;1m` for a bright bold green, `\033[0;0m` to
#: reset), and a pattern that matched a single number left those on screen as
#: "[92;1m" and "[0;0m" around every line it printed. Semicolons, empty codes
#: and the other SGR forms are all matched here and removed; which of the codes
#: inside decides the colour is LEVELS' business, below.
ANSI = re.compile(r"\033\[[0-9;]*m")
#: the numbers inside one of those, so a sequence like 92;1 is read as 92 and 1
ANSI_CODES = re.compile(r"[0-9]+")
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
    # Neither mean nor star_basis is here: the static part is one star spectrum
    # per order parity and the star is one cubic B-spline, both decided in the
    # code (config.DEFAULTS). The other modes are still read from a config or a
    # variant file, for redoing the runs that were made on them.
    ("velocity_term", "twoframe.velocity_term", "bool"),
    ("iters", "twoframe.iters", "int"),
    ("shrink", "correct.shrink", "bool"),
    # two buttons, one of which is always down: the metric is a choice between
    # two spellings of the same thing, and a tick box can only name one of them.
    # The kind carries (value, what the button says), and the DEFAULT is the
    # config's, which is the derivative (pca2d.gui._build_options)
    ("weight", "correct.weight",
     ("radio", ("flux", "F"), ("velocity", "(dF/dv)\u00b2"))),
    # no correct.mask here either: one line set for the whole campaign
    # ("common"), decided in the code
    ("width_kms", "highpass.width_kms", "float"),
    ("dv", "domain.dv", "float"),
    ("nightly_stack", "input.nightly_stack", ("auto", "true", "false")),
    # lbl.run is in the LBL window, not here: the stages already say whether
    # the lbl step happens at all, and two boxes for one step is one too many
]
#: the `lbl:` block, in its own window: what the wrapper would have been asked
OPTIONS_LBL = [
    ("run", "lbl.run", "bool"),
    ("lbl_prepare", "lbl.prepare", "bool"),
    ("lbl_before", "lbl.before", "bool"),
    ("lbl_after", "lbl.after", "bool"),
    ("lbl_star_template", "lbl.star_template", "bool"),
    ("lbl_strpca", "lbl.strpca", "bool"),
    ("lbl_suffix", "lbl.suffix", "text"),
    ("lbl_teff", "lbl.teff", "text"),
    ("lbl_template", "lbl.template", "text"),
    ("lbl_steps", "lbl.steps", "list"),
]
#: shown with the paths, beside the LBL folder it is about: a setting of what
#: goes INTO that folder, and on the LBL page nobody found it (2026-09-16)
OPTIONS_PATHS = [
    ("lbl_link", "lbl.link", ("symlink", "copy")),
]
#: every option the window can change, wherever it is shown
ALL_OPTIONS = OPTIONS + OPTIONS_LBL + OPTIONS_PATHS
#: the config keys the window sets through a PATH field rather than an option:
#: a folder on this machine, never exported to a variant, but a setting of the
#: run all the same, and on the report's front page with the others
PATH_SETTINGS = [("lbl_dir", "lbl.directory")]


def window_settings():
    """Every config key the window sets, in the order the report lists them."""
    return ([path for _key, path, _kind in OPTIONS + OPTIONS_LBL]
            + [path for _key, path in PATH_SETTINGS]
            + [path for _key, path, _kind in OPTIONS_PATHS])


#: the flag each LBL setting travels under (cli.SETTING_FLAGS). Until
#: 2026-09-16 none of them travelled at all: the LBL page was shown, changed,
#: and ignored by the run, which took the configuration's own values
LBL_FLAGS = {
    "run": "--lbl-run", "lbl_prepare": "--lbl-prepare",
    "lbl_before": "--lbl-before", "lbl_after": "--lbl-after",
    "lbl_star_template": "--lbl-star-template", "lbl_strpca": "--lbl-strpca",
    "lbl_suffix": "--lbl-suffix", "lbl_teff": "--lbl-teff",
    "lbl_template": "--lbl-template", "lbl_steps": "--lbl-steps",
    "lbl_link": "--lbl-link",
}
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
    "fits_dir": "products disk (optional)",
    "help_fits_dir":
        "Where a run's products are KEPT: each run folder becomes one link to"
        " this disk, so the corrected spectra, the fit and the figures land"
        " there and not on the internal disk. Empty keeps everything under the"
        " output root, which works. A path is a thing on ONE machine, so this"
        " field is emptied whenever what it names is not there, and a run whose"
        " disk is missing stops rather than quietly fill the internal one.",
    "log_no_disk":
        "the configuration names %s as the disk to keep products on, and it is"
        " not there: the field is empty, so this run would keep everything"
        " under the output root. Fill it in if the disk should be mounted.",
    "out_dir": "output root (optional)",
    "lbl_dir": "LBL output folder",
    "help_lbl_dir":
        "LBL's own tree (its DATA_DIR): the science folders it reads, and the"
        " templates, masks, per-line tables and rdb it writes. Proposed as lbl"
        " under the output root, and it follows that root until you type"
        " another. One tree for every run under a root, since the delivered"
        " object's LBL is the same for all of them and hours to make; point it"
        " at an existing tree and that tree is used as it is.",
    "browse": "Browse", "rescan": "Rescan",
    "objects": "objects", "settings": "settings", "stages": "stages",
    "command": "the command this runs", "output": "output",
    "col_object": "object", "col_files": "files", "col_instrument": "instrument",
    "col_snr": "SNR", "col_exptime": "exp (s)", "col_mag": "mag",
    'run_name': 'reduction name', 'auto': 'auto',
    'help_auto_button':
        'Proposes a name again from the targets and the settings as they'
        ' stand: the targets, then six characters of a hash of everything that'
        ' makes this reduction a different result. The same parameters give the'
        ' same six, one number different gives another six.',
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
    "openpdf": "Open compil PDF",
    "savedefaults": "Save as defaults...",
    "all": "all", "none": "none",
    # the language the window is in, on the button that opens the others
    "idle": "idle", "running": "running", "lang": "\U0001F1EC\U0001F1E7 English",
    "snr_berv": "signal-to-noise against barycentric velocity",
    "snr_none": "tick a target to see where its best nights sit",
    "help_snr_berv":
        "One dot per spectrum, in its star's colour: the signal-to-noise APERO"
        " measured against the barycentric velocity it was taken at. The"
        " histogram above says which velocities a campaign covers; this says"
        " what it covers them WITH. The fit weighs a spectrum by 1/sigma^2, so"
        " a range covered at one end by the worst nights of a campaign is not"
        " the range the fit really sees, and the two frames come apart less"
        " well than the coverage promises.",
    "tab_targets": "  targets  ", "tab_settings": "  settings  ",
    "tab_lbl": "  LBL  ", "tab_run": "  analysis  ",
    "tab_clean": "  cleanup  ",
    "clean_title": "what pca2d has left on these disks",
    "help_clean":
        "Every place this program puts bytes, measured. The cube cache, the"
        " scratch of a fit too big for memory and LBL's own logs can go at any"
        " time: they are made again from what is still here, and one button"
        " empties all of them. Everything else can go too, row by row: pick it"
        " in the list and Delete the selected. Fits, templates, masks and"
        " per-line tables come back only by running the hours that made them,"
        " and the corrected spectra, the reports and the velocities only by"
        " running again from the spectra; the window says which of those you"
        " are about to spend before it deletes anything. A folder moved to"
        " another disk is counted where it really is; a link to somebody"
        " else's spectra weighs nothing, because deleting it frees nothing.",
    "clean_measure": "Measure", "clean_purge": "Delete what can go",
    "clean_selected": "Delete the selected",
    "help_clean_measure":
        "Walks each folder and adds up what is in it. On a network disk this"
        " takes a few seconds and it runs off to one side, so the window stays"
        " usable while it counts.",
    "help_clean_purge":
        "Empties the folders marked as scratch or rebuildable, and nothing"
        " else: the ones that cost time and only time to have back. It asks"
        " first, and names the total it is about to free. The folders"
        " themselves stay: the next run expects to find them. For anything"
        " else, pick the rows and use Delete the selected.",
    "help_clean_selected":
        "Deletes the rows picked in the list, whatever kind they are: a cube"
        " cache, the templates that took an afternoon, the velocities"
        " themselves. Shift or command picks several. It names every folder"
        " and says what having them back would cost before it deletes"
        " anything, and that is the only thing between you and a disk with"
        " room on it.",
    "clean_measuring": "measuring…",
    "clean_totals": "%s in all, of which %s can be freed",
    "clean_none": "nothing measured yet: press Measure",
    "clean_confirm_title": "delete the rebuildable files",
    "clean_confirm":
        "About to free %s from %d places:\n\n%s\n\nNothing here is a result:"
        " the cubes are read again from the spectra, and LBL's logs are"
        " written again by LBL. Go ahead?",
    "clean_freed": "freed %s from %d places",
    "clean_nothing": "nothing to free: there is no scratch or cache here",
    "clean_pick": "pick one row or several in the list first: the button"
                  " deletes what is picked",
    "clean_while_running":
        "A RUN IS GOING, and it reads from these folders. It builds again"
        " whatever it finds missing, so deleting now costs that run the time"
        " to make it a second time, in the middle of the stage it is in.",
    "clean_confirm_pick_title": "delete what is selected",
    "clean_confirm_pick":
        "About to delete %s from %d places:\n\n%s\n\n%s\n\nGo ahead?",
    "clean_cost_rebuildable":
        "All of it is made again from what stays here: minutes, and a run that"
        " reads the spectra once more.",
    "clean_cost_expensive":
        "Some of it comes back only by running what made it: a fit, or an LBL"
        " pass over every exposure. Hours, not minutes.",
    "clean_cost_results":
        "SOME OF IT IS A RESULT: corrected spectra, a report, or the"
        " velocities an rdb holds. Nothing here remakes those. Only running"
        " the whole thing again does, from the spectra, and that is the hours"
        " it took the first time.",
    "col_size": "size", "col_nfiles": "files", "col_kind": "kind",
    "kind_scratch": "scratch", "kind_rebuildable": "rebuildable",
    "kind_expensive": "expensive", "kind_results": "results",
    # what each line of the cleanup list is, said when the row is picked
    "clean_cache":
        "The cubes: every spectrum of a campaign on one wavelength grid, so a"
        " second run does not read them all again. Usually the biggest thing"
        " here, and never a result. Rebuilt in minutes per campaign.",
    "clean_spill":
        "The mapped scratch of a fit too big to hold in memory. Nothing reads"
        " it once the fit has ended, so anything still here belongs to a run"
        " that was interrupted.",
    "clean_results":
        "The corrected spectra, the reports and the fits. What the whole thing"
        " was for; never offered for deletion.",
    "clean_lbl_science":
        "The links LBL measures through, one per exposure. Remade by the lbl"
        " stage. They are symlinks, so this frees almost nothing and the"
        " spectra they point at are not touched.",
    "clean_lbl_plots": "LBL's own figures. Drawn again whenever LBL runs.",
    "clean_lbl_log": "LBL's logs. Nothing reads them but a person.",
    "clean_lbl_lblrv":
        "LBL's per-line velocity tables, one file per exposure, and the"
        " largest thing LBL writes. Remade only by running LBL again, which is"
        " hours.",
    "clean_lbl_templates":
        "The templates LBL built, one per object and run. Remade only by"
        " running LBL's template step again.",
    "clean_lbl_masks": "The line masks LBL built, one per object and run.",
    "clean_lbl_models": "LBL's models.",
    "clean_lbl_calib": "LBL's calibrations.",
    "clean_lbl_lblreftable": "LBL's reference tables, one per object and run.",
    "clean_lbl_lblrdb":
        "The velocities: the rdb every RV page of every report is drawn from."
        " A result, and the smallest thing on this list; deleting it frees"
        " almost nothing and throws away the measurement.",
    "clean_pycache":
        "Compiled Python, remade the next time the package is imported.",
    "quit": "Quit", "quit_title": "quit pca2d-preclean",
    "quit_yes": "Quit anyway", "quit_no": "Stay",
    "quit_running":
        "A run is going, and it is a subprocess of this window: quitting stops"
        " it. What the stages before it wrote stays where it is, the stage it is"
        " in is lost. Quit anyway?",
    "help_quit_button":
        "Closes the window. The settings are written at every change, so nothing"
        " here is lost by leaving; a run that is going is stopped, and it asks"
        " before doing that.",
    # said instead of "idle" for as long as there is nowhere to read from
    "pick_root": "pick a data root: Browse, beside the field at the top",
    "command_pending":
        "tick a target and the command appears here, in full, before it runs",
    "log_pick_root":
        "no data root yet. Browse to the folder that holds ONE FOLDER PER"
        " TARGET of t.fits spectra; nothing is ever written in it. The output"
        " root is proposed beside it once it is chosen, and the config is the"
        " one that came with this installation.",
    "log_other_clone":
        "the configuration you point at lives in ANOTHER copy of this package"
        " (%s), while the code running is %s. The run uses the code running;"
        " the settings shown are that other copy's. Point the config at the"
        " same place unless you mean to mix them.",
    "log_no_config":
        "no config.yaml came with this installation: Browse to one, or run the"
        " window from a checkout.",
    "lbl_title": "LBL settings",
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
    "log_watch_added": "new in the data root: %s. Reading it.",
    "log_watch_gone": "no longer in the data root: %s",
    "log_watch_grew": "more spectra in %s than a moment ago. Reading them.",
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
    "log_no_report":
        "no compilation PDF yet at %s. The figures stage writes it, and the LBL"
        " stage adds the velocity pages to it",
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
    "opt_weight": "correction fit metric",
    "help_weight":
        "The metric the correction's amplitudes are measured in. `flux`, the"
        " nominal: every sample as the fit saw it. `velocity`: each sample"
        " weighted by the star's own derivative there, (dT/dv)^2, because what"
        " a contaminant does to a radial velocity is its overlap with that"
        " derivative, and a contaminant flat where the star has structure moves"
        " no line. It implies a refit of the amplitudes for every exposure, so"
        " the correction is slower; it changes nothing else, not what is"
        " divided out, not which samples are blanked, not the shrinkage."
        " Measured on TOI-2120, where the correction gains a factor three, the"
        " two are indistinguishable: 15.4 +- 1.2 against 15.0 +- 1.2 m/s. The"
        " targets where the correction COSTS are the ones that will decide.",
    "opt_width_kms": "high pass (km/s)", "opt_dv": "grid step (km/s)",
    "opt_nightly_stack": "coadd each night", "opt_run": "run LBL (hours)",
    "opt_lbl_prepare": "write LBL's tree",
    "opt_lbl_before": "measure the delivered spectra",
    "opt_lbl_after": "measure the corrected spectra",
    "opt_lbl_star_template": "use our star as LBL's template",
    "opt_lbl_strpca": "extra star components as RESPROJ",
    "opt_lbl_suffix": "corrected name",
    "opt_lbl_teff": "effective temperature", "opt_lbl_template": "template file",
    "opt_lbl_steps": "steps", "opt_lbl_link": "spectra in it as",
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
        "How the spectra get into LBL's science folders. `symlink`, the"
        " default, puts a link there: nothing is stored twice, and LBL only"
        " reads them. `copy` puts every spectrum there a second time, tens of"
        " gigabytes for a campaign, and a folder that no longer needs the data"
        " disk. A disk that cannot hold a link gets copies whatever this says,"
        " and the run says so at the top.",
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
    "help_openpdf_button":
        "Opens this run's own PDF, the one everything is bound into: the"
        " spectra before and after, the components, the correlations, and the"
        " velocity pages at the end once LBL has measured them. It is"
        " <object>_<tag>.pdf in the run's folder, written by the figures"
        " stage, so it is there once that stage has run.",
    "help_lang":
        "Choose the window's language: click the one you want. English,"
        " French, Spanish and Portuguese.",
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
    "fits_dir": "disque des produits (optionnel)",
    "help_fits_dir":
        "Où les produits d'un passage sont CONSERVÉS : chaque dossier de"
        " passage devient un lien vers ce disque, donc les spectres corrigés,"
        " l'ajustement et les figures y vont et non sur le disque interne."
        " Vide, tout reste sous le dossier de sortie, ce qui fonctionne. Un"
        " chemin est une chose propre à UNE machine, donc ce champ est vidé dès"
        " que ce qu'il nomme n'est pas là, et un passage dont le disque manque"
        " s'arrête plutôt que de remplir le disque interne sans le dire.",
    "log_no_disk":
        "la configuration nomme %s comme disque où conserver les produits, et"
        " il n'est pas là : le champ est vide, donc ce passage garderait tout"
        " sous le dossier de sortie. Remplissez-le si le disque doit être"
        " monté.",
    "out_dir": "dossier de sortie (optionnel)",
    "lbl_dir": "dossier de sortie LBL",
    "help_lbl_dir":
        "L'arbre du LBL (son DATA_DIR) : les dossiers science qu'il lit, et les"
        " gabarits, masques, tables raie par raie et rdb qu'il écrit. Proposé"
        " comme lbl sous le dossier de sortie, il suit ce dossier tant que vous"
        " n'en tapez pas un autre. Un seul arbre pour tous les passages d'un"
        " dossier, puisque le LBL de l'objet livré est le même pour tous et"
        " prend des heures ; pointé sur un arbre existant, cet arbre est"
        " utilisé tel quel.",
    "browse": "Parcourir", "rescan": "Relire",
    "objects": "objets", "settings": "réglages", "stages": "étapes",
    "command": "la commande qui sera lancée",
    "output": "sortie",
    "col_object": "objet", "col_files": "fichiers", "col_instrument": "instrument",
    "col_snr": "SNR", "col_exptime": "pose (s)", "col_mag": "mag",
    'run_name': 'nom de la réduction', 'auto': 'auto',
    'help_auto_button':
        'Repropose un nom à partir des cibles et des réglages tels'
        " qu'ils sont : les cibles, puis six caractères d'une empreinte de tout"
        ' ce qui fait de cette réduction un résultat différent. Les mêmes'
        ' paramètres donnent les mêmes six, un nombre changé en donne six'
        ' autres.',
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
    "openpdf": "Ouvrir le PDF de compilation",
    "savedefaults": "Enregistrer comme défauts...",
    "all": "tout", "none": "rien",
    "idle": "au repos", "running": "en cours",
    "lang": "\U0001F1EB\U0001F1F7 Français",
    "snr_berv": "rapport signal sur bruit en fonction du BERV",
    "snr_none": "cochez une cible pour voir où sont ses meilleures nuits",
    "help_snr_berv":
        "Un point par spectre, dans la couleur de son étoile : le rapport"
        " signal sur bruit mesuré par APERO, en fonction de la vitesse"
        " barycentrique à laquelle il a été pris. L'histogramme du dessus dit"
        " quelles vitesses une campagne couvre ; celui-ci dit AVEC QUOI elle"
        " les couvre. L'ajustement pondère un spectre en 1/sigma^2, donc une"
        " plage couverte à une extrémité par les pires nuits d'une campagne"
        " n'est pas la plage que l'ajustement voit vraiment, et les deux"
        " référentiels se séparent moins bien que la couverture ne le promet.",
    "tab_targets": "  cibles  ", "tab_settings": "  réglages  ",
    "tab_lbl": "  LBL  ", "tab_run": "  analyse  ",
    "tab_clean": "  nettoyage  ",
    "clean_title": "ce que pca2d a laissé sur ces disques",
    "help_clean":
        "Chaque endroit où ce programme met des octets, mesuré. Le cache des"
        " cubes, le brouillon d'un ajustement trop gros pour la mémoire et les"
        " journaux de LBL peuvent partir à tout moment : ils se refont à partir"
        " de ce qui reste ici, et un bouton les vide tous. Le reste peut partir"
        " aussi, ligne par ligne : choisissez-la dans la liste et Effacer la"
        " sélection. Les ajustements, les templates, les masques et les tables"
        " raie par raie ne reviennent qu'en relançant les heures qui les ont"
        " faits ; les spectres corrigés, les rapports et les vitesses qu'en"
        " relançant tout depuis les spectres. La fenêtre dit lesquelles de ces"
        " heures vous êtes sur le point de dépenser avant d'effacer quoi que ce"
        " soit. Un dossier déplacé sur un autre disque est compté là où il est"
        " vraiment ; un lien vers les spectres de quelqu'un d'autre ne pèse"
        " rien, puisque l'effacer ne libère rien.",
    "clean_measure": "Mesurer", "clean_purge": "Effacer ce qui peut partir",
    "clean_selected": "Effacer la sélection",
    "help_clean_measure":
        "Parcourt chaque dossier et additionne ce qu'il contient. Sur un disque"
        " réseau cela prend quelques secondes, et cela tourne à côté : la"
        " fenêtre reste utilisable pendant le compte.",
    "help_clean_purge":
        "Vide les dossiers marqués brouillon ou refaisable, et rien d'autre :"
        " ceux dont le retour ne coûte que du temps. La fenêtre demande"
        " d'abord, et annonce le total qu'elle va libérer. Les dossiers"
        " eux-mêmes restent : le prochain passage s'attend à les trouver. Pour"
        " tout le reste, choisissez les lignes et Effacer la sélection.",
    "help_clean_selected":
        "Efface les lignes choisies dans la liste, quelle que soit leur nature :"
        " un cache de cubes, les templates d'un après-midi, les vitesses"
        " elles-mêmes. Majuscule ou commande en choisit plusieurs. La fenêtre"
        " nomme chaque dossier et dit ce que les récupérer coûterait avant"
        " d'effacer, et c'est tout ce qu'il y a entre vous et un disque qui"
        " respire.",
    "clean_measuring": "mesure en cours…",
    "clean_totals": "%s en tout, dont %s peuvent être libérés",
    "clean_none": "rien de mesuré encore : appuyez sur Mesurer",
    "clean_confirm_title": "effacer les fichiers refaisables",
    "clean_confirm":
        "Sur le point de libérer %s à %d endroits :\n\n%s\n\nRien ici n'est"
        " un résultat : les cubes se relisent depuis les spectres, et les"
        " journaux de LBL sont réécrits par LBL. On y va ?",
    "clean_freed": "%s libérés à %d endroits",
    "clean_nothing": "rien à libérer : ni brouillon ni cache ici",
    "clean_pick": "choisissez d'abord une ligne ou plusieurs dans la liste :"
                  " le bouton efface ce qui est choisi",
    "clean_while_running":
        "UN PASSAGE EST EN COURS, et il lit dans ces dossiers. Il reconstruit"
        " ce qu'il ne trouve plus : effacer maintenant coûte à ce passage le"
        " temps de le refaire, au milieu de l'étape où il en est.",
    "clean_confirm_pick_title": "effacer la sélection",
    "clean_confirm_pick":
        "Sur le point d'effacer %s à %d endroits :\n\n%s\n\n%s\n\nOn y va ?",
    "clean_cost_rebuildable":
        "Tout cela se refait à partir de ce qui reste ici : des minutes, et un"
        " passage qui relit les spectres une fois de plus.",
    "clean_cost_expensive":
        "Une partie ne revient qu'en relançant ce qui l'a faite : un"
        " ajustement, ou un passage de LBL sur chaque exposition. Des heures,"
        " pas des minutes.",
    "clean_cost_results":
        "UNE PARTIE EST UN RÉSULTAT : des spectres corrigés, un rapport, ou"
        " les vitesses d'un rdb. Rien ici ne les refait. Seul un passage"
        " complet, depuis les spectres, les refait, et ce sont les heures que"
        " cela a pris la première fois.",
    "col_size": "taille", "col_nfiles": "fichiers", "col_kind": "nature",
    "kind_scratch": "brouillon", "kind_rebuildable": "refaisable",
    "kind_expensive": "coûteux", "kind_results": "résultats",
    # ce qu'est chaque ligne de la liste, dit quand la ligne est choisie
    "clean_cache":
        "Les cubes : tous les spectres d'une campagne sur une même grille de"
        " longueurs d'onde, pour qu'un second passage ne les relise pas tous."
        " D'habitude le plus gros poste ici, et jamais un résultat. Se refait"
        " en quelques minutes par campagne.",
    "clean_spill":
        "Le brouillon mappé d'un ajustement trop gros pour la mémoire. Plus"
        " rien ne le lit une fois l'ajustement fini : ce qui reste ici"
        " appartient à un passage interrompu.",
    "clean_results":
        "Les spectres corrigés, les rapports et les ajustements. Ce pour quoi"
        " tout cela existe ; jamais proposé à l'effacement.",
    "clean_lbl_science":
        "Les liens que LBL suit pour mesurer, un par pose. Refaits par l'étape"
        " lbl. Ce sont des liens symboliques : les effacer ne libère presque"
        " rien et ne touche pas aux spectres visés.",
    "clean_lbl_plots":
        "Les figures de LBL. Redessinées à chaque passage de LBL.",
    "clean_lbl_log": "Les journaux de LBL. Personne d'autre qu'un humain ne"
                     " les lit.",
    "clean_lbl_lblrv":
        "Les tables de vitesses raie par raie de LBL, un fichier par"
        " exposition, et ce que LBL écrit de plus gros. Cela ne se refait qu'en"
        " relançant LBL, ce qui prend des heures.",
    "clean_lbl_templates":
        "Les templates construits par LBL, un par objet et par passage. Ne se"
        " refont qu'en relançant l'étape template de LBL.",
    "clean_lbl_masks":
        "Les masques de raies construits par LBL, un par objet et par passage.",
    "clean_lbl_models": "Les modèles de LBL.",
    "clean_lbl_calib": "Les calibrations de LBL.",
    "clean_lbl_lblreftable":
        "Les tables de référence de LBL, une par objet et par passage.",
    "clean_lbl_lblrdb":
        "Les vitesses : le rdb dont chaque page RV de chaque rapport est tirée."
        " Un résultat, et la plus petite chose de cette liste : l'effacer ne"
        " libère presque rien et jette la mesure.",
    "clean_pycache":
        "Du Python compilé, refait au prochain import du paquet.",
    "quit": "Quitter", "quit_title": "quitter pca2d-preclean",
    "quit_yes": "Quitter quand même", "quit_no": "Rester",
    "quit_running":
        "Un passage est en cours, et c'est un processus fils de cette fenêtre :"
        " quitter l'arrête. Ce que les étapes précédentes ont écrit reste en"
        " place, l'étape en cours est perdue. Quitter quand même ?",
    "help_quit_button":
        "Ferme la fenêtre. Les réglages sont écrits à chaque changement, rien"
        " n'est donc perdu en partant ; un passage en cours est arrêté, et la"
        " fenêtre le demande avant.",
    "pick_root": "choisissez un dossier de données : Parcourir, en haut",
    "command_pending":
        "cochez une cible : la commande s'écrit ici, en entier, avant de partir",
    "log_pick_root":
        "aucun dossier de données pour l'instant. Choisissez le dossier qui"
        " contient UN DOSSIER PAR CIBLE de spectres t.fits ; rien n'y est jamais"
        " écrit. Le dossier de sortie est proposé à côté une fois celui-ci"
        " choisi, et la configuration est celle livrée avec cette installation.",
    "log_other_clone":
        "la configuration que vous visez est dans UNE AUTRE copie de ce paquet"
        " (%s), alors que le code qui tourne est %s. Le passage utilise le code"
        " qui tourne ; les réglages affichés sont ceux de l'autre copie. Visez"
        " le même endroit, sauf si vous voulez vraiment les mélanger.",
    "log_no_config":
        "aucun config.yaml n'est livré avec cette installation : choisissez-en"
        " un, ou lancez la fenêtre depuis un dépôt cloné.",
    "lbl_title": "réglages LBL",
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
    "log_watch_added": "nouveau dans le dossier de données : %s. Lecture.",
    "log_watch_gone": "n'est plus dans le dossier de données : %s",
    "log_watch_grew": "plus de spectres dans %s qu'il y a un instant. Lecture.",
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
    "log_no_report":
        "pas encore de PDF de compilation à %s. C'est l'étape des figures qui"
        " l'écrit, et l'étape LBL qui y ajoute les pages de vitesses",
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
    "opt_weight": "métrique d'ajustement de la correction",
    "help_weight":
        "La métrique dans laquelle les amplitudes de la correction sont"
        " mesurées. `flux`, le nominal : chaque échantillon tel que"
        " l'ajustement l'a vu. `velocity` : chaque échantillon pondéré par la"
        " dérivée de l'étoile à cet endroit, (dT/dv)^2, parce que ce qu'un"
        " contaminant fait à une vitesse radiale est son recouvrement avec"
        " cette dérivée, et qu'un contaminant plat là où l'étoile a de la"
        " structure ne déplace aucune raie. Cela implique de réajuster les"
        " amplitudes de chaque pose, donc la correction est plus lente ; rien"
        " d'autre ne change, ni ce qui est divisé, ni les échantillons"
        " blanchis, ni le rétrécissement. Mesuré sur TOI-2120, où la correction"
        " gagne un facteur trois, les deux sont indiscernables : 15,4 +- 1,2"
        " contre 15,0 +- 1,2 m/s. Ce sont les cibles où la correction COÛTE qui"
        " trancheront.",
    "opt_width_kms": "passe-haut (km/s)", "opt_dv": "pas de grille (km/s)",
    "opt_nightly_stack": "empiler chaque nuit", "opt_run": "lancer LBL (heures)",
    "opt_lbl_prepare": "écrire l'arbre du LBL",
    "opt_lbl_before": "mesurer les spectres livrés",
    "opt_lbl_after": "mesurer les spectres corrigés",
    "opt_lbl_star_template": "notre étoile comme gabarit du LBL",
    "opt_lbl_strpca": "composantes stellaires en RESPROJ",
    "opt_lbl_suffix": "nom du corrigé",
    "opt_lbl_teff": "température effective", "opt_lbl_template": "fichier gabarit",
    "opt_lbl_steps": "étapes", "opt_lbl_link": "spectres dedans en",
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
        " `symlink`, le défaut, y met un lien : rien n'est stocké deux fois, et"
        " le LBL ne fait que les lire. `copy` y met chaque spectre une seconde"
        " fois, des dizaines de gigaoctets pour une campagne, et un dossier qui"
        " n'a plus besoin du disque de données. Un disque incapable de porter"
        " un lien reçoit des copies quoi que dise ce choix, et le passage le"
        " dit en tête.",
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
    "help_openpdf_button":
        "Ouvre le PDF de ce passage, celui où tout est relié : les spectres"
        " avant et après, les composantes, les corrélations, et les pages de"
        " vitesses à la fin une fois que LBL les a mesurées. C'est"
        " <objet>_<tag>.pdf dans le dossier du passage, écrit par l'étape des"
        " figures : il est là dès que cette étape a tourné.",
    "help_lang":
        "Choisit la langue de la fenêtre : cliquez sur celle que vous voulez."
        " Anglais, français, espagnol et portugais.",
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

from .gui_es import ES  # noqa: E402
from .gui_pt import PT  # noqa: E402

TEXTS = {"en": EN, "fr": FR, "es": ES, "pt": PT}
#: the order the language buttons are in, and switch_language cycles
LANGUAGES = ("en", "fr", "es", "pt")


def text(lang, key, default=None):
    """The label or the explanation, in the window's language."""
    return TEXTS.get(lang, EN).get(key, default if default is not None else key)


#: The stages in the order they happen, which is the order they depend on each
#: other in: the cube feeds the fit, the fit feeds the correction, and LBL
#: measures what the correction wrote. `figures` is deliberately not in it: it
#: draws what the fit left, and nothing waits for a drawing.
CHAIN = ("cube", "fit", "correct", "lbl")


def follow_stages(ticked, changed):
    """The stage boxes after `changed` was just ticked or unticked.

    Unticking one unticks everything AFTER it in the chain: nothing downstream
    has its input any more, and a run that says it will correct spectra it is
    not fitting is a run that will fail twenty minutes in. Ticking `lbl` ticks
    `correct`, since LBL measures what the correction writes.

    Ticking does NOT pull the rest of the chain: the fit is always redone when
    it is asked for, so correcting again with the fit that is already there is
    a real thing to want, and the window has to let it be said.
    """
    out = dict(ticked)
    if changed in CHAIN and not out.get(changed):
        for name in CHAIN[CHAIN.index(changed) + 1:]:
            out[name] = False
    if changed == "lbl" and out.get("lbl"):
        out["correct"] = True
    return out


#: how often the window looks at the data root for folders that appeared
WATCH_MS = 10000


def folder_news(found, shown, seen=None):
    """What a look at the data root found that the window does not have.

    `found` is {campaign: files in its folder}, `shown` what the list holds,
    `seen` the previous look. Returns the campaigns that APPEARED, the ones that
    are GONE, and the ones whose file count MOVED since the last look.

    The count is compared with the previous look and never with the list,
    because the list counts what the index holds: a spectrum the scan could not
    read is missing from it for good, and comparing the two would ask for a
    rescan every ten seconds for ever.
    """
    added = sorted(set(found) - set(shown))
    gone = sorted(set(shown) - set(found))
    grown = sorted(name for name in found
                   if seen and name in seen and found[name] != seen[name])
    return added, gone, grown


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


def run_folder(state, out_root):
    """The folder the run these settings describe writes into, or None.

    The same three pieces cli.resolve joins, in the same order: the output root
    with `_<name>` under it when the run is named, then the object, or
    `joint/<A+B>` when several are fitted together, then the tag. The tag is
    the two counts AND the `v` of the velocity term: a run with it and a run
    without it are two folders on purpose, and a tag without the v points at
    the other one's.
    """
    names = state.get("objects") or []
    if not names:
        return None
    tag = "%s-%s%s" % (state.get("n_star") or 0, state.get("n_earth") or 3,
                       "v" if state.get("velocity_term") else "")
    named = str(state.get("run_name") or "").strip()
    where = os.path.join(absolute(out_root), "_" + named if named else "")
    if len(names) > 1:
        return os.path.join(where, "joint", "+".join(names), tag)
    return os.path.join(where, names[0], tag)


def report_pdf(state, out_root):
    """The run's compilation PDF: <folder>/<object>_<tag>.pdf, or None.

    Named after the first object even when several were fitted together, which
    is how cli.main names it.
    """
    folder = run_folder(state, out_root)
    if not folder:
        return None
    return os.path.join(folder, "%s_%s.pdf" % ((state.get("objects") or [""])[0],
                                               os.path.basename(folder)))


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


def package_home():
    """The folder holding the `pca2d` package this window is running from.

    What it is for: the run is spawned with the config file's folder as its
    working directory, and `python -m` puts that folder first on the import
    path. If it happens to hold another copy of the package, the run is that
    copy. That is not a hypothetical: on 2026-09-15 a window that had just
    written --no-fits-dir launched a pipeline that had never heard of it,
    because the config it was pointed at lived in a second clone.
    """
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def preclean_argv():
    """How to run the pipeline: THIS window's own code, through this python.

    `[sys.executable, "-m", "pca2d.cli"]` and not the `pca2d-preclean` script
    beside the interpreter, which is what this used to return. The two are the
    same thing only when one copy of the package is installed. With two
    editable installs, which happens the moment a second clone is pip-installed
    for a test, WHICH ONE RUNS DEPENDS ON THE DIRECTORY: on 2026-09-15 a window
    that knew --no-fits-dir launched a pipeline that did not, and the run died
    on `unrecognized arguments` after the window had written the flag itself.
    The module is imported the way the window imported it, so the two are one
    version by construction.

    The command SHOWN stays `pca2d-preclean ...`, which is what a person types;
    this is only what gets spawned.
    """
    try:
        import pca2d.cli                                        # noqa: F401
    except Exception:                                           # noqa: BLE001
        beside = os.path.join(os.path.dirname(os.path.abspath(sys.executable)),
                              "pca2d-preclean")
        if os.path.exists(beside) and os.access(beside, os.X_OK):
            return [beside]
        found = shutil.which("pca2d-preclean")
        return [found] if found else None
    # -P, and it is the whole point: `python -m` puts the WORKING DIRECTORY
    # first on the import path, ahead of PYTHONPATH, so pinning this package
    # through the environment is not enough on its own. The run starts in the
    # config file's folder, and a folder holding another pca2d/ then wins.
    # Measured from a second clone: `-m` alone ran that clone, `-P -m` ran this
    # one. Python 3.11 and later; before that, the environment is all there is.
    safe = ["-P"] if sys.version_info >= (3, 11) else []
    return [sys.executable] + safe + ["-m", "pca2d.cli"]


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


def build_command(state, defaults=None):
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
    # An empty field is a DECISION, keep the products under the output root,
    # and it has its own flag: `--fits-dir ""` is not a command anybody can
    # paste, since a shell drops the empty word and argparse then asks for the
    # argument it was promised.
    disk = str(state.get("fits_dir") or "").strip()
    if disk:
        argv += ["--fits-dir", disk]
    elif "fits_dir" in state:
        argv += ["--no-fits-dir"]
    for flag, key in (("--min-rjd", "min_rjd"), ("--max-rjd", "max_rjd")):
        value = str(state.get(key) or "").strip()
        if value:
            argv += [flag, value]
    if state.get("n_star") not in (None, ""):
        argv += ["--n-star", str(state["n_star"])]
    if state.get("n_earth") not in (None, ""):
        argv += ["--n-earth", str(state["n_earth"])]
    # EVERY other setting this window can change, each under its own flag.
    # Until 2026-09-15 only the two counts travelled: the high pass, the
    # shrinkage, the sweeps, the grid step, the coadding and the velocity term
    # were shown, changed, and then silently ignored by the run, which used the
    # configuration's own values. A flag also puts what was asked in the run's
    # own log, where a config edited afterwards would not be.
    for key, flag in (("velocity_term", "--velocity-term"),
                      ("iters", "--iters"), ("shrink", "--shrink"),
                      ("weight", "--weight"), ("width_kms", "--high-pass"),
                      ("dv", "--dv"), ("nightly_stack", "--nightly-stack")):
        if key not in state:
            continue
        value = state[key]
        if isinstance(value, bool):
            value = "true" if value else "false"
        value = str(value).strip()
        if value:
            argv += [flag, value]
    argv += lbl_flags(state, defaults)
    stages = [s for s in STAGES if state.get("stage_" + s)]
    if stages and len(stages) != len(STAGES):
        argv += ["--stages", ",".join(stages)]
    if state.get("dry_run"):
        argv.append("--dry-run")
    return argv


def lbl_flags(state, defaults=None):
    """The LBL folder, and the LBL settings the window changed, as flags.

    The folder always travels when there is one, like the output root: it is a
    path on this machine, and the command should say where LBL will write. The
    settings travel when they differ from the configuration (`defaults`), or
    all of them when there is no configuration to compare with; a command that
    repeats a dozen values the configuration already holds says nothing more
    and is harder to read.
    """
    argv = []
    tree = str(state.get("lbl_dir") or "").strip()
    if tree:
        argv += ["--lbl-dir", tree]
    changed = (variant_yaml(state, defaults).get("lbl") or {}
               if defaults is not None else None)
    for key, path, kind in OPTIONS_LBL + OPTIONS_PATHS:
        flag = LBL_FLAGS.get(key)
        if flag is None or key not in state or state[key] in (None, ""):
            continue
        if changed is not None and path.split(".")[1] not in changed:
            continue
        value = state[key]
        if kind == "bool":
            value = "true" if value else "false"
        elif kind == "list":
            value = ",".join(str(value).replace(",", " ").split())
        argv += [flag, str(value).strip()]
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


#: How much of a campaign has been read, as a slice of a small disc, beside the
#: name in the list. Nothing at all once every file is in: a finished list must
#: not be a column of symbols, so the mark means "still reading" by existing.
PIE = ("\u25cb", "\u25d4", "\u25d1", "\u25d5", "\u25cf")


def pie_glyph(known, total):
    """The slice for `known` of `total` files read, or "" when there is no more.

    Five steps rather than a number: what a row needs to say while a scan runs
    is whether its numbers are nearly settled, not how many headers are left.
    """
    try:
        known, total = int(known or 0), int(total or 0)
    except (TypeError, ValueError):
        return ""
    if total <= 0 or known >= total:
        return ""
    step = int(max(0.0, min(1.0, known / float(total))) * len(PIE))
    return PIE[min(len(PIE) - 1, step)]


def pie_step(known, total, steps=12):
    """Which step of a `steps`-step disc `known` of `total` files is, or None.

    None where there is nothing to say: a folder with no files, or one whose
    files are all read. A disc that never goes away would be decoration.
    """
    try:
        known, total = int(known or 0), int(total or 0)
    except (TypeError, ValueError):
        return None
    if total <= 0 or known >= total:
        return None
    fraction = max(0.0, min(1.0, known / float(total)))
    return min(steps - 1, int(fraction * steps))


def rounded_image(width, height, radius, fill, outline=None, thickness=1,
                  corners=(True, True, True, True)):
    """A rounded rectangle as a PIL image, or None without PIL.

    Supersampled and shrunk, like the row marks: ImageDraw has no antialiasing,
    and a corner drawn at its final size is a staircase. What it is FOR is the
    nine-patch borders ttk needs, since no theme here draws a round corner.
    """
    try:
        from PIL import Image, ImageDraw
    except Exception:                                           # noqa: BLE001
        return None
    k = max(1, int(SUPERSAMPLE))
    image = Image.new("RGBA", (width * k, height * k), (0, 0, 0, 0))
    draw = ImageDraw.Draw(image)
    draw.rounded_rectangle([0, 0, width * k - 1, height * k - 1],
                           radius=radius * k, fill=fill,
                           outline=outline, width=thickness * k if outline else 0,
                           corners=corners)
    return image.resize((width, height), Image.LANCZOS)


def rounded_photo(tk_module, *args, **kwargs):
    """`rounded_image`, as something a ttk element can be made of."""
    drawing = rounded_image(*args, **kwargs)
    if drawing is None:
        return None
    from PIL import ImageTk

    return ImageTk.PhotoImage(drawing)


#: how many times over the marks are drawn before being shrunk to their size.
#: ImageDraw has no antialiasing: a circle drawn at 18 pixels IS a staircase,
#: and the only way to round it is to draw it large and shrink it with a filter
#: that averages.
SUPERSAMPLE = 8


def row_drawing(ticked, step=None, steps=12, size=20, pie=18, gap=6,
                edge=None, tick=None, ground=None, dial=None):
    """The tick box, and the disc of a campaign still being read, as one PIL
    image, or None without PIL.

    One image, because a Treeview row has exactly one slot for one, before its
    text. Drawn rather than written for the same reason in both cases: a
    character is the size of the text beside it, and neither the thing aimed at
    with a mouse nor the fraction read across the room should be that small.

    Everything is drawn SUPERSAMPLE times over and shrunk with Lanczos at the
    end, which is where the smooth edges come from.
    """
    try:
        from PIL import Image, ImageDraw
    except Exception:                                           # noqa: BLE001
        return None
    edge, tick = edge or LINE, tick or ACCENT
    ground, dial = ground or SURFACE, dial or MUTED
    width = size if step is None else size + gap + pie
    k = max(1, int(SUPERSAMPLE))
    image = Image.new("RGBA", (width * k, size * k), (0, 0, 0, 0))
    draw = ImageDraw.Draw(image)
    draw.rounded_rectangle([1 * k, 1 * k, (size - 2) * k, (size - 2) * k],
                           radius=4 * k, outline=tick if ticked else edge,
                           width=max(1, int(1.6 * k)), fill=ground)
    if ticked:
        draw.line([(size * 0.27 * k, size * 0.53 * k),
                   (size * 0.44 * k, size * 0.71 * k),
                   (size * 0.75 * k, size * 0.30 * k)],
                  fill=tick, width=max(1, int(1.8 * k)), joint="curve")
    if step is not None:
        left = (size + gap) * k
        top = ((size - pie) // 2) * k
        box = [left, top, left + pie * k - 1, top + pie * k - 1]
        draw.ellipse(box, outline=dial, width=max(1, int(1.1 * k)), fill=ground)
        if step > 0:
            # from twelve o'clock, clockwise, like any other dial
            draw.pieslice(box, start=-90, end=-90 + 360.0 * step / steps,
                          fill=tick, outline=tick)
    return image.resize((width, size), Image.LANCZOS)


def row_image(tk, ticked, step=None, **kwargs):
    """`row_drawing`, as something a Treeview row can hold."""
    drawing = row_drawing(ticked, step, **kwargs)
    if drawing is None:
        return None
    from PIL import ImageTk

    return ImageTk.PhotoImage(drawing)


def suggested_run_name(state, digits=6):
    """A name for this reduction: the targets, then a short hash of the rest.

    Two things have to be read off a folder name in an output root months
    later: WHAT was reduced, and whether it was reduced the same way as the one
    beside it. The targets give the first, and no list of settings short enough
    to be a folder name gives the second, so the settings are hashed: same
    parameters, same six characters; one number different anywhere, a different
    six. Unique in the only sense that matters here, which is that a clash is
    improbable rather than impossible.

    The hash covers what makes a run a different result: the targets, the
    component counts, the velocity term, the sweeps, the shrinkage, the high
    pass, the grid step, the nightly coadding, and the date window.
    """
    import hashlib

    # sorted, like the command line's own joint name: the same set of targets
    # in another order is the same reduction
    names = sorted(str(n).strip() for n in (state.get("objects") or [])
                   if str(n).strip())
    if not names:
        head = "run"
    elif len(names) <= 3:
        head = "+".join(names)
    else:
        # three names is already a long folder name; past that, say how many
        head = "%s+%d" % ("+".join(names[:2]), len(names) - 2)
    head = re.sub(r"[^0-9A-Za-z._+-]", "_", head)[:40].strip("_+") or "run"

    keys = ("n_star", "n_earth", "velocity_term", "iters", "shrink",
            "width_kms", "dv", "nightly_stack", "min_rjd", "max_rjd")
    payload = "|".join(["+".join(sorted(names))] +
                       ["%s=%s" % (key, state.get(key, "")) for key in keys])
    short = hashlib.blake2b(payload.encode("utf-8"),
                            digest_size=8).hexdigest()[:digits]
    return "%s_%s" % (head, short)


def command_line(state, pending, defaults=None):
    """What the command box shows: the command, or what is still missing.

    `pca2d-preclean --config ... --n-star 1` with no object in it is not a
    command anybody can paste, and showing it while nothing is ticked made a
    window with no data root look like a run that was ready to go.
    """
    return (" ".join(build_command(state, defaults)) if state.get("objects")
            else pending)


def installed_config():
    """The config.yaml that came with THIS installation, or "" if there is none.

    Not `config.yaml` resolved against the working directory: a window started
    from a home folder would then open on a file that does not exist, or worse
    on somebody else's. The nominal to open on is the one belonging to the code
    that is running, which sits beside the package in a checkout and inside it
    in a copy that ships one.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    for path in (os.path.join(os.path.dirname(here), "config.yaml"),
                 os.path.join(here, "config.yaml")):
        if os.path.isfile(path):
            return path
    return ""


def opening_paths(saved, config_value=None):
    """The paths a window opens on: what was kept, else a fresh start.

    A fresh start is the installation's own config.yaml and NO roots. The data
    root and the output root are a choice about somebody's disks, and a window
    that opens with a plausible-looking path in them invites a run against a
    folder nobody picked. Empty, and said in the status line, is the honest
    opening state.

    The products disk (`output.fits_directory`) is the same thing one step
    further: it is a path on ONE machine, so a configuration that carries one
    hands every fresh install a disk it has never heard of, and a run that
    finds it missing stops. It is taken from what this window kept, else from
    the configuration, and **blanked when it is not there**, whichever it came
    from: an empty field is a run that keeps its products under the output
    root, which is a thing that works.
    """
    disk = saved.get("fits_dir")
    if disk is None:
        disk = config_value or ""
    disk = absolute(disk)
    if disk and not os.path.isdir(disk):
        disk = ""
    out_dir = saved.get("out_dir") or ""
    return {"data_dir": absolute(saved.get("data_dir") or ""),
            "config": absolute(saved.get("config") or installed_config()),
            "out_dir": out_dir,
            "lbl_dir": saved.get("lbl_dir") or lbl_proposal(out_dir),
            "fits_dir": disk}


def lbl_proposal(out_dir):
    """LBL's tree as the window proposes it: `lbl` under the output root.

    Empty while there is no output root, since a tree under nothing is the
    `lbl` beside wherever the run happens to start, which is how SMETHELLS_20's
    velocities came to be in another clone's folder (2026-09-15).
    """
    out_dir = str(out_dir or "").strip()
    return os.path.join(absolute(out_dir), "lbl") if out_dir else ""


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
        self._build_header(root)
        # Three pages, because everything at once on one page was a wall: WHAT
        # is being reduced, HOW, and the run itself with its log. Each of them
        # then has room to be read, and the window can be any size.
        book = ttk.Notebook(root)
        book.pack(fill="both", expand=True, padx=12, pady=(2, 8))
        self.tabs = []
        pages = {}
        for key in ("tab_targets", "tab_settings", "tab_lbl", "tab_run",
                    "tab_clean"):
            page = ttk.Frame(book, padding=8)
            book.add(page, text=self.t(key))
            self.tabs.append((book, page, key))
            pages[key] = page
        self.book = book

        # what is being reduced: the roots, the campaigns, and what they cover
        self._build_paths(pages["tab_targets"])
        panes = ttk.Panedwindow(pages["tab_targets"], orient="horizontal")
        panes.pack(fill="both", expand=True)
        left, right = ttk.Frame(panes), ttk.Frame(panes)
        panes.add(left, weight=3)
        panes.add(right, weight=4)
        # the list on one side, what the ticked campaigns cover on the other:
        # the two panels under the list left the right half of the page empty
        self._build_objects(left)
        self._build_berv(right)
        self._build_snr(right)
        self._build_time(right)

        # how: every setting that changes a result, and the buttons that write
        # them down somewhere
        self._build_options(pages["tab_settings"])
        self._build_buttons(pages["tab_settings"], "settings")

        # what LBL is asked, which is a setting like the others and used to be
        # a window of its own
        self._build_lbl(pages["tab_lbl"])

        # and the run: what it will be called, the command in full, the buttons
        # that start and end it, and everything it says
        self._build_name(pages["tab_run"])
        self._build_command(pages["tab_run"])
        self._build_buttons(pages["tab_run"], "run")
        self._build_log(pages["tab_run"])

        # and what all of it has left on the disks, with the one button that
        # takes any of it away
        self._build_clean(pages["tab_clean"])
        self.stop_button.configure(state="disabled")
        self._propose_out()      # on opening, not only when the data root moves
        if not self.vars["config"].get().strip():
            # a copy installed without one: say it here rather than let a run
            # fail on a config.yaml resolved against the working directory
            self._say("log_no_config", level="warn")
        else:
            self._warn_other_clone(self.vars["config"].get())
        self.refresh_objects()
        self._drain_id = self.root.after(80, self._drain)
        self._seen = {}
        self.root.after(WATCH_MS, self._watch_root)
        self.root.protocol("WM_DELETE_WINDOW", self.quit_window)

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
        # the language buttons: a flag and two letters, as small as a button
        # can be and still be clicked
        style.configure("Lang.TButton", font=(body, 10), padding=(3, 0))
        style.configure("TNotebook", background=BG, borderwidth=0,
                        bordercolor=BG, lightcolor=BG, darkcolor=BG,
                        tabmargins=(2, 4, 2, 0))
        style.configure("TNotebook.Tab", background=BG, foreground=MUTED,
                        padding=(16, 8), font=(body, 12, "bold"),
                        borderwidth=0)
        style.map("TNotebook.Tab",
                  background=[("selected", SURFACE), ("active", ACCENT_SOFT)],
                  foreground=[("selected", ACCENT)],
                  # the theme shrinks the selected tab, 16 8 down to 6 4 6 2,
                  # and a row of tabs that changes size as it is clicked reads
                  # as the window flinching. One size, whatever is selected.
                  padding=[("selected", (16, 8)), ("active", (16, 8))])
        self._round(style)
        style.configure("Treeview", background=SURFACE, fieldbackground=SURFACE,
                        foreground=INK, rowheight=26, bordercolor=LINE,
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
        if self.scanning:
            return
        if self.proc:
            self.status.configure(text=self.t("running"))
        elif not self.vars["data_dir"].get().strip():
            # nothing can be run from here, and an empty table does not say why
            self.status.configure(text=self.t("pick_root"))
        else:
            self.status.configure(text=self.t("idle"))


    def _round(self, style):
        """Round the corners ttk draws square, out of images.

        No theme here draws a rounded corner, and the only way to have one is to
        hand ttk a picture of it: a nine-patch, which it stretches along its
        middle and leaves alone at the corners. Every image is kept on the
        instance, or the garbage collector takes it and the widget comes back
        empty. Without PIL none of this happens and the window is the square one
        it was.
        """
        tk = self.tk
        self._art = getattr(self, "_art", [])
        r = 9

        def art(fill, outline=None, corners=(True, True, True, True), h=30):
            image = rounded_photo(tk, 2 * r + 4, h, r, fill, outline=outline,
                                  thickness=1, corners=corners)
            if image is not None:
                self._art.append(image)     # kept, or it is collected and blank
            return image

        # A tab is round on top only, and the SELECTED one is the same ground as
        # the page under it: that is what makes it read as the page's own edge
        # rather than a white card floating over a grey one.
        top = (True, True, False, False)
        tabs = (art("#dfe6f1", "#dfe6f1", top), art(BG, LINE, top),
                art(ACCENT_SOFT, ACCENT_SOFT, top))
        buttons = (art(SURFACE, LINE), art(ACCENT_SOFT, ACCENT),
                   art("#f4f7fb", LINE))
        fields = (art(SURFACE, LINE, h=28), art(SURFACE, ACCENT, h=28))
        # the same box and tick the list draws, so one window speaks one language
        ticks = (row_image(tk, False, size=16), row_image(tk, True, size=16))
        for image in ticks:
            if image is not None:
                self._art.append(image)
        primary = (art(ACCENT, ACCENT), art("#0d7cc2", "#0d7cc2"),
                   art("#9dbdd4", "#9dbdd4"))
        if not all(tabs + buttons + primary + fields + ticks):
            return                          # no PIL: the square window, as before

        border = (r, r, r, 2)
        try:
            style.element_create(
                "round.tab", "image", tabs[0],
                ("selected", tabs[1]), ("active", tabs[2]),
                border=border, sticky="nsew")
            style.element_create(
                "round.button", "image", buttons[0],
                ("pressed", buttons[1]), ("active", buttons[1]),
                ("disabled", buttons[2]),
                border=(r, r, r, r), sticky="nsew")
            style.element_create(
                "primary.button", "image", primary[0],
                ("pressed", primary[1]), ("active", primary[1]),
                ("disabled", primary[2]),
                border=(r, r, r, r), sticky="nsew")
            style.element_create(
                "round.field", "image", fields[0], ("focus", fields[1]),
                border=(r, r, r, r), sticky="nsew")
            style.element_create(
                "round.check", "image", ticks[0], ("selected", ticks[1]),
                border=0, sticky="")
        except tk.TclError:
            return                          # a second window in one process
        style.layout("TNotebook.Tab", [
            ("round.tab", {"sticky": "nsew", "children": [
                ("Notebook.padding", {"side": "top", "sticky": "nsew",
                                      "children": [
                                          ("Notebook.label",
                                           {"side": "top", "sticky": ""})]})]})])
        for name, element in (("TButton", "round.button"),
                              ("Run.TButton", "primary.button")):
            style.layout(name, [
                (element, {"sticky": "nsew", "children": [
                    ("Button.padding", {"sticky": "nsew", "children": [
                        ("Button.label", {"sticky": "nsew"})]})]})])
        try:
            style.layout("TEntry", [
                ("round.field", {"sticky": "nsew", "children": [
                    ("Entry.padding", {"sticky": "nsew", "children": [
                        ("Entry.textarea", {"sticky": "nsew"})]})]})])
            style.layout("TCombobox", [
                ("round.field", {"sticky": "nsew", "children": [
                    ("Combobox.downarrow", {"side": "right", "sticky": "ns"}),
                    ("Combobox.padding", {"expand": "1", "sticky": "nsew",
                                          "children": [
                                              ("Combobox.textarea",
                                               {"sticky": "nsew"})]})]})])
            style.layout("TCheckbutton", [
                ("Checkbutton.padding", {"sticky": "nsew", "children": [
                    ("round.check", {"side": "left", "sticky": ""}),
                    ("Checkbutton.focus", {"side": "left", "sticky": "w",
                                           "children": [
                                               ("Checkbutton.label",
                                                {"sticky": "nsew"})]})]})])
        except tk.TclError:
            pass                            # the square field, then
        style.configure("TButton", padding=(12, 6))
        style.configure("Run.TButton", padding=(18, 6))
        style.configure("TEntry", padding=(8, 5))
        style.configure("TCombobox", padding=(8, 4))
        style.configure("TCheckbutton", padding=(0, 2))
        # What a rounded corner leaves transparent shows the widget's OWN
        # background, and these were white: every corner had a white notch
        # sitting outside it, on a page that is not white. The fill now comes
        # from the picture and the background is the window's ground, so the
        # corner has the page behind it and disappears into it. Only once the
        # images are really in use: without them the background IS the fill.
        for name in ("TButton", "Run.TButton", "TEntry", "TCombobox"):
            style.configure(name, background=BG, lightcolor=BG, darkcolor=BG,
                            bordercolor=BG)
        # the tabs have corners too, and the theme still painted the selected
        # one white behind its picture
        style.configure("TNotebook.Tab", background=BG, lightcolor=BG,
                        darkcolor=BG, bordercolor=BG)
        style.map("TNotebook.Tab",
                  background=[("selected", BG), ("active", BG)],
                  lightcolor=[("selected", BG), ("active", BG)],
                  foreground=[("selected", ACCENT), ("active", ACCENT)],
                  padding=[("selected", (16, 8)), ("active", (16, 8))])
        style.map("TButton", background=[("active", BG), ("pressed", BG),
                                         ("disabled", BG)])
        style.map("Run.TButton", background=[("active", BG), ("pressed", BG),
                                             ("disabled", BG)])

    def _register(self, widget, key, how="text"):
        self.labels.append((widget, key, how))
        return widget

    def _tip(self, widget, key):
        Tip(self, widget, key)
        return widget

    def switch_language(self):
        """The next language in LANGUAGES, round again after the last."""
        order = list(LANGUAGES)
        here = order.index(self.lang) if self.lang in order else -1
        self.set_language(order[(here + 1) % len(order)])

    def _draw_languages(self):
        """A button for each language the window is not in, flag and name.

        Drawn again at every change, since which languages are "the others"
        is what changed.
        """
        bar = getattr(self, "lang_bar", None)
        if bar is None:
            return
        for child in bar.winfo_children():
            child.destroy()
        for code in LANGUAGES:
            if code == self.lang:
                continue
            flag = TEXTS[code]["lang"].split()[0]
            button = self.ttk.Button(bar, text="%s %s" % (flag, code.upper()),
                                     style="Lang.TButton",
                                     command=lambda c=code: self.set_language(c))
            button.pack(side="left", padx=(0, 2))
            self._tip(button, "help_lang")

    def set_language(self, lang):
        """Say everything again in `lang`: every label, tab and heading."""
        self.lang = lang if lang in TEXTS else "en"
        self._draw_languages()
        for book, page, key in getattr(self, "tabs", []):
            try:
                book.tab(page, text=self.t(key))
            except Exception:                                   # noqa: BLE001
                pass
        for widget, key, how in self.labels:
            try:
                if how == "text":
                    widget.configure(text=self.t(key))
                elif how == "heading":
                    pass                  # redrawn together, below
            except Exception:                                 # noqa: BLE001
                pass
        self._draw_headings()
        if hasattr(self, "clean_tree"):
            self._draw_clean_headings()
            self._clean_rows()
        self._state()
        self._sync()

    # ---- widgets ------------------------------------------------------
    def _build_header(self, parent):
        """The banner: who made the spectra, what this is, and the language."""
        ttk = self.ttk
        frame = ttk.Frame(parent)
        frame.pack(fill="x", padx=14, pady=(10, 4))
        # APERO's own logo, then the name: every spectrum this window reads was
        # reduced by APERO, so the pipeline signs the top left corner. The file
        # is in the package, never fetched at run time.
        logo = self._logo("apero_logo.png", height=22)
        if logo is not None:
            self.logo_label = ttk.Label(frame, image=logo, background=BG)
            self.logo_label.image = logo      # or the garbage collector eats it
            self.logo_label.pack(side="left", padx=(0, 8))
            self._tip(self.logo_label, "help_apero")
        ttk.Label(frame, text="pca2d", style="Head.TLabel").pack(side="left")
        # Quit in the banner, not on one page: it ends the window, which is not
        # a property of whichever page happens to be open, and looking for it
        # meant going back to the run page to leave from anywhere else. Packed
        # first, so it sits at the far corner with the language button beside it.
        quit_button = ttk.Button(frame, text=self.t("quit"), width=10,
                                 command=self.quit_window)
        quit_button.pack(side="right")
        self._register(quit_button, "quit")
        self._tip(quit_button, "help_quit_button")
        # every other language, a small flag and its two letters, one click
        # each: with four a button that flips between two reaches none of them
        self.lang_bar = ttk.Frame(frame)
        self.lang_bar.pack(side="right", padx=(0, 6))
        self._draw_languages()
        self._register(ttk.Label(frame, style="Hint.TLabel", wraplength=780,
                                 justify="left", text=self.t("subtitle")),
                       "subtitle").pack(side="left", padx=10)

    def _build_paths(self, parent):
        """The two roots and the configuration, each shown in full."""
        ttk, tk = self.ttk, self.tk
        frame = ttk.Frame(parent)
        frame.pack(fill="x", pady=(0, 6))
        # shown in full, since a relative path means a different folder from a
        # different working directory and these data are reached through a link
        opening = opening_paths(self.saved, self._config_fits_dir())
        # the LBL folder follows the output root for as long as it is the
        # proposal, as the reduction name follows the targets (_follow_lbl_dir)
        proposal = lbl_proposal(opening["out_dir"])
        self._proposed_lbl = (opening["lbl_dir"]
                              if opening["lbl_dir"] in ("", proposal) else None)
        for i, (key, default) in enumerate((
                ("data_dir", opening["data_dir"]),
                ("config", opening["config"]),
                ("out_dir", opening["out_dir"]),
                ("lbl_dir", opening["lbl_dir"]),
                ("fits_dir", opening["fits_dir"]))):
            label = ttk.Label(frame, text=self.t(key))
            label.grid(row=i, column=0, sticky="w", pady=2)
            self._register(label, key)
            var = tk.StringVar(value=default)
            self.vars[key] = var
            entry = ttk.Entry(frame, textvariable=var, width=70)
            entry.grid(row=i, column=1, sticky="we", padx=6)
            self._tip(entry, "help_" + key)
            self._tip(label, "help_" + key)
            var.trace_add("write", lambda *_: self._sync())
            if key == "data_dir":
                # proposed when the data root is settled, not at every keystroke
                entry.bind("<FocusOut>", lambda _e: self._propose_out())
                entry.bind("<Return>", lambda _e: self._propose_out())
            browse = ttk.Button(frame, text=self.t("browse"),
                                command=lambda k=key: self._browse(k))
            browse.grid(row=i, column=2)
            self._register(browse, "browse")
        # said once, on opening, when a disk a configuration names is not there
        if self._lost_disk:
            self._say("log_no_disk", self._lost_disk, level="warn")
        # beside the LBL folder, what goes into it: links or copies
        how = ttk.Frame(frame)
        how.grid(row=3, column=3, sticky="w", padx=4)
        for key, path, kind in OPTIONS_PATHS:
            label = ttk.Label(how, text=self.t("opt_" + key))
            label.pack(side="left")
            self._register(label, "opt_" + key)
            var = tk.StringVar(value=str(self.saved.get(
                key, self._config_default(path, kind) or kind[0])))
            box = ttk.Combobox(how, textvariable=var, values=list(kind),
                               width=8, state="readonly")
            box.pack(side="left", padx=(4, 0))
            var.trace_add("write", lambda *_: self._sync())
            self.vars[key] = var
            self._tip(box, "help_" + key)
            self._tip(label, "help_" + key)
        rescan = ttk.Button(frame, text=self.t("rescan"),
                            command=self.refresh_objects)
        rescan.grid(row=0, column=3, padx=4)
        self._register(rescan, "rescan")
        self._tip(rescan, "help_rescan")
        frame.columnconfigure(1, weight=1)
        # the date window's two bounds live with the timeline that sets them
        for key in ("min_rjd", "max_rjd"):
            self.vars[key] = tk.StringVar(value=self.saved.get(key, ""))
            self.vars[key].trace_add("write", lambda *_: self._sync())

    def _warn_other_clone(self, config_path):
        """Say it when the configuration lives in ANOTHER copy of this package.

        A configuration is read from where it is, and the code is the code this
        window is running: those are two different things, and when the config
        sits in a second clone they are two different VERSIONS. That is how a
        window came to write a flag the run refused (2026-09-15), and the
        settings read on this page were yesterday's while the run was today's.
        """
        folder = os.path.dirname(os.path.abspath(str(config_path)))
        theirs = os.path.join(folder, "pca2d", "__init__.py")
        mine = os.path.join(package_home(), "pca2d", "__init__.py")
        if os.path.exists(theirs) and os.path.abspath(theirs) != mine:
            self._say("log_other_clone", folder, package_home(), level="warn")

    def _config_fits_dir(self):
        """What the configuration says the products disk is, if it says one."""
        try:
            config = self._config()
        except Exception:                                       # noqa: BLE001
            return ""
        value = ((config or {}).get("output") or {}).get("fits_directory") or ""
        self._lost_disk = (str(value) if value and not os.path.isdir(str(value))
                           else "")
        return str(value)

    def _build_name(self, parent):
        """What this reduction is called, which is what its folder is called."""
        ttk, tk = self.ttk, self.tk
        row = ttk.Frame(parent)
        row.pack(fill="x", pady=(2, 6))
        label = ttk.Label(row, text=self.t("run_name"))
        label.pack(side="left")
        self._register(label, "run_name")
        self._tip(label, "help_run_name")
        # A name is PROPOSED rather than left empty: a window that opens with
        # an empty field says a reduction needs no name, and then two of them
        # land in one folder. Only when none was kept: a field somebody emptied
        # stays empty.
        kept = self.saved.get("run_name")
        opening = dict(self.saved)
        opening["objects"] = self.saved.get("checked") or []
        proposal = suggested_run_name(opening)
        if not str(kept or "").strip():
            kept = proposal
        # It goes on following the targets for as long as it IS the proposal,
        # which a name restored from the last session still can be: see
        # _follow_name, which stops the moment somebody types their own.
        self._proposed = kept if kept == proposal else None
        self.vars["run_name"] = tk.StringVar(value=kept)
        entry = ttk.Entry(row, textvariable=self.vars["run_name"], width=34)
        entry.pack(side="left", padx=(6, 4))
        self._tip(entry, "help_run_name")
        auto = ttk.Button(row, text=self.t("auto"), width=6,
                          command=self._auto_name)
        auto.pack(side="left", padx=(0, 14))
        self._register(auto, "auto")
        self._tip(auto, "help_auto_button")
        self.vars["run_name"].trace_add("write", lambda *_: self._sync())
        self.exists = ttk.Label(row, style="Hint.TLabel", text="")
        self.exists.pack(side="left", padx=(6, 0))

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
        # kept on the instance or the garbage collector takes them and the
        # column goes blank. One per (ticked, step), made when first needed:
        # two states and thirteen steps, none of them drawn twice
        self._row_images = {}
        self.tree.column("#0", width=210)
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
        # BOTH directions, and expanding: this panel and the timeline share
        # whatever height the page has beside the list of targets, and both
        # redraw themselves on <Configure>, so dragging the window taller gives
        # the bars the room rather than leaving a grey band under them
        box.pack(fill="both", expand=True, padx=6, pady=(0, 6))
        self._register(box, "berv")
        self.berv_canvas = self.tk.Canvas(box, height=86, highlightthickness=0,
                                          background=SURFACE)
        self.berv_canvas.pack(fill="both", expand=True, padx=6, pady=(4, 2))
        self.berv_note = ttk.Label(box, style="Hint.TLabel", text="")
        self.berv_note.pack(anchor="w", padx=8, pady=(0, 4))
        self._tip(self.berv_canvas, "help_berv")
        self._tip(self.berv_note, "help_berv")
        self.berv_canvas.bind("<Configure>", lambda _e: self._draw_berv())

    def _build_snr(self, parent):
        """Signal-to-noise against barycentric velocity, one point per spectrum.

        The pair the histogram cannot show. A campaign can cover the whole
        range and cover one end of it with its worst nights, and the fit weighs
        a spectrum by 1/sigma^2: where the signal-to-noise SITS along the range
        is what decides how well the two frames come apart.
        """
        ttk = self.ttk
        box = ttk.Labelframe(parent, text=self.t("snr_berv"))
        box.pack(fill="both", expand=True, padx=6, pady=(0, 6))
        self._register(box, "snr_berv")
        self.snr_canvas = self.tk.Canvas(box, height=96, highlightthickness=0,
                                         background=SURFACE)
        self.snr_canvas.pack(fill="both", expand=True, padx=6, pady=(4, 6))
        self._tip(self.snr_canvas, "help_snr_berv")
        self.snr_canvas.bind("<Configure>", lambda _e: self._draw_snr())

    def _draw_snr(self):
        """One dot per spectrum, in its star's colour: SNR against BERV."""
        canvas = getattr(self, "snr_canvas", None)
        if canvas is None:
            return
        canvas.delete("all")
        names = self.picked()
        width = max(int(canvas.winfo_width()), 50)
        height = max(int(canvas.winfo_height()), 40)
        if not names:
            canvas.create_text(width // 2, height // 2,
                               text=self.t("snr_none"), fill="#888",
                               font=("Helvetica", 10))
            return
        series = {}
        for name in names:
            berv, snr = scan.snr_against_berv(self.index, name)
            if berv.size:
                series[name] = (berv, snr)
        if not series:
            canvas.create_text(width // 2, height // 2, text=self.t("berv_wait"),
                               fill="#b26a00", font=("Helvetica", 10))
            return

        pad, foot, gutter = 6, 16, 26
        top = height - foot
        highest = max(float(np.nanmax(s)) for _b, s in series.values())
        highest = max(highest, 1.0)
        span = scan.BERV_AXIS

        def x_of(v):
            return gutter + (v + span) / (2 * span) * (width - gutter - pad)

        def y_of(value):
            return top - (value / highest) * (top - pad)

        for level in (0, highest / 2.0, highest):
            y = y_of(level)
            canvas.create_line(gutter, y, width - pad, y,
                               fill="#e8edf3" if level else "#bbb")
            canvas.create_text(gutter - 4, y, text="%d" % round(level),
                               anchor="e", fill="#555", font=("Helvetica", 8))
        for value in (-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0):
            x = x_of(value)
            canvas.create_line(x, top, x, top + 3, fill="#888")
            canvas.create_text(x, top + 9, text="%+.0f" % value, fill="#555",
                               font=("Helvetica", 8))
        for name, (berv, snr) in series.items():
            colour = self._star_colour(name)
            step = max(1, berv.size // 900)   # a canvas is not a plot library
            for b, sn in zip(berv[::step], snr[::step]):
                if not (np.isfinite(b) and np.isfinite(sn)):
                    continue
                x, y = x_of(b), y_of(sn)
                canvas.create_oval(x - 1.6, y - 1.6, x + 1.6, y + 1.6,
                                   fill=colour, outline="")

    def _build_time(self, parent):
        """When the ticked campaigns were observed, and what to keep of them.

        A reduced Julian date says nothing to anybody, so the dates are drawn
        and read in the calendar, and chosen with two sliders rather than typed.
        """
        ttk = self.ttk
        box = ttk.Labelframe(parent, text=self.t("timeline"))
        box.pack(fill="both", expand=True, padx=6, pady=(0, 6))
        self._register(box, "timeline")
        self.time_canvas = self.tk.Canvas(box, height=80, highlightthickness=0,
                                          background=SURFACE)
        self.time_canvas.pack(fill="both", expand=True, padx=6, pady=(4, 2))
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
        # the marks grow with the band they sit in: a panel dragged twice as
        # tall used to hold the same 8-pixel ticks with twice the white space
        # between them. Clamped at the bottom so ten campaigns stay legible,
        # and at 0.35 of the band so two rows never touch.
        half = max(3.0, min(band * 0.35, 90.0))
        for i, (name, values) in enumerate(sorted(times.items())):
            y = top + band * (i + 0.5)
            colour = self._star_colour(name)
            step = max(1, len(values) // 700)     # a canvas is not a plot library
            for t in values[::step]:
                x = x_of(t)
                inside = keep_lo <= t <= keep_hi
                canvas.create_line(x, y - half, x, y + half,
                                   fill=colour if inside else "#d6dbe4")
            canvas.create_text(pad, y - half - 5, text=name, anchor="w",
                               fill=colour, font=("Helvetica", 8))
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
        if item and self.tree.identify_column(event.x) == "#0" and event.x <= 36:
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

    def _draw_check(self, item, name, known=None, total=None):
        """The tick box, how much of that campaign is read, and the name.

        Both marks are in the row's one image slot, which is before the text:
        a Treeview cell holds no image, so a disc drawn at any size can only go
        there. Written as a character it could sit after the name, but then it
        is the size of the name.
        """
        if known is None or total is None:
            known, total = self._read_counts(name)
        image = self._row_image(name in self.checked, known, total)
        if image is not None:
            self.tree.item(item, image=image, text="  %s" % name)
            return
        glyph = CHECKED if name in self.checked else UNCHECKED
        slice_ = pie_glyph(known, total)
        self.tree.item(item, text=("%s  %s %s"
                                   % (glyph, name, slice_)).rstrip())

    def _row_image(self, ticked, known, total):
        """The image for this row, drawn once per (ticked, step)."""
        key = (bool(ticked), pie_step(known, total))
        if key not in self._row_images:
            made = row_image(self.tk, key[0], key[1])
            if made is None:
                return None
            self._row_images[key] = made
        return self._row_images[key]

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
        self._draw_snr()
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
        box.pack(fill="x", pady=4)
        self._register(box, "settings")
        grid = ttk.Frame(box)
        grid.pack(fill="x", padx=8, pady=6)
        # three columns, as evenly as the list divides: four rows each was
        # written when there were eleven settings, and left the ninth alone in
        # a column of its own when two of them were settled
        per_column = -(-len(OPTIONS) // 3)
        for i, (key, path, kind) in enumerate(OPTIONS):
            row, col = i % per_column, (i // per_column) * 3
            label = ttk.Label(grid, text=self.t("opt_" + key))
            label.grid(row=row, column=col, sticky="w", pady=2)
            self._register(label, "opt_" + key)
            default = self.saved.get(key, self._config_default(path, kind))
            if kind == "bool":
                var = tk.BooleanVar(value=bool(default))
                widget = ttk.Checkbutton(grid, variable=var)
            elif isinstance(kind, tuple) and kind and kind[0] == "radio":
                # One variable, two buttons, so choosing one releases the other
                # by construction rather than by a callback that has to
                # remember. The variable holds the config's OWN word, so the
                # command still says --weight velocity and nothing downstream
                # translates anything back.
                var = tk.StringVar(value=str(default))
                widget = ttk.Frame(grid)
                for value, caption in kind[1:]:
                    button = ttk.Radiobutton(widget, text=caption, value=value,
                                             variable=var)
                    button.pack(side="left", padx=(0, 10))
                    self._tip(button, "help_" + key)
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
        for stage in STAGES:
            self.vars["stage_" + stage] = tk.BooleanVar(
                value=self.saved.get("stage_" + stage, True))
        # the chain, with an arrow between each pair, because they happen in
        # that order and depend on each other in that order
        column = 1
        for i, stage in enumerate(CHAIN):
            box_ = ttk.Checkbutton(run, text=stage,
                                   variable=self.vars["stage_" + stage],
                                   command=lambda s=stage: self._stage(s))
            box_.grid(row=0, column=column, padx=(0, 1))
            self._tip(box_, "help_stage_" + stage)
            column += 1
            if i < len(CHAIN) - 1:
                arrow = ttk.Label(run, text="\u2192", style="Hint.TLabel")
                arrow.grid(row=0, column=column, padx=1)
                column += 1
        # figures apart, further along the line and with no arrow to it: it
        # draws what the fit left, and nothing waits for a drawing
        drawn = ttk.Checkbutton(run, text="figures",
                                variable=self.vars["stage_figures"],
                                command=lambda: self._stage("figures"))
        drawn.grid(row=0, column=column, padx=(26, 0))
        self._tip(drawn, "help_stage_figures")
        # No variant picker: the parameters converged, and a second set of
        # settings offered beside the settings is a window that contradicts
        # itself. `pca2d-preclean --object X --variant NAME` still runs one, and
        # Export YAML still writes one, which is how the runs in
        # variants/README.md stay reproducible.

    #: which buttons belong to which tab: what changes a setting is with the
    #: settings, what starts or ends a run is with the run
    BUTTONS = {
        "settings": (("export", "export", None, "help_export_button"),
                     ("savedefaults", "save_defaults", None,
                      "help_savedefaults_button")),
        "run": (("run", "start", "run_button", "help_run_button"),
                ("stop", "stop", "stop_button", "help_stop_button"),
                ("dry", "dry_run", None, "help_dry_button"),
                ("savelog", "save_log", None, "help_savelog_button"),
                ("openout", "open_outputs", None, "help_openout_button"),
                ("openpdf", "open_report", None, "help_openpdf_button")),
    }

    def dry_run(self):
        """Resolve everything and touch nothing."""
        self.start(dry=True)

    def _build_buttons(self, parent, which, quit_too=False):
        """One row of buttons, the ones that belong on this tab."""
        ttk = self.ttk
        bar = ttk.Frame(parent)
        bar.pack(fill="x", pady=(2, 6))
        for key, method, attr, tip in self.BUTTONS[which]:
            button = ttk.Button(bar, text=self.t(key),
                                command=getattr(self, method),
                                style="Run.TButton" if key == "run"
                                else "TButton")
            button.pack(side="left", padx=(0 if key == "run" else 5, 0))
            self._register(button, key)
            self._tip(button, tip)
            if attr:
                setattr(self, attr, button)
        if quit_too:
            # kept for a bar that wants its own: the window's own Quit lives in
            # the banner now, where every page can reach it
            quit_button = ttk.Button(bar, text=self.t("quit"),
                                     command=self.quit_window)
            quit_button.pack(side="right")
            self._register(quit_button, "quit")
            self._tip(quit_button, "help_quit_button")
        return bar

    def _build_command(self, parent):
        """The command this would run, in full, before it runs."""
        ttk, tk = self.ttk, self.tk
        box = ttk.Labelframe(parent, text=self.t("command"))
        box.pack(fill="x", pady=(2, 4))
        self._register(box, "command")
        self.command = tk.Text(box, height=2, wrap="word",
                               font=self.fonts["mono"], background=SURFACE,
                               foreground=INK, relief="flat", padx=8, pady=6,
                               highlightthickness=1,
                               highlightbackground=LINE, highlightcolor=LINE)
        self.command.pack(fill="x", padx=6, pady=6)
        self._tip(self.command, "help_command")

    def _build_log(self, parent):
        ttk, tk = self.ttk, self.tk
        # What the window is doing, over the box that shows what it says rather
        # than beside the buttons: a target and a count that change every ten
        # spectra, in small grey text at the right edge, is where nobody looks.
        self.status = ttk.Label(parent, text=self.t("idle"), style="Hint.TLabel")
        self.status.pack(anchor="w", padx=4, pady=(2, 0))
        box = ttk.Labelframe(parent, text=self.t("output"))
        box.pack(fill="both", expand=True, pady=(2, 4))
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

    def state(self):
        out = {"objects": self.picked(), "lang": self.lang}
        for key, var in self.vars.items():
            out[key] = var.get()
        return out

    def _sync(self, *_args):
        self._follow_name()
        self._follow_lbl_dir()
        state = self.state()
        self.command.delete("1.0", "end")
        self.command.insert("1.0", command_line(state, self.t("command_pending"),
                                                self._defaults()))
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
        if not state.get("objects"):
            label.configure(text="")
            return
        root = (state.get("out_dir")
                or (self._config().get("output") or {}).get("directory")
                or "outputs")
        where = run_folder(state, root)
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
        if self.scanning:
            # one scan at a time: two threads walking the same index would
            # overwrite each other's answers
            self._say("log_busy", level="warn")
            return
        if not str(root).strip():
            # the opening state of a fresh window: nothing to read yet, and an
            # empty table on its own does not say what is missing
            self._fill([])
            self._say("log_pick_root", level="warn")
            self._state()
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
        # a gutter on the left for the counts: a histogram whose height means
        # "how many spectra" and never says how many is a shape, not a
        # measurement
        pad, foot, gutter = 6, 16, 26
        top = height - foot
        tallest = max(1, max(int(c.max()) for c in counts.values() if c.size))
        n_bins = len(edges) - 1
        step = (width - gutter - pad) / max(n_bins, 1)

        def x_of(v):
            return gutter + (v - edges[0]) / max(edges[-1] - edges[0], 1e-9) * (
                width - gutter - pad)

        def y_of(count):
            return top - (count / float(tallest)) * (top - pad)

        # round numbers, and never more of them than there is room for
        ticks = [0, tallest]
        for fraction in (0.5, 0.25, 0.75):
            if (top - pad) > 60:
                ticks.append(int(round(tallest * fraction)))
        for count in sorted(set(t for t in ticks if 0 <= t <= tallest)):
            y = y_of(count)
            canvas.create_line(gutter, y, width - pad, y,
                               fill="#e8edf3" if count else "#bbb")
            canvas.create_text(gutter - 4, y, text="%d" % count, anchor="e",
                               fill="#555", font=("Helvetica", 8))

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
            x0 = gutter + i * step
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
        canvas.create_line(gutter, top, width - pad, top, fill="#bbb")
        # a legend, or the colours say nothing: one square and one name per star,
        # in the order they are stacked
        x = gutter + 2
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

    def _read_counts(self, name):
        """(files read of this object, files its folder holds)."""
        known = len(((self.index.get("objects") or {}).get(name) or {})
                    .get("files") or {})
        return known, int((self.rows.get(name) or {}).get("files") or 0)

    def _berv_partial(self, names):
        """Whether some spectrum of a ticked object has not been read yet."""
        if self.scanning:
            return True
        for name in names:
            known, total = self._read_counts(name)
            if known < total:
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
                self._draw_check(item, name, known, total)
                break
        # The coverage panel and the timeline are drawn from what has been
        # READ, so a campaign that has just been finished changes both of them,
        # and the banner saying they are partial has to go when it stops being
        # true. Only when the campaign is done and it is ticked: every ten
        # spectra would be a redraw of both panels a few hundred times.
        if known >= total and name in self.picked():
            self._draw_berv()
            self._draw_snr()
            self._draw_time()

    def _watch_root(self):
        """Look at the data root every WATCH_MS: campaigns appear while this is
        open.

        Spectra are copied in while the window sits there, and a list that only
        changes when somebody presses Rescan is a list that is quietly wrong.
        One folder listing per campaign, nothing opened, and off the main thread
        because these data usually live on a disk that can be slow to answer.
        """
        self.root.after(WATCH_MS, self._watch_root)
        root = self.vars["data_dir"].get().strip()
        if not root or self.scanning or getattr(self, "_watching", False):
            return
        self._watching = True

        def look():
            try:
                found = dict(objects_in(root))
            except OSError:
                found = None
            self.lines.put(("watched", root, found))

        threading.Thread(target=look, daemon=True).start()

    def _watched(self, root, found):
        """One look at the data root, applied on the main thread."""
        self._watching = False
        if found is None or self.scanning:
            return
        if root != self.vars["data_dir"].get().strip():
            return                      # the root moved while we were looking
        added, gone, grown = folder_news(found, self.rows,
                                         getattr(self, "_seen", None))
        self._seen = found
        if not (added or gone or grown):
            return
        if added:
            self._say("log_watch_added", ", ".join(added), level="value")
        if gone:
            self._say("log_watch_gone", ", ".join(gone), level="warn")
        if grown and not added and not gone:
            self._say("log_watch_grew", ", ".join(grown), level="value")
        self.refresh_objects()

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
        self._draw_snr()
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
        argv = build_command(state, self._defaults())
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
        # this window's own package, ahead of everything: the run starts in the
        # config file's folder, and whatever copy of pca2d sits there would
        # otherwise be the one that runs
        env["PYTHONPATH"] = os.pathsep.join(
            [package_home()] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH")
                                else []))
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
            elif item[0] == "watched":
                self._watched(*item[1:])
            elif item[0] == "scanned":
                self._scanned(*item[1:])
            elif item[0] == "measured":
                self._measured(item[1])
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
            for sequence in ANSI.findall(line):
                for code in ANSI_CODES.findall(sequence):
                    # 92 is 32's bright twin, 91 is 31's, and so on: a bright
                    # colour is the same colour as far as this window cares
                    code = str(int(code) - 60) if code in ("90", "91", "92",
                                                           "93", "94", "95",
                                                           "96", "97") else code
                    if code in LEVELS:
                        tag = LEVELS[code][0]
                        break
                if tag != "plain":
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

    def _out_root(self):
        """The output root in full, as the run will read it.

        A relative root is relative to the CONFIGURATION, not to wherever the
        window was started from, which is the same rule the run follows.
        """
        root = self.vars["out_dir"].get() or (
            (self._config().get("output") or {}).get("directory") or "outputs")
        return os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(
            self.vars["config"].get())), root))

    def _open(self, path):
        """Hand a file or a folder to whatever opens it here."""
        opener = ("open" if sys.platform == "darwin"
                  else "explorer" if os.name == "nt" else "xdg-open")
        subprocess.Popen([opener, path])

    def open_outputs(self):
        """Show the output root in the file browser: a run leaves one PDF and a
        folder of corrected spectra, and they are easier to find by looking."""
        path = self._out_root()
        if not os.path.isdir(path):
            self._say("log_no_outputs", path, level="warn")
            return
        self._open(path)

    def open_report(self):
        """Open the one PDF this run binds everything into.

        The output root holds every run ever made, and finding this one's
        report in it means knowing that a run is named by its counts and by
        the name typed beside them. The window knows all of that already, so
        it opens the file rather than the folder above it.
        """
        path = report_pdf(self.state(), self._out_root())
        if not path or not os.path.exists(path):
            self._say("log_no_report", path or self._out_root(), level="warn")
            return
        self._open(path)

    def _build_lbl(self, parent):
        """The `lbl:` block, on its own page: the wrapper's own questions.

        A window of its own was one window too many once the settings had
        pages: this is a setting like the others, it is simply LBL's rather
        than the fit's, and it is the step that produces velocities.
        """
        ttk, tk = self.ttk, self.tk
        note = ttk.Label(parent, style="Hint.TLabel", wraplength=900,
                         justify="left", text=self.t("help_lblwin_button"))
        note.pack(anchor="w", pady=(0, 8))
        self._register(note, "help_lblwin_button")
        box = ttk.Labelframe(parent, text=self.t("lbl_title"))
        box.pack(fill="x", pady=4)
        self._register(box, "lbl_title")
        grid = ttk.Frame(box)
        grid.pack(fill="x", padx=8, pady=6)
        per_column = -(-len(OPTIONS_LBL) // 2)
        for i, (key, path, kind) in enumerate(OPTIONS_LBL):
            row, col = i % per_column, (i // per_column) * 2
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
                                      width=14, state="readonly")
            else:
                var = tk.StringVar(value="" if default is None else str(default))
                widget = ttk.Entry(grid, textvariable=var, width=30)
            widget.grid(row=row, column=col + 1, sticky="w", padx=(8, 28))
            var.trace_add("write", lambda *_: self._sync())
            self.vars[key] = var
            self._tip(widget, "help_" + key)
            self._tip(label, "help_" + key)

    # ---- the disks ----------------------------------------------------
    def _build_clean(self, parent):
        """What pca2d has left on the disks, and the one button that removes it.

        The disks fill with cubes: a campaign's cache is a few GB and there is
        one per configuration, so the cache is usually the biggest thing here
        and none of it is a result. It was invisible until this page, and the
        way to find it was to know that `cache/` existed.
        """
        ttk = self.ttk
        note = ttk.Label(parent, style="Hint.TLabel", wraplength=900,
                         justify="left", text=self.t("help_clean"))
        note.pack(anchor="w", pady=(0, 8))
        self._register(note, "help_clean")
        box = ttk.Labelframe(parent, text=self.t("clean_title"))
        box.pack(fill="both", expand=True, pady=4)
        self._register(box, "clean_title")

        self.clean_tree = ttk.Treeview(
            box, columns=("size", "files", "kind"), show="tree headings",
            selectmode="extended", height=12)
        self.clean_headings = (("#0", "col_object"), ("size", "col_size"),
                               ("files", "col_nfiles"), ("kind", "col_kind"))
        self.clean_tree.column("#0", width=240)
        self.clean_tree.column("size", width=90, anchor="e")
        self.clean_tree.column("files", width=80, anchor="e")
        self.clean_tree.column("kind", width=110)
        self._draw_clean_headings()
        self.clean_tree.pack(fill="both", expand=True, padx=6, pady=(6, 2))
        # one line per item is not room for what deleting it would cost, so the
        # selected row says it underneath, in full
        self.clean_tree.bind("<<TreeviewSelect>>", self._clean_selected)

        self.clean_what = ttk.Label(box, style="Hint.TLabel", wraplength=880,
                                    justify="left", text="")
        self.clean_what.pack(anchor="w", padx=8, pady=(0, 4))
        self.clean_totals = ttk.Label(box, text=self.t("clean_none"))
        self.clean_totals.pack(anchor="w", padx=8, pady=(0, 4))
        # not _register'ed: it says a MEASUREMENT, and a registered label is
        # overwritten with its own translation the next time the language
        # changes, which would put "nothing measured yet" over a survey

        bar = ttk.Frame(box)
        bar.pack(fill="x", padx=6, pady=(0, 8))
        self.clean_button = ttk.Button(bar, text=self.t("clean_measure"),
                                       command=self.measure_disks)
        self.clean_button.pack(side="left", padx=(0, 6))
        self._register(self.clean_button, "clean_measure")
        self._tip(self.clean_button, "help_clean_measure")
        # NOT disabled when there is nothing to free: a greyed button says
        # "you cannot", and what was meant was "there is nothing here yet".
        # It says that itself, in the log, when it is pressed
        self.purge_button = ttk.Button(bar, text=self.t("clean_purge"),
                                       command=self.purge_disks)
        self.purge_button.pack(side="left")
        self._register(self.purge_button, "clean_purge")
        self._tip(self.purge_button, "help_clean_purge")
        # and the other half: whatever is picked in the list goes, of whatever
        # kind. An afternoon of LBL on a disk with no room left is still a
        # choice somebody is entitled to make, so it is offered and priced
        self.select_button = ttk.Button(bar, text=self.t("clean_selected"),
                                        command=self.purge_selected)
        self.select_button.pack(side="left", padx=(6, 0))
        self._register(self.select_button, "clean_selected")
        self._tip(self.select_button, "help_clean_selected")
        self.clean_items = []
        self._measuring = False

    def _draw_clean_headings(self):
        for column, key in self.clean_headings:
            self.clean_tree.heading(column, text=self.t(key))

    def _clean_state(self):
        """(config path, output root) as the window has them right now.

        The fields, not the file: the window's own out directory is where this
        run's products would go, and it is the one the reader is looking at.
        """
        return (self.vars["config"].get().strip() or None,
                self.vars["out_dir"].get().strip() or None,
                self.vars["lbl_dir"].get().strip() or None)

    def measure_disks(self):
        """Walk the folders off the main thread and fill the list."""
        if self._measuring:
            return
        self._measuring = True
        self.clean_button.configure(state="disabled")
        self.purge_button.configure(state="disabled")
        self.select_button.configure(state="disabled")
        self.clean_totals.configure(text=self.t("clean_measuring"))
        config_path, out_root, lbl_dir = self._clean_state()

        def look():
            from . import housekeeping
            from .config import load_config
            try:
                config = load_config(config_path) if config_path else {}
            except (Exception, SystemExit):                   # noqa: BLE001
                # a config that will not load is not a reason to refuse to
                # count: the cache and the LBL tree are then the defaults
                config = {}
            try:
                items = housekeeping.survey(
                    config, config_path, out_root,
                    package=os.path.dirname(os.path.abspath(__file__)),
                    lbl_dir=lbl_dir)
            except OSError:
                items = []
            self.lines.put(("measured", items))

        threading.Thread(target=look, daemon=True).start()

    def _measured(self, items):
        """One survey, applied on the main thread."""
        from .housekeeping import totals
        self._measuring = False
        self.clean_items = items
        self.clean_button.configure(state="normal")
        self.purge_button.configure(state="normal")
        self.select_button.configure(state="normal")
        self._clean_rows()

    def _clean_rows(self):
        """The survey in the list, in the window's language.

        Apart from _measured because the language button redraws it and a
        survey may be going at that moment: switching language must not tell
        the window that the measurement it is waiting for has arrived.
        """
        from .housekeeping import human
        for row in self.clean_tree.get_children():
            self.clean_tree.delete(row)
        for i, item in enumerate(self.clean_items):
            self.clean_tree.insert(
                "", "end", iid=str(i), text=item["name"],
                values=(human(item["bytes"]), item["files"],
                        self.t("kind_" + item["kind"])))
        self._clean_totals()

    def _clean_totals(self, items=None):
        """The two numbers under the list, in the window's language."""
        from .housekeeping import human, totals
        items = self.clean_items if items is None else items
        if not items:
            self.clean_totals.configure(text=self.t("clean_none"))
            return
        total, free = totals(items)
        self.clean_totals.configure(
            text=self.t("clean_totals") % (human(total), human(free)))

    def _clean_selected(self, _event=None):
        rows = self.clean_tree.selection()
        if not rows:
            return
        item = self.clean_items[int(rows[0])]
        path = item["path"]
        where = ", ".join(path) if isinstance(path, list) else (path or "")
        key = "clean_" + str(item.get("key") or "")
        what = self.t(key) if key in EN else item["what"]
        self.clean_what.configure(
            text="%s\n%s" % (where, what) if what else where)

    def purge_disks(self):
        """Delete the scratch and the rebuildable, after saying what and how
        much. The one press that needs no choosing: everything it touches comes
        back by itself, at a cost in minutes."""
        if not self.clean_items:
            self._say("clean_none", level="warn")
            return
        going = [it for it in self.clean_items if it.get("removable")
                 and it["bytes"]]
        if not going:
            self._say("clean_nothing", level="warn")
            return
        from .housekeeping import human
        freed = sum(it["bytes"] for it in going)
        if not self._ask_delete(self.t("clean_confirm_title"),
                                self.t("clean_confirm")
                                % (human(freed), len(going),
                                   self._listed(going))):
            return
        self._delete(going)

    def purge_selected(self):
        """Delete whatever is picked in the list, of whatever kind.

        The expensive and the results are on the list too, and a row that can
        be read but never chosen is a row that lies about what the disk holds.
        What separates them from the cache is not whether they may be deleted,
        it is what it costs to have them back, so that is what the question
        says, in the terms of the most expensive kind in the selection.
        """
        from .housekeeping import EXPENSIVE, RESULTS, human
        rows = self.clean_tree.selection()
        going = [self.clean_items[int(r)] for r in rows if r.isdigit()]
        going = [it for it in going if it["bytes"]]
        if not rows or not going:
            self._say("clean_pick" if not rows else "clean_nothing",
                      level="warn")
            return
        kinds = {it["kind"] for it in going}
        cost = ("clean_cost_results" if RESULTS in kinds else
                "clean_cost_expensive" if EXPENSIVE in kinds else
                "clean_cost_rebuildable")
        if not self._ask_delete(
                self.t("clean_confirm_pick_title"),
                self.t("clean_confirm_pick")
                % (human(sum(it["bytes"] for it in going)), len(going),
                   self._listed(going), self.t(cost))):
            return
        self._delete(going)

    @staticmethod
    def _listed(items):
        """What is about to go, one line each, biggest first."""
        from .housekeeping import human
        return "\n".join("  %s   %s" % (human(it["bytes"]), it["name"])
                          for it in sorted(items, key=lambda it: -it["bytes"]))

    def _ask_delete(self, title, question):
        """The question, with a run going on this machine added to it.

        A run reads from the folders on this page. It survives them going,
        since a stage that finds a cube missing builds it again (cli.
        rebuild_missing_cubes), but building it again is twenty minutes that
        run is in the middle of, and that is worth knowing before rather than
        after.
        """
        from tkinter import messagebox
        if getattr(self, "proc", None) is not None:
            question = "%s\n\n%s" % (self.t("clean_while_running"), question)
        return bool(messagebox.askyesno(title, question))

    def _delete(self, going):
        """Remove exactly these items, and measure again so the list is true."""
        from .housekeeping import human, purge
        gained, gone = purge(going, kinds={it["kind"] for it in going})
        self._say("clean_freed", human(gained), len(gone), level="info")
        self.measure_disks()

    def _stage(self, changed):
        """A stage box moved: the ones that depend on it follow."""
        ticked = {name: bool(self.vars["stage_" + name].get())
                  for name in STAGES}
        for name, value in follow_stages(ticked, changed).items():
            if bool(self.vars["stage_" + name].get()) != value:
                self.vars["stage_" + name].set(value)
        self._sync()

    def _auto_name(self):
        """Propose a name again, from the targets and settings as they stand.

        And put the name back under the window's care: Auto is how one asks
        for the proposal after having typed something else.
        """
        self._proposed = suggested_run_name(self.state())
        self.vars["run_name"].set(self._proposed)

    def _defaults(self):
        """The configuration the command is compared with, or None."""
        try:
            return self._config() or None
        except Exception:                                       # noqa: BLE001
            return None

    def _follow_lbl_dir(self):
        """Keep the proposed LBL folder under the output root as it moves.

        Only while the field holds the window's own proposal, or nothing: an
        empty field means the default anyway, which IS the proposal. A folder
        somebody typed or browsed to, an existing LBL tree above all, is a
        decision and is left alone.
        """
        var = self.vars.get("lbl_dir")
        out = self.vars.get("out_dir")
        if var is None or out is None or getattr(self, "_naming_lbl", False):
            return
        current = str(var.get() or "").strip()
        if current and current != getattr(self, "_proposed_lbl", None):
            return
        fresh = lbl_proposal(out.get())
        if fresh == current:
            return
        self._naming_lbl = True
        try:
            self._proposed_lbl = fresh
            var.set(fresh)
        finally:
            self._naming_lbl = False

    def _follow_name(self):
        """Keep the proposed name in step with what it names.

        The name IS the folder, and the folder is what says, months later,
        which targets were reduced and whether the settings were these ones. A
        name still reading GJ1+PROXIMA over a run of TOI-2120 is worse than no
        name at all, so it follows the ticks and the settings that go into it.

        It follows only while it is still the window's own proposal. A name
        somebody typed is a decision and is left alone, and so is a field
        somebody emptied; Auto is how the proposal is asked for back.
        """
        var = self.vars.get("run_name")
        if var is None or getattr(self, "_naming", False):
            return
        current = str(var.get() or "").strip()
        if not current or current != getattr(self, "_proposed", None):
            return
        fresh = suggested_run_name(self.state())
        if fresh == current:
            return
        # the write fires this trace again; the flag is what stops it being a
        # loop, and the second pass would in any case find nothing to change
        self._naming = True
        try:
            self._proposed = fresh
            var.set(fresh)
        finally:
            self._naming = False

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

    def quit_window(self, confirm=None):
        """Leave, and ask first when a run would be stopped by leaving.

        The settings are written at every change, so nothing of the window is
        lost by closing it. A run is another matter: it is a subprocess of this
        one, and it goes when this does.
        """
        if self.proc is not None:
            if not (confirm or self._ask_quit)():
                return
        self._close()

    #: what the quit question wears. A run is an hour of somebody's afternoon,
    #: and the stage it is in goes with the window; the plain grey line of text
    #: tkinter offers for that is not the size of what it is asking.
    STARTLED = "\U0001F633"

    def _ask_quit(self):
        """The quit question, in its own window, with a face on it. True to go.

        messagebox.askyesno takes no image, so this is a Toplevel: the same
        words, the face beside them, and the two answers named rather than
        called Yes and No, since neither of those is the question.
        """
        tk, ttk = self.tk, self.ttk
        win = tk.Toplevel(self.root)
        win.title(self.t("quit_title"))
        win.transient(self.root)
        win.resizable(False, False)
        answer = {"leave": False}
        frame = ttk.Frame(win, padding=18)
        frame.pack(fill="both", expand=True)
        ttk.Label(frame, text=self.STARTLED, font=(EMOJI_FONT, 46)).grid(
            row=0, column=0, rowspan=2, padx=(0, 16), sticky="n")
        ttk.Label(frame, text=self.t("quit_running"), wraplength=380,
                  justify="left").grid(row=0, column=1, sticky="w")
        bar = ttk.Frame(frame)
        bar.grid(row=1, column=1, sticky="e", pady=(14, 0))

        def leave():
            answer["leave"] = True
            win.destroy()

        ttk.Button(bar, text=self.t("quit_no"),
                   command=win.destroy).pack(side="right")
        ttk.Button(bar, text=self.t("quit_yes"),
                   command=leave).pack(side="right", padx=(0, 8))
        win.bind("<Escape>", lambda _event: win.destroy())
        win.protocol("WM_DELETE_WINDOW", win.destroy)
        win.update_idletasks()
        # over the window it interrupts, rather than wherever the system puts it
        x = self.root.winfo_rootx() + (self.root.winfo_width()
                                       - win.winfo_width()) // 2
        y = self.root.winfo_rooty() + (self.root.winfo_height()
                                       - win.winfo_height()) // 3
        win.geometry("+%d+%d" % (max(x, 0), max(y, 0)))
        win.grab_set()
        self.root.wait_window(win)
        return answer["leave"]

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
