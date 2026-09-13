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
import subprocess
import sys
import threading
import time

from . import scan
from .logger import stamp

HOME_STATE = os.path.expanduser("~/.pca2d_gui.json")
ANSI = re.compile(r"\033\[(\d+)m")
#: the logger's colours, and how this window paints them
LEVELS = {"32": ("info", "#1a7f37"), "34": ("value", "#0a58ca"),
          "33": ("warn", "#b26a00"), "31": ("error", "#c1121f")}
STAGES = ("cube", "fit", "figures", "correct", "lbl")
#: what the window can change, and the config key each one is
OPTIONS = [
    ("n_star", "twoframe.n_star", "int"),
    ("n_earth", "twoframe.n_earth", "int"),
    ("mean", "twoframe.mean", ("star", "offset", "full", "iterate")),
    ("star_basis", "twoframe.star_basis", ("spline", "grid")),
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
    "log_scan_start": "reading the data root %s",
    "log_scan_done": "%s: %d objects, %d spectra read, %d already known, %d gone",
    "log_scan_none": "no folder with spectra under %s",
    "log_busy": "still reading the data root: wait for it to finish",
    "log_no_root": "not a folder: %s",
    "log_no_object": "no object ticked: tick at least one in the list",
    "log_mixed": "two instruments ticked (%s): a run is one instrument",
    "log_mixed_refused": "refusing to run: %s are from two instruments (%s)",
    "log_command": "running: %s",
    "log_ended": "the run ended, exit code %d",
    "log_stopping": "asking the run to stop",
    "log_export": "variant written: %s",
    "log_defaults": "defaults written into %s: %s",
    "log_nothing": "nothing to write: every setting is the configuration's own",
    "log_saved_log": "log written: %s",
    "log_no_outputs": "nothing written there yet: %s",
    "log_failed": "could not start it: %s",
    "opt_n_star": "star components", "opt_n_earth": "observer components",
    "opt_mean": "static part", "opt_star_basis": "star basis",
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
        "The target's brightness as ITS OWN pipeline recorded it, with the band"
        " it is in: NIRPS writes the J magnitude, SPIRou writes H, and the two"
        " differ by about a magnitude on an M dwarf, so the band is shown rather"
        " than assumed. Read from the headers, never from a catalogue.",
    "help_check":
        "Tick a target to run it. Several ticked are fitted TOGETHER against one"
        " observer basis. Click the box, double-click the row, or press the space"
        " bar. They must all come from the same instrument.",
    "help_all_button": "Tick every target in the list, or untick every one.",
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
        "Where a run writes, if not the config's own output root. A run puts its"
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
        "The static part of the model. `star`, the nominal, takes one star"
        " spectrum per order parity out before the fit and divides out the"
        " observer block alone. `offset` and `full` put an observer-frame mean"
        " back into the correction, which injected about 47 m/s on Proxima."
        " `iterate` re-estimates means in both frames every sweep and did not"
        " converge on a whole campaign.",
    "help_star_basis":
        "How the star side is carried and updated. `spline`: one cubic B-spline"
        " with a knot per sample, evaluated at each exposure's shifted position"
        " and updated exactly; the fit is about a third faster. `grid`: samples"
        " carried by the Lanczos kernel, the older path. Their velocities agree"
        " within the run-to-run scatter.",
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
    "log_scan_start": "lecture du dossier de données %s",
    "log_scan_done": "%s : %d objets, %d spectres lus, %d déjà connus, %d disparus",
    "log_scan_none": "aucun dossier contenant des spectres sous %s",
    "log_busy": "lecture du dossier de données en cours : attendre la fin",
    "log_no_root": "ce n'est pas un dossier : %s",
    "log_no_object": "aucun objet coché : en cocher au moins un dans la liste",
    "log_mixed": "deux instruments cochés (%s) : un passage, c'est un instrument",
    "log_mixed_refused": "passage refusé : %s viennent de deux instruments (%s)",
    "log_command": "lancement : %s",
    "log_ended": "passage terminé, code de sortie %d",
    "log_stopping": "demande d'arrêt du passage",
    "log_export": "variante écrite : %s",
    "log_defaults": "défauts écrits dans %s : %s",
    "log_nothing": "rien à écrire : tous les réglages sont ceux de la configuration",
    "log_saved_log": "journal écrit : %s",
    "log_no_outputs": "rien n'y est encore écrit : %s",
    "log_failed": "impossible de le lancer : %s",
    "opt_n_star": "composantes stellaires",
    "opt_n_earth": "composantes observateur",
    "opt_mean": "partie statique", "opt_star_basis": "base stellaire",
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
        "Où le passage écrit, si ce n'est pas la racine de sortie de la"
        " configuration. Un passage y dépose sa configuration résolue, son"
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
        "La partie statique du modèle. `star`, le nominal, retire un spectre"
        " stellaire par parité d'ordre avant l'ajustement et ne divise que le"
        " bloc observateur. `offset` et `full` remettent une moyenne en"
        " référentiel observateur dans la correction, ce qui injectait environ"
        " 47 m/s sur Proxima. `iterate` réestime des moyennes dans les deux"
        " référentiels à chaque itération, et n'a pas convergé sur une campagne"
        " entière.",
    "help_star_basis":
        "Comment le côté stellaire est décalé et mis à jour. `spline` : une"
        " spline cubique, un nœud par échantillon, évaluée à la position"
        " décalée de chaque pose et mise à jour exactement ; l'ajustement est"
        " environ un tiers plus rapide. `grid` : des échantillons décalés par le"
        " noyau de Lanczos, l'ancienne voie. Leurs vitesses s'accordent à la"
        " dispersion entre passages près.",
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
                 background="#ffffe0", relief="solid", borderwidth=1,
                 wraplength=460, font=("Helvetica", 11), padx=8, pady=6).pack()

    def leave(self, _event=None):
        if self.after is not None:
            self.widget.after_cancel(self.after)
            self.after = None
        if self.window is not None:
            self.window.destroy()
            self.window = None


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
        self.lbl_window = None
        root.title("pca2d-preclean")
        root.geometry("1200x780")
        style = ttk.Style()
        for theme in ("aqua", "clam", "default"):
            if theme in style.theme_names():
                style.theme_use(theme)
                break
        style.configure("Head.TLabel", font=("Helvetica", 13, "bold"))
        style.configure("Hint.TLabel", foreground="#666")

        self.vars = {}
        self.status = ttk.Label(root, text=self.t("idle"), style="Hint.TLabel")
        self._build_top(root)
        panes = ttk.Panedwindow(root, orient="horizontal")
        panes.pack(fill="both", expand=False, padx=10)
        left, right = ttk.Frame(panes), ttk.Frame(panes)
        panes.add(left, weight=1)
        panes.add(right, weight=2)
        self._build_objects(left)
        self._build_options(right)
        self._build_command(root)
        self._build_log(root)
        self.refresh_objects()
        self.root.after(80, self._drain)
        self.root.protocol("WM_DELETE_WINDOW", self._close)

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
                    self.tree.heading(key[0], text=self.t(key[1]))
            except Exception:                                 # noqa: BLE001
                pass
        self._state()
        if self.lbl_window is not None:
            self.lbl_window.title(self.t("lbl_title"))
        self._sync()

    # ---- widgets ------------------------------------------------------
    def _build_top(self, parent):
        ttk, tk = self.ttk, self.tk
        frame = ttk.Frame(parent)
        frame.pack(fill="x", padx=10, pady=(10, 6))
        ttk.Label(frame, text="pca2d-preclean", style="Head.TLabel").grid(
            row=0, column=0, sticky="w")
        self._register(ttk.Label(frame, style="Hint.TLabel",
                                 text=self.t("subtitle")), "subtitle").grid(
            row=0, column=1, columnspan=2, sticky="w", padx=8)
        button = ttk.Button(frame, text=self.t("lang"),
                            command=self.switch_language, width=10)
        button.grid(row=0, column=3, sticky="e")
        self._register(button, "lang")
        self._tip(button, "help_lang")
        for i, (key, default) in enumerate((
                ("data_dir", self.saved.get("data_dir", "data")),
                ("config", self.saved.get("config", "config.yaml")),
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
            browse = ttk.Button(frame, text=self.t("browse"),
                                command=lambda k=key: self._browse(k))
            browse.grid(row=i + 1, column=2)
            self._register(browse, "browse")
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
                                 height=13)
        for column, key in (("#0", "col_object"), ("files", "col_files"),
                            ("snr", "col_snr"), ("exptime", "col_exptime"),
                            ("mag", "col_mag"), ("instrument", "col_instrument")):
            self.tree.heading(column, text=self.t(key))
            self._register(self.tree, (column, key), how="heading")
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
        self.count = ttk.Label(bar, style="Hint.TLabel", text="")
        self.count.pack(side="right")

    #: what each column of the list is, for the explanation that follows the
    #: pointer: the box, the counts, then the two numbers read from the headers
    COLUMN_HELP = {"#0": "help_check", "#1": "help_objects", "#2": "help_col_snr",
                   "#3": "help_col_exptime", "#4": "help_col_mag",
                   "#5": "help_objects"}

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
        self._sync()

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
            row, col = i % 6, (i // 6) * 3
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
        box.pack(fill="x", padx=10, pady=6)
        self._register(box, "command")
        self.command = tk.Text(box, height=2, wrap="word", font=("Menlo", 11),
                               background="#f6f6f6", relief="flat")
        self.command.pack(fill="x", padx=6, pady=6)
        self._tip(self.command, "help_command")
        bar = ttk.Frame(parent)
        bar.pack(fill="x", padx=10)
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
            button = ttk.Button(bar, text=self.t(key), command=command)
            button.pack(side="left", padx=(0 if key == "run" else 4, 0))
            self._register(button, key)
            self._tip(button, tip)
            if attr:
                setattr(self, attr, button)
        self.stop_button.configure(state="disabled")
        self.status.pack(side="right")

    def _build_log(self, parent):
        ttk, tk = self.ttk, self.tk
        box = ttk.Labelframe(parent, text=self.t("output"))
        box.pack(fill="both", expand=True, padx=10, pady=8)
        self._register(box, "output")
        self.log = tk.Text(box, wrap="none", font=("Menlo", 11),
                           background="white", relief="flat")
        bar = ttk.Scrollbar(box, command=self.log.yview)
        self.log.configure(yscrollcommand=bar.set)
        bar.pack(side="right", fill="y")
        self.log.pack(fill="both", expand=True, padx=6, pady=6)
        for _code, (name, colour) in LEVELS.items():
            self.log.tag_configure(name, foreground=colour)
        self.log.tag_configure("plain", foreground="#222")
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
        keep = {k: v for k, v in state.items() if k != "objects"}
        keep["checked"] = sorted(self.checked)
        _write_state(keep)

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
        """Draw one row per object, keeping whatever was ticked."""
        self.tree.delete(*self.tree.get_children())
        self.names, self.rows = {}, {}
        for row in rows:
            name = row["object"]
            item = self.tree.insert("", "end", values=self._values(row))
            self.names[item] = name
            self.rows[name] = row
            self._draw_check(item, name)
        self.checked &= set(self.rows)     # an object that is gone is not run
        self.count.configure(text="%d / %d" % (len(self.picked()),
                                               len(self.names)))
        self._sync()

    def _values(self, row, approximate=False):
        """One row's columns.

        A number not read yet is blank, never a zero. One read from the first
        few spectra of a campaign carries a tilde: ten of them already give the
        signal-to-noise and the exposure time to the precision anybody picks a
        target with, and the value sharpens as the rest are read.
        """
        tilde = "~" if approximate else ""
        mag = ("" if row.get("mag") is None else
               "%s%s=%.1f" % (tilde, row.get("mag_band") or "?", row["mag"]))
        return (row["files"],
                "" if row.get("snr") is None else "%s%.0f" % (tilde, row["snr"]),
                "" if row.get("exptime") is None
                else "%s%.0f" % (tilde, row["exptime"]),
                mag,
                row.get("instrument") or "?")

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
                return

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
        self._say("log_command", " ".join(argv), level="value")
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
        self.root.after(80, self._drain)

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

    def _browse(self, key):
        from tkinter import filedialog
        if key == "config":
            path = filedialog.askopenfilename(title=self.t("config"),
                                              filetypes=[("YAML", "*.yaml")])
        else:
            path = filedialog.askdirectory(title=self.t(key))
        if path:
            self.vars[key].set(path)
            if key == "data_dir":
                self.refresh_objects()

    def _close(self):
        if self.proc is not None:
            self.proc.terminate()
        self.root.destroy()


def main(argv=None):
    import tkinter as tk

    root = tk.Tk()
    App(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    sys.exit(main())
