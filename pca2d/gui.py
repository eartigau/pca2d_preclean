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

EN = {
    "subtitle": "two-frame precleaning, then LBL. Pick a data root, pick"
                " objects, look at the command, run it.",
    "data_dir": "data root", "config": "config",
    "out_dir": "output root (optional)",
    "browse": "Browse", "rescan": "Rescan",
    "objects": "objects", "settings": "settings", "stages": "stages",
    "variant": "variant", "command": "the command this runs", "output": "output",
    "col_object": "object", "col_files": "files", "col_instrument": "instrument",
    "run": "Run", "stop": "Stop", "dry": "Dry run", "export": "Export YAML...",
    "savelog": "Save log...", "openout": "Open outputs",
    "idle": "idle", "running": "running", "lang": "Français",
    "no_object": "no object selected\n",
    "mixed": "these objects are not from one instrument (%s): a joint fit needs"
             " one domain and one grid, and the run will refuse them\n",
    "nothing_export": "every setting is the configuration's own, so a variant"
                      " file would say nothing.",
    "export_title": "save these settings as a variant",
    "opt_n_star": "star components", "opt_n_earth": "observer components",
    "opt_mean": "static part", "opt_star_basis": "star basis",
    "opt_velocity_term": "fit a velocity per exposure",
    "opt_iters": "sweeps at most",
    "opt_shrink": "divide only what is significant",
    "opt_mask": "samples a corrected file blanks",
    "opt_width_kms": "high pass (km/s)", "opt_dv": "grid step (km/s)",
    "opt_nightly_stack": "coadd each night", "opt_run": "run LBL (hours)",
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
    "help_rescan": "Read the data root again: use it after copying new spectra in.",
    "help_objects":
        "One object: a solo run. SEVERAL: they are fitted together against ONE"
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
    "run": "Lancer", "stop": "Arrêter", "dry": "Essai à blanc",
    "export": "Exporter le YAML...", "savelog": "Enregistrer le journal...",
    "openout": "Ouvrir les sorties",
    "idle": "au repos", "running": "en cours", "lang": "English",
    "no_object": "aucun objet sélectionné\n",
    "mixed": "ces objets ne viennent pas du même instrument (%s) : un ajustement"
             " conjoint exige un seul domaine et une seule grille, et le passage"
             " les refusera\n",
    "nothing_export": "tous les réglages sont ceux de la configuration : un"
                      " fichier de variante ne dirait rien.",
    "export_title": "enregistrer ces réglages comme variante",
    "opt_n_star": "composantes stellaires",
    "opt_n_earth": "composantes observateur",
    "opt_mean": "partie statique", "opt_star_basis": "base stellaire",
    "opt_velocity_term": "ajuster une vitesse par pose",
    "opt_iters": "itérations au plus",
    "opt_shrink": "ne diviser que le significatif",
    "opt_mask": "échantillons blanchis",
    "opt_width_kms": "passe-haut (km/s)", "opt_dv": "pas de grille (km/s)",
    "opt_nightly_stack": "empiler chaque nuit", "opt_run": "lancer LBL (heures)",
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
        "Relire le dossier de données : utile après y avoir copié des spectres.",
    "help_objects":
        "Un objet : passage solo. PLUSIEURS : ils sont ajustés ensemble contre"
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
    """[(name, number of files)] for every object folder under the data root."""
    if not root or not os.path.isdir(root):
        return []
    out = []
    for name in sorted(os.listdir(root)):
        path = os.path.join(root, name)
        if not os.path.isdir(path):
            continue
        n = len(glob.glob(os.path.join(path, pattern)))
        if n:
            out.append((name, n))
    return out


def instrument_of(root, name, pattern="*t.fits"):
    """INSTRUME of the first file that opens, or '?'."""
    from astropy.io import fits
    for path in sorted(glob.glob(os.path.join(root, name, pattern)))[:3]:
        try:
            with fits.open(path) as hdulist:
                value = str(hdulist[0].header.get("INSTRUME", "")).strip()
        except Exception:                                     # noqa: BLE001
            continue
        if value:
            return value
    return "?"


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
    for key, path, kind in OPTIONS:
        if key not in state or state[key] in (None, ""):
            continue
        value = state[key]
        if kind == "int":
            value = int(value)
        elif kind == "float":
            value = float(value)
        elif kind == "bool":
            value = bool(value)
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

    The same text also goes to the status line, so it can be read without
    waiting and without covering anything.
    """

    def __init__(self, app, widget, key):
        self.app, self.widget, self.key = app, widget, key
        self.window = None
        self.after = None
        widget.bind("<Enter>", self.enter, add="+")
        widget.bind("<Leave>", self.leave, add="+")
        widget.bind("<ButtonPress>", self.leave, add="+")

    def enter(self, _event=None):
        self.app.status.configure(text=self.app.t(self.key)[:110])
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
        tk.Label(self.window, text=self.app.t(self.key), justify="left",
                 background="#ffffe0", relief="solid", borderwidth=1,
                 wraplength=460, font=("Helvetica", 11), padx=8, pady=6).pack()

    def leave(self, _event=None):
        if self.after is not None:
            self.widget.after_cancel(self.after)
            self.after = None
        if self.window is not None:
            self.window.destroy()
            self.window = None
        self.app.status.configure(
            text=self.app.t("running" if self.app.proc else "idle"))


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
        self.status.configure(text=self.t("running" if self.proc else "idle"))
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
        self.tree = ttk.Treeview(box, columns=("files", "instrument"),
                                 show="tree headings", selectmode="extended",
                                 height=12)
        for column, key in (("#0", "col_object"), ("files", "col_files"),
                            ("instrument", "col_instrument")):
            self.tree.heading(column, text=self.t(key))
            self._register(self.tree, (column, key), how="heading")
        self.tree.column("#0", width=170)
        self.tree.column("files", width=60, anchor="e")
        self.tree.column("instrument", width=90)
        self.tree.pack(fill="both", expand=True, padx=6, pady=6)
        self.tree.bind("<<TreeviewSelect>>", lambda _e: self._sync())
        self._tip(self.tree, "help_objects")

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
                ("export", self.export, None, "help_export_button"),
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
        out = {"objects": [self.tree.item(i, "text")
                           for i in self.tree.selection()], "lang": self.lang}
        for key, var in self.vars.items():
            out[key] = var.get()
        return out

    def _sync(self, *_args):
        state = self.state()
        self.command.delete("1.0", "end")
        self.command.insert("1.0", " ".join(build_command(state)))
        _write_state({k: v for k, v in state.items() if k != "objects"})

    def refresh_objects(self):
        root = self.vars["data_dir"].get()
        self.tree.delete(*self.tree.get_children())
        found = objects_in(root)
        for name, n in found:
            self.tree.insert("", "end", text=name, values=(n, "..."))
        self._cfg = None
        self.variant_box.configure(values=self._variants())
        self._sync()
        if found:
            threading.Thread(target=self._instruments, args=(root, found),
                             daemon=True).start()

    def _instruments(self, root, found):
        for item, (name, _n) in zip(self.tree.get_children(), found):
            value = instrument_of(root, name)
            self.root.after(0, lambda i=item, v=value:
                            self.tree.set(i, "instrument", v))

    # ---- running ------------------------------------------------------
    def start(self, dry=False):
        if self.proc is not None:
            return
        state = self.state()
        if not state["objects"]:
            self._write(self.t("no_object"), "error")
            return
        state["dry_run"] = dry
        instruments = {self.tree.set(i, "instrument")
                       for i in self.tree.selection()} - {"...", "?", ""}
        if len(instruments) > 1:
            self._write(self.t("mixed") % ", ".join(sorted(instruments)), "warn")
        argv = build_command(state)
        self._write("\n%s\n" % (" ".join(argv)), "value")
        env = dict(os.environ, PCA2D_COLOUR="1", PYTHONUNBUFFERED="1")
        try:
            self.proc = subprocess.Popen(
                argv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1, env=env,
                cwd=os.path.dirname(os.path.abspath(self.vars["config"].get()))
                or None)
        except OSError as exc:
            self._write("%s\n" % exc, "error")
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
        self.lines.put("\n[%d]\n" % code)
        self.proc = None
        self.root.after(0, self._finished)

    def _finished(self):
        self.run_button.configure(state="normal")
        self.stop_button.configure(state="disabled")
        self.status.configure(text=self.t("idle"))

    def stop(self):
        if self.proc is not None:
            self.proc.terminate()

    def _drain(self):
        while True:
            try:
                line = self.lines.get_nowait()
            except queue.Empty:
                break
            self._write(line)
        self.root.after(80, self._drain)

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
        self._write("%s\n" % path, "value")
        self.variant_box.configure(values=self._variants())

    def save_log(self):
        from tkinter import filedialog
        path = filedialog.asksaveasfilename(title=self.t("savelog"),
                                            defaultextension=".log")
        if path:
            with open(path, "w") as handle:
                handle.write(self.log.get("1.0", "end"))
            self._write("%s\n" % path, "value")

    def open_outputs(self):
        """Show the output root in the file browser: a run leaves one PDF and a
        folder of corrected spectra, and they are easier to find by looking."""
        root = self.vars["out_dir"].get() or (
            (self._config().get("output") or {}).get("directory") or "outputs")
        path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(
            self.vars["config"].get())), root))
        if not os.path.isdir(path):
            self._write("%s\n" % path, "warn")
            return
        opener = ("open" if sys.platform == "darwin"
                  else "explorer" if os.name == "nt" else "xdg-open")
        subprocess.Popen([opener, path])

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
