#!/usr/bin/env python
"""Un objet LBL de plus sur des fichiers qui existent deja, avec le template
d'un autre objet ou le sien. Pour la matrice des templates de la nuit du
2026-10-17: {livre, corrige} x {template de LBL, template du fit}, sans
recrire un seul spectre.

A lancer par le python de LBL, jamais celui de pca2d (lbl.environment):

    /opt/anaconda3/envs/lbl-rapide/bin/python tools/lbl_extra.py \
        --run-script outputs/GJ1_N17/1-3/run_lbl.py \
        --name GJ1_N17_DFT --files 'data/GJ1_N17/*t.fits' \
        --template-of GJ1_N17_PCA2D_1-3

Les reglages sont ceux du dictionnaire AFTER du script LBL que le run a ecrit
(meme instrument, meme arbre, meme lbl_config.yaml, meme Teff); seuls les
noms d'objet changent. --template-of own laisse LBL construire le template
des fichiers eux-memes; un nom d'objet lui fait prendre le template et le
masque de cet objet (OBJECT_COMPARISON). Les fichiers sont lies, pas copies,
dans <DATA_DIR>/science/<name>/. Un dossier deja plein n'est pas touche:
relancer reprend la ou LBL s'etait arrete.
"""

from __future__ import annotations

import argparse
import glob
import os
import runpy
import sys
import time

# LE DOSSIER DE DONNEES lbl/ DU DEPOT N'EST PAS LE PAQUET LBL. Le PYTHONPATH
# du profil finit par ':', ce qui met le dossier courant dans sys.path; lance
# depuis la racine du depot, `import lbl` trouve alors l'arbre de donnees
# comme paquet-espace-de-noms et echoue sur lbl.science (vu le 2026-10-16).
# Tout dossier qui contient un lbl/ sans __init__.py est retire du chemin.
sys.path[:] = [entry for entry in sys.path
               if not (os.path.isdir(os.path.join(entry or os.getcwd(), "lbl"))
                       and not os.path.exists(os.path.join(entry or os.getcwd(),
                                                           "lbl", "__init__.py")))]


def stamp(message):
    print("%s | %s" % (time.strftime("%y%m%d %H:%M:%S"), message), flush=True)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-script", required=True,
                   help="le run_lbl.py ecrit par le run de reference")
    p.add_argument("--name", required=True, help="le nouvel objet LBL")
    p.add_argument("--files", required=True,
                   help="un motif glob, entre guillemets, des fichiers a mesurer")
    p.add_argument("--template-of", required=True,
                   help="'own', ou l'objet LBL dont prendre template et masque")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    # executer le script SANS son __main__: on veut ses dictionnaires
    space = runpy.run_path(args.run_script, run_name="lbl_extra")
    params = dict(space["AFTER"])
    files = sorted(glob.glob(args.files))
    if not files:
        raise SystemExit("aucun fichier ne repond a %s" % args.files)
    comparison = args.name if args.template_of == "own" else args.template_of
    params["OBJECT_SCIENCE"] = [args.name]
    params["OBJECT_COMPARISON"] = [comparison]
    # les tables STRPCA decrivent le fit de l'objet corrige, pas celui-ci
    params["RESPROJ_TABLES"] = {k: v for k, v in
                                (params.get("RESPROJ_TABLES") or {}).items()
                                if not k.startswith("STRPCA")}
    params["REPORT"] = False
    science = os.path.join(params["DATA_DIR"], "science", args.name)
    stamp("LBL object %s: %d files from %s, template of %s, tree %s"
          % (args.name, len(files), args.files, comparison, params["DATA_DIR"]))
    if args.dry_run:
        stamp("dry run: nothing linked, nothing run")
        return None
    os.makedirs(science, exist_ok=True)
    present = set(os.listdir(science))
    made = 0
    for path in files:
        base = os.path.basename(path)
        if base in present:
            continue
        os.symlink(os.path.abspath(path), os.path.join(science, base))
        made += 1
    stamp("%d links made, %d already there" % (made, len(files) - made))
    from lbl.recipes import lbl_wrap
    lbl_wrap.main(params)
    stamp("LBL done for %s" % args.name)
    return None


if __name__ == "__main__":
    main()
